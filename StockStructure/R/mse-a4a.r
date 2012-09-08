#' runs an simulation of a stock, stock assessment using a simulated index 
#' and stock recruit function.
#'
#' The stock is projected using the function \code{fwd} from the FLash library
#' the assessment model is the a4a model from the FLa4a package
#' while the management procedure is a Harvest Control Rule (HCR)
#' implemented using the function \code{hcr} from the FLAdvice package.
#' 
#' @param OM an FLStock object
#' @param start the first year of projections 
#' @return an FLStock object
#' @note \code{mse} is based on code written by 
#'       Ernesto Jardim \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @seealso \code{\link{FLStock}}, \code{\link{fwd}}, \code{\link{hcr}}
#' @export
mse <- function (OM, iniyr, sr.model1, sr.model2,
                 survey.q = rep(1e6, dims(OM) $ age),
                 sr.residuals = FLQuant(1, dimnames = dimnames(window(rec(OM), start = iniyr))), 
                 CV = 0.15, Ftar = 1, Btrig = 1, refpt,
                 seed = NULL, ...) {
  
  #--------------------------------------------------------------------
  # set year's info  
  #	current.year - i year in the loop, it coincides with the intermediate year
  #			in the usual ICES settings (loop) 
  #	data.year - last year with data which is the last year for which there 
  #			are estimates (loop)
  #	advice.year - advice year, the year for which advice is being given (loop)     
  #--------------------------------------------------------------------
  
  time0 <- proc.time()
  
  #--------------------------------------------------------------------
  # general settings
  #--------------------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)
  # final year of OM for prediction only
  data.years       <- with(dims(OM), seq(minyear, maxyear - 1))
  assessment.years <- seq(iniyr, dims(OM) $ maxyear - 1)
  estimate.years <- with(dims(OM), seq(minyear, maxyear - 2))

  #--------------------------------------------------------------
  # Get stock assessment data
  #--------------------------------------------------------------
  
  # weights at age are weighted means from each population - 
  # presumably estimated each year from a survey
  assessment.stock  <- collapseUnit(OM)[,ac(data.years)]
  # index = q * N1 + q * N2 = q * (N1 + N2)
  assessment.index  <- FLIndex(index = unitSums( stock.n(OM)[,ac(data.years)] ) * survey.q)
  
  #--------------------------------------------------------------------
  # introduce OEM on historical data (adding to all years but
  # assessment years are overwritten anyway
  #--------------------------------------------------------------------
  index(assessment.index) <- 
    index(assessment.index) * rlnorm( prod(dim(index(assessment.index))), 0, CV)  
  catch.n(assessment.stock) <- 
    catch.n(assessment.stock) * rlnorm( prod(dim(m(assessment.stock))), 0, CV)
  catch(assessment.stock) <- computeCatch(assessment.stock)
  
  #--------------------------------------------------------------------
  # objects to store TACs, catch, fbar and observations of these ...
  # might put them all into one
  #--------------------------------------------------------------------
  OM.summaries <- catch(OM)[,ac(data.years)]
  MP.summaries <- catch(assessment.stock)[,ac(data.years)]
  # hack to overcome exported method by reshape
  OM.summaries <- FLCore::expand(OM.summaries, age=c("all", "catch", "fbar", "forecast.f"))
  MP.summaries <- FLCore::expand(MP.summaries, age=c("all", "catch", "fbar", "forecast.f"))
  dimnames(OM.summaries)[[1]][1] <- "TAC" 
  dimnames(MP.summaries)[[1]][1] <- "TAC"
  
  # put in some values
  OM.summaries["fbar"] <- fbar(OM)[,ac(data.years)] 
  OM.summaries["catch"] <- catch(OM)[,ac(data.years)]
  MP.summaries["catch"] <- catch(assessment.stock)
  
  
  #--------------------------------------------------------------------
  # Go !!
  #--------------------------------------------------------------------
  for (current.year in assessment.years) {
    cat("===============", current.year, "of", max(assessment.years), "===============\n")
    
    #current.year <- assessment.years[1]
    #current.year <- current.year + 1
    data.year    <- ac(current.year - 1)
    advice.year  <- ac(current.year + 1)

    #--------------------------------------------------------------
    # Observation of data year catch and index
    # build up index and catch data year by year
    #--------------------------------------------------------------
  
    assessment.stock[,data.year] <- collapseUnit( OM[,data.year] )
    catch.n(assessment.stock)[,data.year] <- 
        catch.n(assessment.stock)[,data.year] * 
        rlnorm(prod(dim(m(assessment.stock))[-2]), 0, CV)
    
    index(assessment.index)[,data.year] <- 
        unitSums( stock.n(OM)[,data.year] ) * 
        survey.q * 
        rlnorm(prod(dim(index(assessment.index))[-2]), 0, CV)
    
    # the stock assessment data!
    current.stock <- window(assessment.stock, end = an(data.year))
    current.index  <- window(assessment.index, end = an(data.year))
     
    #--------------------------------------------------------------
    # management procedure
    #--------------------------------------------------------------
 
    # run assessments
    year.df <- floor(dims(current.stock) $ year * 0.4)
    for (i in 1:dims(current.stock) $ iter) {
      cat("      Iteration", i, "of", dims(current.stock) $ iter, "\n")
      out <- 
         try(a4aFit(current.stock[,,,,,i], current.index[,,,,,i], 
                    f.model = ~ bs(age, 3) + bs(year, year.df),
                    q.model = ~ bs(survey.age, 3),
                    r.model = ~ factor(cohort),
                    control = list(trace = 0, do.fit = TRUE),
                    ...)
            )
      if (class(out) != "try-error") {
        attr(out, "env") <- NULL
        current.stock[,,,,,i] <- out
        cat("fbar: " ,round(c(fbar(current.stock[,,,,,i])[,data.year]), 3),
         "   SSB: "  ,round(c(ssb(current.stock[,,,,,i])[,data.year]), 0),
            "\n")
      } else {
        # drop the iteration!?
      }
    }
    
    # save estimated F and N
    stock.n(assessment.stock[,data.year]) <- stock.n(current.stock[,data.year])
    harvest(assessment.stock[,data.year]) <- harvest(current.stock[,data.year])
    
    # save estimated fbar
    MP.summaries[c("fbar"), ac(current.year)] <- fbar(current.stock[,data.year]) 

    # Harvest Control Rule (HCR)
    # what is the forecast F in advice.year
    forecast.f <- suppressWarnings( hcr(ssb(current.stock[,data.year]), refpt) )

    # save forecast f
    MP.summaries["forecast.f", ac(current.year)] <- forecast.f
    
    # get recruitments for forecasts - GM of recruitment history
    recruitments <- exp( apply(log(rec(current.stock)), 6, mean) )
    recruitments <- FLQuant(rep(c(recruitments), each = 2), 
                            dimnames = list(age  = 1, 
                                            year = c(ac(current.year), advice.year),
                                            iter = seq(1, dims(current.stock) $ iter)
                            ))
    
    # add space for forecast - try capture.output to stock screen output...
    current.stock <- fwdWindow(current.stock, FLBRP(current.stock), end = an(advice.year))
    
    # set up forecast control - TAC is taken in intermediate year, Fmsy in forecast year
    trgtArray <- array(NA, 
                       dim = c(2, 3, dims(current.stock) $ iter), 
                       dimnames = list(1:2, c("min", "val", "max"), 
                                       iter = seq(1, dims(current.stock) $ iter)))
    trgtArray[1,2,] <- MP.summaries["TAC", data.year] # last years TAC
    trgtArray[2,2,] <- c(forecast.f)
    
    ctrl <- fwdControl( data.frame(year = an(c(current.year, advice.year)), quantity = c("catch", "f")), trgtArray = trgtArray)
    
    # what selection pattern is being used?
    current.stock <- fwd(current.stock, ctrl = ctrl, sr = list(model = "geomean", params = FLPar(a = 1)), sr.residuals = recruitments)

    tac <- catch(current.stock)[,ac(advice.year)]

    # update TAC record, all advYrs get the same TAC, may need more options
    MP.summaries["TAC", ac(current.year)] <- tac
      
    #--------------------------------------------------------------
    # project OM
    #--------------------------------------------------------------
    
    assign("tmp", OM, envir = .GlobalEnv)

    # first we need to find out what the F gives the TAC
    
    sel <- sweep(harvest(OM), 2:6, fbar(OM), "/")[,ac(data.year)]
    m <- m(OM)[, ac(advice.year)]
    n <- stock.n(OM)[, ac(advice.year)]
    wt <- catch.wt(OM)[, ac(advice.year)]
    
    findMult <- function(loga, iter = 1) {
      a <- exp(loga)
      z <- a * sel[,,,,,iter] + m[,,,,,iter]
      abs(tac[,,,,,iter] - 
          sum(a * sel[,,,,,iter] / z * (1 - exp(- z)) * 
                  n[,,,,,iter] * wt[,,,,,iter]
             )
         )
    }
    
    
    OM.f <- rep(NA, dims(OM) $ iter)
    for (i in 1:dims(OM) $ iter) 
      OM.f[i] <- exp( optimize(findMult, c(-10,10), iter = i) $ minimum )
 
    # set up forecast control - TAC is taken in intermediate year, Fmsy in forecast year
    trgtArray <- array(NA, 
        dim = c(1, 3, dims(OM) $ iter), 
        dimnames = list(1, c("min", "val", "max"), 
                        iter = seq(1, dims(OM) $ iter)))
    trgtArray[,2,] <- OM.f
    
    
    ctrl <- fwdControl( data.frame(year = an(advice.year), quantity = "f"),
                        trgtArray = trgtArray)

    OM[,,1] <- fwd(OM[,,1], ctrl = ctrl, sr = sr.model1, sr.residuals = sr.residuals[,advice.year,1])
    OM[,,2] <- fwd(OM[,,2], ctrl = ctrl, sr = sr.model2, sr.residuals = sr.residuals[,advice.year,2])
    
    
    # check this should be pretty much zero so we have found the f that results
    # in the TAC :)
    # unitSums(catch(OM)[,advice.year]) - tac
    OM.summaries["TAC"       , ac(current.year)] <- catch(OM)[,ac(advice.year)]
    OM.summaries["catch"     , ac(current.year)] <- NA # catch(OM)[,ac(current.year)]
    OM.summaries["fbar"      , ac(current.year)] <- fbar(OM)[,ac(advice.year)]
    OM.summaries["forecast.f", ac(current.year)] <- NA # OM.f
    
  } # end year loop
  
  harvest(assessment.stock)[,ac(seq(dims(OM) $ minyear, iniyr - 2))] <- NA # take out non-estimated data
  stock.n(assessment.stock)[,ac(seq(dims(OM) $ minyear, iniyr - 2))] <- NA # take out non-estimated data
  
  catch(assessment.stock) <- computeCatch(assessment.stock)
  attr(OM, "summaries") <- list(OM = OM.summaries, MP = MP.summaries)
  attr(OM, "assessment.data") <- list(stock = assessment.stock[,ac(estimate.years)], 
                                      index = index(assessment.index)[,ac(estimate.years)])
  
  time0 <- c(proc.time() - time0)[3]                                
  cat("\ntotal time:", floor(time0/60/60), "h", floor((time0 - floor(time0/60/60))/60), time0 - floor(time0/60), "s\n\n")
                                  
  return(OM[,ac(data.years)])
}





