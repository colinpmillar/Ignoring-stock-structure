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
#' @seealso \code{\link{FLStock}}, \code{\link{\fwd}}, \code{\link{hcr}}
#' @export
mse <- function (OM, iniyr, sr.model1, sr.model2,
                 survey.q = rep(1e6, dims(OM) $ age),
                 sr.residuals = FLQuant(1, dimnames = dimnames(window(rec(OM), start = iniyr))), 
                 CV = 0.15, Ftar = 0.75, Btrig = 0.75, refpt) {
  
  #--------------------------------------------------------------------
  # set year's info  
  #	lastPyr - operating model maxyear
  #	iniPyr - start - first year for projections (user set)
  #	nPyr - number of year to project (computed)
  #	aYr - assessment year (computed)
  #	iYr - i year in the loop, it coincides with the intermediate year
  #			in the usual ICES settings (loop) 
  #	dtYr - last year with data which is the last year for which there 
  #			are estimates (loop)
  #	advYr - advice year, the year for which advice is being given (loop)     
  #--------------------------------------------------------------------
  
  Ftar <- 0.75
  Btrig <- 0.75
  Fmin <- 0.3
  Blim <- 0.3
  
  
  #--------------------------------------------------------------------
  # general settings
  #--------------------------------------------------------------------
  #set.seed(seed)
  iniPyr <- iniyr
  lastPyr <- OM @ range["maxyear"]
  nPyr <- lastPyr - iniPyr + 1
  # assessment years
  aYrs <- seq(iniPyr, lastPyr)
  
  
  #--------------------------------------------------------------------
  # introduce OEM on historical data
  #--------------------------------------------------------------------
  index <- FLIndex(index = stock.n(OM) * survey.q)
  # abundance index variability
  #index(index) <- index(index) * rlnorm(prod(dim(index(index))), 0, CV)
  # catch observation error ... is this right?
  #filt <- ac( dims(OM) $ minyear:(iniyr - 2) )
  #catch.n(OM)[,filt] <- catch.n(OM)[,filt] * rlnorm(prod(dim(catch.n(OM)[,filt])), 0, CV)
  #catch(OM)[,filt] <- computeCatch(OM[,filt])
  
  #--------------------------------------------------------------------
  # objects to store TACs, catch, fbar and observations of these ...
  # might put them all into one
  #--------------------------------------------------------------------
  true.summary <- catch(OM)
  # hack to overcome exported method by reshape
  true.summary <- FLCore::expand(true.summary, age=c("all", "catch", "fbar"))
  dimnames(true.summary)[[1]][1] <- "TAC"	
  observed.summary <- true.summary
  # put in some values
  true.summary["fbar"] <- fbar(OM) 
  
  #--------------------------------------------------------------------
  # Go !!
  #--------------------------------------------------------------------
  for (initial.year in iniPyr:lastPyr) {
    cat("===================", initial.year, "===================\n")
    
    initial.year <- 2000
    data.year  <- ac(initial.year - 1) 
    
    #--------------------------------------------------------------
    # Observation error
    #--------------------------------------------------------------
  
    input.data  <- collapseUnit( window(OM, end = an(data.year)) )
    input.index  <- qapply(window(index, end = an(data.year)), unitSums)
     
    # the new catch observation
    catch.n(input.data)[, data.year] <- unitSums( catch.n(OM)[, data.year] * rlnorm(prod(dim(m(input.data)[, data.year])), 0, CV) )
    catch(input.data)[, data.year] <- computeCatch(input.data[, data.year])
    
    # the new survey observation
    index(input.index)[, data.year] <- unitSums(stock.n(OM)[, data.year]) * survey.q * rlnorm(prod(dim(m(input.data)[, data.year])), 0, CV)
    
    
    #--------------------------------------------------------------
    # management procedure
    #--------------------------------------------------------------
 
    # run assessments
    for (i in 1:dims(input.data) $ iter) {
      cat("      Iteration", i, "of", dims(input.data) $ iter, "\n")
      out <- a4aFit(input.data[,,,,,i], input.index[,,,,,i], 
                    f.model = ~ bs(age, 3) + bs(year, 8),
                    q.model = ~ bs(survey.age, 3),
                    r.model = ~ 1,
                    control = list(trace = 0, do.fit = TRUE))
      attr(out, "env") <- NULL
      input.data[,,,,,i] <- out
    }
    
    # save estimated fbar
    observed.summary[c("fbar"), ac(initial.year)] <- fbar(input.data[,ac(initial.year - 1)]) 

    # estimate SR relationship for forecast
    input.sr.model <- as.FLSR(input.data, model = "ricker")
    input.sr.model <- fmle(input.sr.model, control = list(trace = 0))
    
    # projections considering TAC was caught in initial.year
    # to deal with a missing feature in fwd the last year of data must also be included
    ct <- unitSums(true.summary["TAC", c(data.year, initial.year)])
    ct[, data.year] <- catch(input.data[, data.year])
    
    # add space for assessment year
    input.data <- fwdWindow(input.data, FLBRP(input.data, input.sr.model), end = initial.year)
    
    # forecast the assessment
    #ctrl <- fwdControl( data.frame(year = data.year:initial.year, quantity = "catch", val = ct) )
    input.data <- fwd(input.data, sr = input.sr.model, catch = ct)
 
    # Harvest Control Rule (HCR)
    # the lag on the hcr method is between the advice year and the data year, so it's 2
    forecast.f <- suppressWarnings( hcr(ssb(input.data[,ac(initial.year)]), refpt) )
 
    # check F is not above maxF and replace if it is
    # hv[hv>maxHR] <- maxHR
    #PARrec["HCR:hr", ac(iYr)] <- hv
    
    # add space for forecast year
    input.data <- fwdWindow(input.data, FLBRP(input.data, input.sr.model), end = initial.year + 1)
    
    # do forecast
    trgtArray <- array(NA, dim = c(1,3,10), dimnames = list(1, c("min", "val", "max"), iter = 1:10))
    trgtArray[,2,] <- c(forecast.f)
    
    ctrl <- fwdControl( data.frame(year = initial.year + 1, quantity = "f"), trgtArray = trgtArray)
    
    input.data <- fwd(input.data, ctrl = ctrl, sr = input.sr.model)

    tac <- catch(input.data)[,ac(initial.year + 1)]
    
    # update TAC record, all advYrs get the same TAC, may need more options
    true.summary["TAC",ac(initial.year + 1)] <- tac
      
    #--------------------------------------------------------------
    # fwd control to project OM
    #--------------------------------------------------------------
    ctrl <- 
      fwdControl(
        data.frame(year     = an(c(data.year, initial.year + rep(0:1, each = 2))),
                   quantity = c("catch", rep(c("catch", "f"), 2))
                   ))
    dms <- dimnames(ctrl @ trgtArray)
    dms $ iter <- 1:nits
    ctrl @ trgtArray                              <- array(NA, lapply(dms, length), dms)
    ctrl @ trgtArray[1, "val", ]                  <- catch(OM)[, ac(data.year)]
    ctrl @ trgtArray[2 * (1:(aLag + 1)), "val", ] <- PARrec["HCR:TAC", c(iYr, advYrs)]

    OM <- fwd(OM, ctrl = ctrl, sr = sr, sr.residuals = srRsdl) 

  } # end year loop
  
  attr(OM, "PARs") <- PARrec
  
  return(OM)
}





