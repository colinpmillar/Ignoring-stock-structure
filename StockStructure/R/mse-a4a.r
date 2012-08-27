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
mse <- function (OM, start, sr, 
                srresidual = FLQuant(1, dimnames = dimnames(window(rec(OM), start = start))), 
                CV = 0.15, Ftar = 0.75, Btrig = 0.75) {
  
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
  
  #--------------------------------------------------------------------
  # general settings
  #--------------------------------------------------------------------
  set.seed(seed)
  iniPyr <- start
  lastPyr <- OM @ range["maxyear"]
  nPyr <- lastPyr - iniPyr + 1
  # assessment years
  aYrs <- seq(iniPyr, lastPyr)
  
  #--------------------------------------------------------------------
  # introduce OEM on historical data
  #--------------------------------------------------------------------
  bd <- as(OM, "FLBioDym")
  # abundance index variability
  index(bd) <- index(bd) * rlnorm(prod(dim(index(bd))), 0, CV)
  
  #--------------------------------------------------------------------
  # Management procedure - add bounds for BioDym	
  #--------------------------------------------------------------------
  if (!is.null(bounds)) bd @ bounds <- bounds
  
  #--------------------------------------------------------------------
  # object to register TACs, start year catch taken from OM 
  #--------------------------------------------------------------------
  PARrec <- catch(OM)
  # hack to overcome exported method by reshape
  PARrec <- FLCore::expand(PARrec, age=c("all", "HCR:hr", paste("BioDym", dimnames(params(bd))$params, sep=":")))
  PARrec["HCR:hr"] <- catch(OM)/stock(OM) 
  dimnames(PARrec)[[1]][1] <- "HCR:TAC"	
  
  #--------------------------------------------------------------------
  # Go !!
  #--------------------------------------------------------------------
  for (iYr in iniPyr:lastPyr) {
    cat("===================", iYr, "===================\n")
    
    dtaYr  <- ac(iYr - 1)
    dtaYrs <- ac((iYr - aLag):(iYr - 1))
    advYrs <- ac((iYr + 1):(iYr + aLag)) 
    
    #--------------------------------------------------------------
    # Observation error
    #--------------------------------------------------------------
    bd <- window(bd, end = an(dtaYr))

    # abundance index variability
    index(bd)[, dtaYrs] <- stock(OM)[, dtaYrs] * rlnorm(prod(dim(index(bd)[, dtaYrs])), 0, CV)
 
    # catch
    catch(bd)[, dtaYrs] <- computeCatch(OM)[, dtaYrs]
    
    #--------------------------------------------------------------
    # management procedure
    #--------------------------------------------------------------
 
  # assessment
    bd <- admbBD(bd)
    PARrec[-c(1, 2), ac(iYr)] <- c(params(bd)) 
  
    # projections considering TAC was caught in iYr
    # to deal with a missing feature in fwd the last year of data must also be included
    ct <- PARrec["HCR:TAC", c(dtaYr, iYr)]
    ct[, dtaYr] <- catch(bd)[, dtaYr]
    bd <- fwd(bd, catch = ct)
 
    # Harvest Control Rule (HCR)
    # the lag on the hcr method is between the advice year and the data year, so it's 2
    hv <- hcr(bd, FLPar(Ftar = Ftar, Btrig = Btrig, Fmin = Fmin, Blim = Blim), lag = 2)
 
    # check F is not above maxF and replace if it is
    # hv[hv>maxHR] <- maxHR
    PARrec["HCR:hr", ac(iYr)] <- hv
    tac <- TAC(bd, hv) # this is just hv*b

    # update TAC record, all advYrs get the same TAC, may need more options
    PARrec["HCR:TAC",advYrs] <- tac[,rep(1, length(advYrs))]
      
      
    #--------------------------------------------------------------
    # fwd control to project OM
    #--------------------------------------------------------------
    ctrl <- 
      fwdControl(
        data.frame(year = an(c(dtaYr, iYr, iYr, rep(advYrs, 2))), 
                   max = NA, 
                   quantity = c("catch", rep(c("catch", "f"), 1))
                   ))
    dms <- dimnames(ctrl @ trgtArray)
    dms $ iter <- 1:nits
    ctrl @ trgtArray                              <- array(NA, lapply(dms, length), dms)
    ctrl @ trgtArray[1, "val", ]                  <- catch(OM)[, ac(dtaYr)]
    ctrl @ trgtArray[2 * (1:(aLag + 1)), "val", ] <- PARrec["HCR:TAC", c(iYr, advYrs)]

    OM <- fwd(OM, ctrl = ctrl, sr = sr, sr.residuals = srRsdl) 

  } # end year loop
  
  attr(OM, "PARs") <- PARrec
  
  return(OM)
}





