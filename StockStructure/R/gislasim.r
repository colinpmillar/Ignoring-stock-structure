#==============================================================================
# gislasim - cleaned version
#==============================================================================

#' Methods to generate lifehistory params
#'
#' @param linf the output from the mse function
#'
#' @return an FLPar
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname gislasim-methods
#'
setGeneric("gislasim", function(linf, ...) standardGeneric("gislasim"))

#' @rdname gislasim-methods
#' @aliases gislasim, numeric
setMethod("gislasim", c("numeric"), function (linf, t0 = -0.1, a = 1e-05, b = 3, ato95 = 1, sl = 2, sr = 5000, s = 0.9, v = 1000, asym=1, bg=b, iter=1, k="missing", M1="missing", M2="missing", a50="missing", a1="missing"){
      if(missing(k))  k <- 3.15 * linf^(-0.64)
      if(missing(M1)) M1 <- 0.55 + 1.44 * log(linf) + log(k) 
      if(missing(M2)) M2 <- -1.61
      if(missing(a50)) a50 <- FLAdvice:::invVonB(FLPar(linf=linf, t0=t0, k=k), 0.72 * linf^0.93)
      if(missing(a1)) a1 <- a50
      par <- FLPar(linf=linf, k=k, t0 = t0, a = a, b = b, asym=asym, bg=bg, sl=sl, sr=sr, s=s, v=v, M1=M1, M2=M2, ato95 = ato95, a50=a50, a1=a1, iter=iter)
      attributes(par)$units = c("cm", "kg", "1000s")
      return(par)
    })

#' @rdname gislasim-methods
#' @aliases gislasim, FLPar
setMethod("gislasim", c("FLPar"), function (linf){
      # Renaming to avoid confusing the argument with the object.
      # linf here is an FLPar object that can contain several parameters 
      object <- linf
      rm(linf)
      # now the real thing
      v0 <- dimnames(object)$params	    
      if(!("linf" %in% v0)) stop("The function requires linf.")
      par <- FLPar(c(linf=NA, t0 = -0.1, a = 1e-05, b = 3, ato95 = 1, sl = 2, sr = 5000, s = 0.9, v = 1000, asym=1, bg=3, k=NA, M1=NA, M2=NA, a50=NA, a1=NA), iter=ncol(object))
      dimnames(par)$iter <- dimnames(object)$iter 
      par[dimnames(object)$params] <- object
      if(!("bg" %in% v0)) par["bg"] = par["b"]
      if(!("k" %in% v0)) par["k"] = 3.15 * par["linf"]^(-0.64)
      if(!("M1" %in% v0)) par["M1"] = 0.55 + 1.44 * log(par["linf"]) + log(par["k"])
      if(!("M2" %in% v0)) par["M2"] = -1.61
      if(!("a50" %in% v0)) par["a50"] = FLAdvice:::invVonB(FLPar(linf=par["linf"], t0=par["t0"], k=par["k"]), c(0.72 * par["linf"]^0.93))
      if(!("a1" %in% v0)) par["a1"] = par["a50"]
      attributes(par)$units = c("cm", "kg", "1000s")
      return(par)
    })
