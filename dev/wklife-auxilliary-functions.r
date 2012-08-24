###############################################################################
# EJ(20120413)
# Auxiliary functions for LH generation of datasets
# Parameters: 
#	a1, sr, sl = 50% selectivity age, right variance, left variance
###############################################################################

iterMedians <- function(x, ...){
	return(apply(x, c(1:5), median, na.rm = FALSE))
}

iterSums <- function(x, ...){
	return(apply(x, c(1:5), sum, na.rm = FALSE))
}

iterCv <- function(object, ...){
	sqrt(iterVars(object))/iterMeans(object)
}


#==============================================================================
# gislasim - cleaned version
#==============================================================================

setGeneric("gislasim", function(linf, ...) standardGeneric("gislasim"))

setMethod("gislasim", signature(linf="numeric"), function (linf, t0 = -0.1, a = 1e-05, b = 3, ato95 = 1, sl = 2, sr = 5000, s = 0.9, v = 1000, asym=1, bg=b, iter=1, k="missing", M1="missing", M2="missing", a50="missing", a1="missing"){
    if(missing(k))  k <- 3.15 * linf^(-0.64)
    if(missing(M1)) M1 <- 0.55 + 1.44 * log(linf) + log(k) 
    if(missing(M2)) M2 <- -1.61
    if(missing(a50)) a50 <- FLAdvice:::invVonB(FLPar(linf=linf, t0=t0, k=k), 0.72 * linf^0.93)
    if(missing(a1)) a1 <- a50
    par <- FLPar(linf=linf, k=k, t0 = t0, a = a, b = b, asym=asym, bg=bg, sl=sl, sr=sr, s=s, v=v, M1=M1, M2=M2, ato95 = ato95, a50=a50, a1=a1, iter=iter)
    attributes(par)$units = c("cm", "kg", "1000s")
    return(par)
})

setMethod("gislasim", signature(linf="FLPar"), function (linf){
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

##==============================================================================
## lh
##==============================================================================

#lh <- function(par, growth=vonB, fnM=function(par,len) exp(par["M1"]+par["M2"]*log(len)), fnMat=logistic, fnSel=dnormal, sr="bevholt", range=c(min=1, max=40, minfbar=1, maxfbar=40, plusgroup=40), spwn=0, fish=0.5, units=if("units" %in% names(attributes(par))) attributes(par)$units else NULL, ...){

#	# Check that m.spwn and harvest.spwn are 0 - 1
#	if (spwn > 1 | spwn < 0 | fish > 1 | fish < 0) 
#		stop("spwn and fish must be in the range 0 to 1\n")
#	if (("m.spwn" %in% names(args)))
#		m.spwn <- args[["m.spwn"]]
#	else 
#		m.spwn <- FLQuant(spwn, dimnames=list(age=range["min"]:range["max"]))
#	if (("harvest.spwn" %in% names(args)))
#		harvest.spwn <- args[["harvest.spwn"]]
#	else
#		harvest.spwn <- FLQuant(spwn, dimnames=list(age=range["min"]:range["max"]))
#	age <- propagate(FLQuant(range["min"]:range["max"],dimnames=list(age=range["min"]:range["max"])),length(dimnames(par)$iter))

#	# Get the lengths through different times of the year
#	stocklen   <- growth(par[c("linf","t0","k")],age+m.spwn)	# stocklen is length at spawning time
#	catchlen   <- growth(par[c("linf","t0","k")],age+fish)	# catchlen is length when fishing happens
#	midyearlen <- growth(par[c("linf","t0","k")],age+0.5)	# midyear length used for natural mortality
#	# Corresponding weights
#	swt <- par["a"]*stocklen^par["b"]
#	cwt <- par["a"]*catchlen^par["b"]
#	if ("bg" %in% dimnames(par)$param) swt <- par["a"]*stocklen^par["bg"]
#  	args <- list(...)
#	m. <- fnM(par=par,len=midyearlen)	# natural mortality is always based on mid year length
#	mat. <- fnMat(par,age + m.spwn)	# maturity is biological therefore + m.spwn
#	sel. <- fnSel(par,age + fish)	# selectivty is fishery  based therefore + fish

#	# create a FLBRP object to   calculate expected equilibrium values and ref pts
#	dms <- dimnames(m.)
#	res=FLBRP(stock.wt	=swt, 
#		landings.wt	=cwt, 
#		discards.wt	=cwt, 
#		bycatch.wt	=cwt, 
#		m		=m., 
#		mat		=mat., 
#		landings.sel	=sel., 
#		discards.sel	=FLQuant(0, dimnames=dms),
#		bycatch.harvest	=FLQuant(0, dimnames=dms), 
#		harvest.spwn	=harvest.spwn, 
#		m.spwn		=m.spwn, 
#		availability	=FLQuant(1, dimnames=dms), 
#		range		=range)


#	# replace any slot passed in as an arg
#	for(slt in names(args)[names(args) %in% names(getSlots("FLBRP"))[names(getSlots("FLBRP"))!="fbar"]])
#		slot(res, slt) <- args[[slt]]
#	params(res) <- propagate(params(res),dims(res)$iter)

#	# Stock recruitment relationship
#	model(res) <- do.call(sr,list())$model
#	if (dims(par)$iter>1) {
#		warning("Scarab, iters dont work for SRR:sv/ab etc")
#		params(res) <- FLPar(c(a=NA,b=NA),iter=dims(par)$iter)
#		for (i in seq(dims(par)$iter)){ 
#			params(res)[,i][] <- unlist(c(ab(par[c("s","v"),i], sr, spr0=iter(spr0(res),i))[c("a","b")]))
#		}
#		warning("iter(params(res),i)=ab(par[c(s,v),i],sr,spr0=iter(spr0(res),i))[c(a,b)] assignment doesnt work")
#		warning("iter(FLBRP,i) doesnt work")
#	} else {
#		params(res) <- ab(par[c("s","v")], sr, spr0=spr0(res))[c("a","b")]
#	}
#	dimnames(refpts(res))$refpt[5] <- "crash"
#	res <- brp(res)
#   
#	if ("fbar" %in% names(args)) 
#		fbar(res) <- args[["fbar"]] else 
#		if (any((!is.nan(refpts(res)["crash","harvest"])))) 
#      			fbar(res) <- FLQuant(seq(0,1,length.out=101))*refpts(res)["crash","harvest"]
#  
#	res <- brp(res)
#	#res <- setUnits(res, par)
#	return(res)
#}

##==============================================================================
## vonB & inVonB 
##==============================================================================

#setGeneric('invVonB', function(params,data, ...) standardGeneric('invVonB'))

#setMethod("invVonB", signature(params="FLPar",data="ANY"), function(params, data="missing",...){
#	if (!missing(data) & "FLQuant" %in% is(data) ||  "FLCohort" %in% is(data)) data <- ages(data)
#	res <- -log(1.0-(data/params["linf"]))/params["k"]+params["t0"]
#	res
#})

#setGeneric('vonB', function(params,data, ...) standardGeneric('vonB'))

#setMethod("vonB", signature(params="FLPar",data="ANY"), function(params, data="missing",...){
#	if (!missing(data) & "FLQuant" %in% is(data) || "FLCohort" %in% is(data)) data <- ages(data)
#	params["linf"]%*%(1.0-exp(-params["k"]%*%(data-params["t0"])))
#})

##==============================================================================
## ages 
##==============================================================================

#setGeneric('ages', function(data, ...) standardGeneric('ages'))

#setMethod("ages", signature(data="FLQuant"),function(data,timing=NULL){
#	res <- FLQuant(dimnames(data)$age,dimnames=dimnames(data))
#	if (is.null(timing))
#		res <- sweep(res,4,(1:dim(res)[4]-1)/dim(res)[4],"+") else
#		res <- sweep(res,4,timing,"+")
#	return(res)
#})

#setMethod("ages", signature(data="FLCohort"), function(data,timing=NULL){
#	res <- FLCohort(dimnames(data)$age,dimnames=dimnames(data))
#	if (is.null(timing))
#		res <- sweep(res,4,(1:dim(res)[4]-1)/dim(res)[4],"+") else
#         	res <- sweep(res,4,timing,"+")
#	return(res)
#})


##==============================================================================
## setUnits 
##==============================================================================

#setUnits <- function(res, par){
#    units <- attributes(par)$units
#    allUnits <- list("params"	= "",          
#               "refpts"		= "",            
#               "fbar"		= "",        
#               "fbar.obs"	= "",    
#               "landings.obs"	= paste(units[2],units[3]),
#               "discards.obs"	= paste(units[2],units[3]),
#               "rec.obs"	= units[3],         
#               "ssb.obs"	= paste(units[2],units[3]),
#               "stock.obs"	= paste(units[2],units[3]),
#               "profit.obs"	= "",     
#               "revenue.obs"	= "",    
#               "landings.sel"	= "",    
#               "discards.sel"	= "", 
#               "bycatch.harvest"="",        
#               "stock.wt"	= units[2],     
#               "landings.wt"	= units[2],     
#               "discards.wt"	= units[2],      
#               "bycatch.wt"	= units[2],               
#               "m"		= "",             
#               "mat"		= "proportion", 
#               "harvest.spwn"	= "proportion",          
#               "m.spwn"		= "proportion",    
#               "availability"	= "proportion",           
#               "price"		= "",           
#               "vcost"		= "",           
#               "fcost"		= "")            
#    units(res)[names(allUnits)] <- allUnits
#    return(res)
#}

##==============================================================================
## logistic, dnormal 
##==============================================================================

#logistic <- function(params,data) { #x, a50, ato95){
#	if ((params["a50"]-data)/params["ato95"] > 5) res <- 0 
#	if ((params["a50"]-data)/params["ato95"] < -5) res <- 0 else 
#		res <- params["asym"]/(1.0+19.0^((params["a50"]-data)/params["ato95"]))
#	return(res)
#}

#dnormal <- function(params,data){    
#    a1 <- FLQuant(1,dimnames=dimnames(data))%*%params["a1"]
#    s <- FLQuant(1,dimnames=dimnames(data))%*%params["sl"]
#    sr <- FLQuant(1,dimnames=dimnames(data))%*%params["sr"]
#    if (dims(data)$iter==1 &  dims(a1)$iter>1) data <- propagate(data,dims(a1)$iter)
#    s[data>=a1] <- sr[data>=a1]
#    2.0^(-((data-a1)/s*(data-a1)/s))
#}


