###############################################################################
# EJ&CO(20120416)
# Generate data from WKLIFE life history parameters list.
# NOTE #1: The fishery is simulated to follow 3 distinct patterns
# ToDo: run with other selectivities, fix maximum age, run with other S/R,
#	introduce variability, 
###############################################################################

library(FLAdvice)
library(FLAssess)
library(Matrix)
library(numDeriv)
source("funs.R")
source("../../dev/a4a9.linear.R")

#==============================================================================
# Read data and massage
#==============================================================================

# read data
wklifeLst <- read.table("allStockslifeHistoryParam.txt", sep="\t", head=TRUE)
wklifeLst <-transform(wklifeLst, value=as.numeric(as.character(value)))
wklifeLst <- wklifeLst[!is.na(wklifeLst$value),]
# remove "fle-2232" the parameters are duplicated and seem inconsistent
wklifeLst <- subset(wklifeLst, stock != "fle-2232")

#==============================================================================
# simulate
#==============================================================================
#------------------------------------------------------------------------------
# get LH pars 
#------------------------------------------------------------------------------
set.seed(123)
wklife.brp <- lapply(split(wklifeLst, wklifeLst$stock), function(x){
    cat(as.character(x$stock)[1], "\n")
    # get parameters
    par <- FLPar(x$value, tolower(x$param))
    if(!("linf" %in% dimnames(par)$params) & "lmax" %in% dimnames(par)$params){
        dimnames(par)$params[dimnames(par)$params == "lmax"] <- "linf"
    } 
    if("linf" %in% dimnames(par)$params){
      # complete with gislasim
      dnms <- dimnames(par)$params
      par <- par[dnms %in% dimnames(gislasim(0))$params]
      par <- gislasim(par)
      # run LH
	  if("tmax" %in% x$param){
		  tmax <- x[x$param=="tmax", "value"]
	      res <- lh(par, range=c(min=1, max=tmax, minfbar=1, maxfbar=tmax, plusgroup=tmax))    
	  } else {
	      res <- lh(par)    
	  } 
      res@desc <- as.character(x$stock[1])
      res
    } else {
      NULL
    }
})

wklife.brp <- wklife.brp[!unlist(lapply(wklife.brp, is.null))]

#------------------------------------------------------------------------------
# simulate F trajectory
#------------------------------------------------------------------------------

wklife.stk <- lapply(wklife.brp, function(x){
  cat(x@desc, "\n")
  Fc <- c(refpts(x)["crash","harvest"]*0.7)
  if(!is.na(Fc)){	
    Fmsy <- c(refpts(x)["msy","harvest"])
    Ftrg <- c(seq(0, Fc, len=19), rep(Fc, 20), seq(Fc, Fmsy, len=10))
    trg <- fwdControl(data.frame(year=c(2:50), quantity=rep('f', 49), val=Ftrg))
    ex.stk <- as(x, "FLStock")[,1:50]
    #ex.sr <- as.FLSR(ex.stk, model=x@model, params=x@params)
    ex.sr <- fmle(as.FLSR(ex.stk, model="bevholt"), control=list(trace=0))
    fwd(ex.stk, ctrl=trg, sr=ex.sr)
  } else {
    NULL
  }
})

wklife.stk <- wklife.stk[!unlist(lapply(wklife.stk, is.null))]
wklife.stk <- FLStocks(wklife.stk)
plot(wklife.stk)

#------------------------------------------------------------------------------
# run model
#------------------------------------------------------------------------------

lst <- list()
length(lst) <- length(wklife.stk)

for(i in 1:length(wklife.stk)){
	try(obj <- a4aLinearFit(x, as(x[1:5], "FLIndex")), TRUE)
	
	save(lst[[i]], file=deparse(i))
}





#------------------------------------------------------------------------------
# tests
#------------------------------------------------------------------------------

#ex.brp <- wklife.brp[[1]]
#ex.stk <- as(ex.brp, "FLStock")[,1:50]
#ex.sr <- as.FLSR(ex.stk, model="bevholt", params=ex.brp@params)

#Fc.8 <- c(refpts(ex.brp)["crash","harvest"]*0.8)
#Fmsy <- c(refpts(ex.brp)["msy","harvest"])
#Ftrg <- c(seq(0, Fc.8, len=19), rep(Fc.8, 20), seq(Fc.8, Fmsy, len=10))

#trg <- fwdControl(data.frame(year=c(2:50), quantity=rep('f', 49), val=Ftrg))

#stk <- fwd(ex.stk, ctrl=trg, sr=ex.sr)

#trg# developing fishery
#Bvirgin <- 1000
#iniF <- 0.1

#nextF <- rlnorm(1, 1, 1)

#library(FLAdvice)
#ex.brp <- lh(gislasim(FLPar(linf=100)))
#Fc.8 <- c(refpts(ex.brp)["crash","harvest"]*0.8)
#Fmsy <- c(refpts(ex.brp)["msy","harvest"])
#Ftrg <- c(seq(0, Fc.8, len=19), rep(Fc.8, 20), seq(Fc.8, Fmsy, len=10))
#trg <- fwdControl(data.frame(year=c(2:50), quantity=rep('f', 49), val=Ftrg))
#ex.stk <- as(ex.brp, "FLStock")[,1:50]
#ex.sr <- as.FLSR(ex.stk, model="bevholt")
#ex.sr@params <- ex.brp@params
#stk <- fwd(ex.stk, ctrl=trg, sr=ex.sr)
#plot(stk)






