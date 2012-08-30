
library(devtools)
library(testthat)
library(roxygen2)

if (0) {
roxygenize("../../a4a/packages/FLa4a")

pkg <- as.package("../../a4a/packages/FLa4a")
build(pkg)
#check(pkg)
install(pkg)
}


roxygenize("../StockStructure")
pkg <- as.package("../StockStructure")
build(pkg)
#check(pkg)
install(pkg)
rm(list = ls())

#####################################################################
# first runs for testing
# using a4a assessment model
#--------------------------------------------------------------------
# year definitions: 
#	aLag - assessment lag in years (user set)
#	lastPyr - last year for projections (user set)
#	iniPyr - first year for projections (user set)
#	nPyr - number of year to project (computed)
#	aYr - assessment year (computed)
#	iYr - i year in the loop, it coincides with the intermediate year
#			in the usual ICES settings (loop) 
#	dtYr - last year with data which is the last year for which there 
#			are estimates (loop)
#	advYr - advice year, the year for which advice is being given (loop)     
#--------------------------------------------------------------------
#####################################################################
options(width = 150)
library(StockStructure)
data(wklife.stk)

library(FLa4a)
data(wklife.brp)

#====================================================================
# Simulation settings
#====================================================================

nits     <- 10                   # number of iterations
iniyr    <- 2000                 # first year in projections
npyr     <- 10                    # number of years to project
lastyr   <- iniyr + npyr         # last year in projections - note need one 
                                 # extra year of data for predictions
srsd     <- 0.3 			           # sd for S/R
units    <- 2                    # number of stock units
nhyr     <- 18                   # number of historical years
max.age  <- 8                    # plus group age
survey.q <- 1e-6 * exp(-2 * 0:7) # survey catchability at age
CV       <- 0.15                 # variability of catch.n and index observations
stock.id <- 7                    # which wklife stock to use

#====================================================================
# Use the stock history from one of the WKLife stocks
#====================================================================


#--------------------------------------------------------------------
# True stock msy reference points
#--------------------------------------------------------------------
refpt <- wklife.brp[[stock.id]] @ refpts["msy", c("ssb", "harvest")]

#--------------------------------------------------------------------
# True stock history
#--------------------------------------------------------------------
start.yr <- 30
pop1 <- window(wklife.stk[[stock.id]], start = start.yr, end = start.yr + nhyr - 1)[1:max.age,]
pop2 <- window(wklife.stk[[stock.id]], start = start.yr, end = start.yr + nhyr - 1)[1:max.age,]
dimnames(pop1) <- dimnames(pop2) <- list(year = 2001 - nhyr:1)

true.stock <- pop1
true.stock <- FLCore::expand(true.stock, unit = c("unique", "pop1", "pop2"))[,,2:3]

true.stock[,,1] <- pop1
true.stock[,,2] <- pop2

# change fbar here also ...
range(true.stock)[c("minfbar", "maxfbar")] <- c(2,5)
range(true.stock)["plusgroup"] <- max.age + 1 # to deal sepVPA and a4aFit


#====================================================================
# Stock recruit model
#====================================================================

# assume a stock recruitment model for each population and estimate 
# the parameters from it - this serves as the underlying
# stock recruitment model
sr.model1 <- as.FLSR(true.stock[,,1], model = "ricker")
sr.model1 <- fmle(sr.model1, control = list(trace = 0))

sr.model2 <- as.FLSR(true.stock[,,2], model = "ricker")
sr.model2 <- fmle(sr.model2, control = list(trace = 0))

# Residuals - simulate residuals as lognormal with sd=srsd
sr.residuals <- FLQuant(rlnorm(npyr * nits, sd = srsd), 
		                    dimnames = list(age  = 1, 
				                            year = dims(sr.model1)$minyear:lastyr, 
										                unit = c("pop1", "pop2"),
										                iter = 1:nits
						                   )) 
								
#--------------------------------------------------------------------
# create OM object
# Note: stocks are projected at Fsq and the values for the first
#	intermediate year are used in the projections
# TODO need to check this for my implentation... 
#--------------------------------------------------------------------

# fwdWindow extends filling with the contents of FLBRP object
# default in FLBRP is to use the mean of the last three years
OM <- fwdWindow(true.stock, FLBRP(true.stock), end = lastyr)

# propagate
OM <- propagate(OM, nits)

# project to the end of projections at last year F level
ctrl <- fwdControl( data.frame(year = iniyr:lastyr, quantity = "f", val = .15) )
# if no sr residuals the no noise!
OM[,,1] <- fwd(OM[,,1], ctrl = ctrl, sr = sr.model1)
OM[,,2] <- fwd(OM[,,2], ctrl = ctrl, sr = sr.model2)


#====================================================================
# first simulation
#====================================================================

base <- mse(OM = OM, iniyr = iniyr, 
		        sr.model1 = sr.model1, sr.model2 = sr.model2,
            sr.residuals = sr.residuals, 
		        Btrig = 0.5, Ftar = 0.75, refpt = refpt,
			      seed = 12345)

    
myrec <- function (object) stock.n(object)[1,]
        
plotComp(base[,,1], 
         fn = list(SSB = ssb, Recruits = myrec, Yield = catch, F = fbar), 
         probs = c(0.75, 0.5, 0.25), 
         size= c(0.5, 1, 0.5), 
         lty = c(2,1,2), 
         facet = facet_wrap(~qname, scale = "free"))


dev.new()
assessment <- attr(base, "assessment.data") $ stock
plotComp(assessment, 
    fn = list(SSB = ssb, Recruits = myrec, Yield = catch, F = fbar), 
    probs = c(0.75, 0.5, 0.25), 
    size= c(0.5, 1, 0.5), 
    lty = c(2,1,2), 
    facet = facet_wrap(~qname, scale = "free"))


#attr(base, "summaries")
