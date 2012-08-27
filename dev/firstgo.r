
library(devtools)
library(testthat)
library(roxygen2)

roxygenize("../StockStructure")
pkg <- as.package("../StockStructure")
build(pkg)
#check(pkg)
install(pkg)
rm(pkg)

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

library(StockStructure)
data(wklife.stk)

#====================================================================
# Simulation settings
#====================================================================

nits   <- 10                # number of iterations
iniyr  <- 2000              # first year in projections
npyr   <- 50                # number of years to project
lastyr <- iniyr + npyr - 1  # last year in projections
srsd   <- 0.3 			    # sd for S/R

#====================================================================
# Use the stock history from one of the WKLife stocks
#====================================================================

#--------------------------------------------------------------------
# True stock history, nits of them
#--------------------------------------------------------------------
true.stock <- window(wklife.stk[[1]], start = 30, end = 50)
true.stock <- setPlusGroup(true.stock, plusgroup = 8)
dimnames(true.stock) <- list(year = 1980 + 0:20)
# change fbar here also ...

#====================================================================
# Stock recruit model
#====================================================================

sr.model <- as.FLSR(true.stock, model = "ricker")
sr.model <- fmle(sr.model, control = list(trace = 0))

# Residuals - simulate residuals as lognormal with sd=srsd
sr.residuals <- FLQuant(rnorm(npyr * nits, sd = srsd), 
		                dimnames = list(age  = with(dims(true.stock), min:max), 
				                        year = dims(sr.model) $ minyear:lastyr, 
										iter = 1:nits
						                )) 
								
#--------------------------------------------------------------------
# create OM object
# Note: this object is projected at Fsq and the values for the first
#	intermediate year are used in the projections
#--------------------------------------------------------------------

# window with FLBRP expands the object including weights, etc, the 
# brp doesn't seem to do anything except dispatching. it replaces "stf". 
OM <- window(true.stock, FLBRP = FLBRP(true.stock, sr.model), end = lastyr)

# trick to get iterations, start with M and fwd will add to other slots
m(OM) <- propagate(m(OM), nits)

# project to the end of projections at last year F level
#ctrl <- fwdControl( data.frame(year = iniyr:lastyr, quantity = "f", val = 2) )
OM <- fwd(OM, sr = sr.model, sr.residuals = sr.residuals)

#====================================================================
# first simulation need to find optimum Btrig and Ftar
#====================================================================

base <- mse(OM=OM, start = iniyr, 
		    sr = sr.model, sr.residuals = srRsdl, 
		    Btrig = 0.5, Ftar = 0.75)

