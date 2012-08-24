#####################################################################
# first runs for testing
# pelat: B[t+1] = B[t] + r/p*B[t]*(1-(B[t]/k)^p)) - C[t]
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

library(FLAdvice)
library(FLBioDym)
library(plyr)
source("mse.biomass.dynamic.r")

#====================================================================
# Simulation settings
#====================================================================

nits <- 150 			# number of iterations
iniyr <- 2011 			# first year in projections
lastyr <- 2061 			# last year in projections
npyr <- lastyr-iniyr+1 	# number of years to project
srsd <- 0.3 			# sd for S/R

#====================================================================
# Rebuild the stock history with simulated data
#====================================================================

#--------------------------------------------------------------------
# Sardinella data
#--------------------------------------------------------------------
nc <- read.table('../data/nc_sar.dat', header=TRUE)
ia <- read.table('../data/cpue_sar.dat', header=TRUE)
ca <- FLQuant(nc$catch, dimnames=list(age='all', year=nc$year))
cp <- FLQuant(c(rep(NA, 7), ia$cpue), dimnames=list(age='all', year=nc$year))

# CECAF SA results
sa <- read.table('../data/sa.dat', sep='\t', header=T)
saq <- FLQuants(HR=FLQuant(sa$hr, dimnames=list(age='all', year=sa$year)),
  SSB=FLQuant(sa$ssb, dimnames=list(age='all', year=sa$year)))

#--------------------------------------------------------------------
# gislasim and projection
#--------------------------------------------------------------------
# Use CECAF SA B0 as starting point
par <- as(data.frame(linf=35, k=0.4, t0=-0.1, s=0.8, v=1750, a50=1), 'FLPar')
sar <- lh(gislasim(par), range=c(min=1, max=8, minfbar=1, maxfbar=6))

# stk with initial F closest to estimated HR
saa <- as(sar, 'FLStock')
saa <- saa[,7]
dimnames(saa) <- list(year=1989)

# prepare for projection
saa <- window(saa, FLBRP=sar, end=iniyr-1)

trg <- fwdControl(data.frame(year=1990:2010, val=c(ca), quantity="catch"))

# optim for exceptional recruitment scaler 

# foo
foo <- function(m) {
  sr.residuals <- FLQuant(c(rep(1, 15), m[1], 1, m[2], 1, 1, 1), 
    dimnames=list(year=1990:2010))
  saar <- fwd(saa, trg, sr=list(model="bevholt", params=params(sar)),
    sr.residuals=sr.residuals)
  return(sum((ssb(saar)[,-1] - saq$SSB)^ 2))
}

# get scalers
m <- optim(c(1.5,3), foo)

# apply w/ iters
sr.residuals <- rlnorm(nits, log(FLQuant(c(rep(1, 15),
  m$par[1], 1, m$par[2], 1, 1, 1), dimnames=list(year=1990:2010))), srsd)

bwplot(data~year, sr.residuals)

# project
saa <- fwd(saa, trg, sr=list(model="bevholt", params=params(sar)),
  sr.residuals=sr.residuals)

xyplot(data~year|qname, as.data.frame(FLQuants(SA=saq$SSB, OM=iterMeans(ssb(saa)[,-1]))), type='b')

par(mfrow=c(2,1))
plot(sa$hr~sa$year, type="l", ylim=c(0,1), main="Harvest Rate")
lines(c(iter(catch(saa)/stock(saa), 1))[-1]~sa$year, col=2)
plot(sa$ssb~sa$year, type="l", ylim=c(0,2000), main="SSB")
lines(c(iter(ssb(saa),1))[-1]~sa$year, col=2)

saa <- propagate(saa, nits)

#====================================================================
# Conditioning
#====================================================================

om <- saa
srBH <- as.FLSR(om,model="bevholt")
params(srBH) <- params(sar)

# Residuals - simulate residuals as lognormal with sd=srsd
set.seed(123)
#srRsdl <- FLQuant(rlnorm(npyr*nits, 0, srsd), dimnames=list(year=srBH@range["minyear"]:lastyr, iter=1:nits)) 
srRsdl <- FLQuant(rlnorm(npyr*nits, 0, srsd) * (rpois(npyr*nits, 2/20) * 0.5 + 1), dimnames=list(year=2011:2061, iter=1:nits))

#--------------------------------------------------------------------
# create OM object
# Note: this object is projected at Fsq and the values for the first
#	intermediate year are used in the projections
#--------------------------------------------------------------------

# window with FLBRP expands the object including weights, etc, the 
# brp doesn't seem to do anything except dispatching. it replaces "stf". 
OM <- window(om, FLBRP=sar, end=lastyr)

# trick to get iterations, start with M and fwd will add to other slots
m(OM) <- propagate(m(OM), nits)

# project to the end of projections at last year F level
ctrl <- fwdControl(data.frame(year=iniyr:lastyr, quantity="f", val=c(2)))
OM <- fwd(OM, ctrl=ctrl, sr=srBH, sr.residuals=srRsdl)

#====================================================================
# MP settings
#====================================================================

# these bounds are related with ADMB adjusted to CECAF usage
bounds <- bounds(FLBioDym())
bounds["sigma","start"]=0.50
bounds["q",    "start"]=1.0
bounds["q",    1]      =1.0
bounds["b0",   c("phase","start")]=c(-1,0.5)
bounds["p",   c("phase","start")]=c(-1,1)
bounds[,"lower"]=bounds[,"start"]*0.1
bounds[,"upper"]=bounds[,"start"]*10.0
bounds["r",    c("phase","lower","upper","start")] = c(1, 0.32, 1.28, 0.64)
bounds["K",    c("phase","lower","upper","start")] = c(2, 1750, 7000, 2500)

#====================================================================
# Scenarios and simulations
#====================================================================

scn <- mkScn(ref="base") 

base <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=1, cthBias=1, IEM="linear")

#--------------------------------------------------------------------
# assessment lag effect
#--------------------------------------------------------------------

scn[2,] <- mkScn(ref="aLag3", aLag=3, runid=2) 
scn[3,] <- mkScn(ref="aLag5", aLag=5, runid=3) 

base.aLag3 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=3, srvBias=1, cthBias=1, IEM="linear")

base.aLag5 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=5, srvBias=1, cthBias=1, IEM="linear")

#--------------------------------------------------------------------
# underreporting effect
#--------------------------------------------------------------------

scn[4,] <- mkScn(ref="cthBias0.5", runid=4, cthBias=0.5) 

bounds["r",    c("phase","lower","upper","start")] = c(1, 0.1, 2, 0.64)
bounds["K",    c("phase","lower","upper","start")] = c(2, 1000, 7000, 2500)

base.cthBias0.5 <- mseBD2(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=1, cthBias=0.5, IEM="linear")

bounds["r",    c("phase","lower","upper","start")] = c(1, 0.32, 1.28, 0.64)
bounds["K",    c("phase","lower","upper","start")] = c(2, 1750, 7000, 2500)

#--------------------------------------------------------------------
# survey coverage
#--------------------------------------------------------------------
scn[5,] <- mkScn(ref="srvBias0.5", runid=5, srvBias=0.5) 
scn[6,] <- mkScn(ref="srvBias1.5", runid=6, srvBias=1.5) 

base.srvBias0.5 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=0.5, cthBias=1, IEM="linear")

base.srvBias1.5 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=1.5, cthBias=1, IEM="linear")

#--------------------------------------------------------------------
# mp effect
#--------------------------------------------------------------------

scn[7,] <- mkScn(ref="srv", runid=7, Btrig=NA, Ftar=NA, maxHR=NA, slag=5, clag=5, am="srv", b0=NA) 
scn[8,] <- mkScn(ref="cth", runid=8, Btrig=NA, Ftar=NA, maxHR=NA, clag=5, am="cth", b0=NA) 

base.srv <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=1, cthBias=1, IEM="linear", slag=5, clag=5, am="srv")

base.cth <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=1, cthBias=1, IEM="linear", clag=5, am="cth")

#--------------------------------------------------------------------
# fox model
#--------------------------------------------------------------------

scn[9,] <- mkScn(ref="fox", runid=9, p=0.05) 

bounds["p",   c("phase","start")]=c(-1,0.05)

base.fox <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=1, cthBias=1, IEM="linear")

bounds["p",   c("phase","start")]=c(-1,1)

#--------------------------------------------------------------------
# underreporting effect 0.8
#--------------------------------------------------------------------

scn[10,] <- mkScn(ref="cthBias0.8", runid=9, cthBias=0.8) 

base.cthBias0.8 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=1, cthBias=0.8, IEM="linear")

#--------------------------------------------------------------------
# MP based on survey (need mseBD2)
#--------------------------------------------------------------------

scn[11,] <- mkScn(ref="srvBias0.75", runid=10, Btrig=NA, Ftar=NA, maxHR=NA, slag=5, clag=5, aLag=5, am="srv", b0=NA, srvBias=0.75) 

base.srvBias0.75 <- mseBD2(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=5, srvBias=0.75, cthBias=1, IEM="linear", slag=5, clag=5, am="srv")

#--------------------------------------------------------------------
# worst case
#--------------------------------------------------------------------

scn[12,] <- mkScn(ref="worstCase", runid=11, Btrig=0.5, Ftar=1, maxHR=1, aLag=5, am="bd", b0=0.5, srvBias=0.5, cthBias=0.5) 

bounds["r",    c("phase","lower","upper","start")] = c(1, 0.1, 2, 0.64)
bounds["K",    c("phase","lower","upper","start")] = c(2, 1000, 7000, 2500)

base.worstCase <- mseBD2(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=1, aLag=5, srvBias=0.5, cthBias=0.5, IEM="linear", am="bd")

bounds["r",    c("phase","lower","upper","start")] = c(1, 0.32, 1.28, 0.64)
bounds["K",    c("phase","lower","upper","start")] = c(2, 1750, 7000, 2500)

#====================================================================
# Analysis
#====================================================================

lst <- list(base, base.aLag3, base.aLag5, base.cthBias0.5, base.srvBias0.5, base.srvBias1.5, base.srv, base.cth, base.fox, base.cthBias0.8, base.srvBias0.75, base.worstCase)

base.summ <- mseSumm(lst, scn)
# trick to get the panels on a better order
base.summ$par <- factor(base.summ$par, levels=c("BioDym:b0", "BioDym:K", "HCR:TAC", "HCR:hr", "BioDym:p", "BioDym:r", "catch", "fbar",  "BioDym:sigma", "BioDym:q", "rec", "ssb"))

save(base, base.aLag3, base.aLag5, base.cthBias0.5, base.srvBias0.5, base.srvBias1.5, base.srv, base.cth, base.fox, base.cthBias0.8, base.srvBias0.75, base.worstCase, base.summ, scn, file="../report/RData.mse5")

xyplot(data~factor(year)|par, groups=qtl, data=subset(base.summ, ref=="base"), scales=list(y=list(relation="free")), par.settings=list(superpose.line=list(lty=c(2,1,2), col=c(2,1,2))), type="l", main="base case")

xyplot(data~factor(year), groups=par, data=subset(base.summ, ref=="base" & par %in% c("HCR:TAC", "catch") & qtl==0.5), type="l")

plot(FLStocks(base=iter(base, 1), lag3=iter(base.aLag3, 1), lag5=iter(base.aLag5, 1), cthBias=iter(base.cthBias0.5, 1), srvBias=iter(base.srvBias0.5, 1), srvBias1.5=iter(base.srvBias1.5, 1), srvBias=iter(base.srvBias0.75, 1), srv=iter(base.srv, 1), cth=iter(base.cth,1), wc=iter(base.worstCase,1), fox=iter(base.fox,1), cthBias0.8=iter(base.cthBias0.8, 1)))

q("yes")

