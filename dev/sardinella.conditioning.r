# explore.R - DESC
# explore.R

# Copyright 2003-2012 FLR Team. Distributed under the GPL 2 or later
# Maintainer: Iago Mosqueira, JRC
# $Id: $

# http://www.fao.org/fishery/rfb/cecaf/en

library(FLAdvice)
library(FLBioDym)
library(FLAssess)
library(plyr)
library(reshape)

# Sardinella data
nc <- read.table('data/sar/2010/nc_sar.dat', header=TRUE)
ia <- read.table('data/sar/2010/cpue_sar.dat', header=TRUE)

ca <- FLQuant(nc$catch, dimnames=list(age='all', year=nc$year))
cp <- FLQuant(c(rep(NA, 7), ia$cpue), dimnames=list(age='all', year=nc$year))

# CECAF SA results
sa <- read.table('data/sar/2010/sa.dat', sep='\t', header=T)
saq <- FLQuants(HR=FLQuant(sa$hr, dimnames=list(age='all', year=sa$year)),
  SSB=FLQuant(sa$ssb, dimnames=list(age='all', year=sa$year)))


# SA using FLBioDym {{{

bd <- FLBioDym(catch=ca, index=cp)

# bounds & start
# r
bd@bounds[1,] <- c(1, 0.32, 1.28, 0.64)
# K
bd@bounds[2,] <- c(2, 1750, 7000, 2500)
# p
bd@bounds[3,] <- c(-1, 1, 1, 1)
# b0
bd@bounds[4,] <- c(-1, 0.25, 0.95, 0.5)
# q
bd@bounds[5,] <- c(1, 0.1, 10, 1)
# sigma
bd@bounds[6,] <- c(1, 0.01, 1, 0.1)

# admbBD assessment
res <- admbBD(bd)
# }}}


# gislasim-based OM {{{
# LH parameters from FishBase & FLBioDym
par <- as(data.frame(linf=28.5, k=0.4, t0=-0.1, s=0.8, v=4500, a50=1), 'FLPar')

# Use CECAF SA B0 as starting point
par <- as(data.frame(linf=28.5, k=0.4, t0=-0.1, s=0.8, v=1750, a50=1), 'FLPar')

# gislasim
sar <- lh(gislasim(par), range=c(min=1, max=8, minfbar=1, maxfbar=6))

# stk with initial F closest to estimated HR
stk <- as(sar, 'FLStock')
stk <- stk[,12]
dimnames(stk) <- list(year=1989)
# }}}


# standard OM projections {{{

# prepare for projection
saa <- stf(stk, 21, 1)

# projection control
trg <- fwdControl(data.frame(year=1990:2010, val=c(ca), quantity="catch"))

# catch projection
saa <- fwd(saa, trg, sr=list(model="bevholt", params=params(sar)))

name(saa) <- "SAA"
desc(saa) <- "Simulated Sardinella aurita, partly conditioned on CECAF SA 2010"

save(saa, file=paste("SAA_", format(Sys.time(), "%Y%m%d_%H%M"), ".RData", sep=""))
# }}}


# Recruitment patterns in conditioning {{{

# prepare for projection
saar <- stf(stk, 21, 1)

trg <- fwdControl(data.frame(year=1990:2010, val=c(ca), quantity="catch"))

# optim for exceptional recruitment scaler 

# foo
foo <- function(m) {
  sr.residuals <- FLQuant(c(rep(1, 15), m*3, 1, m*2, 1, 1, 1), 
    dimnames=list(year=1990:2010))
  saar <- fwd(saar, trg, sr=list(model="bevholt", params=params(sar)),
    sr.residuals=sr.residuals)
  return(sum((ssb(saar)[,-1] - saq$SSB)^ 2))
}

# get scaler
m <- optimize(foo, interval=c(0.1, 5))$minimum

# apply w/ iters
sr.residuals <- rlnorm(5, log(FLQuant(c(rep(1, 15),
  m*3, 1, m*2, 1, 1, 1), dimnames=list(year=1990:2010))), 0.3)

bwplot(data~year, sr.residuals)


# project
saar <- fwd(saar, trg, sr=list(model="bevholt", params=params(sar)),
  sr.residuals=sr.residuals)

xyplot(data~year|qname, as.data.frame(FLQuants(SA=saq$SSB, OM=iterMeans(ssb(saar)[,-1]))), type='b')

bwplot(data~year, ssb(saar)[,-1],
  scales = list(x = list(tck=1, at = seq(1, 20, length=5),
  labels = seq(1990, 2010, by=5))))

saar <- propagate(saar, 5)

save(saar, file=paste("SAA_", format(Sys.time(), "%Y%m%d_%H%M"), ".RData", sep=""))
# }}}


# Use fwd on both ssb and catch {{{
# prepare for projection
saa <- stf(stk, 21, 1)

# projection control
df <- rbind(data.frame(year=1990:2010, val=c(ca), quantity="catch"),
  data.frame(year=1990:2010, val=sa$ssb, quantity="ssb"))

trg <- fwdControl(df[order(df$year),])

# catch projection
saa <- fwd(saa, trg, sr=list(model="bevholt", params=params(sar)))

xyplot(data~year|qname, as.data.frame(FLQuants(SA=saq$SSB, OM=ssb(saa)[,-1])), type='b')

# }}}


# 2011 SA inputs {{{
library(ROpenOffice)                                                                             
ncf <- read.ods('data/sar/2011/Sardinella_CECAF_2011.ods')
ncf <- melt(ncf, 1:2)
names(ncf)[3:4] <- c('year', 'catch')

ggplot(ncf, aes(year, catch)) + geom_line(aes(group=fleet, colour=fleet)) + facet_wrap(~country, scales='free')

ggplot(ncf, aes(year, catch)) + geom_line(aes(group=country, colour=country)) + facet_wrap(~fleet, scales='free')

# }}}


# fwd recruitment residuals
resid <- FLQuant(rlnorm(50*5, 0, 0.3) * (rpois(50*5, 2/20) * 0.5 + 1), dimnames=list(year=2011:2061, iter=1:5))
