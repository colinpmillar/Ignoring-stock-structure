
rm(list = ls(all = TRUE))

mybuild <- function() {
  library(devtools)
  library(testthat)
  library(roxygen2)

  roxygenize("../../a4a/packages/FLa4a")

  pkg <- as.package("../../a4a/packages/FLa4a")
  build(pkg)
  #check(pkg)
  install(pkg)

  roxygenize("../StockStructure")
  pkg <- as.package("../StockStructure")
  build(pkg)
  #check(pkg)
  install(pkg)
}

mybuild()

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
library(FLa4a)


#====================================================================
# Choose which stocks to play with
#====================================================================

# look through ICES reports for the values of a few stocks and rerun the simulations
# use
# cod like
# plaice like

#==============================================================================
# get ICES recruitment data
#==============================================================================

ices <- read.csv("../data/icesdata.csv", stringsAsFactors = FALSE)
stocks <- c("nop-34", "cod-347d", "had-34", "ple-eche", "ple-nsea", "sol-eche", "sol-nsea", "whg-47d")
ices <- subset(ices, FishStock %in% stocks)[,c("FishStock", "Year", "Recruitment", "SSB","MeanF")]
rm(stocks)
ices $ species <- substring(ices $ FishStock, 1, 3)
ices $ SSB <- ices $ SSB * 1000 # kg
ices $ Recruitment <- ices $ Recruitment # 1000s

srr <- 
    sapply(
        split(ices, ices $ FishStock),
        function(x) {
          x <- subset(x, Year > 1990)
		  a <- exp(mean(log(x $ SSB), na.rm = TRUE))
		  x $ a <- a
          b <- coef(glm(Recruitment ~ I(1/ (SSB * a)) + offset(1/a) - 1, Gamma, x)  )
		  slope <- unname(a/b)
          c(a = a, b = b, spr0(FLQuant(x $ SSB), FLQuant(x $ Recruitment), FLQuant(x $ MeanF)))
        })

srr <- data.frame(FishStock = colnames(srr), bha = srr[1,], bhb = srr[2,], spr0 = srr[3,])
rownames(srr) <- NULL
srr $ species <- substring(srr $ FishStock, 1, 3)

# get LH params for these stocks ... 
# need - cod, had, plaice, sol, whg
LH <- read.csv("../data/lh.csv", stringsAsFactors = FALSE)
LH $ a <- LH $ a / 1000  # convert output to kg

LH <- merge(srr, LH, all = TRUE)
rm(srr)

LH <- subset(LHlist, !is.na(bha))

LH[c("s", "v", "spr0")] <- svPars("bevholt", spr0 = LH $ spr0, a = LH $ bha, b = LH $ bhb)

LH

#==============================================================================
# simulate
#==============================================================================
#------------------------------------------------------------------------------
# get LH pars 
#------------------------------------------------------------------------------

sim.design <- expand.grid(
                v90    = 20000, # biomass at 90% recruitment
		            vrange = c(1.1, 0.9),
                linf   = c(60, 80, 100),
                rmax   = 300000)

ASC.brp <- 
  lapply(1:nrow(sim.design), 
    function(i) {
      x <- sim.design[i,]
      nam <- with(x, paste0("linf:", linf, "-v90:", v90 * vrange))
      cat(nam, "\n")

      # complete M, mat and k with gislasim
      par <- gislasim(linf = x $ linf,
			          v = 5000, s = 0.9, # change these later
					  a = 5e-6, b = 2.9, # weight length
					  sl = 2, sr = 5000, a1 = 4, # selectivity params
					  ato95 = 1,  # slope of maturity ogive
					  t0 = -1)    # fixed vonB param

      # run LH
      res <- lh(par, 
                range = c(min = 1, max = 20, 
                          minfbar =4, maxfbar = 8, 
                          plusgroup = 20), 
                spwn = 0, fish = 0.5,
                fbar = seq(0, 2, length = 101))
	  
      # adjust recruitment params
      params(res) <- FLPar(a = x $ rmax, 
                           b = x $ v90 * x $ vrange / 0.9)
			res <- brp(res)
                       
      res @ desc <- nam
      res
    })

rm(sim.design)

# refpts
refs <- sapply(ASC.brp, function(x) drop(refpts(x)[c("msy", "crash"),"harvest"] @ .Data))
colnames(refs) <- sapply(ASC.brp, function(x) x @ desc)
refs

msyrefs <- round(sapply(ASC.brp, function(x) drop(refpts(x)["msy",c("harvest","yield","rec","ssb")] @ .Data)), 2)
colnames(msyrefs) <- sapply(ASC.brp, function(x) x @ desc)
msyrefs


# plot a summary of the BRP
funs <- c("m", "mat", "stock.wt", "landings.sel")
sum <- do.call(rbind, lapply(ASC.brp, function(brp) 
          do.call(rbind, lapply(funs, function(f)
             cbind(as.data.frame(do.call(f, list(brp)))[c("age","data")],
                   type = f, desc = brp @ desc)
       ))))
rm(funs)

xyplot(data ~ age | type, data = sum, groups = desc, type = "l",
       scales = list(relation = "free"))

   
   
#------------------------------------------------------------------------------
# simulate F trajectory
#------------------------------------------------------------------------------

ASC.stk <- 
  lapply(ASC.brp, 
    function(x) {
      cat(x @ desc, "\n")
      Fc <- c( refpts(x)["crash", "harvest"] * 0.7)
      Fmsy <- c(refpts(x)["msy", "harvest"])
      nFc <- 60
      Ftrg <- c(exp( seq(log(Fmsy), log(Fc), len = nFc) ), 
                seq(Fc, Fmsy, len = 5), 
                rep(Fmsy, 95 - nFc))
      trg <- fwdControl(data.frame(year = 2:101, quantity = rep('f', 100), val = Ftrg))
      ex.stk <- as(x, "FLStock")
      out <- fwd(ex.stk, ctrl = trg, sr = list(model = "bevholt", params = params(x)))[,-(1)]
      plot(out)
      out
    })

ASC.stk <- FLStocks(ASC.stk)




do.one <- function(stock.id, Ftar = 3, Btar = 60000) {
  
#--------------------------------------------------------------------
# True stock msy reference points
#--------------------------------------------------------------------
  refpt1 <- ASC.brp[[stock.id[1]]] @ refpts["msy", c("ssb", "harvest")]
  refpt2 <- ASC.brp[[stock.id[2]]] @ refpts["msy", c("ssb", "harvest")]
  
  
#--------------------------------------------------------------------
# True stock history
#--------------------------------------------------------------------
  pop1 <- window(ASC.stk[[stock.id[1]]], 
      start = start.yr, end = start.yr + nhyr - 1)[1:max.age,]
  pop2 <- window(ASC.stk[[stock.id[2]]], 
      start = start.yr, end = start.yr + nhyr - 1)[1:max.age,]
  dimnames(pop1) <- dimnames(pop2) <- list(year = 2001 - nhyr:1)
  
  true.stock <- pop1
  true.stock <- FLCore::expand(true.stock, unit = c("unique", "pop1", "pop2"))[,,2:3]
  
  true.stock[,,1] <- pop1
  true.stock[,,2] <- pop2
  rm(pop1, pop2)
  
# change fbar here also ...
  range(true.stock)[c("minfbar", "maxfbar")] <- c(4,8)
  range(true.stock)["plusgroup"] <- max.age + 1 # to deal sepVPA and a4aFit
  
  
#====================================================================
# Stock recruit model
#====================================================================
  
# assume a stock recruitment model for each population and estimate 
# the parameters from it - this serves as the underlying
# stock recruitment model
  sr.model1 <- list(model = "bevholt", params = params(ASC.brp[[stock.id[1]]]))
  sr.model2 <- list(model = "bevholt", params = params(ASC.brp[[stock.id[2]]]))
  
# stock recruit residuals - simulate residuals as lognormal with sd=srsd
  sr.residuals <- FLQuant(rlnorm(npyr * nits, sd = srsd), 
      dimnames = list(age  = 1, 
          year = dims(true.stock)$minyear:lastyr, 
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
# this line fails validFLBRP check on unadultarated FLCore
  OM <- fwdWindow(true.stock, FLBRP(true.stock), end = lastyr)
  
# propagate
  OM <- propagate(OM, nits)
  
# project to the end of projections at last year F level
  Fsq <- c(fbar(OM[, ac(iniyr)]))[1]
  ctrl <- fwdControl( data.frame(year = iniyr:lastyr, quantity = "f", val =  Fsq))
# if no sr residuals then no noise!
  OM[,,1] <- fwd(OM[,,1], ctrl = ctrl, sr = sr.model1)
  OM[,,2] <- fwd(OM[,,2], ctrl = ctrl, sr = sr.model2)
  rm(ctrl)
  
#====================================================================
# first simulation
#====================================================================
  
  
  refpt <- refpt1
  refpt[] <- c(Btar, Ftar) # estimates of SSB and F msy
  
  base <- mse(OM = OM, iniyr = iniyr, 
      sr.model1 = sr.model1, sr.model2 = sr.model2,
      sr.residuals = sr.residuals, 
      refpt = refpt)
  
  base        
}        


#====================================================================
# Simulation settings
#====================================================================

start.yr <- 40                   # where on the stock sims to start
nits     <- 100                  # number of iterations
iniyr    <- 2000                 # first year in projections
npyr     <- 50                   # number of years to project
lastyr   <- iniyr + npyr         # last year in projections - note need one 
                                 # extra year of data for predictions
srsd     <- 0.3 			           # sd for S/R
units    <- 2                    # number of stock units
nhyr     <- 15                   # number of historical years
max.age  <- 10                   # plus group age
#survey.q <- 1e-6 * exp(-2 * 0:7) # survey catchability at age
CV       <- 0.01                 # variability of catch.n and index observations


#====================================================================
# Use the stock history from one of the WKLife stocks
#====================================================================

choices <- cbind( rep(1:5, 5:1), unlist(lapply(2:6, function(i) i:6)))

Ftar <- .5
Btar <- 1000

#mybuild()

for (i in 1:nrow(choices)) {
  out <- do.one(choices[i,], Ftar = Ftar, Btar = Btar)

  dev.new()
  plotOM(out)
  dev.new()
  plotAssessment(out)

  save(out, file = paste0("comp:",i,"Syr:", start.yr,"-Ftar:", Ftar, "-Btar:", Btar,".rda"))
}     


#attr(base, "summaries")
