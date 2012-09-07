
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
library(FLa4a)

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


#====================================================================
# Choose which stocks to play with
#====================================================================

# look through ICES reports for the values of a few stocks and rerun the simulations
# use
# cod like
# plaice like

#==============================================================================
# get NS demersal data for conditioning
#==============================================================================

ices <- read.csv("../data/icesdata.csv", stringsAsFactors = FALSE)
stocks <- c("nop-34", "cod-347d", "had-34", "ple-eche", "ple-nsea", "sol-eche", "sol-nsea", "whg-47d")
ices <- subset(ices, FishStock %in% stocks)[,c("FishStock", "Year", "Recruitment", "SSB","MeanF")]
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
          b <- coef(glm(Recruitment ~ I(1/(SSB * a)) + offset(1/a) - 1, Gamma, x)  )
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

LHlist <- merge(srr, LH, all = TRUE)

LHlist $ t0 [is.na(LHlist $ t0)] <- -0.1

LHlist <- subset(LHlist, !is.na(bha))

LHlist[c("s", "v", "spr0")] <- svPars("bevholt", spr0 = LHlist $ spr0, a = LHlist $ bha, b = LHlist $ bhb)

LHlist

#==============================================================================
# simulate
#==============================================================================
#------------------------------------------------------------------------------
# get LH pars 
#------------------------------------------------------------------------------

B90 <- 15000

sim.design <- mapply(function(v, linf) list(v = v, linf = linf), 
                     v    = rep(c(1.2, 0.8) * B90, each = 3), 
                     linf = rep(c(60 , 80, 100), 2), 
                     SIMPLIFY = FALSE)
                 
ASC.brp <- 
  lapply(sim.design, 
    function(x) {
      nam <- paste0("linf:", x $ linf, "-v:", x $ v) 
      cat(nam, "\n")

      # complete with gislasim
      par <- gislasim(linf = x $ linf, v = 5000, 
                      s = 1, # slope at 20% of virgin biomass - if asymtotic R is 
                      #a50 = 2, ato95 = 1, # a50 and age to a95 for maturity ogive
                      sl = 2, sr = 5000, a1 = 4, # left and right slope and fully selected age for fishery selectivity
                      a = 5e-6, b = 2.9, # weight length relationship - similar to demersal types
                      t0 = -1) # t0 in vonB growth model - could be better... 
                      # all other params - k, M1, M2 from gislason et al. 
                      # bg set to b, bg is b in WLR for stock weight
      
      # run LH          
      res <- lh(par, sr = "bevholt", # choose beverton holt SRR 
                range = c(min = 1, max = 40, # full age range in simulations
                          minfbar = 4, maxfbar = 8, # fbar over fully expoited ages 
                          plusgroup = 40), # plus group at 30cm
                fish = 0.5, # fishing occurs half way through the year.  
                fbar = seq(0, 1, length = 101)) 
            
      # adjust recruitment params for fixed assymptote and biomass at 90% of assymptote
      params(res) <- FLPar(a = 300000, b = x $ v / .9)
      res <- brp(res)
      
      res @ desc <- nam
      res
    })

# stock summary comparisons
funs <- c("catch.sel", "mat", "m", "stock.wt")
sum <- do.call(rbind,
         lapply(ASC.brp, 
           function(brp) 
             do.call(rbind,
                 lapply(funs, function(x) cbind(as.data.frame(do.call(x, list(brp)))[c("age","data")], fun = x, desc = brp @ desc))
             )))
  
# stock recruit comparisons
sr.sum <- do.call(rbind,
            lapply(ASC.brp, function(brp) { 
              ssb <- seq(1, virginB, length = 100)
              rec <- eval(model(brp)[[3]], 
                          c(as.list(params(brp) @ .Data[,1]), list(ssb = ssb)))
              data.frame(desc = brp @ desc, data = rec, age = ssb, fun = "srr")
            }))


xyplot(I(data * exp(rnorm(1, 0, .5))) ~ age | fun, groups = desc, data = rbind(sum, sr.sum), type = "l", scales = list(relation = "free"))


# ref points
refpt <- sapply(ASC.brp, function(x) drop(refpts(x)[c("msy", "crash"),c("harvest")] @ .Data))
colnames(refpt) <- sapply(ASC.brp, function(x) x @ desc)
refpt

#------------------------------------------------------------------------------
# simulate F trajectory
#------------------------------------------------------------------------------

ASC.stk <- 
  lapply(ASC.brp, 
    function(x) {
      cat(x @ desc, "\n")
      Fc <- c( refpts(x)["crash", "harvest"] * 0.9)
      Fmsy <- c(refpts(x)["msy", "harvest"])
      nFc <- 60
      Ftrg <- c(exp( seq(log(Fmsy), log(Fc), len = nFc)), seq(Fc, Fmsy, len = 5), rep(Fmsy, 95 - nFc))
      trg <- fwdControl(data.frame(year = 2:101, quantity = rep('f', 100), val = Ftrg))
      ex.stk <- as(x, "FLStock")
      out <- fwd(ex.stk, ctrl = trg, sr = list(model = "bevholt", params = params(x)))[,-1]
      plot(out)
      out
    })


ASC.stk <- FLStocks(ASC.stk)

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
nhyr     <- 15                   # number of historical years
max.age  <- 10                    # plus group age
survey.q <- 1e-6 * exp(-2 * 0:7) # survey catchability at age
CV       <- 0.15                 # variability of catch.n and index observations
stock.ids <- c(1,2)                    # which stocks to use

#====================================================================
# Use the stock history from one of the WKLife stocks
#====================================================================


#--------------------------------------------------------------------
# True stock msy reference points
#--------------------------------------------------------------------
refpt1 <- refpts( ASC.brp[[stock.ids[1]]] )["msy", c("ssb", "harvest")]
refpt2 <- refpts( ASC.brp[[stock.ids[2]]] )["msy", c("ssb", "harvest")]


#--------------------------------------------------------------------
# True stock history
#--------------------------------------------------------------------
start.yr <- 10
pop1 <- window(ASC.stk[[stock.ids[1]]], 
               start = start.yr, end = start.yr + nhyr - 1)[1:max.age,]
pop2 <- window(ASC.stk[[stock.ids[2]]], 
               start = start.yr, end = start.yr + nhyr - 1)[1:max.age,]
dimnames(pop1) <- dimnames(pop2) <- list(year = 2001 - nhyr:1)

true.stock <- pop1
true.stock <- FLCore::expand(true.stock, unit = c("unique", "pop1", "pop2"))[,,2:3]

true.stock[,,1] <- pop1
true.stock[,,2] <- pop2

# change fbar here also ...
range(true.stock)[c("minfbar", "maxfbar")] <- c(4,8)
range(true.stock)["plusgroup"] <- max.age + 1 # to deal sepVPA and a4aFit


#====================================================================
# Stock recruit model
#====================================================================

# assume a stock recruitment model for each population and estimate 
# the parameters from it - this serves as the underlying
# stock recruitment model
sr.model1 <- list(model = "bevholt", params = params(ASC.brp[[stock.ids[1]]]))
sr.model2 <- list(model = "bevholt", params = params(ASC.brp[[stock.ids[2]]]))

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
OM <- fwdWindow(true.stock, FLBRP(true.stock), end = lastyr)

# propagate
OM <- propagate(OM, nits)

# project to the end of projections at last year F level
ctrl <- fwdControl( data.frame(year = iniyr:lastyr, quantity = "f", val = .15) )
# if no sr residuals then no noise!
OM[,,1] <- fwd(OM[,,1], ctrl = ctrl, sr = sr.model1)
OM[,,2] <- fwd(OM[,,2], ctrl = ctrl, sr = sr.model2)


#====================================================================
# first simulation
#====================================================================

base <- mse(OM = OM, iniyr = iniyr, 
		        sr.model1 = sr.model1, sr.model2 = sr.model2,
            sr.residuals = sr.residuals, 
		        Btrig = 0.5, Ftar = 0.75, refpt = refpt1,
			      seed = NULL)


        
# some plotting functions:        
plotOM <- function(OM, probs = c(0.025, .5, .975)) {
  assessment <- attr(OM, "assessment.data") $ stock
  
  myrec <- function (object) stock.n(object)[1,]               
  fn    <- list(SSB = ssb, Recruits = myrec, Yield = catch, F = fbar)
  res1  <- ggplotFL:::whooow(OM[,,1], fn, probs)
  res2  <- ggplotFL:::whooow(OM[,,2], fn, probs)
  res <- rbind(cbind(res1, population = "1"), cbind(res2, population = "2"))
  res $ group <- paste(res $ pop, res $ iter)
  
  p1 <- ggplot(res) + 
        geom_line(aes(x = year, y = data, group = group, 
                      size = iter, lty = iter, col = population)) + 
        scale_size_manual(values = c(0.5, 1, 0.5), name = "Quantile") + 
        scale_linetype_manual(values = c(2, 1, 2),  name = "Quantile") +
        expand_limits(y = 0) + xlab("Year") + ylab("") +
        facet_wrap(~ qname, scale = "free")

  print(p1)
  invisible(p1)  
}        

plotAssessment <- function(OM) {
  assessment <- attr(OM, "assessment.data") $ stock
  plotComp(assessment, 
      fn = list(SSB = ssb, Recruits = myrec, Yield = catch, F = fbar), 
      probs = c(0.75, 0.5, 0.25), 
      size= c(0.5, 1, 0.5), 
      lty = c(2,1,2), 
      facet = facet_wrap(~qname, scale = "free"))  
}

plotOM(OM)
plotOM(base)
plotAssessment(base)


     


#attr(base, "summaries")
