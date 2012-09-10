
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

if (0) { ## not run - only for exploration of LH params
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

LH <- subset(LH, !is.na(bha))

LH[c("s", "v", "spr0")] <- svPars("bevholt", spr0 = LH $ spr0, a = LH $ bha, b = LH $ bhb)

LH
}

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

p <- 
  xyplot(data ~ age | type, data = sum, groups = desc, type = "l",
    scales = list(relation = "free"))

plot(p)


#------------------------------------------------------------------------------
# simulate F trajectory
#------------------------------------------------------------------------------

ASC.stk <- 
    lapply(ASC.brp, 
        function(x) {
          cat(x @ desc, "\n")
          Fc <- c(refpts(x)["crash", "harvest"]) * .8
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

#ASC.stk <- FLStocks(ASC.stk)




#====================================================================
# Simulation settings
#====================================================================

start.yr <- 40                   # where on the stock sims to start
nits     <- 20                  # number of iterations
iniyr    <- 2000                 # first year in projections
npyr     <- 20                   # number of years to project
lastyr   <- iniyr + npyr         # last year in projections - note need one 
                                 # extra year of data for predictions
srsd     <- 0.3 			           # sd for S/R
units    <- 2                    # number of stock units
nhyr     <- 15                   # number of historical years
max.age  <- 10                   # plus group age
#survey.q <- 1e-6 * exp(-2 * 0:7) # survey catchability at age
CV       <- 0.1                 # variability of catch.n and index observations


#====================================================================
# Use the stock history from one of the WKLife stocks
#====================================================================

choices <- data.frame(s1 = rep(1:5, 5:1), s2 = unlist(lapply(2:6, function(i) i:6)))

# start both stocks from an expoited state
choices $ start.yr1 <- 40
choices $ start.yr2 <- 40

choices $ Btrig <- 0.75
choices <- rbind(choices, choices)
choices $ Ftarg <- rep(c(0.75, 1), each = nrow(choices)/2)

fname <- 
    function(i) 
      paste0("comp:", i, 
             paste0(c(rbind("-", names(choices), ":", unlist(choices[i,]))), 
                    collapse = ""), 
             ".rda")

#mybuild()

set.seed(845863298)
time0 <- proc.time()
mc.out <-
mclapply(1:4, # nrow(choices),    
  function(i) {
    x <- choices[i,]
    out <- doOne(stock.id = c(x $ s1, x $ s2), 
                  Ftarg = x $ Ftarg, Btrig = x $ Btrig,
                  start.yr = c(x $ start.yr1, x $ start.yr2),
                  which.ref = "f0.1")
  
    save(out, file = fname(i))

    return(fname(i))
  }
)     
time1 <- proc.time()
print(time1 - time0)


#attr(base, "summaries")