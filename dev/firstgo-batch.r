#####################################################################
# first runs for testing
# using a4a assessment model
#--------------------------------------------------------------------
# year definitions: 
# TODO out of date!!  update!!
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

# add cran mirror to this line...
install.packages(c("numDeriv", "Matrix", "multicore", 
                   "ggplot2", "plyr", "akima", "lattice", "MASS"),
                 repos = "http://star-www.st-andrews.ac.uk/cran/")

# install.packages should install packages in correct order - FLCore first etc...
pkgs <- c("FLCore_2.5.0", "FLAssess_2.5.0", "FLXSA_2.5", "FLash_2.5.0",
          "ggplotFL_0.1", "FLAdvice_1.0", "FLa4a_0.1.1", "StockStructure_0.1")
pkgs <- paste0("software/", pkgs, ".tar.gz")    

install.packages(pkgs, repos = NULL)

library(StockStructure)

#==============================================================================
# simulate
#==============================================================================

sim.design <- expand.grid(
    v90     = 200000, # biomass at 90 % recruitment
    vrange = c(1),
    linf   = c(60, 80, 100),
    rmax   = 350000,
    rrange = c(1, 1.5),
    sr     = c("bevholt", "ricker"),
    stringsAsFactors = FALSE) 

ASC.brp <- 
    mclapply(1:nrow(sim.design), 
        function(i) {
          x <- sim.design[i,]
          nam <- with(x, paste0("linf:", linf, "-v90:", v90 * vrange,"-rmax:", rmax * rrange, "-sr:",sr))
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
              fbar = seq(0, 2, length = 101),
              sr = x$sr)
          
          # adjust recruitment params
          if (x $ sr == "bevholt") {
            params(res) <- FLPar(a = x$rmax * x$rrange, b = with(x, (1-.9) / .9 * v90 * vrange))
          } else {
            b <- 1 / (x$v90 * x$vrange / .9 * .4)
            a <- x$rmax * x$rrange * b * exp(1)
            params(res) <- FLPar(a = a, b = b)
          }  
          res <- brp(res)
          
          res @ desc <- nam
          res
        })

#------------------------------------------------------------------------------
# simulate F trajectory
#------------------------------------------------------------------------------

ASC.stk <- 
    mclapply(ASC.brp, 
        function(x) {
          cat(x @ desc, "\n")
          Fc <- c(refpts(x)["crash", "harvest"])
          Fmsy <- c(refpts(x)["msy", "harvest"])
          nFc <- 60
          Ftrg <- c(exp( seq(log(Fmsy), log(Fc), len = nFc) ), 
              seq(Fc, Fmsy, len = 5), 
              rep(Fmsy, 95 - nFc))
          trg <- fwdControl(data.frame(year = 2:101, quantity = rep('f', 100), val = Ftrg))
          ex.stk <- as(x, "FLStock")
          fwd(ex.stk, ctrl = trg, sr = list(model = model(x), params = params(x)))[,-(1)]
        })



#====================================================================
# Simulation settings
#====================================================================
               
nits     <- 30                  # number of iterations
iniyr    <- 2000                 # first year in projections
npyr     <- 30                   # number of years to project
lastyr   <- iniyr + npyr         # last year in projections
srsd     <- 0.3 			           # sd for S/R
units    <- 2                    # number of stock units
nhyr     <- 15                   # number of historical years
max.age  <- 10                   # plus group age
CV       <- 0.1                 # variability of catch.n and index observations


#====================================================================
# Use the stock history from one of the WKLife stocks
#====================================================================

#the following is hard coded for 6 stock unit types per SRR
n.stks <- length(ASC.brp) / 2
choices <- data.frame(s1 = rep(1:n.stks, n.stks:1), s2 = unlist(lapply(1:n.stks, function(i) i:n.stks)))

# where on the stock sims to start
#  -- start both stocks from an expoited state
choices $ syr1 <- 40
choices $ syr2 <- 40

# bevholt comes fist! make sure!! use 0.75 of f0.1 and 1 of fmsy
choices $ Btrig <- 0.75
choices <- rbind(choices, choices)
choices $ Ftarg <- rep(c(0.75, 1), each = nrow(choices)/2)

choices $ s1 <- choices $ s1 + rep(c(0, n.stks), each = nrow(choices)/2)
choices $ s2 <- choices $ s2 + rep(c(0, n.stks), each = nrow(choices)/2)
choices $ ref <- rep(c("f0.1", "msy"), each = nrow(choices)/2)
choices $ srr <- rep(c("bevh", "rick"), each = nrow(choices)/2)


fname <- 
    function(i) 
      paste0("comp:", i, 
             paste0(c(rbind("-", names(choices), ":", unlist(choices[i,]))), 
                    collapse = ""), 
             ".rda")

#mybuild()
set.seed(1734876)
seeds <- sample(1e7, 42)
# for set.seed(1734876) we get:
#4238224 8105965  380861 7140042 1197589 8984019 8242492 7539883 1822141 9835998 6339942 4149751 2488340 1205709 2121215 7025869  909545 2923516
#6670903 9087699  864331 8126975 3416187 2130411 8218030 1577699 9410389 9167451 7350541 1061413 6208053 7427966 9170855 2296642 3436897 1593588
#1812319 4647062 7237014 5040458 8799765 2213020

time0 <- proc.time()
mc.out <-
mclapply(nrow(choices):1,    
  function(i) {
    set.seed( seeds[i] )
    x <- choices[i,]
    out <- doOne(stock.id = c(x $ s1, x $ s2), 
                  Ftarg = x $ Ftarg, Btrig = x $ Btrig,
                  start.yr = c(x $ syr1, x $ syr2),
                  which.ref = "f0.1")
    attr(out, "choices") <- x          
    assign(paste0("comp", i), out, env = environment())
    save(list = paste0("comp", i), file = fname(i))

    return(fname(i))
  }
)     
time1 <- proc.time()
print(time1 - time0)


