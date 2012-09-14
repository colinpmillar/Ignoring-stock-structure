
#####################################################################
# Collate the rda files from runs on different pcs
#####################################################################

library(StockStructure)

setwd("../presentation")
load("ASC.rda")

ibind <- function(obj1, obj2) {
  # check dims
  # TODO check dimensions
  
  iters1 <- dims(obj1) $ iter
  iters2 <- dims(obj2) $ iter
  
  obj <- propagate(obj1, iters1 + iters2)
  
  obj[,,,,,1:iters1] <- c(obj1)
  obj[,,,,,iters1 + 1:iters2] <- c(obj2)
  
  obj
}

ibind.stk <- function(obj1, obj2) {
  iters1 <- dims(obj1) $ iter
  iters2 <- dims(obj2) $ iter
  
  obj <- propagate(obj1, iters1 + iters2)
  
  funs <- c("harvest", 
            "catch"   , "catch.n"   , "catch.wt",
            "discards", "discards.n", "discards.wt",
            "landings", "landings.n", "landings.wt",
            "stock.n" , "stock.wt"  , "stock", 
            "m", "mat", "harvest", "harvest.spwn", "m.spwn")
  for (i in funs) {
    fun <- match.fun(i)
    slot(obj, i) <- ibind(fun(obj1), fun(obj2))
  } 

  obj
}


# first job - combine the output
files <- data.frame(comp = 1:nrow(choices))
files[paste0("out", 1:2)] <- 
    lapply(paste0("out", 3:4), 
      function(x) 
        1:42 %in% sort(as.numeric(gsub("-", "", substring(dir(x), 6, 7)))))

comps <- lapply(1:42,
    function(i) {
      out1 <- out2 <- NULL
      if (files $ out1[i]) {
        load(paste0("out3/", fname(i)))
        out1 <- get(paste0("comp", i))
      } 

      if (files $ out2[i]) {
        load(paste0("out4/", fname(i)))
        out2 <- get(paste0("comp", i))
      } 

      out <- NULL
      if (is.null(out1)) {
        out <- out2
      } else if (is.null(out2)) {
        out <- out1
      } else if (!is.null(out1) & !is.null(out2)){
        out <- ibind.stk(out1, out2)
      }

      out
    })
  
# okay - collate the summaries

my.summary <- function(x, p = c(.1, .9)) {
  out <- c(quantile(x, min(p)), mean(x), quantile(x, max(p)))
  names(out)[2] <- "mean"
  out
}

# F2029 / Ftarg
f.ratio <-
  mclapply(1:42, function(i) {
      tmp <- comps[[i]]

      f29 <- c(fbar(tmp[,"2029"])[,,1])
      fref <- choices $ Ftarg[i] *
              c(c(refpts( ASC.brp[[ choices $ s1[i] ]] )[choices $ ref[i], "harvest"]),
              c(refpts( ASC.brp[[ choices $ s2[i] ]] )[choices $ ref[i], "harvest"])) 

      fratio <- outer(f29, fref, "/")
      apply(fratio, 2, my.summary)
    })


propSummary <- function(obj, year = "2029", fun = "tsb") {
  fun <- match.fun(fun)
  out <- drop(fun(obj[,ac(year)]) @ .Data)
  apply(sweep(out, 2, colSums(out), "/"), 1, my.summary)
}

ratioSummary <- function(obj, year = "2029", fun = "tsb") {
  fun <- match.fun(fun)
  out <- drop(fun(obj[,ac(year)]) @ .Data)
  out <- out[2,] / out[1,]
  my.summary(out)
}

# TSB and catch ratio pop2 / pop1 in 2029
tsb.ratio <- lapply(1:42, function(i) ratioSummary(comps[[i]], fun = tsb))
catch.ratio <- lapply(1:42, function(i) ratioSummary(comps[[i]], fun = catch))

# TSB and catch proprtions in 2029
tsb.prop <- lapply(1:42, function(i) propSummary(comps[[i]], fun = tsb))
catch.prop <- lapply(1:42, function(i) propSummary(comps[[i]], fun = catch))

save.image("presentation-data.rda")

pdf("tex/full-result1.pdf")
plotOM(comps[[1]])
dev.off()

pdf("tex/full-result2.pdf")
plotOM(comps[[3]])
dev.off()



## results summary


res <- choices

res $ tsb.prop <- sapply(tsb.prop, function(x) x["mean", "pop1"])
res $ tsb.ratio <- sapply(tsb.ratio, function(x) x["mean"])
res $ catch.prop <- sapply(catch.prop, function(x) x["mean", "pop1"])
res $ catch.ratio <- sapply(catch.ratio, function(x) x["mean"])

    
subset(res, s1 %in% 1)

xyplot(tsb.prop ~ I(s2 %% 6 + 1) | factor(s1), data = res)




