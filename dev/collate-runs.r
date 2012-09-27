#####################################################################
# Collate the rda files from runs on different pcs
#####################################################################
options(width=80)
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
  #out <- c(quantile(x, min(p)), mean(x), quantile(x, max(p)))
  #names(out)[2] <- "mean"
  #out
  mean(x)
}

# F2029 / Ftarg
f.ratio <-
  sapply(1:42, function(i) {
      tmp <- comps[[i]]

      f29 <- drop(fbar(tmp[,"2029"]) @ .Data)
#      fref <- choices $ Ftarg[i] *
#              c(c(refpts( ASC.brp[[ choices $ s1[i] ]] )[choices $ ref[i], "harvest"]),
#              c(refpts( ASC.brp[[ choices $ s2[i] ]] )[choices $ ref[i], "harvest"])) 
      fref <-
          c(c(refpts( ASC.brp[[ choices $ s1[i] ]] )["msy", "harvest"]),
              c(refpts( ASC.brp[[ choices $ s2[i] ]] )["msy", "harvest"])) 
      
      fratio <- sweep(f29, 1, fref, "/")
      apply(fratio, 1, my.summary)
    })

b.ratio <-
    sapply(1:42, function(i) {
          tmp <- comps[[i]]
          
          b29 <- drop(ssb(tmp[,"2029"]) @ .Data)
          bref <- 
              c(c(refpts( ASC.brp[[ choices $ s1[i] ]] )["msy", "ssb"]),
                  c(refpts( ASC.brp[[ choices $ s2[i] ]] )["msy", "ssb"])) 
          
          bratio <- sweep(b29, 1, bref, "/")
          apply(bratio, 1, my.summary)
        })

choices $ type <- ""
choices $ type[c(1, 7, 12, 16, 19, 21)] <- "control"
choices $ type[c(2, 3, 8, 17, 18, 20)] <- "growth"
choices $ type[c(4, 10, 15)] <- "recruitment"
choices $ type[c(9, 13, 14, 5, 6, 11)] <- "interaction"

choices $ type2 <- choices $ type
choices $ type2[c(9, 13, 14)] <- "unclear"

choices $ type3 <- ""
choices $ type3[c(1, 2, 3, 7, 8, 12)] <- "growth at low"
choices $ type3[16:21] <- "growth at high"
choices $ type3[c(4, 5, 6, 10, 11, 15)] <- "interaction"


choices $ type[22:42] <- choices $ type[1:21]
choices $ type2[22:42] <- choices $ type2[1:21]
choices $ type3[22:42] <- choices $ type3[1:21]


choices $ linf1 <- sim.design $ linf[choices $ s1]
choices $ linf2 <- sim.design $ linf[choices $ s2]

choices $ rec1 <- sim.design $ rrange[choices $ s1]
choices $ rec2 <- sim.design $ rrange[choices $ s2]


# order comps - which is the most productive
choices $ high <- 2
choices $ high[c(9, 13, 14)] <- 1
choices $ low <- 3 - choices $ high 

#choices $ frat.diff <- f.ratio[cbind(choices $ high, 1:42)] - f.ratio[cbind(choices $ low, 1:42)]    
choices $ brat.diff <- b.ratio[cbind(choices $ low, 1:42)] / b.ratio[cbind(choices $ high, 1:42)]    


propSummary <- function(obj, year = "2029", fun = "tsb") {
  fun <- match.fun(fun)
  out <- drop(fun(obj[,ac(year)]) @ .Data)
  apply(sweep(out, 2, colSums(out), "/"), 1, my.summary)
}

# TSB and catch proprtions in 2029
tsb.str <- sapply(1:42, function(i) propSummary(comps[[i]], fun = tsb, year = 2000))
tsb.end <- sapply(1:42, function(i) propSummary(comps[[i]], fun = tsb))

choices $ tsbrat.str <- tsb.str[cbind(choices $ low, 1:42)] / tsb.str[cbind(choices $ high, 1:42)]    
choices $ tsbrat.end <- tsb.end[cbind(choices $ low, 1:42)] / tsb.end[cbind(choices $ high, 1:42)]    


choices2 <- subset(choices, type2 != "unclear")
choices2 $ srr <- ifelse(choices2 $ srr == "rick", "Ricker", "Beverton Holt")
choices2 $ type3 <- ordered(choices2 $ type3, levels = c("growth at low", "growth at high", "interaction"))


p <- with(choices2,
xyplot(factor(linf1) ~ factor(linf2) | type3 * srr,
       panel = function(x, y, subscripts, ...) {
         cex <- choices2 $ brat.diff[subscripts]
         cex <- cex * 8
         lpoints(x, y, cex = cex, col = grey(0.5), pch = 16)
         lpoints(x, y, cex = cex, col = 1)
       }, ylab = "Population 1 Linf", xlab = "Population 2 Linf")
)

pdf("tex/Bmsy-plot.pdf", width = 9, height = 7)
print(p)
dev.off()




choices3 <- rbind(choices, choices)
choices3 $ year <- rep(c("2000", "2029"), each = nrow(choices))
choices3 $ tsbrat <- with(choices, c(tsbrat.str, tsbrat.end))
choices3 $ comp <- rep(1:42, 2)
    
choices3 <- subset(choices3, type2 != "unclear")
choices3 $ srr <- ifelse(choices3 $ srr == "rick", "Ricker", "Beverton Holt")
choices3 $ type3 <- ordered(choices3 $ type3, levels = c("growth at low", "growth at high", "interaction"))


p <- with(choices3,
    xyplot(tsbrat ~ factor(year) | type3 * srr, 
        group = paste(linf1, linf2, sep=":"), lwd = 2,
        ylab = "TSB ratio", xlab = "", type = "l", ylim = c(0, 1.3),
        auto.key = list(points = FALSE, lines = TRUE, columns = 3),
        #lty = c(1,1,1,1,2,1), 
        #col = c("grey70", "grey10", "grey70", "grey50", "grey50", "grey70"),
        par.settings = list(superpose.line = list(lty = c(1,1,1,1,2,1), 
                col = c("grey70", "grey10", "grey70", "grey50", "grey50", "grey70")
                )))
)

pdf("tex/Bratio-plot.pdf", width = 9, height = 7)
print(p)
dev.off()






choices $ group <- paste0(choices $ srr, choices $ type2)
unique(choices $ group)

p1 <- xyplot(frat.diff ~ brat.diff | type, group = group, data = choices,
       pch = c(1, 2, 3, 4, 5), col = rep(c(1,2), each = 5),
       panel = function(x, y, ...) {
         panel.superpose(x, y, ...)
         panel.abline(h = 0, v = 0, col = "black")
         panel.superpose(x, y, ...)
       }, ylab = "F/Fmsy ratio", xlab = "B/Bmsy ratio")
   
pdf("tex/stock-plot.pdf")
print(p1)
dev.off()


state <- rbind(choices, choices)
state $ pop <- rep(1:2, each = nrow(choices))
state $ frat <- c(t(f.ratio))
state $ brat <- c(t(b.ratio))
state $ pop[ state $ s1 == state $ s2 ] <- 0
revpop <- c(17, 18, 20)
state $ pop[revpop] <- 2
state $ pop[revpop + 21] <- 2
state $ pop[revpop + 42] <- 1
state $ pop[revpop + 21 + 42] <- 1


par(oma = c(0,0,0,0), mar = c(0,0,0,0))
plot(state $ brat, state $ frat, type="n")
abline(h=1, v = 1)
for (i in which(!(choices $ type %in% c("base-case", "odd-case")) &
        choices $ srr == "bevh")) {
  with(state[i + c(0, 42),], {
        lines(brat, frat, col = grey(0.5))
        #points(brat, frat, col = pop + 1, pch = 16)
        text(brat, frat, i %% 21, col = pop + 1, cex = 0.8)
      })
}

par(oma = c(0,0,0,0), mar = c(0,0,0,0))
plot(state $ brat, state $ frat, type="n")
abline(h=1, v = 1)
for (i in which(!(choices $ type %in% c("base-case", "odd-case")) &
                 choices $ srr == "rick")) {
  with(state[i + c(0, 42),], {        
        lines(brat, frat, col = grey(0.5))
        #points(brat, frat, col = grey(0.5), pch = 16, cex = 2)
        text(brat, frat, i %% 21, col = pop + 1, cex = 0.8)
      })
}



stat <- choices
stat $ f <- f.ratio[1,]/f.ratio[2,]
stat $ b <- b.ratio[1,]/b.ratio[2,]


par(oma = c(0,0,0,0), mar = c(0,0,0,0))
plot(state $ brat, state $ frat, type="n")
abline(h=1, v = 1)

make.mat <- function(x) {
  mat <- matrix(NA, 7, 7)
  mat[lower.tri(mat)] <- x[1:21]
  mat <- t(mat)
  mat[lower.tri(mat)] <- x[22:42]
  mat
}

plot.mat <- function(x) 
  image(Matrix(x), sub ="", xlab = "", ylab = "", 
      col.regions = heat.colors(100),
      useRaster = TRUE)

fmat <- make.mat(stat $ f)
bmat <- make.mat(stat $ b)

colnames(fmat) <- rownames(fmat) <- colnames(bmat) <- rownames(bmat) <- 
  c("", paste0(paste0("Linf:", c(60, 80, 100)), rep(c("low", "hi"), each = 3)))

round(fmat - 1, 2)
pdf("tex/fmat.pdf")
plot.mat(fmat)
dev.off()

round(bmat - 1, 2)
pdf("tex/bmat.pdf")
plot.mat(bmat)
dev.off()

plot.mat <- function(col = matrix(NA, 6, 6)) { 
  op <- par()
  par(mar=c(1,1,5,7))
  plot(0, 0, xlim = c(1, 7), ylim = c(1, 7), ann = FALSE, type = "n", axes = FALSE)
  for (i in 1:6)
    for (j in i:6)
      polygon(c(0,1,1,0) + j, c(0,0,1,1) + (7-i), col = col[i,j])
  mtext(rep(c(60, 80, 100), 2), at = 1:6 + .5, line = -.5)
  mtext(c("Low", "High"), at = c(2.5, 5.5), line = 1)
  mtext("Population 2", at = 4, line = 2)
  mtext(rep(c(60, 80, 100), 2), at = rev(1:6 + .5), line = -.5, side = 4, las = 1)
  mtext(c("Low", "High"), at = rev(c(2.5, 5.5)), line = 1.2, side = 4, las = 1)  
  mtext("Population 1", at = 7, line = 1, side = 4, las = 1)
  par(op)
}
mat <- matrix(NA, 6, 6)
diag(mat) <- 1

pdf("tex/scenario0.pdf")
plot.mat()
dev.off()

col <- matrix(NA, 6, 6)
diag(col) <- "seagreen"
pdf("tex/scenario1.pdf")
plot.mat(col = col)
dev.off()

col <- matrix(NA, 6, 6)
col[c(7,13:14)] <- "blue"
col[c(28, 34:35)] <- "red"
pdf("tex/scenario2.pdf")
plot.mat(col = col)
dev.off()

col <- matrix(NA, 6, 6)
col[c(19, 26, 33)] <- "brown"
pdf("tex/scenario3.pdf")
plot.mat(col = col)
dev.off()

col <- matrix(NA, 6, 6)
col[c(25, 31:32)] <- "orange"
pdf("tex/scenario4.pdf")
plot.mat(col = col)
dev.off()

col <- matrix(NA, 6, 6)
diag(col) <- "seagreen"
col[c(7,13:14)] <- "blue"
col[c(28, 34:35)] <- "red"
col[c(19, 26, 33)] <- "brown"
col[c(25, 31:32)] <- "orange"
pdf("tex/scenario5.pdf")
plot.mat(col = col)
dev.off()

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
tsb.ratio <- sapply(1:42, function(i) ratioSummary(comps[[i]], fun = tsb))
catch.ratio <- sapply(1:42, function(i) ratioSummary(comps[[i]], fun = catch))

# TSB and catch proprtions in 2029
tsb.prop <- sapply(1:42, function(i) propSummary(comps[[i]], fun = tsb))
catch.prop <- sapply(1:42, function(i) propSummary(comps[[i]], fun = catch))

choices $ tsb.str <- 1/sapply(1:42, function(i) ratioSummary(comps[[i]], fun = tsb, year = "2000"))
choices $ tsb.end <- 1/sapply(1:42, function(i) ratioSummary(comps[[i]], fun = tsb))

plot(rep(1:2, each = 42), c(choices $ tsb.str, choices $ tsb.end), type = "n")
segments(1, choices $ tsb.str, 2, choices $ tsb.end)


at <- c(.5, 1, 1.5, 2, 3, 4, 6, 8, 10)
p1 <- xyplot(log(tsb) ~ factor(type) | srr, groups = type2, data = choices,
       panel = function(x, y, ...) {
         panel.superpose(x, y, ...)
         panel.abline(h = 0)
       },
       scales = list(y = list(at = log(at), label = at)))

pdf("biomass-ratio.pdf")
print(p1)
dev.off()


save.image("presentation-data.rda")

th <- theme_grey()
th $ panel.background $ fill <- "black"
th $ panel.grid.major $ colour <- grey(0.1)
th $ panel.grid.minor $ colour <- "black"

theme_set( th )


## plot a summary of the BRP
funs <- list(M = "m", Mat = "mat", weight = "stock.wt", selectivity = "landings.sel")
sum <- do.call(rbind, lapply(ASC.brp[1:3], function(brp) 
          do.call(rbind, lapply(funs, function(f)
                    cbind(as.data.frame(do.call(f, list(brp)))[c("age","data")],
                        type = f, Linf = brp @ desc)
              ))))
sum $ Linf <- ordered(gsub("-", "", substring(sum $ Linf, 6, 8)), levels = c(60, 80, 100))
sum <- subset(sum, age <= 10)

p <- ggplot(sum) + geom_line(aes(x = age, y = data, group = Linf, colour = Linf)) + facet_wrap(~type, scale = "free")
p <- p + expand_limits(y = 0) + xlab("age") + ylab("")

pdf("tex/LH-choices1.pdf")
print(p)
dev.off()

ssb <- seq(1, max(sim.design $ v90/.9), length = 100)
rec.df <- expand.grid(ssb = ssb, group = 1:length(ASC.brp))
rec.df $ rec <- unlist( lapply(ASC.brp, function(x) eval( model(x)[[3]], c(list(ssb = ssb), as.list(drop(params(x) @ .Data))))))

rec.df <- subset(rec.df, group %in% c(1, 4, 7, 10))
rec.df $ group <- factor(rec.df $ group)

p <- ggplot(rec.df) + geom_line(aes(x = ssb, y = rec, colour = group, group = group))
p <- p + expand_limits(y = 0, x = 0) + xlab("SSB") + ylab("Recruitment")

pdf("tex/LH-choices2.pdf")
print(p)
dev.off()


plotComp <-
function (x, fn = list(SSB = ssb, Recruits = rec, Yield = catch, F = fbar), 
          facet = facet_wrap(~qname, scale = "free")) 
{
  res = ggplotFL:::whooow(x, fn, 1)
  names(res)[names(res)=="iter"] <- "stock"
  p1 = ggplot(res) + geom_line(aes(x = year, y = data))
  p1 = p1 + expand_limits(y = 0) + xlab("Year") + ylab("") + 
      facet
  p1 = p1 + geom_line(aes(year, data, group = stock, colour = stock), 
        data = transform(as.data.frame(FLQuants(lapply(fn, 
                            function(f, x) f(x), x = x))), 
                stock = factor(iter)))
  print(p1)
  invisible(p1)
}

pdf("tex/stock-unit1.pdf")
tmp <- ASC.stk[[1]]
tmp <- propagate(tmp, 3)
tmp[,,,,,2] <- ASC.stk[[2]]
tmp[,,,,,3] <- ASC.stk[[3]]
plotComp(window(tmp, end = 40))
dev.off()

pdf("tex/stock-unit2.pdf")
tmp <- ASC.stk[[1]]
tmp <- propagate(tmp, 4)
tmp[,,,,,2] <- ASC.stk[[4]]
tmp[,,,,,3] <- ASC.stk[[7]]
tmp[,,,,,4] <- ASC.stk[[10]]
plotComp(window(tmp, end = 70))
dev.off()



pdf("tex/full-result1.pdf")
plotOM(comps[[1]])
dev.off()

pdf("tex/full-result2.pdf")
plotOM(comps[[24]])
dev.off()



## results summary


res <- choices

res $ tsb.prop <- sapply(tsb.prop, function(x) x["mean", "pop1"])
res $ tsb.ratio <- sapply(tsb.ratio, function(x) x["mean"])
res $ catch.prop <- sapply(catch.prop, function(x) x["mean", "pop1"])
res $ catch.ratio <- sapply(catch.ratio, function(x) x["mean"])

    
subset(res, s1 %in% 1)

xyplot(tsb.prop ~ I(s2 %% 6 + 1) | factor(s1), data = res)




