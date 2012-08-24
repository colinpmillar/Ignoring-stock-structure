library(FLAdvice)
library(plyr)

tiff("fig722.tiff", 6.5, 6.5, units="in", res=250, compression="zip")
xyplot(data~year|par, groups=qtl, data=subset(base.summ, ref=="base"), scales=list(y=list(relation="free")), par.settings=list(superpose.line=list(lty=c(2,1,2), col=c(2,1,2))), type="l", par.strip.text=list(cex=0.5), ylab="")
dev.off()

stks <- FLStocks(i1=iter(base,1),i2=iter(base,2), i3=iter(base,3))
tiff("fig723.tiff", 6.5, 6.5, units="in", res=250, compression="zip")
plot(stks)
dev.off()

tiff("fig724.tiff", 6.5, 6.5, units="in", res=250, compression="zip")
xyplot(data~year|par, groups=ref, data=subset(base.summ, ref %in% c("base", "aLag3", "aLag5") & qtl==0.5), scales=list(y=list(relation="free")), par.settings=list(superpose.line=list(col=c("red", "green", "grey60")), superpose.symbol=list(col=c("red", "green", "grey40"), pch=19, cex=0.2)), type="b", par.strip.text=list(cex=0.5), auto.key=list(lines=TRUE, points=FALSE, columns=3), ylab="")
dev.off()

tiff("fig725.tiff", 6.5, 6.5, units="in", res=250, compression="zip")
xyplot(data~year|par, groups=ref, data=subset(base.summ, ref %in% c("base", "cthBias0.5", "cthBias0.8") & qtl==0.5), scales=list(y=list(relation="free")), par.settings=list(superpose.line=list(col=c("grey60", "red", "green"))), type="l", par.strip.text=list(cex=0.5), auto.key=list(lines=TRUE, points=FALSE, columns=3), ylab="")
dev.off()

tiff("fig726.tiff", 6.5, 6.5, units="in", res=250, compression="zip")
xyplot(data~year|par, groups=ref, data=subset(base.summ, ref %in% c("base", "srvBias0.5", "srvBias1.5") & qtl==0.5), scales=list(y=list(relation="free")), par.settings=list(superpose.line=list(col=c("grey60", "red", "green"))), type="l", par.strip.text=list(cex=0.5), auto.key=list(lines=TRUE, points=FALSE, columns=3), ylab="")
dev.off()

tiff("fig727.tiff", 6.5, 6.5, units="in", res=250, compression="zip")
xyplot(data~year|par, groups=ref, data=subset(base.summ, ref %in% c("base", "fox") & qtl==0.5), scales=list(y=list(relation="free")), par.settings=list(superpose.line=list(col=c("grey60", "red"))), type="l", par.strip.text=list(cex=0.5), ylab="", auto.key=list(lines=TRUE, points=FALSE, columns=2))
dev.off()

tiff("fig728.tiff", 6.5, 6.5, units="in", res=250, compression="zip")
xyplot(data~year|par, groups=ref, data=subset(base.summ, ref %in% c("base", "worstCase") & qtl==0.5), scales=list(y=list(relation="free")), par.settings=list(superpose.line=list(col=c("grey60", "red")), superpose.symbol=list(pch=19, col=c("grey60", "red"), cex=0.2)), type="b", par.strip.text=list(cex=0.5), ylab="", auto.key=list(lines=TRUE, points=FALSE, columns=2))
dev.off()

tiff("fig729.tiff", 6.5, 6.5, units="in", res=250, compression="zip")
xyplot(data~year|par, groups=ref, data=subset(base.summ, ref %in% c("base", "srv", "cth") & qtl==0.5), scales=list(y=list(relation="free")), par.settings=list(superpose.line=list(col=c("grey60", "red", "green"))), type="l", par.strip.text=list(cex=0.5), ylab="", auto.key=list(lines=TRUE, points=FALSE, columns=3))
dev.off()




