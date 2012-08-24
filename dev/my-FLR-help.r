
library(FLCore)
library(FLAssess)
library(FLash)
library(FLAdvice)
library(FLBioDym)
library(plyr)
library(reshape)

slotNames(FLQuant())
slotNames(FLStock())
slotNames(FLIndex())
slotNames(FLBRP())
slotNames(FLModel())
slotNames(FLPar())
slotNames(FLSR())

slotNames(fwdControl(data.frame(year = 1:10)))