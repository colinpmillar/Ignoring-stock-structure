
library(devtools)
library(testthat)
library(roxygen2)

roxygenize("../StockStructure")
pkg <- as.package("../StockStructure")
build(pkg)
#check(pkg)
install(pkg)

library(StockStructure)

helloWorld(882346)

data(wklife.stk)



# Jobs