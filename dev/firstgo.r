
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

getFunctions <- function(x, show = TRUE) 
{
  funs <- lapply(findMethods(x) @ .Data, slot, ".Data")
  names(funs) <- sapply(findMethods(x) @ .Data, function(x) paste(slot(x, "target") @ .Data, collapse = ", "))
  if (show) funs else names(funs)
}

getFunctions("fwd")

# Jobs