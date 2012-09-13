#' Hello world
#'
#' 
#' @param x random seed
#' @return Some text
#' @note \code{helloWorld} is just a test
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
#' @examples
#' helloWorld(3435)
helloWorld <- function(x = 1) {
	set.seed(x)
  
  cat("hello there ", sample(100, 1), "!\n", sep = "")
  invisible(NULL)
}

#' get methods defined in a generic function
#'
#' 
#' @param x the name of a generic function as a character
#' @return a names list of functions.  The names are the method dispatch types
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
#' @examples
#' getFunctions("show")
getFunctions <- function(x, show = TRUE) 
{
	funs <- lapply(findMethods(x) @ .Data, slot, ".Data")
	names(funs) <- sapply(findMethods(x) @ .Data, function(x) paste(slot(x, "target") @ .Data, collapse = ", "))
	if (show) funs else names(funs)
}

#' seconds to hrs, mins, secs
#'
#' 
#' @param x seconds
#' @return Some text
#' @note \code{helloWorld} is just a test
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
#' @examples
#' helloWorld(3435)
frmtSeconds <- function(x = 1) {

  hrs <- floor(x/60^2)
  mins <- floor(x/60 - hrs * 60)
  secs <- x %% 60
  paste(hrs, "hrs ", mins, "mins ", secs, "s", sep = "")
}



#' build a4a and stock structure package
#'
#' 
#' @return NULL
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
#' @examples
#' x <- 1
mybuild <- function() {
  require(devtools)
  require(testthat)
  require(roxygen2)
  
  #roxygenize("../../a4a/packages/FLa4a")
  
  #pkg <- as.package("../../a4a/packages/FLa4a")
  #build(pkg)
  #check(pkg)
  #install(pkg)
  
  roxygenize("../StockStructure")
  pkg <- as.package("../StockStructure")
  build(pkg)
  #check(pkg)
  install(pkg)
  
  invisible(NULL)
}