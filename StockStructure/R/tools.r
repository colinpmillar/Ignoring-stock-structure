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