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
