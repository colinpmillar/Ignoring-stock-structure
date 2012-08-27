#' Methods to convert objects to numeric
#'
#' currently implemented are character and factor to numeric
#' the factor to numeric assumes that the factor levels are
#' numerics stored as characters
#'
#' @param object the thing you want to convert
#'
#' @return something the same class and dimension as \code{object} 
#'         but with numeric contents
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname an-methods
#'
#' @examples
#' an(paste(1:10))
#' an(factor(paste(1990:2010))
setGeneric("an",              function(object) standardGeneric("an")            )

#' @rdname an-methods
#' @aliases an, character
setMethod ("an", "character", function(object) as.numeric(object)               )

#' @rdname an-methods
#' @aliases an, factor
setMethod ("an", "factor"   , function(object) as.numeric(as.character(object)) )
