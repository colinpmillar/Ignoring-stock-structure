#' Methods to collpase the units of an object
#'
#' currently implemented are collapsing accross
#' units of an FLStock
#'
#' @param object the thing you want to collapse
#'
#' @return something the same class and dimension as \code{object} 
#'         but with the unit dimension equal to 1
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname collapseUnit-methods
setGeneric("collapseUnit", function(object) standardGeneric("collapseUnit")            )

#' @rdname collapseUnit-methods
#' @aliases collapseUnit, FLStock
setMethod ("collapseUnit", "FLStock", function(object) {
  # simply go through each one doing the right thing!
  out <- qapply(object, unitSums)
  
  catch.wt(out) <- unitSums( catch.wt(object) * catch.n(object) ) / catch.n(out)
  discards.wt(out) <- unitSums( discards.wt(object) * discards.n(object) ) / discards.n(out)
  landings.wt(out) <- unitSums( landings.wt(object) * landings.n(object) ) / landings.n(out)
  stock.wt(out) <- unitSums( stock.wt(object) * stock.n(object) ) / stock.n(out)

  # these are not right - but i am only using one M and one mat
  # and the rest is not used
  m(out) <- unitMeans(m(object))
  mat(out) <- unitMeans(mat(object))
  harvest(out) <- unitMeans(harvest(object))
  FLCore::units( harvest(out) ) <- "f"
  harvest.spwn(out) <- unitMeans(harvest.spwn(object))
  m.spwn(out) <- unitMeans(m.spwn(object))

  out
})
