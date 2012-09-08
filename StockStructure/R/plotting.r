#' Hello world
#'
#' 
#' @param OM random seed
#' @return null
#' @note \code{helloWorld} is just a test
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
#' @examples
#' helloWorld(3435)
plotOM <- function(OM, probs = c(0.025, .5, .975)) {
  assessment <- attr(OM, "assessment.data") $ stock
  
  myrec <- function (object) stock.n(object)[1,]               
  fn    <- list(SSB = ssb, Recruits = myrec, Yield = catch, F = fbar)
  res1  <- ggplotFL:::whooow(OM[,,1], fn, probs)
  res2  <- ggplotFL:::whooow(OM[,,2], fn, probs)
  res <- rbind(cbind(res1, population = "1"), cbind(res2, population = "2"))
  res $ group <- paste(res $ pop, res $ iter)
  
  p1 <- ggplot(res) + 
      geom_line(aes(x = year, y = data, group = group, 
              size = iter, lty = iter, col = population)) + 
      scale_size_manual(values = c(0.5, 1, 0.5), name = "Quantile") + 
      scale_linetype_manual(values = c(2, 1, 2),  name = "Quantile") +
      expand_limits(y = 0) + xlab("Year") + ylab("") +
      facet_wrap(~ qname, scale = "free")
  
  print(p1)
  invisible(p1)  
}        

#' Hello world
#' 
#' @param OM random seed
#' @return NULL
#' @note \code{helloWorld} is just a test
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
#' @examples
#' helloWorld(3435)
plotAssessment <- function(OM) {
  assessment <- attr(OM, "assessment.data") $ stock
  myrec <- function (object) stock.n(object)[1,]               
  plotComp(assessment, 
      fn = list(SSB = ssb, Recruits = myrec, Yield = catch, F = fbar), 
      probs = c(0.75, 0.5, 0.25), 
      size= c(0.5, 1, 0.5), 
      lty = c(2,1,2), 
      facet = facet_wrap(~qname, scale = "free"))  
}
