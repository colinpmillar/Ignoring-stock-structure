\name{getFunctions}
\alias{getFunctions}
\title{get methods defined in a generic function}
\usage{
  getFunctions(x, show = TRUE)
}
\arguments{
  \item{x}{the name of a generic function as a character}
}
\value{
  a names list of functions.  The names are the method
  dispatch types
}
\description{
  get methods defined in a generic function
}
\examples{
getFunctions("show")
}
\author{
  Colin Millar \email{colin.millar@jrc.ec.europa.eu}
}

