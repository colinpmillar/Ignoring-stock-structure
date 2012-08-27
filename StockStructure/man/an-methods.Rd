\docType{methods}
\name{an}
\alias{an}
\alias{an,}
\alias{character}
\alias{factor}
\title{Methods to convert objects to numeric}
\arguments{
  \item{object}{the thing you want to convert}
}
\value{
  something the same class and dimension as \code{object}
  but with numeric contents
}
\description{
  currently implemented are character and factor to numeric
  the factor to numeric assumes that the factor levels are
  numerics stored as characters
}
\examples{
an(paste(1:10))
an(factor(paste(1990:2010))
}
\seealso{
  \code{\link{print}} and \code{\link{cat}}
}

