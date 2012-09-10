\name{doOne}
\alias{doOne}
\title{wrapper for mse function - specific for this task!}
\usage{
  doOne(stock.id, Ftarg = 0.75, Btrig = 0.75,
    start.yr = c(1, 1), which.ref = "f0.1")
}
\arguments{
  \item{stock.id}{an FLStock object}

  \item{Ftarg}{the first year of projections}

  \item{Btrig}{}

  \item{start.yr}{}
}
\value{
  an FLStock object
}
\description{
  wrapper for mse function - specific for this task!
}
\note{
  \code{mse} is based on code written by Ernesto Jardim
  \email{ernesto.jardim@jrc.ec.europa.eu}
}
\author{
  Colin Millar \email{colin.millar@jrc.ec.europa.eu}
}
\seealso{
  \code{\link{FLStock}}, \code{\link{fwd}},
  \code{\link{hcr}}
}

