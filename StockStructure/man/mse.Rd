\name{mse}
\alias{mse}
\title{runs an simulation of a stock, stock assessment using a simulated index
and stock recruit function.}
\usage{
  mse(OM, iniyr, sr.model1, sr.model2,
    survey.q = rep(1e+06, dims(OM)$age),
    sr.residuals = FLQuant(1, dimnames = dimnames(window(rec(OM), start = iniyr))),
    CV = 0.15, Ftar = 1, Btrig = 1, refpt, seed = NULL,
    ...)
}
\arguments{
  \item{OM}{an FLStock object}

  \item{start}{the first year of projections}
}
\value{
  an FLStock object
}
\description{
  The stock is projected using the function \code{fwd} from
  the FLash library the assessment model is the a4a model
  from the FLa4a package while the management procedure is
  a Harvest Control Rule (HCR) implemented using the
  function \code{hcr} from the FLAdvice package.
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

