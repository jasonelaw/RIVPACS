\name{calculateNullExpected}
\alias{calculateNullExpected}
\title{Calculate the expected number of species from a null model.}
\usage{
  calculateNullExpected(occur)
}
\arguments{
  \item{occur}{the site * species binary matrix of
  calibration data used to fit the RIVPACS model}
}
\description{
  Calculates the expected number of species based on a null
  model (all calibration data is in 1 group).  These are
  just the probability of occurrence of each species in the
  calibration data.
}
\examples{
x <- matrix(rbinom(100, size = 1, prob = 0.5), ncol = 10)
calculateNullExpected(x)
}

