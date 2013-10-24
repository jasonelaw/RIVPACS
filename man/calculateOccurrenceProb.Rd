\name{calculateOccurrenceProb}
\alias{calculateOccurrenceProb}
\title{Probability of occurrence for each species in each group}
\usage{
  calculateOccurrenceProb(occur, group)
}
\arguments{
  \item{occur}{the site * species binary matrix of
  calibration data used to fit the RIVPACS model}

  \item{group}{a vector of group indicators for the
  calibration data}
}
\value{
  a group * species matrix of probabilities
}
\description{
  This function takes a binary presence-absence matrix and
  a vector of group indicators and returns the occurrence
  probability of each species in each group.
}
\examples{
x <- matrix(rbinom(100, size = 1, prob = 0.5), ncol = 10)
calculateOccurrenceProb(x, rep(1:2, each = 5))
}

