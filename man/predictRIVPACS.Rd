\name{predictRIVPACS}
\alias{predictRIVPACS}
\title{Use a RIVPACS model to calculate metrics}
\usage{
  predictRIVPACS(mod, newdata, cutoff = cutoff,
    outlier.probs = c(0.01, 0.05))
}
\arguments{
  \item{x}{a biodiversity object}

  \item{mod}{a RIVPACS model}

  \item{cutoff}{a probability; taxa with a group occurrence
  probability < cutoff will not be used to calculate
  metrics.}

  \item{outlier.probs}{p-value thresholds to report which
  data points are outliers.}
}
\description{
  This function predicts the expected number of taxa for a
  site or sample based on a RIVPACS model.
}

