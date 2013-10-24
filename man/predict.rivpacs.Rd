\name{predict.rivpacs}
\alias{predict.rivpacs}
\title{Use a RIVPACS model to calculate metrics}
\usage{
  \method{predict}{rivpacs} (object, newdata,
    cutoff = cutoff, outlier.probs = c(0.01, 0.05), ...)
}
\arguments{
  \item{x}{a RIVPACS model}

  \item{newdata}{a \link{biodiversity} object}

  \item{cutoff}{a probability; taxa with a group occurrence
  probability < cutoff will not be used to calculate
  metrics.}

  \item{outlier.probs}{p-value thresholds to report which
  data points are outliers. Only used if method used to fit
  the RIVPACS model was 'lda'}
}
\description{
  This function predicts the expected number of taxa for a
  site or sample based on a RIVPACS model.
}
\details{
  If newdata is missing, fitted values for the data used to
  fit the model will be returned.
}

