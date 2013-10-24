\name{buildRIVPACS}
\alias{buildRIVPACS}
\alias{getCalibrationData}
\alias{getCalibrationData.rivpacs}
\alias{getCalibrationGroups}
\alias{getCalibrationGroups.rivpacs}
\alias{getCalibrationVars}
\alias{getCalibrationVars.rivpacs}
\alias{getFormula}
\alias{getFormula.rivpacs}
\alias{getGroupInvCov}
\alias{getGroupInvCov.rivpacs}
\alias{getGroupMeans}
\alias{getGroupMeans.rivpacs}
\alias{getModel}
\alias{getModel.rivpacs}
\alias{getSamples.rivpacs}
\alias{getSites.rivpacs}
\title{Build a RIVPACS model}
\usage{
  buildRIVPACS(formula, data,
    type = c("lda", "multinom", "rf"), ...)

  \method{getSamples}{rivpacs} (x, ...)

  \method{getSites}{rivpacs} (x, ...)

  getCalibrationData(x, ...)

  getGroupMeans(x, ...)

  getGroupInvCov(x, ...)

  getCalibrationGroups(x, ...)

  getCalibrationVars(x, ...)

  getFormula(x, ...)

  getModel(x, ...)

  \method{getCalibrationData}{rivpacs} (mod)

  \method{getGroupMeans}{rivpacs} (mod)

  \method{getGroupInvCov}{rivpacs} (mod)

  \method{getCalibrationGroups}{rivpacs} (mod)

  \method{getCalibrationVars}{rivpacs} (mod)

  \method{getFormula}{rivpacs} (mod)

  \method{getModel}{rivpacs} (mod)
}
\arguments{
  \item{formula}{a formula with the groups on the left and
  the predictors on the right}

  \item{data}{the data frame of site/sample data}

  \item{type}{the type of model to build}

  \item{\ldots}{additional arguments to be passed to the
  modeling functions}
}
\value{
  a RIVPACS model object (a function)
}
\description{
  This function builds a RIVPACS model using either linear
  discriminant analysis, random forests, or a multinomial
  log linear model.  The get* functions can be used to
  extract model components.
}

