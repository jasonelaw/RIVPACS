#' Build a RIVPACS model
#' 
#' This function builds a RIVPACS model using either linear discriminant analysis,
#' random forests, or a multinomial log linear model.  The get* functions can be
#' used to extract model components.
#' 
#' @param formula a formula with the groups on the left and the predictors on the right
#' @param data the data frame of site/sample data
#' @param type the type of model to build
#' @param \ldots additional arguments to be passed to the modeling functions
#' @return a RIVPACS model object (a function)
#' @importFrom MASS lda
#' @importFrom nnet multinom
#' @importFrom randomForest randomForest
#' @export
buildRIVPACS <- function(formula, data, type = c('lda', 'multinom', 'rf'), ...){
  type <- match.arg(type)
  site   <- getSites(data)
  sample <- getSamples(data)
  model <- switch(type,
                  lda      = lda(formula, site, ...),
                  multinom = multinom(formula, site, ...),
                  rf       = randomForest(formula, site, ...))
  mf     <- model.frame(formula, site)
  datmat <- model.matrix(removeIntercept(formula), site)
  groups <- model.extract(mf, 'response')
  names(groups) <- rownames(sample)
  
  group.means  <- rowsum(datmat, groups) / as.numeric(table(groups))
  covpool.inv  <- pooledCovariance(datmat, groups, inverse = T)
  
  model.info <- list(data           = data, 
                     groups         = groups,
                     group.means    = group.means, 
                     covpool.inv    = covpool.inv,
                     variables.used = colnames(group.means),
                     formula        = formula, 
                     model          = model)
#   f <- function(){
#     return(model.info)
#   }
  class(model.info) <- 'rivpacs'
  return(model.info)
}

#' @rdname buildRIVPACS
#' @S3method getSamples rivpacs
#' @method getSamples rivpacs
#' @export
getSamples.rivpacs <- function(x, ...){
  x <- getCalibrationData(x)
  getSamples(x, ...)
}

#' @rdname buildRIVPACS
#' @S3method getSites rivpacs
#' @method getSites rivpacs
#' @export
getSites.rivpacs <- function(x, ...){
  x <- getCalibrationData(x)
  getSites(x, ...)
}

#' @rdname buildRIVPACS
#' @export
getCalibrationData <- function(x, ...) UseMethod("getCalibrationData")

#' @rdname buildRIVPACS
#' @export
getGroupMeans <- function(x, ...) UseMethod("getGroupMeans")

#' @rdname buildRIVPACS
#' @export
getGroupInvCov <- function(x, ...) UseMethod("getGroupInvCov")

#' @rdname buildRIVPACS
#' @export
getCalibrationGroups <- function(x, ...) UseMethod("getCalibrationGroups")

#' @rdname buildRIVPACS
#' @export
getCalibrationVars <- function(x, ...) UseMethod("getCalibrationVars")

#' @rdname buildRIVPACS
#' @export
getFormula <- function(x, ...) UseMethod("getFormula")

#' @rdname buildRIVPACS
#' @export
getModel <- function(x, ...) UseMethod("getModel")

#' @rdname buildRIVPACS
#' @S3method getCalibrationData rivpacs
#' @method getCalibrationData rivpacs
#' @export
getCalibrationData.rivpacs <- function(mod){
  mod$data
}

#' @rdname buildRIVPACS
#' @S3method getGroupMeans rivpacs
#' @method getGroupMeans rivpacs
#' @export
getGroupMeans.rivpacs <- function(mod){
  mod$group.means
}

#' @rdname buildRIVPACS
#' @S3method getGroupInvCov rivpacs
#' @method getGroupInvCov rivpacs
#' @export
getGroupInvCov.rivpacs <- function(mod){
  mod$covpool.inv
}

#' @rdname buildRIVPACS
#' @S3method getCalibrationGroups rivpacs
#' @method getCalibrationGroups rivpacs
#' @export
getCalibrationGroups.rivpacs <- function(mod){
  mod$groups
}

#' @rdname buildRIVPACS
#' @S3method getCalibrationVars rivpacs
#' @method getCalibrationVars rivpacs
#' @export
getCalibrationVars.rivpacs <- function(mod){
  mod$variables.used
}

#' @rdname buildRIVPACS
#' @S3method getFormula rivpacs
#' @method getFormula rivpacs
#' @export
getFormula.rivpacs <- function(mod){
  mod$formula
}

#' @rdname buildRIVPACS
#' @S3method getModel rivpacs
#' @method getModel rivpacs
#' @export
getModel.rivpacs <- function(mod){
  mod$model
}

#' @importFrom plyr dlply '.'
pooledCovariance <- function(x, group, inverse = FALSE){
  group.size <- table(group)
  df         <- as.data.frame(x)
  df$group   <- group
  covlist <- dlply(df, .(group), function(x){
    x <- x[, setdiff(names(x), 'group')]
    cov(x) * (nrow(x) - 1)
    })
  pooled.cov <- Reduce('+', covlist) / (sum(group.size) - length(group.size))
  if (inverse) {
    pooled.cov <- solve(pooled.cov)
  }
  return(pooled.cov)
}
