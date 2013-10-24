#' Use a RIVPACS model to calculate metrics
#' 
#' This function predicts the expected number of taxa for a site or sample
#' based on a RIVPACS model.
#' 
#' If newdata is missing, fitted values for the data used to fit the model
#' will be returned.
#' 
#' @param x a RIVPACS model
#' @param newdata a \link{biodiversity} object
#' @param cutoff a probability; taxa with a group occurrence probability < cutoff
#' will not be used to calculate metrics.
#' @param outlier.probs p-value thresholds to report which data points are outliers.
#' Only used if method used to fit the RIVPACS model was 'lda'
#' @output a matrix of RIVPACS metrics
#' @S3method predict rivpacs
#' @method predict rivpacs
#' @export
predict.rivpacs <- function(object, newdata, cutoff = cutoff, outlier.probs = c(0.01, 0.05), ...){
  if (missing(newdata)) {
    bug.pa    <- getSamples(object, pa = T)
    site.data <- getSites(object)
  } else {
    bug.pa    <- getSamples(newdata, pa = T)
    site.data <- getSites(newdata)
  }
  nvars    <- length(getCalibrationVars(object))
  calbug   <- getSamples(object, pa = T)
  calgroup <- getCalibrationGroups(object)
  model    <- getModel(object)
  ngroups  <- length(unique(calgroup))
  bug.pa   <- matchTaxa(bug.pa, calbug)
  
  formula   <- removeIntercept(getFormula(object))
  Terms     <- delete.response(terms(formula))
  site.data <- model.matrix(Terms, site.data)
  if (inherits(model, 'lda')) {
    group.prob <- predictSiteGroup(site.data, getGroupMeans(object),
                                   getGroupInvCov(object), calgroup)
    df <- min(c(nvars, (ngroups - 1)))
    outlier <- outer(group.prob$min.dist, qchisq(1 - outlier.probs, df = df), '>')
    colnames(outlier) <- paste0('outlier.', outlier.probs)
    group.prob <- group.prob$posterior
  } else if (inherits(model, 'randomForest')) {
    group.prob <- predict(model, site.data, type = 'prob')
    outlier    <- NULL
  } else if (inherits(model, 'multinom')) {
    group.prob <- predict(model, site.data, type = 'probs')
    outlier    <- NULL
  }
  E   <- calculateExpected(group.prob, calbug, calgroup)  
  N   <- calculateNullExpected(calbug)
  ans <- calculateRIVPACSMetrics(bug.pa, E, N, cutoff = cutoff)
  cbind(ans, outlier)
}

#' Posterior probability of group membership
#' 
#' Uses results of lda to predict which group the site belongs to based on the 
#' data for the site/sample. Tested! Works! nsite * ngroup
#' @param data a matrix of site predictor variables
#' @param center a matrix of group means for the site calibration data
#' @param invcov a pooled inverse coveraince matrix for the site calibration
#'   data
#' @param group a vector of grouping indicators for the calibration data
#' @return a list:  1) a matrix of posterior probabilities of group
#'   membership 2) the minimum distance of each site from the group centroids
#'   (used to determine whether the site is an outlier)
#' @importFrom plyr aaply
predictSiteGroup <- function(data, center, invcov, group){
  #x <- model.frame(formula, data)
  group.size <- as.numeric(table(group))
  distance <- function(center, x, invcov){
    mahalanobis(x, center, invcov, inverted = T)
  }
  dist <- aaply(center, 1, distance, x = data, invcov = invcov)
  group.prob <- group.size * exp(-0.5 * dist)
  min.dist   <- apply(dist, 2, min)
  ans        <- t(group.prob) / colSums(group.prob)
  return(list(posterior = ans, min.dist = min.dist))
}

#' Probability of occurrence for each species in each group
#' 
#' This function takes a binary presence-absence matrix and a vector
#' of group indicators and returns the occurrence probability of each
#' species in each group.
#' 
#' @inheritParams calculateExpected
#' @return a group * species matrix of probabilities
#' @examples
#' x <- matrix(rbinom(100, size = 1, prob = 0.5), ncol = 10)
#' calculateOccurrenceProb(x, rep(1:2, each = 5))
calculateOccurrenceProb <- function(occur, group){
  if (!is.matrix(occur)) {
    occur <- as.matrix(occur)
  }
  group.sum <- rowsum(occur, group)
  row.sum   <- as.numeric(table(group))
  group.sum / row.sum
}

#' Calculate the expected number of species at a site based on a RIVPACS model
#' 
#' Calculate the expected number of species at a site based on a RIVPACS model.
#' The function uses a matrix of group membership probabilities, the calibration
#' matrix of species occurrence data, and the calibration data group indicators.
#' @param x a site * group matrix of predicted group membership probabilities from a RIVPACS model (for new data)
#' @param occur the site * species binary matrix of calibration data used to fit the RIVPACS model
#' @param group a vector of group indicators for the calibration data
#' @return a nsite (new data) * nspecies (in calibration data) matrix of expected number
#' of species
#' @examples
#' x <- matrix(rbinom(100, size = 1, prob = 0.5), ncol = 10)
#' y <- matrix(c(0.1,0.8,0.9,0.2), ncol = 2)
#' group <- rep(1:2, each = 5)
#' calculateExpected(y, x, group)
calculateExpected <- function(x, occur, group) {
  prob <- calculateOccurrenceProb(occur, group)
  expect <- x %*% prob
  expect
}

#' Calculate the expected number of species from a null model.
#' 
#' Calculates the expected number of species based on a null model (all calibration data
#' is in 1 group).  These are just the probability of occurrence of each species
#' in the calibration data.
#' @inheritParams calculateExpected
#' @examples
#' x <- matrix(rbinom(100, size = 1, prob = 0.5), ncol = 10)
#' calculateNullExpected(x)
calculateNullExpected <- function(occur){
  colMeans(occur)
}
 
#' Calculate RIVPACS metrics
#' 
#' This function calculates the observed-expected ratio for both the RIVPACS
#' model and a null model.  In addition, the Bray-Curtis metric of taxonomic
#' dissimilarity are calculated.
#' 
#' @param O the observed occurrence matrix
#' @param E a matrix of expected probabilities for the RIVPACS model
#' @param N a matrix of expected probabilities for the null model
#' @param cutoff a cutoff value to remove rare speces; species with an
#' expected model probability of less than this value are removed.
#' @examples
#' x <- matrix(rbinom(100, size = 1, prob = 0.5), ncol = 10)
#' y <- matrix(c(0.1,0.8,0.9,0.2), ncol = 2)
#' group <- rep(1:2, each = 5)
#' E <- calculateExpected(y, x, group)
#' N <- calculateNullExpected(x)
#' z <- matrix(rbinom(20, size = 1, prob = 0.5), ncol = 10)
#' calculateRIVPACSMetrics(z, E, N)
calculateRIVPACSMetrics <- function(O, E, N, cutoff = 0.5){
  N <- matrix(N, ncol = ncol(O), nrow = nrow(O), byrow = T)
  null.nonrare <- N >= cutoff
  On  <- O * null.nonrare
  En  <- N * null.nonrare
  Ons <- rowSums(On)
  Ens <- rowSums(En)
  BCn <- rowSums(abs(On - En)) / (Ons + Ens)
  
  nonrare <- E >= cutoff
  O  <- O * nonrare
  E  <- E * nonrare
  Os <- rowSums(O)
  Es <- rowSums(E)
  BC <- rowSums(abs(O - E)) / (Os + Es)
  
  cbind(O = Os, E = Es, OE = Os/Es, BC = BC, 
        Onull = Ons, Enull = Ens, OEnull = Ons / Ens, BCnull = BCn)
}

#' Make the taxa in two different sample * taxa matrices match
#' 
#' Data in the from matrix is expanded to match the taxa occurring in the \code{to} matrix
#' The matrices may (and usually will) have different numbers of rows.  Taxa
#' are matched using \link{match} on the column names, so column names must be set
#' and names must match exactly in order to be matched.  Taxa that are in \code{to}
#' but not in \code{from} will be 0 in the output matrix.
#' @param from a nsample * ntaxa matrix
#' @param to a nsample * ntaxa matrix
#' @return a matrix with the data in the from matrix expanded or dropped so that
#' the taxa are the same as the taxa in the to matrix (i.e., they will have the same
#' number of columns and the same column names)
matchTaxa <- function(from, to){
  stopifnot(is.matrix(from), is.matrix(to))
  i <- match(colnames(to), colnames(from))
  ans             <- from[, i]
  colnames(ans)   <- colnames(to)
  ans[is.na(ans)] <- 0
  return(ans)
}
