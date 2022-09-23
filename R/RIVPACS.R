
# For reading tree dumped from NCBI taxonomy
# use readNewick from phylobase

removeIntercept <- function(x){
  update.formula(x, . ~ . + 0)
}

responseToRownames <- function(formula, data){
  formula <- removeIntercept(formula)
  mf   <- model.frame(formula, data)
  resp <- model.extract(mf, 'response')
  mat  <- model.matrix(formula, mf)
  o   <- order(resp)
  mat <- mat[o, ]
  rownames(mat) <- resp[o]
  return(as.data.frame(mat))
}

clusterTaxa <- function(occur) {
  
}

#'@importFrom plyr alply
#'@importFrom MASS lda
#'@importFrom subselect ldaHmat eleaps
selectLDA <- function(formula, data, k, criterion = 'tau2', nsol = 1, ...) {
  frm <- model.frame(formula, data)
  groups <- model.extract(frm, 'response')
  groups <- as.factor(groups)
  mat    <- model.matrix(removeIntercept(formula), data)
  Hmat   <- ldaHmat(mat, groups)
  vsubset <- eleaps(Hmat$mat, kmin = min(k), kmax = max(k), H = Hmat$H, 
                    r = Hmat$r, criterion = criterion, nsol = nsol)
  fitLDA <- function(x, ...){
    form <- sprintf(". ~ %s", paste0(colnames(mat)[as.numeric(x)], collapse = ' + '))
    form <- update(formula, form)
    model  <- lda(form, data = data, CV = F, ...)
    cv     <- lda(form, data = data, CV = T, ...)
    list(model = model, cv = cv, predict = predict(model, data))
  }
  mods <- alply(vsubset$subsets, c(1,3), fitLDA, .dims = T)
  attr(mods, 'select') <- vsubset
  mods
}

#' Resample a vector
#' 
#' This function adapted from the resample function in \code{?sample}. It silently
#' returns the whole vector in random order when size > length(x).
#' @param x a vector
#' @param n sample size
#' @return a vector of length n with the sample
resample <- function(x, size = length(x), ...){
  n <- length(x)
  x[sample.int(n, size = pmin.int(n, size), ...)]
} 

#'Subsample samples of organisms to a fixed count
#'
#'Subsample samples of organizms to a fixed count.
#'@param sample a vector of sample identifiers
#'@param taxon a vector of taxon identifiers
#'@param count a vector of counts for each taxa in each sample
#'@param n the size of the subsample to select
#'@examples
#'d <- expand.grid(sample = 1:2, 
#'                 taxon = c('Ephemeroptera', 'Plecoptera', 'Trichoptera'))
#'d$count <- runif(6, 10, 20)
#'rarify(d$sample, d$taxon, d$count, 30)
#'@importFrom plyr ddply
#'@export
rarify <- function(sample, taxon, count, n){
  doRarify <- function(x, n){
    individuals <- rep(x$taxon, x$count)
    ans <- resample(individuals, size = n, replace = F)
    ans <- as.data.frame(table(taxon = ans), responseName = 'count')
  }
  x <- data.frame(sample, taxon, count)
  ddply(x, .(sample), doRarify, n = n)
}