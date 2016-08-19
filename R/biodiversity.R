#' Create a biodiversity data object
#' 
#' This function creates a simple data structure to hold biodiversity data, 
#' including taxa occurrence data, data about the sites or samples that the taxa
#' were observed at, data on the traits of the taxa observed, and a \link{phylo4}
#' object (from the \link{phylobase} package) containing the taxa 
#' phylogeny.
#' 
#' Many R packages and functions for phylogenetic or diversity analyses require 
#' the data to be in formats something like these.  This function ensures that
#' the taxa and sites match among all the elements and provides warnings (or
#' errors in some cases) if they don't.
#' 
#' @param sample a site/sample identifier describing the unit of data collection
#'   for the study.
#' @param taxon  a taxon identfier that uniquely describes each taxon.
#' @param count
#' @param site a two sided formula describing the site data to use.  The left 
#'   hand side should contain the site/sample identifiers that match the
#'   \code{sample} argument.  The right hand side contains the variables to use.
#' @param site.data a data.frame containing site data
#' @param trait a two sided formula describing the traits to store from the
#'   \code{trait.data} data.frame. The left hand side contains taxa identifiers that match
#'   the \code{taxon} argument.
#' @param trait.data a data.frame containing trait data
#' @param phylo a phylo4 object containing a phylogeny for the the taxa present 
#'   in the taxon variable
#' @importFrom phylobase subset tipLabels
#' @export
biodiversity <- function(sample, taxon, count, 
                         site = NULL, site.data = NULL,
                         trait = NULL, trait.data = NULL,
                         phylo = NULL){
  
  sample  <- castSampleData(sample, taxon, count)
  samples <- rownames(sample)
  taxa    <- colnames(sample)
  if (is.null(trait)) {
    trait <- data.frame(row.names = taxa)
  } else {
    trait <- responseToRownames(trait, trait.data)
  }
  if (is.null(site)) {
    site <- data.frame(row.names = samples)
  } else {
    site  <- responseToRownames(site, site.data)
  }
  if (!is.null(phylo)){
    any.missing <- any(!taxa %in% tipLabels(phylo))
    any.extra   <- any(!tipLabels(phylo) %in% taxa)
    if (any.missing) {
      warning('Some taxa are missing from the provided phylogeny')
    }
    if (any.extra) {
      phylo <- subset(phylo, tips.include = taxa)
    }
  }
  
  stopifnot(identical(nrow(sample), nrow(site)),
            identical(ncol(sample), nrow(trait)))
  
  if (!identical(samples, rownames(site))) {
    warning('Sites in the site.data data.frame without corresponding occurrence data are being removed')
    i <- match(rownames(site), samples)
    site <- site[i,]
  }
  if (!identical(taxa, rownames(trait))) {
    warning('Taxa in the taxa.data data.frame without corresponding occurrence data are being removed')
    i <- match(rownames(trait), taxa)
    taxa <- taxa[i,]
  }
  obj <- list(sample = sample, site = site, trait = trait)
  class(obj) <- 'biodiversity'
  return(obj)
}

#'@export
getSamples <- function(x, ...){
  UseMethod('getSamples')
}

#'@export
getSites <- function(x, ...){
  UseMethod('getSites')
}

#'@export
getPhylo4 <- function(x, ...){
  UseMethod("getPhylo4")
}

#'@export
getTraits <- function(x, ...){
  UseMethod('getTraits')
}

#'@rdname biodiversity
#' @param x a \code{biodiversity} object
#' @param pa logical; if T presence absence data is returned rather than counts
#' @S3method getSamples biodiversity
#' @method getSamples biodiversity
getSamples.biodiversity <- function(x, pa = F, ...){
  retval <- x$sample
  if (pa) {
    retval <- (retval > 0) * 1
  }
  return(retval)
}

#'@rdname biodiversity
#'@S3method getSites biodiversity
#'@method getSites biodiversity
#'@export
getSites.biodiversity <- function(x, ...){
  retval <- x$site
  if (is.null(retval)) {
    stop("biodiversity object contains no site data!")
  }
  return(retval)
}

#'@rdname biodiversity
#'@S3method getSites biodiversity
#'@method getSites biodiversity
getPhylo4.biodiversity <- function(x, ...){
  retval <- x$phylo
  if (is.null(retval)){
    stop("biodiversity object contains no phylogeny data!")
  }
}

#'@rdname biodiversity
#'@S3method getTraits biodiversity
getTraits.biodiversity <- function(x, ...){
  retval <- x$trait
  if (is.null(retval)) {
    stop("biodiversity object contains no site data!")
  }
  return(retval)
}

#'@rdname biodiversity
#'@S3method print biodiversity
#'@method print biodiversity
print.biodiversity <- function(x, ...){
  cat('A biodiversity object:\n')
  samp <- getSamples(x)
  cat(ncol(samp), 'taxa:', toString(head(colnames(samp))), '...\n')
  cat(nrow(samp), 'samples:', toString(head(rownames(samp))), '...\n')
}

#' Create a sample by taxa matrix from data in long format
#' 
#' Creates a sample by taxa matrix from data in long format
#' @param sample a vector of sample/site indicatros
#' @param taxon a vector of taxa indicators
#' @param count a vector of counts (or 0/1 for presence absence data)
#' @param transform.count logical, if T the counts are transformed to presence absence data
#' @importFrom reshape2 acast
#' @examples
#' i <- sample(1:nrow(calibration), 10)
#' pred.data <- calibration[calibration$sample %in% c(1:5),]
#' setUpTaxonData(pred.data$sample, pred.data$taxon, pred.data$is.present,
#' calibration$group, calibration$sample, calibration$taxon, calibration$is.present)
castSampleData <- function(sample, taxon, count, transform.count = F){
  if (transform.count) {
    count <- count > 0
  }
  x  <- data.frame(sample, taxon, count, stringsAsFactors = T)
  ans <- acast(x, sample ~ taxon, value.var = 'count', fill = 0, drop = F, fun.aggregate = sum)
  return(ans)
}


#' Collapse sample data up to a list of higher level taxa.
#' 
#' Collapses the sample data in a biodiversity object up to a list of higher level taxa. For example,
#' this allows sample data of multiple taxonomic resolutions to be integrated.
#' @param sample a vector of sample/site indicatros
#' @param taxon a vector of taxa indicators
#'@importFrom igraph shortest.paths
#'@importFrom ape as.igraph
#'@export
collapseTaxa <- function(x, taxa){
  graph   <- as.igraph(as(getPhylo4(x), 'phylo'))
  samples <- getSamples(x)
  taxa    <- colnames(samples)
  in.path <- shortest.paths(graph, v = taxa, to = nodes, mode = 'in') < Inf
  x$samples <- apply(in.path, 2, function(i) rowSums(samples[,i, drop = F]))
  x
}
