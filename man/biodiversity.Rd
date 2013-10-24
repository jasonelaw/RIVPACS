\name{biodiversity}
\alias{biodiversity}
\alias{getPhylo4.biodiversity}
\alias{getSamples.biodiversity}
\alias{getSites.biodiversity}
\alias{getTraits.biodiversity}
\alias{print.biodiversity}
\title{Create a biodiversity data object}
\usage{
  biodiversity(sample, taxon, count, site = NULL,
    site.data = NULL, trait = NULL, trait.data = NULL,
    phylo = NULL)

  \method{getSamples}{biodiversity} (x, pa = F, ...)

  \method{getSites}{biodiversity} (x, ...)

  \method{getSites}{biodiversity} (x, ...)

  getTraits.biodiversity(x, ...)

  \method{print}{biodiversity} (x, ...)
}
\arguments{
  \item{sample}{a site/sample identifier describing the
  unit of data collection for the study.}

  \item{taxon}{a taxon identfier that uniquely describes
  each taxon.}

  \item{count}{}

  \item{site}{a two sided formula describing the site data
  to use.  The left hand side should contain the
  site/sample identifiers that match the \code{sample}
  argument.  The right hand side contains the variables to
  use.}

  \item{site.data}{a data.frame containing site data}

  \item{trait}{a two sided formula describing the traits to
  store from the \code{trait.data} data.frame. The left
  hand side contains taxa identifiers that match the
  \code{taxon} argument.}

  \item{trait.data}{a data.frame containing trait data}

  \item{phylo}{a phylo4 object containing a phylogeny for
  the the taxa present in the taxon variable}

  \item{x}{a \code{biodiversity} object}

  \item{pa}{logical; if T presence absence data is returned
  rather than counts}
}
\description{
  This function creates a simple data structure to hold
  biodiversity data, including taxa occurrence data, data
  about the sites or samples that the taxa were observed
  at, data on the traits of the taxa observed, and a
  \link{phylo4} object (from the \link{phylobase} package)
  containing the taxa phylogeny.
}
\details{
  Many R packages and functions for phylogenetic or
  diversity analyses require the data to be in formats
  something like these.  This function ensures that the
  taxa and sites match among all the elements and provides
  warnings (or errors in some cases) if they don't.
}

