% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/beta_diversity_fun.R
\name{abundance_calc}
\alias{abundance_calc}
\title{Calculate Abundance}
\usage{
abundance_calc(phylo, taxlev, ntaxa)
}
\arguments{
\item{phylo}{phyloseq object}

\item{taxlev}{character. taxonomic level}

\item{ntaxa}{numeric. Number of most abundant taxa, which should be
displayed. Remaining taxa will be summarised in "Others".}
}
\value{
data.frame
}
\description{
Calculate abundance in regards to taxonomic level.
}
\examples{
data("mice_B6_N")
abundance_calc(mice_B6_N, "Family", 14)
}
