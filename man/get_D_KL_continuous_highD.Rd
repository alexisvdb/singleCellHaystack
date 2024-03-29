% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/haystack_continuous.R
\name{get_D_KL_continuous_highD}
\alias{get_D_KL_continuous_highD}
\title{Calculates the Kullback-Leibler divergence between distributions for the high-dimensional continuous version of haystack.}
\usage{
get_D_KL_continuous_highD(
  weights,
  density.contributions,
  reference.prob,
  pseudo = 0
)
}
\arguments{
\item{weights}{A numerical vector with expression values of a gene.}

\item{density.contributions}{A matrix of density contributions of each cell (rows) to each center point (columns).}

\item{reference.prob}{A reference distribution to calculate the divergence against.}

\item{pseudo}{A pseudocount, used to avoid log(0) problems.}
}
\value{
A numerical value, the Kullback-Leibler divergence
}
\description{
Calculates the Kullback-Leibler divergence between distributions for the high-dimensional continuous version of haystack.
}
