% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/haystack.R
\name{haystack_2D}
\alias{haystack_2D}
\title{The main Haystack function, for 2-dimensional spaces.}
\usage{
haystack_2D(
  x,
  y,
  detection,
  use.advanced.sampling = NULL,
  dir.randomization = NULL
)
}
\arguments{
\item{x}{x-axis coordinates of cells in a 2D representation (e.g. resulting from PCA or t-SNE)}

\item{y}{y-axis coordinates of cells in a 2D representation}

\item{detection}{A logical matrix showing which genes (rows) are detected in which cells (columns)}

\item{use.advanced.sampling}{If NULL naive sampling is used. If a vector is given (of length = no. of cells) sampling is done according to the values in the vector.}

\item{dir.randomization}{If NULL, no output is made about the random sampling step. If not NULL, files related to the randomizations are printed to this directory.}
}
\value{
An object of class "haystack"
}
\description{
The main Haystack function, for 2-dimensional spaces.
}
