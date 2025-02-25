% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGeneric.R
\name{evalDiffMat}
\alias{evalDiffMat}
\title{create a normalized height function between two HiC matrices}
\usage{
evalDiffMat(object, useLog2 = TRUE, include.zeros = TRUE,
  return_excluded_mat = FALSE)
}
\arguments{
\item{object}{: a \code{treeHiCDataSet} object}

\item{useLog2}{: Logical, default to TRUE to use log2 of two matrices}

\item{include.zeros}{: should partial zero interactions be included? Defaul to TRUE}
}
\value{
a vector with range from 0 to 1 storing height difference
}
\description{
create a normalized height function between two HiC matrices
}
