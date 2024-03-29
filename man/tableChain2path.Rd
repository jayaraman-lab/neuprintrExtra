% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/typesPathways.R
\name{tableChain2path}
\alias{tableChain2path}
\title{Link an arbitrary number of connection tables into a pathway table}
\usage{
tableChain2path(
  ...,
  n_steps = NULL,
  stat = NULL,
  excludeLoops = TRUE,
  type.to = NULL,
  excludeIntermediate = NULL,
  addRecursive = FALSE,
  chunkPath = FALSE,
  progress = TRUE
)
}
\arguments{
\item{...}{An arbitrary number of connection tables or a list of connection tables. In desired downstream order}

\item{stat}{which stat to use for pathway weights. Can be a vector. If NULL (the default), all the available meaningful stats are used: 
(\code{"weightRelative","outputContribution","knownWeightRelative","knownOutputContribution","knownOutputContribution_perType","knownWROutputContribution_perType","knownWeightRelative_perType"})}

\item{excludeLoops}{Logical, whether to exclude paths containing duplicates (TRUE by default)}

\item{type.to}{a metadata data.frame for a set of neurons (typically generated by a search or a getTypesTable call) to arrive to}

\item{addRecursive}{Whether or not to add the recursive paths (of length greater than \code{n_steps}) present in the connection tables generated.}

\item{chunkPath}{Chunking argument to be passed to the function chaining the connection tables together. Useful for deep pathways. Can either be logical,
or an integer specifying the number of connections in the starting table to process at once.}
}
\description{
Link an arbitrary number of connection tables into a pathway table
}
