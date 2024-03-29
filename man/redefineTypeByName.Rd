% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/retyping.R
\name{redefineTypeByName}
\alias{redefineTypeByName}
\title{Retype neurons according to a grep pattern to be run on names}
\usage{
redefineTypeByName(
  connections,
  typeList = NULL,
  pattern = NULL,
  sets = NULL,
  nameModifiers,
  postfix = c("raw", "to", "from"),
  redefinePartners = FALSE,
  perl = FALSE
)
}
\arguments{
\item{connections}{Connectivity table  or neuronBag to modify}

\item{typeList}{: which types to lateralize (by default all the neurons which names
contains L or R)}

\item{pattern}{A grep pattern to be run to match a "name" column}

\item{sets}{List of neuron name sets. To be used in place of pattern.}

\item{nameModifiers}{A vector of strings of lenght 2 if \code{pattern} is used : first string is appended for matched types, second to the absence of matches
(for types in \code{typeList}), or of the same lenght as sets: postfixes
to be appended to the type name of all neurons whose name is in the corresponding set}

\item{postfix}{One of "raw, "to" or "from". Specify if type (and name) columns in table to be modified are postfixed with
to and from or nothing}

\item{redefinePartners}{If table is a neuronBag, should the partners also be retyped?}

\item{perl}{Should the grep match use perl rules?}
}
\description{
Retype neurons according to a grep pattern to be run on names
}
\details{
\code{pattern}, \code{typeList} and \code{perl} are used if \code{pattern} is used. Alternatively \code{sets} allow to retype according to an arbitrary number
of subtypes
}
\examples{
\dontrun{
PFLNames <- getTypesTable(c("PFL1","PFL2","PFL3"))

## Rename only PFL3 according to their strange output pattern
PFLNames3 <- redefineTypeByName(PFLNames,typeList = c("PFL3"),pattern = "(^.*_L(?!.*irreg))|(^.*_R.*irreg)",perl=TRUE,nameModifiers = c("_L*","_R*"))
}
}
