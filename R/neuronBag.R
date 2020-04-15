## Define an S3 class, neuronBag, to hold the connectivity information of a bunch of neurons

#' neuronBag objects contain information about the connectivity of a bunch of neurons
#' @field names : a table of metadata associated with the neurons in the bag -- with an extra column, 'databaseType'
#' which keeps the type name used in the database
#' @field outputs : a type to type connectivity table of the outputs of the set of neurons
#' @field inputs : a type to type connectivity table of the inputs of the set of neurons
#' @field outputs_raw : a neuron to neuron connection table of the outputs of the set of neurons
#' @field inputs_raw : a neuron to neuron connection table of the inputs of the set of neurons
#' @field outputsTableRef : a table holding all instances of all the output types and their associated metadata
neuronBag <- function(outputs,inputs,names,outputs_raw,inputs_raw,outputsTableRef){
  res <- list(outputs=outputs,
              inputs=inputs,
              names=names,
              outputs_raw=outputs_raw,
              inputs_raw=inputs_raw,
              outputsTableRef=outputsTableRef)
  attr(res,"class") <- "neuronBag"
  return(res)
}

#'@export
is.neuronBag <- function(x) inherits(x,"neuronBag")

#' Builds a neuronBag object either from a vector of query strings or a metadata data.frame.
#' @param typeQuery : either a vector of queries (similar to neuprint_search queries) or a
#' metadata data.frame for a set of neurons
#' @param fixed : if typeQuery is a query string, is it fixed?
#' @param by.roi : return results by ROI or just global weights?
#' @param verbose : Inform about progress if TRUE
#' @param selfRef : Should the input data.frame be used as the type reference (use if you already renamed
#  'neurons/types in that data frame)
#' @param ... : to be passed to getConnectionTable
#' @export
create_neuronBag <- function(typeQuery,fixed=FALSE,by.roi=TRUE,...){
  UseMethod("create_neuronBag")}

create_neuronBag.character <- function(typeQuery,fixed=FALSE,by.roi=TRUE,verbose=FALSE,...){
  TypeNames <- distinct(bind_rows(lapply(typeQuery,neuprint_search,field="type",fixed=fixed))) %>%
    mutate(databaseType = type)
  create_neuronBag(TypeNames,fixed=FALSE,by.roi=by.roi,verbose=verbose,...)
}

create_neuronBag.data.frame <- function(typeQuery,fixed=FALSE,selfRef=FALSE,by.roi=TRUE,verbose=FALSE,...){

  if (verbose) message("Calculate raw outputs")
  outputsR <- getConnectionTable(typeQuery,synapseType = "POST",by.roi=by.roi,verbose=verbose,...)

  if (verbose) message("Calculate raw inputs")
  inputsR <- getConnectionTable(typeQuery,synapseType = "PRE",by.roi=by.roi,verbose=verbose,...)

  if (verbose) message("Calculate type to type outputs")
  if (length(outputsR)==0){OUTByTypes <- NULL
  outputsTableRef <- NULL
  unknowns <- NULL
  }else{
    OUTByTypes <- getTypeToTypeTable(outputsR)
    outputsR <- retype.na(outputsR)
    outputsTableRef <- getTypesTable(unique(outputsR$type.to))
    unknowns <- retype.na_meta(neuprint_get_meta(unique(outputsR$to[!(outputsR$to %in% outputsTableRef$bodyid)])))
  }

  if (verbose) message("Calculate type to type inputs")
  if (nrow(inputsR)==0){INByTypes <- NULL}else{
    if (selfRef){
      INByTypes <- getTypeToTypeTable(inputsR,typesTable = typeQuery)
    }else{
      INByTypes <- getTypeToTypeTable(inputsR)
    }
    inputsR <- retype.na(inputsR)}

  return(neuronBag(outputs = OUTByTypes,
                   inputs = INByTypes,
                   names = typeQuery,
                   outputs_raw = outputsR,
                   inputs_raw = inputsR,
                   outputsTableRef = bind_rows(outputsTableRef,unknowns)
  ))
}
