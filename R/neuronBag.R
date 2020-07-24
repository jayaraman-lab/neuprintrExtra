## neuronBags !

## Creator
new_neuronBag <- function(outputs,
                          inputs,
                          names,
                          outputs_raw,
                          inputs_raw,
                          outputsTableRef){
  stopifnot(is.data.frame(outputs))
  stopifnot(is.data.frame(inputs))
  stopifnot(is.data.frame(names))
  stopifnot(is.data.frame(outputs_raw))
  stopifnot(is.data.frame(inputs_raw))
  
  res <- list(outputs=outputs,
              inputs=inputs,
              names=names,
              outputs_raw=outputs_raw,
              inputs_raw=inputs_raw,
              outputsTableRef=outputsTableRef)
  attr(res,"class") <- "neuronBag"
  return(res)
}

## Validator
validate_neuronBag <- function(nBag){
  ## Check essential column names
  nameNames <- c("bodyid","upstream","downstream","name","type","databaseType")
  typeTableNames <- c( "type.from","type.to","databaseType.to","databaseType.from","outputContributionTotal","previous.type.to","previous.type.from","outputContribution","supertype.to1","supertype.from1","supertype.to2",
                       "supertype.from2","supertype.to3","supertype.from3","n_targets","n_type","absoluteWeight","weight","weightRelative","weightRelativeTotal","varWeight","sdWeight","tt")   
  connTableNames <- c("weight","from","to","name.from","type.from","type.to","databaseType.to","databaseType.from","outputContributionTotal","previous.type.to","previous.type.from","outputContribution","supertype.to1","supertype.from1","supertype.to2",
                      "supertype.from2","supertype.to3","supertype.from3","weightRelative","weightRelativeTotal")
  stopifnot(all(nameNames %in% names(nBag$names)))
  stopifnot(all(nameNames %in% names(nBag$outputsTableRef)))
  stopifnot(all(typeTableNames %in% names(nBag$inputs)))
  stopifnot(all(typeTableNames %in% names(nBag$outputs)))
  stopifnot(all(connTableNames %in% names(nBag$inputs_raw)))
  stopifnot(all(connTableNames %in% names(nBag$outputs_raw)))
  nBag
}

#' Test if x is a neuronBag
#' @param x An object to be tested
#' @return TRUE if x is a neuronBag
#'@export
is.neuronBag <- function(x) inherits(x,"neuronBag")

#### Helper -------------------------------------------------------

create_neuronBag <- function(typeQuery,fixed=FALSE,by.roi=TRUE,selfRef=FALSE,verbose=FALSE,omitInputs=FALSE,omitOutputs=FALSE,...){
  .Deprecated("neuronBag")
  UseMethod("neuronBag")
}

#' Create a neuronBag
#' 
#' Builds a neuronBag object either from a vector of query strings or a metadata \code{data.frame}.
#' @param typeQuery Either a vector of queries (similar to neuprint_search queries) or a
#' metadata \code{data.frame} for a set of neurons
#' @param fixed If typeQuery is a query string, is it fixed?
#' @param by.roi Return results by ROI or just global weights?
#' @param verbose Inform about progress if TRUE
#' @param selfRef Should the input data.frame be used as the type reference (use if you already renamed
#' neurons/types in that data frame)
#' @param omitInputs Skip calculation of inputs if TRUE
#' @param omitOutputs Skip calculation of outputs if TRUE
#' @param ... To be passed to getConnectionTable
#' @return An object of class \strong{neuronBag}. The object is a list with fields:
#' \describe{
#'  \item{names}{A table of metadata associated with the neurons in the bag -- with an extra column, 'databaseType'
#' which keeps the type name used in the database.}
#'  \item{outputs}{A type to type connectivity table of the outputs of the set of neurons.}
#'  \item{inputs}{A type to type connectivity table of the inputs of the set of neurons.}
#'  \item{outputs_raw}{A neuron to neuron connection table of the outputs of the set of neurons.}
#'  \item{inputs_raw}{A neuron to neuron connection table of the inputs of the set of neurons.}
#'  \item{outputsTableRef}{A table holding all instances of all the output types and their associated metadata}
#' }
#' @details A \strong{neuronBag} carries information about input and output connectivity of a group of neurons,
#'  at the neuron and type level, while keeping track of eventual type changes that occured. Particularly useful in combination 
#'  with retyping functions. Methods exist for filtering (\code{filter}), concatenating (\code{c}) and all retyping utilities.
#' @seealso \code{\link{lateralize_types}}, \code{\link{cxRetyping}}, \code{\link{redefine_types}} for retyping a bag.  
#' @export
neuronBag <- function(typeQuery,fixed=FALSE,by.roi=TRUE,selfRef=FALSE,verbose=FALSE,omitInputs=FALSE,omitOutputs=FALSE,...){
  UseMethod("neuronBag")}

#' @export
neuronBag.character <- function(typeQuery,fixed=FALSE,by.roi=TRUE,selfRef=FALSE,verbose=FALSE,omitInputs=FALSE,omitOutputs=FALSE,...){
  TypeNames <- distinct(bind_rows(lapply(typeQuery,neuprint_search,field="type",fixed=fixed))) %>%
    mutate(databaseType = type)
  neuronBag(TypeNames,fixed=FALSE,by.roi=by.roi,verbose=verbose,omitInputs=omitInputs,omitOutputs=omitOutputs,...)
}

#' @export
neuronBag.data.frame <- function(typeQuery,fixed=FALSE,selfRef=FALSE,by.roi=TRUE,verbose=FALSE,omitInputs=FALSE,omitOutputs=FALSE,...){
 
  if(!("databaseType" %in% names(typeQuery))){typeQuery <- mutate(typeQuery,databaseType=type)}
  if (!omitOutputs){
    if (verbose) message("Calculate raw outputs")
    outputsR <- getConnectionTable(typeQuery,synapseType = "POST",by.roi=by.roi,verbose=verbose,...)
  }else{
    outputsR <- getConnectionTable(character(),"POST",by.roi=by.roi,verbose=verbose,...)}
  
  if (!omitInputs){
    if (verbose) message("Calculate raw inputs")
    inputsR <- getConnectionTable(typeQuery,synapseType = "PRE",by.roi=by.roi,verbose=verbose,...)
  }else{
    inputsR <- getConnectionTable(character(),"PRE",by.roi=by.roi,verbose=verbose,...)}
  

  if (verbose) message("Calculate type to type outputs")
  
  OUTByTypes <- getTypeToTypeTable(outputsR)
  outputsR <- retype.na(outputsR)
  outputsTableRef <- getTypesTable(unique(outputsR$type.to))
  unknowns <- retype.na_meta(neuprint_get_meta(unique(outputsR$to[!(outputsR$to %in% outputsTableRef$bodyid)])) %>% mutate(databaseType=NA))
  
  if (verbose) message("Calculate type to type inputs")
  if (selfRef){
    INByTypes <- getTypeToTypeTable(inputsR,typesTable = typeQuery)
  }else{
    INByTypes <- getTypeToTypeTable(inputsR)
  }
  inputsR <- retype.na(inputsR)

  validate_neuronBag(new_neuronBag(outputs = OUTByTypes,
                                   inputs = INByTypes,
                                   names = typeQuery,
                                   outputs_raw = outputsR,
                                   inputs_raw = inputsR,
                                   outputsTableRef = rbind(outputsTableRef,unknowns)
  ))
}

#### Methods -------------------------------------------------

## Concatenate neuronBags
c.neuronBag <- function(...){
  full <- list(...)
  out <- new_neuronBag(outputs = distinct(do.call(rbind,lapply(full,function(i){i$outputs}))),
                   inputs = distinct(do.call(rbind,lapply(full,function(i){i$inputs}))),
                   names = distinct(do.call(rbind,lapply(full,function(i) i$names))),
                   outputs_raw = distinct(do.call(rbind,lapply(full,function(i){i$outputs_raw}))),
                   inputs_raw = distinct(do.call(rbind,lapply(full,function(i){i$inputs_raw}))),
                   outputsTableRef = distinct(do.call(rbind,lapply(full,function(i){i$outputsTableRef})))
  )
  validate_neuronBag(out)
}



#' Filter a neuronBag based on its \code{names} field
#' @param .nbag a neuronBag object to filter
#' @param filterPartners : Whether to apply the filter to input/output neurons to
#' @param ... to be passed to a filtering function applied to the \code{names} field
#'
filter.neuronBag <- function(.nbag,filterPartners=FALSE,...){

  .nbag$names <- filter(.nbag$names,...)

  .nbag$outputs <- filter(.nbag$outputs,type.from %in% .nbag$names$type)
  .nbag$outputs_raw <- filter(.nbag$outputs_raw,type.from %in% .nbag$names$type)

  .nbag$inputs <- filter(.nbag$inputs,type.to %in% .nbag$names$type)
  .nbag$inputs_raw <- filter(.nbag$inputs_raw,type.to %in% .nbag$names$type)

  if (filterPartners == TRUE){
    .nbag$outputs <- filter(.nbag$outputs,type.to %in% .nbag$names$type)
    .nbag$outputs_raw <- filter(.nbag$outputs_raw,type.to %in% .nbag$names$type)

    .nbag$inputs <- filter(.nbag$inputs,type.from %in% .nbag$names$type)
    .nbag$inputs_raw <- filter(.nbag$inputs_raw,type.from %in% .nbag$names$type)
  }

  .nbag$outputsTableRef <- filter(.nbag$outputsTableRef,type %in% .nbag$outputs$type.to)
  validate_neuronBag(.nbag)
}

#' Convenience to filter a neuronBag to just connections between the central neurons of a bag
#' @export
getIntraBag <- function(nBag){
  filter(nBag,filterPartners = TRUE,type %in% unique(nBag$names$type))
}

