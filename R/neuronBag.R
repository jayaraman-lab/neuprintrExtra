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

#### Create neuronBags -------------------------------------------------------

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

#' @export
create_neuronBag.character <- function(typeQuery,fixed=FALSE,by.roi=TRUE,verbose=FALSE,...){
  TypeNames <- distinct(bind_rows(lapply(typeQuery,neuprint_search,field="type",fixed=fixed))) %>%
    mutate(databaseType = type)
  create_neuronBag(TypeNames,fixed=FALSE,by.roi=by.roi,verbose=verbose,...)
}

#' @export
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

#### Methods -------------------------------------------------

## Concatenate neuronBags
#' @export
c.neuronBag <- function(...){
  full <- list(...)
  out <- neuronBag(outputs = distinct(bind_rows(lapply(full,function(i) i$outputs))),
                   inputs = distinct(bind_rows(lapply(full,function(i) i$inputs))),
                   names = distinct(bind_rows(lapply(full,function(i) i$names))),
                   outputs_raw = distinct(bind_rows(lapply(full,function(i) i$outputs_raw))),
                   inputs_raw = distinct(bind_rows(lapply(full,function(i) i$inputs_raw))),
                   outputsTableRef = distinct(bind_rows(lapply(full,function(i) i$outputsTableRef)))
  )
  return(out)
}



#' Filter a neuronBag based on its \code{names} field
#' @param nbag a neuronBag object to filter
#' @param filterPartners : Whether to apply the filter to input/output neurons to
#' @param ... to be passed to a filtering function applied to the \code{names} field
#'
#' @export
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
  .nbag
}

#' Convenience to filter a neuronBag to just connections between the central neurons of a bag
#' @export
getIntraBag <- function(nBag){
  filter(nBag,filterPartners = TRUE,type %in% unique(nBag$names$type))
}

#' Build a per roi summary of innervation for neurons in a neuronBag
#' @param neurons : a dataframe as returned by a search or a neuronBag object
#' @param threshold : the minimal average number of synapses for a type in a ROI
#' for it to be included
#' @param rois : a roiset to consider (if NULL consider all rois)
#'
#'@export
getROISummary <- function(neurons,threshold=0,rois=NULL){UseMethod("getROISummary")}

#

#'@export
getROISummary.neuronBag <- function(neurons,threshold=0,rois = NULL){
  getROISummary(neurons$names,rois=rois,threshold=threshold)
}

#'@export
getROISummary.data.frame <- function(neurons,threshold=0,rois = NULL){

  roiSummary <- getRoiInfo(neurons)

  if (!(is.null(rois))){
    roiSummary <- roiSummary %>%
      filter(roi %in% rois$roi)
  }

  countInstances <- group_by(neurons,type) %>% summarize(n=n())

  roiSummary <-
    left_join(roiSummary,
              select(neurons,bodyid,type,databaseType),
              by=c("bodyid")) %>% tidyr::replace_na(list(downstream=0,upstream=0)) %>%
    group_by(roi,type,databaseType) %>%
    summarize(downstream=mean(downstream),
              upstream=mean(upstream)) %>%
    ungroup() %>%
    mutate(fullWeight = downstream+upstream,
           deltaWeight = (downstream - upstream)/fullWeight) %>%
    filter(fullWeight>threshold)

  roiSummary <- left_join(roiSummary,countInstances,by="type")

  if (!is.null(rois)){roiSummary <- left_join(roiSummary,rois,by=roi)}
  return(roiSummary)
}
