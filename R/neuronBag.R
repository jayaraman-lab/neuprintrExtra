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
  typeTableNames <- c( "type.from","type.to","databaseType.to","databaseType.from","outputContributionTotal","outputContribution","supertype1.to","supertype1.from","supertype2.to",
                       "supertype2.from","supertype3.to","supertype3.from","n_targets","n_type","absoluteWeight","weight","weightRelative","weightRelativeTotal","sdWeight")   
  connTableNames <- c("weight","from","to","name.from","type.from","type.to","databaseType.to","databaseType.from","outputContributionTotal","outputContribution","supertype1.to","supertype1.from","supertype2.to",
                      "supertype2.from","supertype3.to","supertype3.from","weightRelative","weightRelativeTotal")
  stopifnot(all(nameNames %in% names(nBag$names)))
  stopifnot(all(nameNames %in% names(nBag$outputsTableRef)))
  stopifnot(all(typeTableNames %in% names(nBag$inputs)))
  stopifnot(all(typeTableNames %in% names(nBag$outputs)))
  stopifnot(all(connTableNames %in% names(nBag$inputs_raw)))
  stopifnot(all(connTableNames %in% names(nBag$outputs_raw)))
  if("ref" %in% names(nBag)){
    refFields <- c("inputs_ref","outputs_ref","allOutsFromIns","allInsToOuts","inputs_outputsTableRef")
    stopifnot(all(refFields %in% names(nBag$ref)))
  }
  nBag
}

#' Test if x is a neuronBag
#' @param x An object to be tested
#' @return TRUE if x is a neuronBag
#'@export
is.neuronBag <- function(x) inherits(x,"neuronBag")

#### Helper -------------------------------------------------------

create_neuronBag <- function(typeQuery,fixed=FALSE,by.roi=TRUE,selfRef=FALSE,verbose=FALSE,omitInputs=FALSE,omitOutputs=FALSE,computeKnownRatio=FALSE,renaming=NULL,...){
  .Deprecated("neuronBag")
  UseMethod("neuronBag")
}

#' Create a neuronBag
#' 
#' Builds a neuronBag object either from a vector of query strings or a metadata \code{data.frame}.
#' @param typeQuery Either a vector of queries (similar to a \code{getTypesTable} queries) or a
#' metadata \code{data.frame} for a set of neurons. If a \code{data.frame} it needs to contain all the
#' instances of the types. 
#' @param fixed If typeQuery is a query string, is it fixed?
#' @param by.roi Return results by ROI or just global weights?
#' @param verbose Inform about progress if TRUE
#' @param selfRef Deprecated argument
#' @param omitInputs Skip calculation of inputs if TRUE
#' @param omitOutputs Skip calculation of outputs if TRUE
#' @param computeKnownRatio Computes relative weights and output contributions that sum to 1 (for types and for neurons).
#' This requires keeping track of all the outputs of the inputs, and all the inputs of the outputs, in a field called \code{ref}
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
neuronBag <- function(typeQuery,fixed=FALSE,by.roi=TRUE,selfRef=FALSE,verbose=FALSE,omitInputs=FALSE,omitOutputs=FALSE,computeKnownRatio=FALSE,renaming=NULL,...){
  UseMethod("neuronBag")}

#' @export
neuronBag.character <- function(typeQuery,fixed=FALSE,by.roi=TRUE,selfRef=FALSE,verbose=FALSE,omitInputs=FALSE,omitOutputs=FALSE,computeKnownRatio=FALSE,renaming=NULL,...){
  TypeNames <- getTypesTable(typeQuery)
  neuronBag(TypeNames,fixed=FALSE,by.roi=by.roi,verbose=verbose,omitInputs=omitInputs,omitOutputs=omitOutputs,computeKnownRatio=computeKnownRatio,renaming=renaming,...)
}

#' @export
neuronBag.data.frame <- function(typeQuery,fixed=FALSE,selfRef=FALSE,by.roi=TRUE,verbose=FALSE,omitInputs=FALSE,omitOutputs=FALSE,computeKnownRatio=FALSE,renaming=NULL,...){
 
  if(is.null(renaming)){renaming <- function(x,postfix="raw",...){
    redefine_types(x,idemtyper,postfix=postfix,...)}
  }
  #typeQuery <- renaming(typeQuery,postfix="raw")
  
  if(!("databaseType" %in% names(typeQuery))){
    warning("No 'databaseType' field. Assuming the 'type' column contains database types.")
    typeQuery <- mutate(typeQuery,databaseType=as.character(type))}
  
  ## NEED TO ADD WARNING ABOUT THAT
  typeQuery <- mutate(typeQuery,ifelse(type!=databaseType |!is.na(databaseType),databaseType,type))
  
  if (!omitOutputs){
    if (verbose) message("Calculate raw outputs")
    outputsR <- getConnectionTable(typeQuery,synapseType = "POST",by.roi=by.roi,verbose=verbose,...)
    outputsTableRef <- getTypesTable(unique(outputsR$type.to))
    outputsR <- retype.na(outputsR)
    unknowns <- getMeta(unique(outputsR$to[!(outputsR$to %in% outputsTableRef$bodyid)]))
    outputsTableRef <- rbind(outputsTableRef,retype.na_meta(unknowns))
    outputsR <- renaming(outputsR,postfix="to")
    outputsR <- renaming(outputsR,postfix="from")
    if (computeKnownRatio){
      if (verbose) message("Calculate full raw inputs to outputs")
      
      allInsToOuts <- getConnectionTable(outputsTableRef,synapseType="PRE",by.roi=by.roi,verbose=verbose,...) %>%
        retype.na() %>%
        renaming(postfix="to") %>%
        renaming(postfix="from") %>%
        group_by(to) %>% mutate(knownTotalWeight=sum(weight[match(from,from)])) %>%
        group_by(to,roi) %>% mutate(knownTotalROIweight=sum(ROIweight)) %>% ungroup() %>%
        mutate(knownWeightRelativeTotal = weight/knownTotalWeight,
               knownWeightRelative = ROIweight/knownTotalROIweight,
               input_completedness = knownTotalROIweight/totalROIweight
               )
      
      outputsR <- group_by(outputsR,from)  %>% mutate(knownTotalPreWeight=sum(weight[match(to,to)])) %>%
        group_by(from,roi) %>% mutate(knownTotalPreROIweight=sum(ROIweight),
                                      knownOutputContribution = ROIweight/knownTotalPreROIweight,
                                      knownOutputContributionTotal = weight/knownTotalPreWeight) %>% ungroup()
      
      outputsR <- mutate(outputsR,
                         knownTotalWeight = allInsToOuts$knownTotalWeight[match(outputsR$to,allInsToOuts$to)],
                         knownWeightRelativeTotal = weight/knownTotalWeight,
                         knownTotalROIweight = allInsToOuts$knownTotalROIweight[match(paste0(outputsR$to,outputsR$roi),paste0(allInsToOuts$to,allInsToOuts$roi))],
                         knownWeightRelative = ROIweight/knownTotalROIweight,
                         output_completedness = knownTotalPreROIweight/totalPreROIweight,
                         input_completedness = knownTotalROIweight/totalROIweight)
      
      allInsToOuts <- right_join(outputsR,allInsToOuts)
    }
    outputsTableRef <- renaming(outputsTableRef)
  }else{
    outputsR <- getConnectionTable(character(),"POST",by.roi=by.roi,verbose=verbose,...)
    outputsTableRef <- getTypesTable(character())
    if(computeKnownRatio) allInsToOuts <- getConnectionTable(character(),"POST",by.roi=by.roi,verbose=verbose,...)
  }
  
  if (!omitInputs){
    if (verbose) message("Calculate raw inputs")
    inputsR <- getConnectionTable(typeQuery,synapseType = "PRE",by.roi=by.roi,verbose=verbose,...)
    inputsTableRef <- getTypesTable(unique(inputsR$type.from))
    inputsR <- retype.na(inputsR)
    unknownsIn <- getMeta(unique(inputsR$from[!(inputsR$from %in% inputsTableRef$bodyid)]))
    inputsFullQuery <- rbind(inputsTableRef,retype.na_meta(unknownsIn))
    inputsR <- renaming(inputsR,postfix="from")
    inputsR <- renaming(inputsR,postfix="to")
   
    if (computeKnownRatio){
      if (verbose) message("Calculate full raw outputs of inputs")
      
      allOutsFromIns <- getConnectionTable(inputsFullQuery,synapseType="POST",by.roi=by.roi,verbose=verbose,...) %>% 
        retype.na() %>%
        renaming(postfix = "from") %>%
        renaming(postfix = "to") %>%
        group_by(from) %>% mutate(knownTotalPreWeight=sum(weight[match(to,to)])) %>%
        group_by(from,roi) %>% mutate(knownTotalPreROIweight=sum(ROIweight)) %>% ungroup() %>%
        mutate(knownOutputContribution = ROIweight/knownTotalPreROIweight,
               knownOutputContributionTotal = weight/knownTotalPreWeight,
               output_completedness = knownTotalPreROIweight/totalPreROIweight
               ) 
        
      inputsR <- group_by(inputsR,to) %>% mutate(knownTotalWeight=sum(weight[match(from,from)])) %>%
        group_by(to,roi) %>% mutate(knownTotalROIweight=sum(ROIweight)) %>% ungroup()
      
      inputsR <- mutate(inputsR,
                        knownTotalPreWeight = allOutsFromIns$knownTotalPreWeight[match(inputsR$from,allOutsFromIns$from)],
                        knownTotalPreROIweight = allOutsFromIns$knownTotalPreROIweight[match(paste0(inputsR$from,inputsR$roi),paste0(allOutsFromIns$from,allOutsFromIns$roi))],
                        knownOutputContribution = ROIweight/knownTotalPreROIweight,
                        knownOutputContributionTotal = weight/knownTotalPreWeight,
                        knownWeightRelativeTotal = weight/knownTotalWeight,
                        knownWeightRelative = ROIweight/knownTotalROIweight,
                        output_completedness = knownTotalPreROIweight/totalPreROIweight,
                        input_completedness = knownTotalROIweight/totalROIweight)
      allOutsFromIns <- right_join(inputsR,allOutsFromIns)
      
      inputsOutTableRef <- getTypesTable(unique(allOutsFromIns$databaseType.to))
      unknownsInOut <- getMeta(unique(allOutsFromIns$to[!(allOutsFromIns$to %in% inputsOutTableRef$bodyid)]))
      inputsOutTableRef <- renaming(rbind(inputsOutTableRef,retype.na_meta(unknownsInOut)))
    }
  }else{
    inputsR <- getConnectionTable(character(),"PRE",by.roi=by.roi,verbose=verbose,...)
    if(computeKnownRatio){allOutsFromIns <- getConnectionTable(character(),"PRE",by.roi=by.roi,verbose=verbose,...)
                          inputsOutTableRef <- getTypesTable(character())
        }
    }
  typeQuery <- renaming(retype.na_meta(typeQuery))   ## Revert typeQuery to the "original"

  if (verbose) message("Calculate type to type outputs")
  if (computeKnownRatio & nrow(outputsR)>0){
    OUTByTypesRef <- getTypeToTypeTable(allInsToOuts,typesTable = outputsTableRef) 
    OUTByTypes <- processTypeToTypeFullOutputs(OUTByTypesRef,outputsR)
    
  }else{
    OUTByTypes <- getTypeToTypeTable(outputsR,typesTable = outputsTableRef)
    if(computeKnownRatio) OUTByTypesRef <- getTypeToTypeTable(allInsToOuts)
  }
  
  
  if (verbose) message("Calculate type to type inputs")
  if (computeKnownRatio & nrow(inputsR)>0){
      INByTypesRef <- getTypeToTypeTable(allOutsFromIns,typesTable = inputsOutTableRef)
      INByTypes <- processTypeToTypeFullInputs(INByTypesRef,inputsR)
    }else{
      INByTypes <- getTypeToTypeTable(inputsR,typesTable = typeQuery)
      if(computeKnownRatio) INByTypesRef <-getTypeToTypeTable(allOutsFromIns)
    }

  nBag <- new_neuronBag(outputs = OUTByTypes,
                                   inputs = INByTypes,
                                   names = typeQuery,
                                   outputs_raw = outputsR,
                                   inputs_raw = inputsR,
                                   outputsTableRef = outputsTableRef
  )
  if (computeKnownRatio){
                        nBag[["ref"]][["allOutsFromIns"]] <- allOutsFromIns
                        nBag[["ref"]][["inputs_ref"]] <- INByTypesRef
                        nBag[["ref"]][["inputs_outputsTableRef"]] <- inputsOutTableRef 
                        nBag[["ref"]][["allInsToOuts"]] <- allInsToOuts
                        nBag[["ref"]][["outputs_ref"]] <- OUTByTypesRef
  }
  
  
  validate_neuronBag(nBag)
  
}

processTypeToTypeFullInputs <- function(INByTypes,
                                    inputsR){
  
  INByTypes <- INByTypes %>% 
    group_by(type.from,roi) %>% 
    mutate(knownOutputContribution_perType=knownOutputContribution/sum(knownOutputContribution),
           knownWROutputContribution_perType=weightRelative/sum(weightRelative)) %>% 
    ungroup() %>% 
    filter(type.to %in% inputsR$type.to) %>%
    group_by(type.to,roi) %>% 
    mutate(knownWeightRelative_perType=knownWeightRelative/sum(knownWeightRelative))
  INByTypes
}

processTypeToTypeFullOutputs <- function(OUTByTypes,outputsR){
  OUTByTypes <- OUTByTypes %>% 
    group_by(type.to,roi) %>% 
    mutate(knownWeightRelative_perType=knownWeightRelative/sum(knownWeightRelative)) %>% 
    ungroup() %>% 
    filter(type.from %in% outputsR$type.from) %>%
    group_by(type.from,roi) %>% 
    mutate(knownOutputContribution_perType=knownOutputContribution/sum(knownOutputContribution),
           knownWROutputContribution_perType=weightRelative/sum(weightRelative))
  OUTByTypes
}
#### Methods -------------------------------------------------

## Concatenate neuronBags
#' @export
c.neuronBag <- function(...){
  full <- rlang::list2(...)
  out <- new_neuronBag(outputs = distinct(do.call(rbind,lapply(full,function(i){i$outputs}))),
                   inputs = distinct(do.call(rbind,lapply(full,function(i){i$inputs}))),
                   names = distinct(do.call(rbind,lapply(full,function(i) i$names))),
                   outputs_raw = distinct(do.call(rbind,lapply(full,function(i){i$outputs_raw}))),
                   inputs_raw = distinct(do.call(rbind,lapply(full,function(i){i$inputs_raw}))),
                   outputsTableRef = distinct(do.call(rbind,lapply(full,function(i){i$outputsTableRef})))
  )
  
  if (all(sapply(full,function(x) "ref" %in% names(x)))){
    out$ref$outputs_ref <- distinct(do.call(rbind,lapply(full,function(i){i$ref$outputs_ref})))
    out$ref$inputs_ref <- distinct(do.call(rbind,lapply(full,function(i){i$ref$inputs_ref})))
    out$ref$allInsToOuts  <- distinct(do.call(rbind,lapply(full,function(i){i$ref$allInsToOuts})))
    out$ref$allOutsFromIns  <- distinct(do.call(rbind,lapply(full,function(i){i$ref$allOutsFromIns})))
    out$ref$inputs_outputsTableRef  <- distinct(do.call(rbind,lapply(full,function(i){i$ref$inputs_outputsTableRef})))
  }
  validate_neuronBag(out)
}



#' Filter a neuronBag based on its \code{names} field
#' @param .nbag a neuronBag object to filter
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
  validate_neuronBag(.nbag)
}

#' Convenience to filter a neuronBag to just connections between the central neurons of a bag
#' @export
getIntraBag <- function(nBag){
  filter(nBag,filterPartners = TRUE,type %in% unique(nBag$names$type))
}

