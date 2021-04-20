#'A cached version of \code{neuprint_get_meta}
#'@export
getMeta <- function(bodyids,...){
  knownMeta <- get("storedMeta",envir=cacheEnv)
  bodyids <- neuprint_ids(bodyids,unique=FALSE,mustWork = FALSE)
  newIDs <- bodyids[!(bodyids %in% knownMeta$bodyid)]
  reusable <- filter(knownMeta,bodyid %in% bodyids)
  
  metad <- neuprint_get_meta(newIDs,...) %>% mutate(databaseType = type)
  if(length(newIDs)>0){assign("storedMeta",rbind(knownMeta,metad),envir=cacheEnv)}
  full <- rbind(metad,reusable)
  full[match(bodyids,full$bodyid),]
}


#' Returns a table with all neurons belonging to types exactly specified
#' @param types A vector of type names to be matched exactly
#' @return A data frame of instances of those types, including a databaseType column (used internally by other functions)
#' @export
getTypesTable <- function(types){
  if(length(types)==0) return(neuprint_get_meta(1) %>% mutate(databaseType=character()))
  
  knownTypes <- get("storedTypes",envir=cacheEnv)
  
  newTypes <- types[!(types %in% knownTypes$type)]
  reusable <- filter(knownTypes,type %in% types)
  
  
  typeQ <- gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\\\\\\\\\\\\\1",newTypes)
  typeQ <- paste0(typeQ,collapse="|")
  typeQ <- paste0("/type:",typeQ)
  typesTable <- neuprint_search(typeQ)
    
  if(is.null(typesTable)){typesTable <- neuprint_get_meta(1)}
  typesTable <- typesTable %>%
      mutate(databaseType = as.character(type))
  
  if(length(newTypes)>0){assign("storedTypes",rbind(knownTypes,typesTable),envir=cacheEnv)
                         assign("storedMeta",distinct(rbind(get("storedMeta",envir=cacheEnv),typesTable)))}
  
  return(rbind(typesTable,reusable))
}

#' Formats the results of a roiInfo query in a nice table with roi/upstream/downstreams columns
#' @param bodyids A collection of bodyids
#' @param ... To be passed to \code{neuprint_get_roiInfo}
#' @export
getRoiInfo <- function(bodyids,...){
  if (length(bodyids)==0){return(data.frame(bodyid=double(),roi=character(),pre=integer(),post=integer(),downstream=integer(),upstream=integer(),stringsAsFactors = FALSE))}
  knownInfo <- get("storedRoiInfo",envir=cacheEnv)
  bodyids <- neuprint_ids(bodyids,unique=FALSE,mustWork = FALSE)
  newIDs <- bodyids[!(bodyids %in% knownInfo$bodyid)]
  reusable <- filter(knownInfo,bodyid %in% bodyids)
  roiInfo <- neuprint_get_roiInfo(newIDs,...)
  
  if (nrow(roiInfo)==0 | ncol(roiInfo) == 1){roiInfo <- data.frame(bodyid=double(),roi=character(),pre=integer(),post=integer(),downstream=integer(),upstream=integer(),stringsAsFactors = FALSE)}else{
    roiInfo <-  tidyr::pivot_longer(roiInfo,cols=-bodyid,names_to=c("roi","field"),names_sep="\\.",values_to="count")
    roiInfo <- tidyr::pivot_wider(roiInfo,names_from = "field",values_from="count")
    ## Ensure all columns are here (even if values are missing)
    roiInfo <- data.frame(bodyid=roiInfo$bodyid,
                          roi=roiInfo$roi) %>% mutate(
                                       pre=if ("pre" %in% names(roiInfo)) roiInfo$pre else NA,
                                       post=if ("post" %in% names(roiInfo)) roiInfo$post else NA,
                                       downstream=if ("downstream" %in% names(roiInfo)) roiInfo$downstream else NA,
                                       upstream=if ("upstream" %in% names(roiInfo)) roiInfo$upstream else NA)
    assign("storedRoiInfo",rbind(knownInfo,roiInfo),envir=cacheEnv)
  }
  rbind(roiInfo,reusable)
}

empty_connTable <- function(has.roi = FALSE){
  res <- data.frame(bodyid=numeric(),partner=numeric(),prepost=integer(),weight=integer())
  if (has.roi){res <- mutate(res,roi=character(),ROIweight=integer())}
  res
}