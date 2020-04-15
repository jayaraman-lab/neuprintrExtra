

#' Returns a table with all neurons belonging to types exactly specified
#' @param types: A vector of type names to be matched exactly
#' @return A data frame of instances of those types, including a databaseType column (used internally by other functions)
#' @export
getTypesTable <- function(types){
  typesTable <- bind_rows(lapply(types,function(t) neuprint_search(t,field="type",fixed=TRUE,exact=TRUE)))

  if (length(typesTable)>0){typesTable <- typesTable %>%
    mutate(databaseType = type)}
  return(typesTable)
}

#' Formats the results of a roiInfo query in a nice table with roi/upstream/downstreams columns
#' @param bodyids A collection of bodyids
#' @param ... To be passed to \code{neuprint_get_roiInfo}
#' @export
getRoiInfo <- function(bodyids,...){
  if (length(bodyids)==0){return(data.frame(bodyid=double(),roi=character(),pre=integer(),post=integer(),downstream=integer(),stringsAsFactors = FALSE))}
  roiInfo <- neuprint_get_roiInfo(bodyids,...)
  roiInfo <-  tidyr::pivot_longer(roiInfo,cols=-bodyid,names_to=c("roi","field"),names_sep="\\.",values_to="count")
  roiInfo <- tidyr::pivot_wider(roiInfo,names_from = "field",values_from="count")
  roiInfo
}