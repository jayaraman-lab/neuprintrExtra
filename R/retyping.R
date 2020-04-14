#' Fill in the type and name field in case they are NAs, using the name field if it exists
#' (removing the _L/_R) or the neuron id. By default expects a table in to/from format.
#'
retype.na <- function(connectionTable){
  connectionTable <- connectionTable %>%
    mutate(name.from = ifelse(is.na(name.from),as.character(from),name.from),
           name.to = ifelse(is.na(name.to),as.character(to),name.to),
           type.from = ifelse(is.na(type.from),gsub("_L$|_R$","",name.from),type.from),
           type.to = ifelse(is.na(type.to),gsub("_L$|_R$","",name.to),type.to)
    )

  return(connectionTable)
}

retype.na_meta <- function(metaTable){
  metaTable <- metaTable %>%
    mutate(
      name = ifelse(is.na(name),as.character(bodyid),name),
      type = ifelse(is.na(type),gsub("_L$|_R$","",name),type)
    )

  return(metaTable)
}
