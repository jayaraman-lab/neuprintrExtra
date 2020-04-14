

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
