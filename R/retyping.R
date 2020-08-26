#' Fill in the type and name field in case they are NAs
#' @param connectionTable A table of connections in the to/from format
#' @return A connection table with name and type fields filled in.
#' @details If types are missing, fill them with a cleaned up (removing the L/R postfixes)
#' version of the name. If the name is missing, use the bodyid for everything. If the name is in parenthesis, append with the bodyid.
#' @export
retype.na <- function(connectionTable){
  ## All missing types become instance + bodyid if there's an instance
  connectionTable$type.from[is.na(connectionTable$type.from) & !is.na(connectionTable$name.from)] <- 
    paste0(connectionTable$name.from[is.na(connectionTable$type.from) & !is.na(connectionTable$name.from)],"_",
           connectionTable$from[is.na(connectionTable$type.from) & !is.na(connectionTable$name.from)])
  connectionTable$type.to[is.na(connectionTable$type.to) & !is.na(connectionTable$name.to)] <- 
    paste0(connectionTable$name.to[is.na(connectionTable$type.to) & !is.na(connectionTable$name.to)],"_",
           connectionTable$to[is.na(connectionTable$type.to) & !is.na(connectionTable$name.to)])
  
  ## Missing names get the bodyid
  connectionTable$name.from[is.na(connectionTable$name.from)] <- as.character(connectionTable$from[is.na(connectionTable$name.from)])
  connectionTable$name.to[is.na(connectionTable$name.to)] <- as.character(connectionTable$to[is.na(connectionTable$name.to)])
  
  ## Same for remaining types
  connectionTable$type.from[is.na(connectionTable$type.from)] <- connectionTable$name.from[is.na(connectionTable$type.from)]
  connectionTable$type.to[is.na(connectionTable$type.to)] <- connectionTable$name.to[is.na(connectionTable$type.to)]
  
  return(connectionTable)
}

retype.na_meta <- function(metaTable){
  metaTable$type[is.na(metaTable$type) & !is.na(metaTable$name)] <- 
    paste0(metaTable$name[is.na(metaTable$type) & !is.na(metaTable$name)],"_",
           metaTable$bodyid[is.na(metaTable$type) & !is.na(metaTable$name)])
  
  metaTable$name[is.na(metaTable$name)] <- as.character(metaTable$bodyid[is.na(metaTable$name)])
  metaTable$type[is.na(metaTable$type)] <- metaTable$name[is.na(metaTable$type)]

  return(metaTable)
}

#' Redefine types in a table/neuronBag
#' @param connections A connection table or a neuronBag object
#' @param retype_func A function taking \code{connections} and \code{postfix} as argument
#' generating new type names from the table
#' @param postfix One of "raw, "to" or "from". Specify if type (and name) columns in table to be modified are postfixed with
#' to and from or nothing
#' @param redefinePartners If table is a neuronBag, should the partners also be retyped?
#' @param ... Extra parameters to be passed to retype_func
#' @return A data frame or neuronBag with the columns \code{type_col} updated.
#' @details This is a low level functions. In most use cases, you can use \code{lateralize_types} or \code{redefineTypeByName} instead
#' @seealso \code{\link{lateralize_types}}, \code{\link{redefineTypeByName}}
#' @export
redefine_types <- function(connections,retype_func,postfix=c("raw","from","to"),redefinePartners=TRUE,...){UseMethod("redefine_types")}

#' @export
redefine_types.data.frame <- function(connections,retype_func,postfix=c("raw","from","to"),redefinePartners=TRUE,...){
  if (length(connections)==0){return(connections)}
  postfix <- match.arg(postfix)
  type_col <- get_col_name("type",postfix)
  newTypes <- retype_func(connections,postfix,...)
  connections[[type_col]] <-  newTypes
  return(connections)
}

#' @export
redefine_types.neuronBag <- function(connections,retype_func,postfix="raw",redefinePartners,...){
  connections$inputs_raw <- redefine_types(connections$inputs_raw,retype_func,postfix="to",...)
  connections$outputs_raw <- redefine_types(connections$outputs_raw,retype_func,postfix="from",...)
  connections$names <- redefine_types(connections$names,retype_func,postfix="raw",...)

  if ("allOutsFromIns" %in% names(connections[["ref"]])){
    connections$ref$allOutsFromIns <- redefine_types(connections$ref$allOutsFromIns,retype_func,postfix="to",...)
  }
  if ("allInsToOuts" %in% names(connections[["ref"]])){
    connections$ref$allInsToOuts <- redefine_types(connections$ref$allInsToOuts,retype_func,postfix="from",...)
  }
  
  if (redefinePartners){
    connections$inputs_raw <- redefine_types(connections$inputs_raw,retype_func,postfix="from",...)
    connections$outputs_raw <- redefine_types(connections$outputs_raw,retype_func,postfix="to",...)
    connections$outputsTableRef <- redefine_types(connections$outputsTableRef,retype_func,postfix="raw",...)
    
    if("allOutsFromIns" %in% names(connections[["ref"]])){
        connections$ref$allOutsFromIns <- redefine_types(connections$ref$allOutsFromIns,retype_func,postfix="from",...)
        connections$ref$inputs_outputsTableRef <- redefine_types(connections$ref$inputs_outputsTableRef,retype_func,postfix="raw",...)
      }
    if ("allInsToOuts" %in% names(connections[["ref"]])){
        connections$ref$allInsToOuts <- redefine_types(connections$ref$allInsToOuts,retype_func,postfix="to",...)
      }
  }
  if ("allOutsFromIns" %in% names(connections[["ref"]])){
    connections$ref$inputs_ref <- getTypeToTypeTable(connections$ref$allOutsFromIns,typesTable=connections$ref$inputs_outputsTableRef,oldTable=connections$ref$inputs_ref)
    connections$inputs <- processTypeToTypeFullInputs(connections$ref$inputs_ref,connections$inputs_raw)
  }else{
    connections$inputs <- getTypeToTypeTable(connections$inputs_raw,typesTable = connections$names, oldTable=connections$inputs)
  }
  
  if ("allInsToOuts" %in% names(connections[["ref"]])){
    connections$ref$outputs_ref <-getTypeToTypeTable(connections$ref$allInsToOuts,typesTable = connections$outputsTableRef,oldTable=connections$ref$outputs_ref)
    connections$outputs <- processTypeToTypeFullOutputs(connections$ref$outputs_ref,connections$outputs_raw)
  }else{
    connections$outputs <- getTypeToTypeTable(connections$outputs_raw,typesTable = connections$outputsTableRef, oldTable=connections$outputs)
  }
  return(connections)
}

#' Small utility to generate "type.from" kind of names
#' @export
get_col_name <- function(col=c("bodyid","name","type","databaseType","supertype",paste0("supertype",1:3)),post=c("raw","from","to")){
  col <- match.arg(col)
  post <- match.arg(post)
  if (post=="raw") return(col) else{
    if (col=="bodyid") return(post) else return(paste0(col,".",post))}
}

#' Retype neurons in a table according to L/R
#' @param connections Connectivity table  or neuronBag to modify
#' @param postfix One of "raw, "to" or "from". Specify if type (and name) columns in table to be modified are postfixed with
#' to and from or nothing
#' @param typeList : which types to lateralize (by default all the neurons which names
#' contains L or R)
#' @param redefinePartners If table is a neuronBag, should the partners also be retyped?
#' @examples
#' \dontrun{
#' PFLNames <- getTypesTable(c("PFL1","PFL2","PFL3"))
#'
#' ## Rename only PFL2
#' PFLNames2 <- lateralize_types(PFLNames,postfix="raw",typeList=c("PFL2"))
#'
#' ##Rename all PFLs
#' PFLNames3 <- lateralize_types(PFLNames,postFix="raw")
#' }
#' @export
lateralize_types <- function(connections,
                             postfix = c("raw","to","from"),
                             typeList=NULL,
                             redefinePartners=TRUE){
  postfix <- match.arg(postfix)
  redefine_types(connections,retype_func = lateralizer,postfix = postfix,redefinePartners=redefinePartners,typeList=typeList)
}

idemtyper <- function(connections,postfix){
  typeCol <- get_col_name(col="type",postfix)
  types <- connections[[typeCol]]
  return(types)
}

lateralizer <- function(connections,postfix,typeList){
    databaseCol <- get_col_name(col="databaseType",postfix)
    typeCol <- get_col_name(col="type",postfix)
    nameCol <- get_col_name(col="name",postfix)

    types <- connections[[typeCol]]

    if (is.null(typeList)) typeList <-  types

    typeList <- filter(connections,((!!as.name(typeCol)) %in% typeList) & (!!as.name(typeCol)) == (!!as.name(databaseCol)))

    typeList <- distinct(typeList,(!!as.name(nameCol)),(!!as.name(typeCol)))
    typeList <- filter(typeList,grepl("_R$|_L$|_R[1-9]$|_L[1-9]$|_R[1-9]/[1-9]$|_L[1-9]/[1-9]$|_L[1-9]_C[1-9]$|_R[1-9]_C[1-9]$|_L[1-9]_C[1-9]_irreg$|_R[1-9]_C[1-9]_irreg$|_L_C[1-9]_irreg$|_R_C[1-9]_irreg$|_L_small$|_R_small$",
                                      (!!as.name(nameCol)))) %>% na.omit()
    typeList <- typeList[[typeCol]]
    condition <- grepl("_L$|_L[1-9]$|_L[1-9]/[1-9]$|_L[1-9]_C[1-9]$|_L[1-9]_C[1-9]_irreg$|_L_C[1-9]_irreg$|_L_small$",connections[[nameCol]])

    for (t in typeList){
      rightType <- paste0(t,"_R")
      leftType <- paste0(t,"_L")
      types[types == t] <-  rightType
      types[types == rightType & condition] <-  leftType
    }
  return(types)
}

#' Retype neurons according to a grep pattern to be run on names
#' @inheritParams lateralize_types
#' @param pattern A grep pattern to be run to match a "name" column
#' @param sets List of neuron name sets. To be used in place of pattern.
#' @param nameModifiers A vector of strings of lenght 2 if \code{pattern} is used : first string is appended for matched types, second to the absence of matches
#' (for types in \code{typeList}), or of the same lenght as sets: postfixes
#' to be appended to the type name of all neurons whose name is in the corresponding set
#' @param perl Should the grep match use perl rules?
#' @details \code{pattern}, \code{typeList} and \code{perl} are used if \code{pattern} is used. Alternatively \code{sets} allow to retype according to an arbitrary number
#' of subtypes
#' @examples
#' \dontrun{
#' PFLNames <- getTypesTable(c("PFL1","PFL2","PFL3"))
#'
#' ## Rename only PFL3 according to their strange output pattern
#' PFLNames3 <- redefineTypeByName(PFLNames,typeList = c("PFL3"),pattern = "(^.*_L(?!.*irreg))|(^.*_R.*irreg)",perl=TRUE,nameModifiers = c("_L*","_R*"))
#' }
#' @export
redefineTypeByName <- function(connections,typeList=NULL,pattern=NULL,sets=NULL,nameModifiers,postfix=c("raw","to","from"),redefinePartners=FALSE,perl=FALSE){
  stopifnot(!is.null(pattern) | !is.null(sets))
  postfix <- match.arg(postfix)
  if (!is.null(pattern)){
    return(redefine_types(connections,retype_func = pattern_renamer,postfix=postfix,redefinePartners = redefinePartners,typeList=typeList,pattern=pattern,newPostFixes=nameModifiers,perl=perl))
  }
  if (!is.null(sets)){
    return(redefine_types(connections,retype_func = sets_retyper,postfix=postfix,sets=sets,nameModifiers=nameModifiers,redefinePartners = redefinePartners,kind="name"))
  }
}

pattern_renamer <- function(connections,postfix,typeList,pattern,newPostFixes,perl){
  name_col <- get_col_name(col="name",postfix)
  typeCol <- get_col_name(col="type",postfix)
  condition <- grepl(pattern,connections[[name_col]],perl=perl)
  types <- connections[[typeCol]]
  if (is.null(typeList)){
    typeList <- unique(types)
  }
  for (t in typeList){
    newNames <- paste0(t,newPostFixes)
    types[types == t] <- newNames[1]
    types[types == newNames[1] & !condition] <-  newNames[2]
  }
  return(types)
}

conditional_renamer <- function(connections,postfix,type,condition,newNames){
  typeCol <- get_col_name(col="type",postfix)
  types <- connections[[typeCol]]
  types[types == type] <- newNames[1]
  types[types == newNames[1] & !condition] <-  newNames[2]
  return(types)
}

sets_retyper <- function(conn,postfix,sets,nameModifiers,kind=c("bodyid","name")){
  kind <- match.arg(kind)
  typeCol <- get_col_name(col="type",postfix)
  nameCol <- get_col_name(col=kind,postfix)
  #nameCol <- ifelse(kind=="name",get_col_name(col="name",postfix),
  #                               ifelse(postfix=="raw",kind,postfix))
  types <- conn[[typeCol]]
  for (i in 1:length(sets)){
    types[conn[[nameCol]] %in% sets[[i]]] <- paste0(types[conn[[nameCol]] %in% sets[[i]]],nameModifiers[i])
  }
  return(types)
}


#' Retype according to sets of bodyids
#' @param connections Connectivity table  or neuronBag to modify
#' @param postfix One of "raw, "to" or "from". Specify if type (and name) columns in table to be modified are postfixed with
#' to and from or nothing
#' @param redefinePartners If table is a neuronBag, should the partners also be retyped?
#' @param sets List of neuron bodyid sets.
#' @param nameModifiers A vector of strings of lenght the same length as sets: postfixes
#' to be appended to the type name of all neurons whose name is in the corresponding set
#' @param nameModifiers Character vector of length the same length as sets containing
#'
#' @export
redefineTypeByBodyId <- function(connections,sets,nameModifiers,postfix=c("raw","to","from"),redefinePartners=FALSE){
  postfix <- match.arg(postfix)
  redefine_types(connections,retype_func = sets_retyper,postfix=postfix,sets=sets,nameModifiers=nameModifiers,redefinePartners = redefinePartners,kind="bodyid")
}

