#' Fill in the type and name field in case they are NAs
#' @param connectionTable A table of connections in the to/from format
#' @return A connection table with name and type fields filled in.
#' @details If types are missing, fill them with a cleaned up (removing the L/R postfixes)
#' version of the name. If the name is missing, use the bodyid for everything.
#' @export
retype.na <- function(connectionTable){
  connectionTable <- connectionTable %>%
    mutate(name.from = ifelse(is.na(name.from),as.character(from),name.from),
           name.to = ifelse(is.na(name.to),as.character(to),name.to),
           type.from = ifelse(is.na(type.from),gsub("_L$|_R$","",name.from),type.from),
           type.to = ifelse(is.na(type.to),gsub("_L$|_R$","",name.to),type.to)
    )

  return(connectionTable)
}

#' Redefine types in a table/neuronBag
#' @param connections A connection table or a neuronBag object
#' @param retype_func A function taking \code{connections} and \code{postfix} as argument
#' generating new type names from the table
#' @param postfix One of "raw, "to" or "from". Specify if type (and name) columns in table to be modified are postfixed with
#' to and from or nothing
#' @param redefinePartners If table is a neuronBag, should the partners also be retyped?
#' @param ... Extra parameters to be passed to retype_func
#' @return A data frame or neuronBag with the columns \code{type_col} updated, and previous
#' types stored in \code{previous.} columns
#'
redefine_types <- function(connections,retype_func,postfix=c("raw","from","to"),redefinePartners=TRUE,...){UseMethod("redefine_types")}

redefine_types.data.frame <- function(table,retype_func,postfix=c("raw","from","to"),redefinePartners=TRUE,...){
  postfix <- match.arg(postfix)
  type_col <- get_col_name("type",postfix)
  table[[paste0("previous.",type_col)]] <- table[[type_col]] ## keeping track of the last named types
  newTypes <- retype_func(table,postfix,...)
  table[[type_col]] <-  newTypes
  return(table)
}

redefine_types.neuronBag <- function(neuronBag,retype_func,postfix="raw",redefinePartners,...){
  neuronBag$inputs_raw <- redefines_types(neuronBag$inputs_raw,retype_func,postfix="to",...)
  neuronBag$outputs_raw <- redefines_types(neuronBag$outputs_raw,retype_func,postfix="from",...)
  neuronBag$names <- redefines_types(neuronBag$names,retype_func,postfix="raw",...)

  if (redefinePartners){
    neuronBag$inputs_raw <- redefines_types(neuronBag$inputs_raw,retype_func,postfix="from",...)
    neuronBag$outputs_raw <- redefines_types(neuronBag$outputs_raw,retype_func,postfix="to",...)
    neuronBag$outputsTableRef <- redefines_types(neuronBag$outputsTableRef,retype_func,postfix="raw",...)
  }

  neuronBag$inputs <- getTypeToTypeTable(neuronBag$inputs_raw,typesTable = neuronBag$names, oldTable=neuronBag$inputs)
  neuronBag$outputs <- getTypeToTypeTable(neuronBag$outputs_raw,typesTable = neuronBag$outputsTableRef, oldTable=neuronBag$outputs)
  return(neuronBag)
}

#' Small utility to generate "type.from" kind of names
get_col_name <- function(col="type",post=c("raw","from","to")){
  post <- match.arg(post)
  return(ifelse(post=="raw",col,paste0(col,".",post)))
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


lateralizer <- function(connections,postfix,typeList){
    databaseCol <- get_col_name(col="databaseType",postfix)
    typeCol <- get_col_name(col="type",postfix)
    nameCol <- get_col_name(col="name",postfix)

    types <- connections[[typeCol]]

    if (is.null(typeList)) typeList <-  types

    typeList <- filter(connections,((!!as.name(typeCol)) %in% typeList) & (!!as.name(typeCol)) == (!!as.name(databaseCol)))

    typeList <- distinct(typeList,(!!as.name(nameCol)),(!!as.name(typeCol)))
    typeList <- filter(typeList,grepl("_R|_L",(!!as.name(nameCol)))) %>% na.omit()
    typeList <- typeList[[typeCol]]
    condition <- grepl("_L",connections[[nameCol]])

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
#' @param nameModifiers A vector of strings of lenght 2 : first string is appended for matched types, second to the absence of matches
#' (for types in \code{typeList})
#' @param perl Should the grep match use perl rules?
#' @examples
#' \dontrun{
#' PFLNames <- getTypesTable(c("PFL1","PFL2","PFL3"))
#'
#' ## Rename only PFL3 according to their strange output pattern
#' PFLNames3 <- redefineTypeByName(PFLNames,typeList = c("PFL3"),pattern = "(^.*_L(?!.*irreg))|(^.*_R.*irreg)",perl=TRUE,nameModifiers = c("_L*","_R*"))
#' }
#' @export
redefineTypeByName <- function(connections,typeList,pattern,nameModifiers,postfix=c("raw","to","from"),redefinePartners=FALSE,perl=FALSE){
  postfix <- match.arg(postfix)
  redefine_types(connections,retype_func = pattern_renamer,postfix=postfix,redefinePartners = redefinePartners,typeList=typeList,pattern=pattern,newPostFixes=nameModifiers,perl=perl)
}

pattern_renamer <- function(connections,postfix,typeList,pattern,newPostFixes,perl){
  name_col <- get_col_name(col="name",postfix)
  condition <- grepl(pattern,connections[[name_col]],perl=perl)
  for (t in typeList){
    newNames <- paste0(t,newPostFixes)
    types <- conditional_renamer(connections,postfix,t,condition,newNames)
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
