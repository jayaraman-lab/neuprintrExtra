
#' Convenience function to do all the CX neurons renaming
#' @param connections Connectivity table  or neuronBag to modify
#' @param postfix One of "raw, "to" or "from". Specify if type (and name) columns in table to be modified are postfixed with
#' to and from or nothing
#' @param redefinePartners If table is a neuronBag, should the partners also be retyped?
#' @param verbose Should messages about progress be displayed
#' @details PFL1,3 and PFR_a and _b are renamed so that the outputs of a subtype go to the same region. All others are renamed L/R
#'
#' @export
cxRetyping <- function(connections,redefinePartners=TRUE,postfix=c("raw","to","from"),verbose=TRUE){
  postfix <- match.arg(postfix)
  if (verbose) message("Renaming PFL3")
  connections <- redefineTypeByName(connections,typeList = c("PFL3"),pattern = "(^.*_L(?!.*irreg))|(^.*_R.*irreg)",perl=TRUE,nameModifiers = c("_L*","_R*"),redefinePartners = redefinePartners,postfix = postfix)
  if (verbose) message("Renaming PFL1/PFR_a")
  connections <- redefineTypeByName(connections,typeList = c("PFR_a","PFL1"),pattern = "_L[2-7]|_R1",nameModifiers = c("_L*","_R*"),redefinePartners = redefinePartners,postfix = postfix)
  if (verbose) message("Renaming PFR_b")
  connections <- redefineTypeByName(connections,typeList = c("PFR_b"),pattern = "(^.*_L(?!.*C9))|(^.*C1.*)",perl=TRUE,nameModifiers = c("_L*","_R*"),redefinePartners = redefinePartners,postfix = postfix)
  if (verbose) message("All other L/R retyping")
  connections <- lateralize_types(connections,redefinePartners = redefinePartners,postfix = postfix)
  return(connections)
}
