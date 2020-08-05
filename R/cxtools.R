
#' Convenience function to do all the CX neurons renaming
#' @param connections Connectivity table  or neuronBag to modify
#' @param postfix One of "raw, "to" or "from". Specify if type (and name) columns in table to be modified are postfixed with
#' to and from or nothing
#' @param redefinePartners If table is a neuronBag, should the partners also be retyped?
#' @param verbose Should messages about progress be displayed
#' @details PFL1,3 and PFR_a and _b are renamed so that the outputs of a subtype go to the same region. All others are renamed L/R
#'
#' @export
cxRetyping <- function(connections,redefinePartners=TRUE,postfix=c("raw","to","from"),verbose=FALSE){
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

#' A color scale for Central complex supertypes
#' @return A list of colors (and breaks) mapping CX supertypes to palette36 colors
#'
#' @export
supertype2Palette <- function(){
  s2 <- c("vDelta","v\u0394","hDelta","h\u0394","Delta7","\u03947","EL","EPG","EPGt","ExR","FBt","FC","FR","FS","LNO","SPS-PB","LPsP","P","PEG","PEN","PFGs","PFL","PFN","PFR","ER","SA")
  pal <- paletteer::paletteer_d("Polychrome::palette36")[c(35,35,32,32,28,28,8,12,33,6,10,9,3,25,18,21,30,31,34,16,27,7,26,1,15,36)]
  names(pal) <- s2
  list(pal=pal,breaks=s2)
}

#' A color scale for Central complex supertypes
#' @return A ggplot scale to be used in plots were color (or fill) map to level 2 supertypes
#'
#' @export
scale_color_CX_supertype <- function(...){
  pal <- supertype2Palette()
  scale_color_manual(values=pal$pal,breaks=pal$breaks,...)
}

#' A fill scale for Central complex supertypes
#' @return A ggplot scale to be used in plots were color (or fill) map to level 2 supertypes
#'
#' @export
scale_fill_CX_supertype <- function(...){
  pal <- supertype2Palette()
  scale_fill_manual(values=pal$pal,breaks=pal$breaks,...)
}
