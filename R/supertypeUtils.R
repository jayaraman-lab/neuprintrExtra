#' Gets a supertype from a type name, with various levels of depth
#' @param types A type name or vector of type names, or a data frame containing databaseType
#'  (or databaseType.to and from) columns, or a neuronBag object
#' @param level Depth of the supertype. Possible values are 1,2 or 3 (default 2), 1 being the finest
#' subdivision and 3 the coarsest
#' @param unicodeDelta Whether or not to use unicode greek characters for Delta (Delta7 and v/h Delta) supertypes. 
#' Unicode by default
#' @details For example, at level 1 vDelta neurons are just divided in vDeltaA to K, at level 2 they are
#' vDelta, and at level 3 they are FB Interneurons
#' @return an object of the same type as the \code{types} inputed, possibly with extra columns
#'
#' @export
supertype <- function(types,level=2,unicodeDelta=TRUE){UseMethod("supertype")}

#' @export
supertype.character <- function(types,level=2,unicodeDelta=TRUE){
  supertype <- types
  supertype[is.na(types)] <- "Unassigned"
  
  supertype <- stringr::str_extract(types,"[A-Z]+")
  supertype[is.na(supertype)] <- "Unassigned"
  supertype[grepl("^FB[1-9].*",types)] <- stringr::str_extract(types,"FB[1-9]")[grepl("^FB[1-9].*",types)]
  supertype[grepl("^vDelta.*",types)] <- stringr::str_extract(types,"vDelta[A-O]")[grepl("^vDelta.*",types)]
  supertype[grepl("^hDelta.*",types)] <- stringr::str_extract(types,"hDelta[A-M]")[grepl("^hDelta.*",types)]
  supertype[grepl("^FC[1-9].*",types)] <- stringr::str_extract(types,"FC[1-9]")[grepl("^FC[1-9].*",types)]
  supertype[grepl("^FS[1-9].*",types)] <- stringr::str_extract(types,"FS[1-9]")[grepl("^FS[1-9].*",types)]
  supertype[grepl("^FR[1-9].*",types)] <- stringr::str_extract(types,"FR[1-9]")[grepl("^FR[1-9].*",types)]
  supertype[grepl("^SA[1-9].*",types)] <- stringr::str_extract(types,"SA[1-9]")[grepl("^SA[1-9].*",types)]
  supertype[grepl("^ER[1-6].*",types)] <- stringr::str_extract(types,"ER[1-6]")[grepl("^ER[1-6].*",types)]
  supertype[grepl("^PFN[a|d|m|p|v].*",types)] <- stringr::str_extract(types,"PFN[a|d|m|p|v]")[grepl("^PFN[a|d|m|p|v].*",types)]
  supertype[grepl("^MBON[0-5][0-9].*",types)] <- stringr::str_extract(types,"MBON[0-5][0-9]")[grepl("^MBON[0-5][0-9].*",types)]
  supertype[grepl("^ExR.*",types)] <-  stringr::str_extract(types,"ExR[1-9]")[grepl("^ExR.*",types)]
  supertype[grepl("SpsP.*",types)] <- "SPS-PB"
  supertype[grepl("^OA_V.*",types)] <- "OA"
  supertype[grepl("^TuBu.*",types)] <- "TuBu"
  supertype[grepl("^Delta7",types)] <- "Delta7"
  supertype[grepl("^LPsP.*",types)] <- "LPsPB"
  supertype[grepl("^KC.*",types)] <- "KC"
  if (unicodeDelta){supertype <- stringr::str_replace(supertype,"Delta","\u0394")}
  
  if (level == 1){return(supertype)}

  supertype[grepl("^FB[1-9].*",types)] <- "FBt"
  supertype[grepl("^vDelta.*",types)] <- "vDelta"
  supertype[grepl("^hDelta.*",types)] <- "hDelta"
  supertype[grepl("^PFN[a|d|m|p|v].*",types)] <- "PFN"
  supertype[grepl("^LNO.*|^LCNO.*|^GLNO.*",types)] <- "LNO"
  supertype[grepl("^FC[1-9].*",types)] <- "FC"
  supertype[grepl("^FR[1-9].*",types)] <- "FR"
  supertype[grepl("^FS[1-9].*",types)] <- "FS"
  supertype[grepl("^ER[1-6].*",types)] <- "ER"
  supertype[grepl("^SA[1-9].*",types)] <- "SA"
  supertype[grepl("^MBON[0-5][0-9].*",types)] <- "MBON"
  supertype[grepl("^ExR.*",types)] <- "ExR"
  supertype[grepl("^Giant.*|^MDN.*|^DN[a-z].*|^DN\\\\\\\\?.*",types)] <- "DN"

  if (unicodeDelta){supertype <- stringr::str_replace(supertype,"Delta","\u0394")}
  if (level == 2){return(supertype)}

  supertype[grepl("^Delta7|^P[1|6].*",types)] <- "PB Interneurons"
  supertype[grepl("^PF.*",types)] <- "FB Columnar"
  supertype[grepl("^EPG.*|^PEG.*|^PEN.*|^EL.*",types)] <- "EB Columnar"
  supertype[grepl("^FC[1-9].*|^FR[1-9].*|^FS[1-9].*",types)] <- "FB Output"
  supertype[grepl("^FB[1-9].*",types)] <- "FB Tangential"
  supertype[grepl("^[h|v]Delta.*",types)] <- "FB Interneuron"
  supertype[grepl("^PPL.*|^PAM.*|^PPM.*",types)] <- "DAN"
  supertype[grepl(".*5-HT.*|.*5HT.*",types)] <- "5HT"
  supertype[grepl("^OA-.*",types)] <- "OA"
  supertype[grepl("SpsP.*|^LPsP.*",types)] <- "PB Input"
  supertype[grepl("^AstA.*|^CRZ.*|^DSKMP.*|^NPFL.*|^PI1.*|^SIF.*",types)] <- "Peptidergic"            
  supertype[grepl("^DN1.*|l-LNv.*|LNd.*|^LPN.*|s-LNv.*",types)] <- "Clock"
  supertype[grepl("^ovi.*|^aSP.*|^aDT.*|^aIP.*|^pC1.*|^vpo.*",types)] <- "Fru"
  supertype[grepl("^HRN.*|^JO.*|^OGC.*|^ORN.*|^TRN.*",types)] <- "Other Sensory"
  supertype[grepl("^LH.*",types)] <- "LH"
  supertype[grepl("^aMe.*|^CTX.*|^DCH.*|^H[1-3].*|^HBeyelet.*|^HS.*|^LC[1-9].*|^Li[1-9].*|^LLPC.*|^LPC.*|^LT[1-9].*|^MC[1-9].*|^VCH.*|^VS.*",types)] <- "Visual PNs"
  supertype[grepl("^AL-.*|^AL[B|I]N.*|^D_ad.*|^D[A|C|L|M|P][1-9].*|^il3.*|^l2LN.*|^lLN[1-9].*|^M_.*|^mAL.*|^MZ_.*|^v2LN[1-9].*|^V[A|C|L][1-9].*|^vLN.*|^V[M|P][1-9].*|^Z_v.*|^Z_lv.*",types)] <- "Antennal lobe"
  
  supertype[is.na(supertype)] <- "Unassigned"
  supertype[! supertype %in% c("PB Interneurons","PB Input","FB Columnar","EB Columnar","ER","ExR","FB Output","FB Tangential","FB Interneuron","LNO","DN","KC","MBON","DAN","5HT","OA","Peptidergic","Clock","Fru","Other Sensory","LH","Visual PNs","Antennal lobe","Unassigned")] <- "Terra incognita"
  
  if (unicodeDelta){supertype <- stringr::str_replace(supertype,"Delta","\u0394")}
  supertype
}

#' @export
supertype.neuronBag <- function(types,unicodeDelta=TRUE){
  for (lev in 1:3){
    for (ty in c(".from",".to")){
      for (tab in c("inputs","outputs","inputs_raw","outputs_raw")){
        types[[tab]][[paste0("supertype",lev,ty)]] <- supertype(types[[tab]][[paste0("databaseType",ty)]],level=lev,unicodeDelta=unicodeDelta)
      }
    }
    types$names[[paste0("supertype",lev)]] <-  supertype(types$names[[paste0("databaseType")]],level=lev,unicodeDelta=unicodeDelta)
    types$outputsTableRef[[paste0("supertype",lev)]] <-  supertype(types$outputsTableRef[[paste0("databaseType")]],level=lev,unicodeDelta=unicodeDelta)
  }
  types
}

#' @export
supertype.data.frame <- function(types,level=1:3,unicodeDelta=TRUE){
  renamable <- names(types)[names(types) %in% c("databaseType","databaseType.from","databaseType.to")]
  for (lev in level){
    for (ty in renamable){
        types[[sub("databaseType",paste0("supertype",lev),ty)]] <- supertype(types[[ty]],level=lev,unicodeDelta=unicodeDelta)
    }
  }
  types
}

supertype.logical <- function(types,level=NULL,unicodeDelta=FALSE){
  stopifnot(length(types)==0)
  return(character(0))
}

#' Select a set of Supertypes from a Supertypes table as the one generated by supertype (when applied to a data frame)
#' @param supertypeTable A supertypes table as returned by \code{supertype}
#' @param default_level An integer specifying which Supertype level to use by default
#' @param exceptions A list of the form list(nameOfType = level) where nameOfType is the name of
#' Supertypes at level \code{exceptionLevelMatch} for which one want to use a different level of description
#' @param exceptionLevelMatch What level to use in the \code{exceptions} list
#' @return a data.frame with the same columns as those returned by \code{supertype} with an extra \code{type} column
#' containing the desired set
#' @details this is to supertypes what \code{selectRoiSet} is to ROIs
#' @seealso \code{supertype}, \code{selectRoiSet}
#' @export

selectSupertypeSet <- function(supertypeTable,default_level=2,exceptions=NULL,exceptionLevelMatch = default_level){
  supertypeTable$supertype0 <- supertypeTable$databaseType
  if (!is.null(exceptions)){
    levelEx <- paste0("supertype",exceptionLevelMatch)
    normalTypes <- supertypeTable %>% filter(!((!!as.name(levelEx)) %in% names(exceptions))) %>%
      mutate(supertype = (!!as.name(paste0("supertype",default_level))))

    exceptionsTypes <- supertypeTable %>% filter(((!!as.name(levelEx)) %in% names(exceptions)))

    typesEx <- as.character(exceptionsTypes[[levelEx]])
    customLev <- sapply(typesEx,function(r) paste0("supertype",exceptions[[r]]))
    exceptionsTypes$supertype <- sapply(1:length(customLev),function(i) as.character(exceptionsTypes[[customLev[i]]][i]))

    types <- bind_rows(normalTypes,exceptionsTypes)
  }else{
    types <- supertypeTable %>% mutate(supertype = (!!(as.name(paste0("supertype",default_level)))))
  }

  types <- types %>% arrange(supertype3) %>%
    mutate(supertype = factor(supertype,levels=unique(supertype)))

  return(distinct(types))
}

#' Create a consistent palette for supertypes for plotting
#' @param typeSelection A supertype selection as returned by \code{selectSupertypeSet}
#' @param my_palette A discrete color palette to use
#' @return A named vector of colors (for example to be used in ggplot2's \code{scale_color_manual})
#' @seealso \code{selectSupertypeSet}, \code{roisPalette}
#' @export
typesPalette <- function(typeSelection,my_palette=paletteer::paletteer_d("Polychrome::palette36")){
  typeSelection <- unique(c(as.character(typeSelection),"Other"))
  if (length(typeSelection)>36) warning(paste0(length(typeSelection)," levels in your palette, this is likely too many."))
  pal <- my_palette[1:length(typeSelection)]
  oPos <- which(typeSelection=="Other")
  typeSelection[c(2,oPos)] <- typeSelection[c(oPos,2)]  ## Exchange to get "other" at a neutral color
  names(pal) <- typeSelection
  pal
}
