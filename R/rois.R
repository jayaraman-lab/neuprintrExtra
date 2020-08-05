#' Get all instances of neurons of types with non zero pre/post counts in a ROI.
#' @param ROI The ROI to look for
#' @param minTypePercentage The minimum proportion of the instances of a type that should be innervating the ROI for
#' it to be considered
#' @param retyping a retyping function to be applied to the types found. \code{identity} by default (so it does nothing)
#' @return  a data frame of metadata for all neurons in the ROI, as returned by \code{neuprint_get_meta}, with extra columns \code{ROI_pre}
#'  and \code{ROI_post}, the counts in the queried ROI.
#' @details  If a type is selected because at least \code{minTypePercentage} of its instances touch the ROI, all instances of the type are returned.
#' This is used internally by \code{getTypesInRoiTable}
#' @seealso  \code{getTypesInRoiTable}
#' @export
getNeuronsInRoiTable <- function(ROI,minTypePercentage=0.5,retyping=identity) {
  roi_Innervate <- neuprint_bodies_in_ROI(ROI) %>%
    mutate(originalInstance = TRUE)
  metaRoi <- neuprint_get_meta(roi_Innervate) %>% tidyr::drop_na(type)

  ## Get all instances of the types touching the ROI
  all_neurons <- getTypesTable(unique(metaRoi$type))

  ## Join to the full table of type members (and fill with zero where the extra instances do not
  ## innervate)
  roi_Innervate <- left_join(all_neurons,roi_Innervate,by=c("bodyid","pre","post","voxels"))
  roi_Innervate <- roi_Innervate %>% select(-c(voxels,cellBodyFiber)) %>%
    tidyr::replace_na(list(ROI_pre = 0,ROI_post = 0,originalInstance=FALSE)) %>%
    mutate(databaseType = as.character(type)) ## Convenience column for when types are changed

  roi_Innervate <- retyping(roi_Innervate)
  roi_Innervate <-roi_Innervate %>% group_by(type) %>%
    mutate(typePercentage = sum(originalInstance)/n()) %>%
    ungroup() %>%
    filter(typePercentage > minTypePercentage)

  return(roi_Innervate)
}

#' Returns a neuronBag object of all the neurons forming significant connections in a ROI.
#' @param ROI The ROI to consider
#' @param retyping a retyping function to be applied to the types found. \code{cxRetyping} by default (for no retyping, set it to \code{identity})
#' @param bagROIs Which ROIs to include in the bag created (by default only the ROI one wants neurons in). If NULL returns all ROIs.
#' @param minTypePercentage The minimum proportion of the instances of a type that should be innervating the ROI for
#' it to be considered (0.5 by default)
#' @param ... Parameters to be passed to \code{neuronBag}
#' @details calls \code{getNeuronsInRoiTable} internally, with \code{minTypePercentage} divided by 2 if \code{lateralize} is TRUE.
#' @seealso  \code{getNeuronsInRoiTable}
#' @export
getTypesInRoiTable <- function(ROI,
                               retyping=cxRetyping,
                               bagROIs=ROI,
                               minTypePercentage=0.5,
                               ...){
  neuronTable <- getNeuronsInRoiTable(ROI,minTypePercentage=minTypePercentage,retyping=retyping) ## Remove types if less than
  ## 25% of the instances touch (l/R)
  roiConnections <- neuronBag(neuronTable,slctROI=bagROIs,selfRef=TRUE,...)
  roiConnections <- retyping(roiConnections)
  roiConnections
}

#' Extract types participating in a significant connection in a ROI
#' @param roiConnections A \code{neuronBag} object
#' @param ROI The ROI to consider
#' @return A dataframe with columns "type" and "databaseType"
#'
#' @export
typesInROI <- function(roiConnections,ROI){
  typesUnfiltered <- roiConnections$names$type
  inputs <- roiConnections$inputs %>% filter((roi == ROI) & (type.to %in% typesUnfiltered) &
                                               (type.from %in% typesUnfiltered))
  outputs <- roiConnections$outputs %>% filter((roi == ROI) & (type.to %in% typesUnfiltered) &
                                                 (type.from %in% typesUnfiltered))

  roiTypes <- data.frame(type = unique(c(inputs$type.to,outputs$type.from))) %>%
    mutate(databaseType = roiConnections$names$databaseType[match(type,roiConnections$names$type)])

  return(roiTypes)
}

#' A wrapper around \code{neuprint_ROI_hierarchy} returning a cleaned up, sorted table with various ROI levels.
#' @return  A table with columns `level` 0 through 4, 0 being the more general and 4 the more detailed, omitting tracts. ROI are somewhat sorted,
#'  from right to central to left (and from periphery to center to periphery). level0 is similar to level1 except that Left/Right distinctions are omitted.
#' @seealso  \code{selectRoiSet} to make a selection of ROIs (for example for a figure) from such a hierarchy
#' @export
getRoiTree <- function(){
  roiH <- neuprint_ROI_hierarchy() %>% mutate_all(as.character)
  roiT <- data.frame(level1 = roiH$roi[roiH$parent == "hemibrain"],stringsAsFactors = F) %>% filter(!(level1 %in% c("hemibrain","AOT(R)","GC","GF(R)","mALT(R)","POC","mALT(L)")))
  roiT <- left_join(roiT,roiH,by=c("level1"="parent")) %>% rename(level2 = roi) %>% mutate(level2 = ifelse(is.na(level2),level1,level2))
  roiT <- left_join(roiT,roiH,by=c("level2"="parent")) %>% rename(level3 = roi) %>% mutate(level3 = ifelse(is.na(level3),level2,level3))
  roiT <- left_join(roiT,roiH,by=c("level3"="parent")) %>% rename(level4 = roi) %>% mutate(level4 = ifelse(is.na(level4),level3,level4))
  roiT <- roiT %>% mutate(side4 = "Central",
                          side2 = "Central")
  roiT$side4[grepl("(L",roiT$level4,fixed=T)] <- "Left"
  roiT$side4[grepl("(R",roiT$level4,fixed=T)] <- "Right"
  roiT$side2[grepl("(L",roiT$level2,fixed=T)] <- "Left"
  roiT$side2[grepl("(R",roiT$level2,fixed=T)] <- "Right"

  roiT$side4 <- factor(roiT$side4,levels=(c("Right","Central","Left")))
  roiT$side2 <- factor(roiT$side2,levels=(c("Right","Central","Left")))
  roiT$level1 <- factor(roiT$level1,levels= c("OL(R)","AL(R)","MB(+ACA)(R)","LH(R)","PENP","GNG","VLNP(R)","SNP(R)","VMNP","INP","LX(R)","CX","LX(L)","SNP(L)","MB(L)","AL(L)"))

  roiT <- arrange(roiT,side2,level1)
  roiT$level0 <- delateralize(roiT$level1)
  roiT <- roiT %>% mutate_at(c("level2","level3","level4"),function(a) factor(a,levels=unique(a)))
  roiT
}

delateralize <- function(roiName){
  gsub("(L)","",gsub("(R)","",roiName,fixed=TRUE),fixed=TRUE)
}

#' Select a set of ROIs from a ROI hierarchy as the one generated by getRoiTree
#' @param roiTree A ROI hierarchy table as returned by \code{getRoiTree}
#' @param default_level An integer specifying which ROI level to use by default
#' @param exceptions A list of the form list(nameOfRoi = level) where nameOfRoi is the name of
#' ROIs at level \code{exceptionLevelMatch} for which one want to use a different level of description
#' @param exceptionLevelMatch What level to use in the \code{exceptions} list. Can be either a scalar or a vector the same
#' length as \code{exceptions}
#' @return a data.frame with the same columns as those returned by \code{getRoiTree} with an extra \code{roi} column
#' containing the desired set
#' @seealso \code{getRoiTree}
#' @export
selectRoiSet <- function(roiTree=getRoiTree(),default_level=2,exceptions=NULL,exceptionLevelMatch = default_level){
  if (!is.null(exceptions)){

    if (length(exceptionLevelMatch) == 1 & length(exceptions)>=1) exceptionLevelMatch <- rep(exceptionLevelMatch,length(exceptions))
    levelEx <- paste0("level",exceptionLevelMatch)

    normalRois <- roiTree
    exceptionsRois <- bind_rows(lapply(1:length(exceptions),function(i) filter(roiTree,
                                                                               ((!!as.name(levelEx[i])) %in% names(exceptions[i]))) %>%
                                         mutate(roi=!!as.name(paste0("level",exceptions[[i]])))))

    for (i in 1:length(exceptions)){
        normalRois <- filter(normalRois,!((!!as.name(levelEx[i])) %in% names(exceptions)[i]))
    }

    normalRois <- mutate(normalRois,roi = (!!as.name(paste0("level",default_level))))

    rois <- bind_rows(normalRois,exceptionsRois)
  }else{
    rois <- roiTree %>% mutate(roi = (!!(as.name(paste0("level",default_level)))))
  }

  rois <- rois %>% arrange(side2,level1) %>%
    mutate(roi = factor(roi,levels=unique(roi)))

  return(distinct(rois))
}

#' Create a consistent palette for ROIs for plotting
#' @param rois A ROI hierarchy as returned by \code{getRoiTree}
#' @param favoriteRegion A brain region (defined at level 1) for which one wants colors
#' defined down to level 2
#' @param my_palette A discrete color palette to use
#' @return A named vector of colors (for example to be used in ggplot2's \code{scale_color_manual})
#' @seealso \code{getRoiTree}
#' @export
roisPalette <- function(rois=getRoiTree(),favoriteRegion="CX",my_palette=paletteer::paletteer_d("Polychrome::palette36")){
  roiL <- unique(delateralize(c(as.character(rois$level1),as.character(rois$level2[rois$level1==favoriteRegion]))))
  pal <- my_palette[1:length(roiL)]
  names(pal) <- roiL
  pal
}

#' Get 2D outlines for a ROI
#' @param roi Either a mesh3d object (as returned by \code{neuprint_ROI_mesh}) or the name of an
#' existing ROI
#' @param alpha Alpha parameter to be passed to the \code{alphahull::ahull} function
#' @param roiName A name to label the outline with
#' @return a data frame with columns x, y (outlines) proj (either xy or xz, the projection considered) and roi (the roi name)
#' @details both xy and xz projections are returned, which can be selected with the \code{proj} column. This is meant to
#' be handy for ggplot2 plotting
#' @export
roiOutline <- function(roi,alpha=100,roiName){UseMethod("roiOutline")}

#' @export
roiOutline.mesh3d <- function(roi,alpha=100,roiName =deparse(substitute(roi))){
  roiPts <-  data.frame(nat::dotprops(roi)$points)
  names(roiPts) <- c("x","y","z")
  roiHullxy <- alphahull::ahull(x=roiPts$x,y=roiPts$y,alpha=alpha)
  roiHullxz <- alphahull::ahull(x=roiPts$x,y=roiPts$z,alpha=alpha)

  roiOutxy <- data.frame(roiHullxy$arcs) %>% mutate(x=c1,y=c2,proj="xy",roi=roiName) %>% select(x,y,proj,roi)
  roiOutxy <-  rbind(roiOutxy,roiOutxy[1,])

  roiOutxz <- data.frame(roiHullxz$arcs) %>% mutate(x=c1,y=c2,proj="xz",roi=roiName) %>% select(x,y,proj,roi)
  roiOutxz <- rbind(roiOutxz,roiOutxz[1,])

  rbind(roiOutxy,roiOutxz)
}

#' @export
roiOutline.character <- function(roi,alpha=100){
  roiMesh <- neuprint_ROI_mesh(roi)
  roiOutline(roiMesh,alpha=alpha,roiName=roi)
}

#' Sum a set of ROIs in a connection table or neuronBag
#' @param connections A raw connection table or a neuronBag
#' @param rois The ROIs to combined
#' @param newRoi The name of the newly formed ROI
#' @param ... Parameters to be passed to getTypeToTypeTable
#'
#' @export
combineRois <- function(connections,rois,newRoi,...){UseMethod("combineRois")}

#' @export
combineRois.data.frame <- function(connections,rois,newRoi){
  ## CHECK IT'S A RAW CONNECTION TABLE
  newRegionTable <- connections %>%
    filter(roi %in% rois) %>%
    group_by(to) %>%
    mutate_at(vars(any_of(c("totalROIweight","knownTotalROIweight"))),~sum(.[match(rois,roi)],na.rm=TRUE)) %>%
    group_by(from) %>%
    mutate_at(vars(any_of(c("totalPreROIweight","knownTotalPreROIweight"))),~sum(.[match(rois,roi)],na.rm=TRUE)) %>%
    group_by_if(names(.) %in% c(paste0(c("","name.","type.","databaseType.","previous.type."),"to"),
                                paste0(c("","name.","type.","databaseType.","previous.type."),"from"),
                                paste0("supertype",1:3,".to"),paste0("supertype",1:3,".from")))

  if ("knownWeightRelative" %in% names(newRegionTable)){
    newRegionTable <- summarize(newRegionTable,roi=newRoi,
                                weight=weight[1],
                                totalPreWeight=totalPreWeight[1],
                                knownTotalPreWeight=knownTotalPreWeight[1],
                                weightRelativeTotal=weightRelativeTotal[1],
                                knownWeightRelativeTotal=knownWeightRelativeTotal[1],
                                outputContributionTotal=outputContributionTotal[1],
                                knownOutputContributionTotal = knownOutputContributionTotal[1],
                                totalROIweight = totalROIweight[1],
                                knownTotalROIweight = knownTotalROIweight[1],
                                totalPreROIweight=totalPreROIweight[1],
                                knownTotalPreROIweight=knownTotalPreROIweight[1],
                                ROIweight=sum(ROIweight),
                                weightROIRelativeTotal=sum(weightROIRelativeTotal)) %>%
      ungroup() %>%
      mutate(weightRelative=ROIweight/totalROIweight,
             knownWeightRelative=ROIweight/knownTotalROIweight,
             outputContribution=ROIweight/totalPreROIweight,
             knownOutputContribution=ROIweight/knownTotalPreROIweight)
  }else{
    newRegionTable <- summarize(newRegionTable,roi=newRoi,
              weight=weight[1],
              totalPreWeight=totalPreWeight[1],
              weightRelativeTotal=weightRelativeTotal[1],
              outputContributionTotal=outputContributionTotal[1],
              totalROIweight = totalROIweight[1],
              totalPreROIweight=totalPreROIweight[1],
              ROIweight=sum(ROIweight),
              weightROIRelativeTotal=sum(weightROIRelativeTotal)) %>%
    ungroup() %>%
    mutate(weightRelative=ROIweight/totalROIweight,
           outputContribution=ROIweight/totalPreROIweight)}

  newRegionTable
}


#' @export
combineRois.neuronBag <- function(connections,rois,newRoi,...){
  new_inputsR <- combineRois(connections$inputs_raw,rois,newRoi)
  new_outputsR <- combineRois(connections$outputs_raw,rois,newRoi)
  
  if("allInsToOuts" %in% names(connections[["ref"]])){
    connections$ref$allInsToOuts <- combineRois(connections$ref$allInsToOuts,rois,newRoi)
    connections$ref$outputs_ref <-getTypeToTypeTable(connections$ref$allInsToOuts,typesTable = connections$ref$outputTableRefFull,...)
    new_outputs <- processTypeToTypeFullOutputs(connections$ref$outputs_ref,new_outputsR)
  }else{
    new_outputs <- getTypeToTypeTable(new_outputsR,typesTable = connections$outputsTableRef,...)
  }
  if("allOutsFromIns" %in% names(connections[["ref"]])){
    connections$ref$allOutsFromIns <- combineRois(connections$ref$allOutsFromIns,rois,newRoi)
    connections$ref$inputs_ref <- getTypeToTypeTable(connections$ref$allOutsFromIns,typesTable=connections$ref$inputsTableRefFull,...)
    new_inputs <- processTypeToTypeFullInputs(connections$ref$inputs_ref,new_inputsR)
  }else{
    new_inputs <- getTypeToTypeTable(new_inputsR,typesTable = connections$names,...)
  }
  nBag <- new_neuronBag(outputs = new_outputs,
                inputs = new_inputs,
                names = connections$names,
                inputs_raw = new_inputsR,
                outputs_raw = new_outputsR,
                outputsTableRef = connections$outputsTableRef)
  if("ref" %in% names(connections)){nBag$ref <- connections$ref}
  nBag
}

#' Build a per roi summary of innervation for neurons in a neuronBag
#' @param neurons : a dataframe as returned by a search or a neuronBag object
#' @param threshold : the minimal average number of synapses for a type in a ROI
#' for it to be included
#' @param rois : a roiset to consider (if NULL consider all rois)
#'
#'@export
getROISummary <- function(neurons,threshold=0,rois=NULL){UseMethod("getROISummary")}

#

#'@export
getROISummary.neuronBag <- function(neurons,threshold=0,rois = NULL){
  getROISummary(neurons$names,rois=rois,threshold=threshold)
}

#'@export
getROISummary.data.frame <- function(neurons,threshold=0,rois = NULL){

  roiSummary <- getRoiInfo(neurons)

  if (!(is.null(rois))){
    roiSummary <- roiSummary %>%
      filter(roi %in% rois$roi)
  }

  countInstances <- group_by(neurons,type) %>% summarize(n=n())

  roiSummary <-
    left_join(roiSummary,
              select(neurons,bodyid,type,databaseType),
              by=c("bodyid")) %>% tidyr::replace_na(list(downstream=0,upstream=0)) %>%
    group_by(roi,type,databaseType) %>%
    summarize(downstream=mean(downstream),
              upstream=mean(upstream)) %>%
    ungroup() %>%
    mutate(fullWeight = downstream+upstream,
           deltaWeight = (downstream - upstream)/fullWeight) %>%
    filter(fullWeight>threshold)

  roiSummary <- left_join(roiSummary,countInstances,by="type")

  if (!is.null(rois)){roiSummary <- left_join(roiSummary,rois,by=roi)}
  return(supertype(roiSummary))
}

