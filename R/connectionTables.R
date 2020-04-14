#' Get connection table of inputs and add meta
#' @return Returns a connection table as data frame. Added columns are \code{weightRelativeTotal} which is
#' the relative weight considering all the synapses (irrespective of the ROI), and if ROI are used (either if
#' \code{slctROI} has a value or \code{by.roi} is \code{TRUE}), \code{weightRelative} is the relative weight in the
#' ROI and \code{totalROIweight} is the absolute number of inputs this neuron receives in that region and
#' \code{weightROIRelativeTotal} is the weight in the ROI normalized by the total number of inputs (in all ROIs)
#' @param bodyIDs: The bodyids of neurons who's connections should be queried or a metadata data.frame
#' @param synapseType: Choose "PRE" or "POST" to get inputs or outputs of the neurons in bodyIDs, respectivly. "PRE" is
#' usually slower as it requires listing all the outputs of the input neurons to get the contribution of the outputs.
#' @param slctROI: String specifying the ROI where connections should be queried. By default all the ROIs.
#' @param by.roi: Passed to neuprint_connection_table. If returning all ROIs, should results be broken down by ROI?
#' @param synThresh: Minimum number of synapses to consider a connection (default 3)
#' @param chunk_meta : to be passed to metadata and roiInfo queries (default TRUE)
#' @param chunk_connection : to be passed to neuprint_connection_table (default TRUE)
#' @param verbose : should the function report on its progress?
#' @param computeKnownRatio : should the function also compute ratios (weightRelative and outputContribution) relative to known synaptic partners rather
#' than relative to the total number of synapses
#' @param ...: Other arguments to be passed to neuprint queries (neuprint_connection_table, getRoiInfo and neuprint_get_meta)
#' @export
getConnectionTable <- function(bodyIDs,synapseType, slctROI,by.roi, synThresh = 3,chunk_connections=TRUE,chunk_meta=TRUE,verbose=FALSE,computeKnownRatio=FALSE,...){
  UseMethod("getConnectionTable")}


getConnectionTable.default = function(bodyIDs, synapseType, slctROI=NULL,by.roi=FALSE, synThresh=3,chunk_connections=TRUE,chunk_meta=TRUE,verbose=FALSE,computeKnownRatio=FALSE,...){
  refMeta <- neuprint_get_meta(bodyIDs,chunk=chunk_meta,...)
  return(getConnectionTable(refMeta,synapseType,slctROI,by.roi,synThresh,chunk_connections=chunk_connections,chunk_meta=chunk_meta,computeKnownRatio=computeKnownRatio,...))
}

getConnectionTable.data.frame <- function(bodyIDs,synapseType, slctROI=NULL,by.roi=FALSE,synThresh=3,chunk_connections=TRUE,chunk_meta=TRUE,verbose=FALSE,computeKnownRatio=FALSE,...){
  refMeta <- bodyIDs
  bodyIDs <- neuprint_ids(bodyIDs$bodyid)

  if (verbose) message("Pull connections")
  myConnections_raw <- neuprint_connection_table(bodyIDs, synapseType, slctROI,by.roi=by.roi,chunk=chunk_connections,...)

  if (by.roi | !is.null(slctROI)){
    myConnections_raw <- myConnections_raw %>% drop_na(ROIweight)
    myConnections <- myConnections_raw %>% filter(ROIweight>synThresh)
  }else{
    myConnections <- myConnections_raw %>% filter(weight>synThresh)
  }

  if (nrow(myConnections)==0){return(NULL)}

  if (verbose) message("Pull metadata")
  partnerMeta <- neuprint_get_meta(myConnections$partner,chunk=chunk_meta,...)
  refMetaOrig <- neuprint_get_meta(myConnections$bodyid,chunk=chunk_meta,...)  ## To get the database type name

  myConnections <- filter(myConnections,partnerMeta$status =="Traced")
  partnerMeta <- filter(partnerMeta,status == "Traced")

  processConnectionTable(myConnections,myConnections_raw,refMeta,partnerMeta,refMetaOrig,synapseType,by.roi,slctROI,verbose,chunk_meta,chunk_connections,computeKnownRatio,...)
}

processConnectionTable <- function(myConnections,myConnections_raw,refMeta,partnerMeta,refMetaOrig,synapseType,by.roi,slctROI,verbose,chunk_meta,chunk_connections,computeKnownRatio,...){


  refMeta <- slice(refMeta,as.integer(sapply(myConnections$bodyid,function(b) match(b,refMeta$bodyid))))
  refMetaOrig <- slice(refMetaOrig,as.integer(sapply(myConnections$bodyid,function(b) match(b,refMetaOrig$bodyid))))

  myConnections <-mutate(myConnections,
                         partnerName = partnerMeta[["name"]],
                         name = refMeta[["name"]],
                         partnerType = partnerMeta[["type"]],
                         type = refMeta[["type"]])

  myConnections <- simplifyConnectionTable(myConnections)
  ## Normalization is always from the perspective of the output (fraction of inputs to the output neuron)
  if (synapseType == "PRE"){
    if (computeKnownRatio){
      knownTablePost <- myConnections_raw  %>% mutate(from = ifelse(prepost==1,bodyid,partner),
                                                      to = ifelse(prepost==1,partner,bodyid)) %>% select(-bodyid,-partner)
      knownTablePre <- neuprint_connection_table(unique(myConnections$from),"POST",slctROI,by.roi=by.roi,chunk=chunk_connections,...) %>% drop_na(ROIweight) %>% mutate(from = ifelse(prepost==1,bodyid,partner),
                                                                                                                                                                        to = ifelse(prepost==1,partner,bodyid)) %>% select(-bodyid,-partner)
    }
    outMeta <- refMeta
    inMeta <- partnerMeta
    myConnections <- mutate(myConnections,databaseType.to = refMetaOrig$type,
                            databaseType.from = type.from)
  } else {
    if (computeKnownRatio){
      knownTablePre <- myConnections_raw  %>% mutate(from = ifelse(prepost==1,bodyid,partner),
                                                     to = ifelse(prepost==1,partner,bodyid)) %>% select(-bodyid,-partner)
      knownTablePost <- neuprint_connection_table(unique(myConnections$to),"PRE",slctROI,by.roi=by.roi,chunk=chunk_connections,...) %>% drop_na(ROIweight) %>% mutate(from = ifelse(prepost==1,bodyid,partner),
                                                                                                                                                                      to = ifelse(prepost==1,partner,bodyid)) %>% select(-bodyid,-partner)
    }
    inMeta <- refMeta
    outMeta <- partnerMeta
    myConnections <- mutate(myConnections,databaseType.to = type.to,
                            databaseType.from = refMetaOrig$type)
  }

  if (computeKnownRatio){
    knownTablePostTotal <- knownTablePost %>% group_by(to) %>% distinct(from,weight) %>% mutate(knownPostWeight = sum(weight)) %>% ungroup()
    knownTablePreTotal <- knownTablePre %>% group_by(from) %>% distinct(to,weight) %>% mutate(knownPreWeight = sum(weight)) %>% ungroup()
  }

  myConnections <-mutate(myConnections,weightRelativeTotal = weight/outMeta[["post"]],
                         totalPreWeight = inMeta[["downstream"]][match(myConnections$from,inMeta$bodyid)],
                         outputContributionTotal = weight/totalPreWeight,
                         previous.type.to = databaseType.to,
                         previous.type.from = databaseType.from
  )
  if  (computeKnownRatio){
    myConnections <-mutate(myConnections,
                           knownWeightRelativeTotal = weight/knownTablePostTotal$knownPostWeight[match(myConnections$to,knownTablePostTotal$to)],
                           knownTotalPreWeight = knownTablePreTotal$knownPreWeight[match(myConnections$from,knownTablePreTotal$from)],
                           knownOutputContributionTotal = weight/knownTotalPreWeight
    )

  }

  if (by.roi | !is.null(slctROI)){
    if (verbose) message("Pull roiInfo")
    myConnections[["weightROIRelativeTotal"]] <- myConnections[["ROIweight"]]/outMeta[["post"]]
    outInfo <- getRoiInfo(unique(myConnections$to),chunk=chunk_meta,...) %>% select(bodyid,roi,post)
    inInfo <- getRoiInfo(unique(myConnections$from),chunk=chunk_meta,...) %>% select(bodyid,roi,downstream)

    myConnections <- left_join(myConnections,inInfo,by=c("from" = "bodyid","roi"="roi")) %>% rename(totalPreROIweight = downstream)

    myConnections <- left_join(myConnections,outInfo,by=c("to" = "bodyid","roi"="roi")) %>% rename(totalROIweight = post) %>%
      mutate(weightRelative=ROIweight/totalROIweight,
             outputContribution=ROIweight/totalPreROIweight)

    if (computeKnownRatio){
      knownTablePostROI <- knownTablePost %>% group_by(to,roi)  %>% summarize(knownPostWeight = sum(ROIweight)) %>% ungroup()
      knownTablePreROI <- knownTablePre %>% group_by(from,roi)  %>% summarize(knownPreWeight = sum(ROIweight)) %>% ungroup()
      ## This is how much this connection accounts for the outputs of the input neuron (not the standard measure)
      myConnections <- myConnections %>% mutate(knownWeightRelative = ROIweight/knownTablePostROI$knownPostWeight[match(paste0(myConnections$to,myConnections$roi),paste0(knownTablePostROI$to,knownTablePostROI$roi))],
                                                knownTotalPreROIweight = knownTablePreROI$knownPreWeight[match(paste0(myConnections$from,myConnections$roi),paste0(knownTablePreROI$from,knownTablePreROI$roi))],
                                                knownOutputContribution = ROIweight/knownTotalPreROIweight) %>% drop_na(weightRelative)  ## NA values can occur in rare cases where
    }
    ## synapse (pre/post) is split between ROIs
  }else{
    myConnections <- mutate(myConnections,roi = "All brain",
                            outputContribution = outputContributionTotal,
                            weightRelative = weightRelativeTotal,
                            ROIweight = weight)
    if (computeKnownRatio){
      myConnections <- mutate(myConnections,
                              knownOutputContribution = knownOutputContributionTotal,
                              knownWeightRelative = knownWeightRelativeTotal)

    }
  }
  return(supertype(myConnections))
}

#' Change from a bodyid/partner/prepost to a from/to format
#' @param connectionTable: A data frame in the bodyid/partner/prepost format (as returned by getConnectionTable)
#' @return A data frame in the from/to/name.from/name.to... format
#'
simplifyConnectionTable <- function(connectionTable){

  if ("from" %in% names(connectionTable)){return(connectionTable)}else{
    connectionTable <- connectionTable %>% mutate(from = ifelse(prepost==1,bodyid,partner),
                                                  to = ifelse(prepost==1,partner,bodyid),
                                                  name.from = as.character(ifelse(prepost==1,name,partnerName)),
                                                  name.to = as.character(ifelse(prepost==1,partnerName,name)),
                                                  type.from = as.character(ifelse(prepost==1,type,partnerType)),
                                                  type.to = as.character(ifelse(prepost==1,partnerType,type))
    ) %>%
      select(-bodyid,-partner,-name,-partnerName,-partnerType,-type,-prepost)
    return(connectionTable)
  }
}
