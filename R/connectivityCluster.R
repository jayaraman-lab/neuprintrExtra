new_connectivityCluster <- function(inputsTable,outputsTable,names,distance,hc,grouping){
  res <- list(inputsTable=inputsTable,
              outputsTable=outputsTable,
              names=names,
              distance=distance,
              hc=hc,
              grouping=grouping)
attr(res,"class") <- "connectivityCluster"
return(res)
}

#' Clustering neurons based on connectivity tables
#' 
#' @param inputsTable The table of input connections to the neurons to cluster
#' @param outputsTable The table of input connections to the neurons to cluster
#' @param ROIs Which ROIs to consider. If NULL (the default) uses all the ROI present
#' @param grouping grouping on which the clustering is to be made (neuron or type, or various supertypes)
#' @param distance Which distance measure to use (defaults to \code{\link{cos_dist}})
#' @param knownStats Whether or not to use \code{knownOutputContribution} and \code{knownWeightRelative} instead of the standard measures
#' @param ... Arguments to be passed to \code{\link{hclust}}
#' @details \itemize{
#'  \item If both \code{inputsTable} and \code{outputsTable} are present, the clustering will be run on the combined connectivity vector. This is
#'  controlled by the \code{clusterOn} argument in \code{clusterBag}.
#'  \item \code{weightRelative} is used for working on input vectors, and \code{outputContribution} to work on output vectors
#'  }
#' @return an object of class \code{connectivityCluster}. The object is a list with components:
#' \describe{
#'   \item{inputsTable}{The table of input connections to the neurons being clustered (can be NULL if clustered on outputs)}
#'   \item{outputsTable}{The table of output connections to the neurons being clustered (can be NULL if clustered on inputs)}
#'   \item{distance}{object containing a dissimilarity measure between the neurons being considered}
#'   \item{hc}{\code{\link{hclust}} object, result of clustering}
#'   \item{grouping}{grouping on which the clustering was made (neuron or type, or various supertypes)}
#' }
#' @seealso \code{\link{clusterBag}}
#' @export
connectivityCluster <- function(inputsTable=NULL,
                                outputsTable=NULL,
                                ROIs=NULL,
                                grouping="type",
                                distance=cos_dist,
                                knownStats=FALSE,
                                ...){
  groupingOld <- grouping
  if (grepl("neuron",grouping) | grepl("bodyid",grouping)){grouping <- ""}else{
    grouping=paste0(grouping,".")}
  from <- paste0(grouping,"from")
  to <- paste0(grouping,"to")
  
  if (!is.null(inputsTable)){
    stat <- ifelse(knownStats,"knownWeightRelative","weightRelative")
    if (!is.null(ROIs)){inputsTable <- filter(inputsTable,roi %in% ROIs)}
    inputsMat <- connectivityMatrix(inputsTable,slctROIs = ROIs,allToAll = FALSE,from = from,to=to,value = stat,ref="outputs")
    connMat <- inputsMat
  }
  
  if (!is.null(outputsTable)){
    stat <- ifelse(knownStats,"knownOutputContribution","outputContribution")
    if (!is.null(ROIs)){outputsTable <- filter(outputsTable,roi %in% ROIs)}
    outputsMat <- connectivityMatrix(outputsTable,slctROIs = ROIs,allToAll = FALSE,from = from,to=to,value = stat,ref="inputs")
    connMat <- outputsMat
  }
  
  if (!is.null(inputsTable) && !is.null(outputsTable)){
    outputsMat <- outputsMat[rownames(inputsMat),,drop=FALSE]
    colnames(outputsMat) <- paste0(colnames(outputsMat),".out")
    colnames(inputsMat) <- paste0(colnames(inputsMat),".in")
    connMat <- cbind(inputsMat,outputsMat)
  }
  
  connDist <- distance(connMat)
  connHC <- hclust(connDist)
  new_connectivityCluster(inputsTable,outputsTable,rownames(connMat),connDist,connHC,groupingOld)
}

#' @describeIn connectivityCluster Cluster neurons in a bag 
#' @param nBag a neuron bag
#' @param clusterOn whether to cluster on the input vector, the output vector, or both
#' @export
clusterBag <- function(nBag,
                       ROIs=NULL,
                       clusterOn = c("inputs","outputs","both"),
                       grouping=c("type","neuron","bodyid","supertype1","supertype2","supertype3","databaseType"),
                       distance=cos_dist,
                       knownStats=FALSE,
                       ...){
  grouping <- match.arg(grouping)
  clusterOn <- match.arg(clusterOn)
  if (grouping=="neuron"){if (clusterOn=="outputs"){inp <- NULL}else{inp <- nBag$inputs_raw}
                      if (clusterOn=="inputs"){out <- NULL}else{out <- nBag$outputs_raw}
                      
                      }else{
                        if (clusterOn=="outputs"){inp <- NULL}else{inp <- nBag$inputs}
                        if (clusterOn=="inputs"){out <- NULL}else{out <- nBag$outputs}
                      }
  connectivityCluster(inp,out,ROIs,grouping,distance,knownStats)
}

#' Test if x is a connectivityCluster object
#' @param x An object to be tested
#' @return TRUE if x is a neuronBag
#'@export
is.connectivityCluster <- function(x) inherits(x,"connectivityCluster")

#' Plot a connectivityCluster object
#' 
#' @param connCl A \code{\link{connectivityCluster}} object
#' @param orderX Optional ordering of the groupings (by default the HClust order)
#' @param orderY Optional ordering of the groupings (by default the HClust order)
#' @param colorAxis Whether or not to color the axis labels per cluster (FALSE by default)
#' @param axesPalette If \code{colorAxis} is TRUE, the palette to use
#' @param theme A ggplot theme to use for the plot
#' @param replaceIdWithNames In case of neuron to neuron clustering, whether or not to replace the bodyids by neuron names on the axis labels
#' @param h To be passed to \code{\link{cutree}} if \code{colorAxis} is TRUE. Either \code{h} or \code{k} must be set in that case. 
#' @param k To be passed to \code{\link{cutree}} if \code{colorAxis} is TRUE. Either \code{h} or \code{k} must be set in that case. 
#' @export
plotClusters <- function(connCl,
                         orderX=connCl$hc$order,
                         orderY=connCl$hc$order,
                         colorAxis=FALSE,
                         axesPalette=paletteer::paletteer_d("Polychrome::palette36")[3:36],
                         theme=theme_minimal(),
                         replaceIdWithNames=TRUE,
                         h=NULL,
                         k=NULL){
    ddM <- as.matrix(connCl$distance)
    
    ddM <- ddM[orderX,orderY]
    
    distDF <- reshape2::melt(ddM,as.is=TRUE) %>% mutate(Var1=factor(Var1,levels=rownames(ddM)),Var2=factor(Var2,levels=rownames(ddM)))
    if (colorAxis){
      clusters <- cutree(connCl$hc,h=h,k=k)
      connCols <- axesPalette[clusters[rownames(ddM)]]
      names(connCols) <- rownames(ddM)
    }
    
    p <- ggplot(distDF) + geom_tile(aes(x=Var1,y=Var2,fill=value)) + coord_fixed() +
      theme + scale_fill_distiller(palette = "Greys",name = "Distance") +
      theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust=0.5))+xlab("")+ylab("")
    if (colorAxis){p <- p +
      theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust=0.5,colour = connCols),axis.text.y = element_text(colour = connCols))
    }
    if (replaceIdWithNames & connCl$grouping=="neuron"){
      if (is.null(connCl$inputsTable)){newNames <- connCl$outputsTable$name.from[match(rownames(ddM),connCl$outputsTable$from)]}else{
        newNames <- connCl$inputsTable$name.to[match(rownames(ddM),connCl$inputsTable$to)]
      }
      p <- p + scale_x_discrete(labels=newNames)+scale_y_discrete(labels=newNames)
    }
    p
} 

#plotClusterChain : comparing clustering of neurons connected to each other

plotClusterChain <- function(connCluster1,
                             connCluster2,
                             chainConnectivity){
  
}