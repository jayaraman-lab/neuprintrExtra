#' Clustering neurons based on connectivity tables
#' 
#' @param inputsTable The table of input connections to the neurons to cluster
#' @param outputsTable The table of input connections to the neurons to cluster
#' @param ROIs Which ROIs to consider. If NULL (the default) uses all the ROI present
#' @param unit Unit on which the clustering is to be made (neuron or type, or various supertypes)
#' @param distance Which distance measure to use (defaults to \code{\link{cos_dist}})
#' @param knownStats Whether or not to use \code{knownOutputContribution} and \code{knownWeightRelative} instead of the standard measures
#' @details \itemize{
#'  \item If both \code{inputsTable} and \code{outputsTable} are present, the clustering will be run on the combined connectivity vector
#'  \item \code{weightRelative} is used for working on input vectors, and \code{outputContribution} to work on output vectors
#'  }
#' @return an object of class \code{connectivityCluster}. The object is a list with components:
#' \describe{
#'   \item{inputsTable}{The table of input connections to the neurons being clustered (can be NULL if clustered on outputs)}
#'   \item{outputsTable}{The table of output connections to the neurons being clustered (can be NULL if clustered on inputs)}
#'   \item{distance}{object containing a dissimilarity measure between the neurons being considered}
#'   \item{hc}{\code{\link{hclust}} object, result of clustering}
#'   \item{unit}{Unit on which the clustering was made (neuron or type, or various supertypes)}
#' }
#' @seealso \code{\link{clusterBag}}
#' @export
connectivityCluster <- function(inputsTable=NULL,
                                outputsTable=NULL,
                                ROIs=NULL,
                                unit=c("type","neuron","supertype1","supertype2","supertype3","databaseType"),
                                distance=cos_dist,
                                knownStats=FALSE){
  unit <-match.arg(unit)
  unitOld <- unit
  if (unit=="neuron") {unit <- ""}else{
    unit=paste0(unit,".")}
  from <- paste0(unit,"from")
  to <- paste0(unit,"to")
  
  if (!is.null(inputsTable)){
    stat <- ifelse(knownStats,"knownWeightRelative","weightRelative")
    inputsMat <- connectivityMatrix(inputsTable,slctROIs = ROIs,allToAll = FALSE,from = from,to=to,value = stat,ref="outputs")
    connMat <- inputsMat
  }
  
  if (!is.null(outputsTable)){
    stat <- ifelse(knownStats,"knownOutputContribution","outputContribution")
    outputsMat <- connectivityMatrix(outputsTable,slctROIs = ROIs,allToAll = FALSE,from = from,to=to,value = stat,ref="inputs")
    connMat <- outputsMat
  }
  
  if (!is.null(inputsTable) && !is.null(outputsTable)){
    outputsMat <- outputsMat[rownames(inputsMat),]
    colnames(outputsMat) <- paste0(colnames(outputsMat),".out")
    colnames(inputsMat) <- paste0(colnames(inputsMat),".in")
    connMat <- cbind(inputsMat,outputsMat)
  }
  
  connDist <- distance(connMat)
  connHC <- hclust(connDist)
  res <- list(inputsTable=inputsTable,
              outputsTable=outputsTable,
              names=rownames(connMat),
              distance=connDist,
              hc=connHC,
              unit=unitOld)
  attr(res,"class") <- "connectivityCluster"
  return(res)
}

#' Cluster neurons in a bag 
#' 
#' @param nBag a neuron bag
#' @param clusterOn whether to cluster on the input vector, the output vector, or both
#' @inheritParams connectivityCluster
#' @return an object of class \code{connectivityCluster}. The object is a list with components:
#' \describe{
#'   \item{inputsTable}{The table of input connections to the neurons being clustered (can be NULL if clustered on outputs)}
#'   \item{outputsTable}{The table of output connections to the neurons being clustered (can be NULL if clustered on inputs)}
#'   \item{distance}{object containing a dissimilarity measure between the neurons being considered}
#'   \item{hc}{\code{\link{hclust}} object, result of clustering}
#'   \item{unit}{Unit on which the clustering was made (neuron or type, or various supertypes)}
#' }
#' @seealso the lower level \code{\link{connectivityCluster}}
#' @export
clusterBag <- function(nBag,
                       ROIs=NULL,
                       clusterOn = c("inputs","outputs","both"),
                       unit=c("type","neuron","supertype1","supertype2","supertype3","databaseType"),
                       distance=cos_dist,
                       knownStats=FALSE){
  unit <- match.arg(unit)
  clusterOn <- match.arg(clusterOn)
  if (unit=="neuron"){if (clusterOn=="outputs"){inp <- NULL}else{inp <- nBag$inputs_raw}
                      if (clusterOn=="inputs"){out <- NULL}else{out <- nBag$outputs_raw}
                      
                      }else{
                        if (clusterOn=="outputs"){inp <- NULL}else{inp <- nBag$inputs}
                        if (clusterOn=="inputs"){out <- NULL}else{out <- nBag$outputs}
                      }
  connectivityCluster(inp,out,ROIs,unit,distance,knownStats)
}

#' Test if x is a connectivityCluster object
#' @param x An object to be tested
#' @return TRUE if x is a neuronBag
#'@export
is.connectivityCluster <- function(x) inherits(x,"connectivityCluster")

#' Plot a connectivityCluster object
#' 
#' @param connCl A \code{\link{connectivityCluster}} object
#' @param order Optional ordering of the units (by default the HClust order)
#' @param colorAxis Whether or not to color the axis labels per cluster (FALSE by default)
#' @param axesPalette If \code{colorAxis} is TRUE, the palette to use
#' @param theming A ggplot theme to use for the plot
#' @param replaceIdWithNames In case of neuron to neuron clustering, whether or not to replace the bodyids by neuron names on the axis labels
#' @param h To be passed to \code{\link{cutree}} if \code{colorAxis} is TRUE. Either \code{h} or \code{k} must be set in that case. 
#' @param k To be passed to \code{\link{cutree}} if \code{colorAxis} is TRUE. Either \code{h} or \code{k} must be set in that case. 
#' @export
plotClusters <- function(connCl,
                         order=NULL,
                         colorAxis=FALSE,
                         axesPalette=paletteer::paletteer_d("Polychrome::palette36")[3:36],
                         theming=theme_minimal,
                         replaceIdWithNames=TRUE,
                         h=NULL,
                         k=NULL){
    ddM <- as.matrix(connCl$distance)
    
    if (is.null(order)){
        ddM <- ddM[connCl$hc$order,connCl$hc$order]
    }else{
      ddM <- ddM[order,order]
    }
    
    
    distDF <- reshape2::melt(ddM,as.is=TRUE) %>% mutate(cluster1=connCl$clusters[Var1],cluster2=connCl$clusters[Var2],Var1=factor(Var1,levels=rownames(ddM)),Var2=factor(Var2,levels=rownames(ddM)))
    if (colorAxis){
      clusters <- cutree(connCl$hc,h=h,k=k)
      connCols <- axesPalette[clusters[rownames(ddM)]]
      names(connCols) <- rownames(ddM)
    }
    
    p <- ggplot(distDF) + geom_tile(aes(x=Var1,y=Var2,fill=value)) + coord_fixed() +
      theming() + scale_fill_distiller(palette = "Greys",name = "Distance") +
      theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust=0.5))+xlab("")+ylab("")
    if (colorAxis){p <- p +
      theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust=0.5,colour = connCols),axis.text.y = element_text(colour = connCols))
    }
    if (replaceIdWithNames & connCl$unit=="neuron"){
      if (is.null(connCl$inputsTable)){newNames <- connCl$outputsTable$name.from[match(rownames(ddM),connCl$outputsTable$from)]}else{
        newNames <- connCl$inputsTable$name.to[match(rownames(ddM),connCl$inputsTable$to)]
      }
      p <- p + scale_x_discrete(labels=newNames)+scale_y_discrete(labels=newNames)
    }
    p
} 
