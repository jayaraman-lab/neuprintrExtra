#' Builds a connectivity matrix from a connection object
#' @param connObj  a connection table
#' @param slctROIs  which ROIs to consider. If NULL (the default), all the ROIs present in the table are used (columns of the matrix are then appended with .roi for each ROI)
#' @param allToAll  whether to build a square matrix or just a from -> to matrix
#' @param from  which field to use as a "source" (default "name.from")
#' @param to  which field to use as a "target" (default "name.to")
#' @param value  which field to use to fill the matrix (default "weightRelative")
#' @param ref  which channel will be used as the "reference" (to be the columns of the output). The
#' other channel gets .roi affixed to their names in case several ROIs are considered
#'
#' @export
connectivityMatrix <- function(connObj,
                               slctROIs=NULL,
                               allToAll=FALSE,
                               from="name.from",
                               to="name.to",
                               value="weightRelative",
                               ref=c("inputs","outputs")){

  ref <- match.arg(ref)
  if (is.null(slctROIs)){slctROIs <- unique(connObj$roi)}
  connObj <- filter(connObj,roi %in% slctROIs)
  if (any(is.na(c(connObj[[to]],connObj[[from]])))){
    warning("Some to/from entries are NA, using retype.na function.")
    connObj <- retype.na(connObj)
  }

  connObj[[to]] <- as.character(connObj[[to]])
  connObj[[from]] <- as.character(connObj[[from]])
  if (nrow(distinct_at(connObj,c(from,to,"roi"))) != nrow(connObj)){
    stop("Multiple entries for some of the from/to/roi combinations. You need to either
         use different from/to or summarize your data.frame beforehand.")}

  if (allToAll){
    bare <- unique(c(connObj[[from]],connObj[[to]]))
    if (length(slctROIs)==1){
      ins <- bare
      outs <- bare
    }else{
      if (ref=="inputs"){
        ins <- unique(c(paste(connObj[[from]],connObj[["roi"]],sep="."),c(paste(connObj[[to]],connObj[["roi"]],sep="."))))
        outs <- bare
      }else{
        ins <- bare
        outs <- unique(c(paste(connObj[[from]],connObj[["roi"]],sep="."),c(paste(connObj[[to]],connObj[["roi"]],sep="."))))
      }
    }
  }else{
    ins <- unique(connObj[[from]])
    outs <- unique(connObj[[to]])
  }

  if (length(slctROIs)>1){
    if (ref=="inputs"){
      ins <- unique(paste(connObj[[from]],connObj[["roi"]],sep="."))
      outMat <- matrix(0,nrow=length(ins),ncol=length(outs),dimnames=list("Inputs"=ins,"Outputs"=outs))

      for (l in 1:nrow(connObj)){
        outMat[paste0(connObj[[from]][l],".",connObj[["roi"]][l]),connObj[[to]][l]] <- connObj[[value]][l]
      }
    }else{
      outs <- unique(paste(connObj[[to]],connObj[["roi"]],sep="."))
      outMat <- matrix(0,nrow=length(ins),ncol=length(outs),dimnames=list("Inputs"=ins,"Outputs"=outs))

      for (l in 1:nrow(connObj)){
        outMat[connObj[[from]][l],paste0(connObj[[to]][l],".",connObj[["roi"]][l])] <- connObj[[value]][l]
      }
      outMat <- t(outMat)
    }
  }else{
    outMat <- matrix(0,nrow=length(ins),ncol=length(outs),dimnames=list("Inputs"=ins,"Outputs"=outs))
    for (l in 1:nrow(connObj)){
      outMat[connObj[[from]][l],connObj[[to]][l]] <- connObj[[value]][l]
    }
    if (ref=="outputs") outMat <-  t(outMat)
  }
  outMat
}

#'Distance measurements
#'
#'@param mat A matrix
#'@return A distance object containing distances between the
#'rows of \code{mat}
#'@details \code{cos_dist} returns the cosine distance, \code{sqrt_cos_dist}
#'the squared cosine distance, \code{cor_dist} one minus the spearman correlation 
#'between vectors, and \code{bin_dist} the binary distance after thresholding
#'@export
cos_dist <- function(mat){
  sim <- mat / sqrt(rowSums(mat * mat))
  sim <- sim %*% t(sim)
  D_sim <- as.dist(1 - sim)
  attr(D_sim,"method") <- "cosine"
  D_sim
}

#' @describeIn cos_dist Squared cosine distance
#' @export
sqrt_cos_dist <- function(mat){
  sim <- sqrt(mat) / sqrt(rowSums(mat))
  sim <- sim %*% t(sim)
  D_sim <- as.dist(1 - sim)
  attr(D_sim,"method") <- "sqcosine"
  D_sim
}

#' @describeIn cos_dist Correlation distance matrix
#' @export
cor_dist <- function(mat){
  connCor <- cor(t(mat),method="spearman")
  D_sim <- as.dist((1-connCor)/2)
  attr(D_sim,"method") <- "correlation"
  D_sim
}

#' @describeIn cos_dist Binary distance
#' @param threshold Threshold to use to binarize the matrix
#' @export
bin_dist <- function(mat,threshold=0.01){
  dist(mat>threshold,method="binary")
}

#' Turn a connectivity matrix back into a data.frame
#conn_mat2df <- function(connMat,
#                        orderIn=NULL,
#                        orderOut=NULL,
#                        thresh=0){
#  if (is.null(orderIn)){orderIn=1:length(dimnames(connMat)$Inputs)}
#  if (is.null(orderOut)){orderOut=1:length(dimnames(connMat)$Outputs)}
#  connMat[connMat<=thresh] <- NA
#  reshape2::melt(connMat,na.rm=TRUE) %>% 
#    mutate(Inputs=factor(Inputs,levels=dimnames(connMat)$Inputs[orderIn]),
#           Outputs=factor(Outputs,levels=dimnames(connMat)$Outputs[orderOut]))
#}