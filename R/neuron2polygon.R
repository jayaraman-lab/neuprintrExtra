#' Transforms a nat neuron object into a polygon that can be plotted in 2D
#' @param neur A neuron or neuronList object
#' @param axis Which axis do we want to project onto. Has to be a length 2 vector of "x", "y" or "z"
#' @param progress When applied to neuronLists, whether to display a progress bar
#' @param cl A cluster parameter, can be either an integer or a cluster object (see \code{pbapply::pbapply})
#' @return a dataframe with columns x, y, proj and bodyid. If \code{neur} is  a neuronlist there will be several bodyids.
#' @export
neuron2polygon <- function(neur,axis=c("x","y"),progress=TRUE,cl=10L){
  UseMethod("neuron2polygon")}

#' @export
neuron2polygon.neuron <- function(neur,axis=c("x","y")){
  refT <- mutate(neur$d,x=!!as.name(toupper(axis[1])),y=!!as.name(toupper(axis[2])))
  refT <- addOutlines(refT)
  roots <- nat::rootpoints(neur)
  ptsOrd <- tryCatch(unlist(lapply(roots,treeOrder,neur$SegList,branchpoints(neur),refTable = refT)),error=function(cond){segOrder(neur,toupper(axis[2]))})
  refT <- left_join(ptsOrd,refT,by = c("idx" = "PointNo"))
  refT <- mutate(refT,x=ifelse(side=="upper",upperX,lowerX),y=ifelse(side=="upper",upperY,lowerY)) %>% mutate(bodyid=neur$bodyid,proj=paste0(axis,collapse=""))
  refT %>% select(x,y,bodyid,proj)
}

#' @export
neuron2polygon.neuronlist <- function(neurL,axis=c("x","y"),progress=TRUE,cl=10L){
  if (progress){
      nD <- pbapply::pblapply(neurL,neuron2polygon,axis=axis,cl=cl)
    }else{
      nD <- lapply(neurL,neuron2polygon,axis=axis)
    }

  outD <- bind_rows(nD)
}

angle <- function(x, y) {
  atan2(y[2] - y[1], x[2] - x[1])
}

addOutlines <- function(neuronD){
  out <- neuronD %>%
    mutate(xPar = ifelse(Parent %in% PointNo,neuronD$x[match(Parent,neuronD$PointNo)],neuronD$x[match(PointNo,neuronD$Parent)]),
                  yPar = ifelse(Parent %in% PointNo,neuronD$y[match(Parent,neuronD$PointNo)],neuronD$y[match(PointNo,neuronD$Parent)]),
                  angle = atan2(y-yPar,x-xPar),
                  dx = W/2 * (cos(angle + pi/2)),
                  dy = W/2 * (sin(angle + pi/2)),
                  lowerX = xPar - dx,
                  lowerY = yPar - dy,
                  upperX = xPar +dx,
                  upperY =yPar + dy)
}

treeOrder <- function(subTreeIdx,segmentList,branchPts,refTable){
  subTree <- segmentList[[subTreeIdx]]
  segmentList <- segmentList[-subTreeIdx]
  if (tail(subTree,1) %in% branchPts){
    bP <- tail(subTree,1)
    childrenIdx <- which(sapply(segmentList,function(l) l[1]==bP))
    childrenTrees <- segmentList[childrenIdx]
    childrenSec <- sapply(childrenTrees,function(s) s[2])
    childrenTab <- arrange(filter(refTable,PointNo %in% childrenSec),desc(upperY))$PointNo
    childrenTrees <- childrenIdx[match(childrenTab,childrenSec)]
    return(bind_rows(list(data.frame(idx=subTree,side="upper",stringsAsFactors = FALSE),
                          bind_rows(lapply(childrenTrees,treeOrder,segmentList,branchPts,refTable)),
                          data.frame(idx=rev(subTree),side="lower",stringsAsFactors = FALSE))))

  }else{
    return(rbind(data.frame(idx=subTree,side="upper",stringsAsFactors = FALSE),data.frame(idx=rev(subTree),side="lower",stringsAsFactors = FALSE)))
  }
}

segOrder <- function(neuronOb,ax="Y"){
  nGraph <- nat::segmentgraph(neuronOb,weights=FALSE,segids=TRUE)
  igraph::E(nGraph)$weight <- neuronOb$d[[ax]][match(igraph::E(nGraph)$segid,neuronOb$d$PointNo)]
  igraph::V(nGraph)$name <-  igraph::V(nGraph)$label
  tree <- igraph::as_data_frame(nGraph)

  roots <- (filter(tree,nat::rootpoints(neuronOb)==from) %>% arrange(desc(weight)))$segid
  baseTree <- unlist(lapply(roots,iterativePreorder,tree))
  startVertices <- match(tree$to,igraph::V(nGraph)$name)

  allSubTrees <- lapply(startVertices,function(r){na.omit(igraph::dfs(nGraph,r,unreachable=FALSE)$order)$name})
  allSubTrees<- lapply(allSubTrees,function(s) filter(tree,to %in% s)$segid)
  allSubTrees <- lapply(1:length(baseTree),function(i) baseTree[baseTree %in% c(i,allSubTrees[[i]])])

  for (i in (baseTree)){
    baseTree <- append(baseTree,-i,after=tail(which(baseTree == tail(allSubTrees[[i]],1),1)))
  }
  idx <- integer()
  side <- character()
  for (s in baseTree){
    if (s>0){
      idx <- append(idx,neuronOb$SegList[[s]])
      side <- append(side,rep("upper",length(neuronOb$SegList[[s]])))
    }else{
      idx <- append(idx,rev(neuronOb$SegList[[-s]]))
      side <- append(side,rep("lower",length(neuronOb$SegList[[-s]])))
    }

  }
  data.frame(idx=idx,side=side)
}

# all argument to return all the subtrees?
iterativePreorder <- function(node,tree){
  if (is.null(node)) return(NULL)
  todo <- c()
  orderedTree <- c()
  todo <- append(todo,node)
  while (length(todo)>0){
    node <- todo[1]
    todo <- todo[-1]

    orderedTree <- append(orderedTree,node)

    children <- filter(tree,from == tree$to[tree$segid == node] & !(segid %in% c(todo,orderedTree))) %>% arrange(desc(weight)) ## Avoid circular stuff
    if (length(children)>0){
      todo <- append(todo,children$segid,0)
    }
  }
  return(orderedTree)
}
