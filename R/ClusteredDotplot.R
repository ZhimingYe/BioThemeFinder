##' @rdname PathwayStatsPlot
##' @exportMethod PathwayStatsPlot
##' @author Zhiming Ye

setGeneric("PathwayStatsPlot",function(x,clusterType=NULL,showStart=NULL,showEnd=NULL,orderBy,...) standardGeneric("PathwayStatsPlot"))


#' @rdname PathwayStatsPlot
#' @title plot dotplot showing the statstical result of BioThemeFinder
#' @param x BioThemeFinder object
#' @param clusterType can be one of "MatrixResult", "NetworkResult"
#' @param showStart showing from <showStart> to <showEnd> in ordered table (ordering is based on the orderBy parameter)
#' @param showEnd showing from <showStart> to <showEnd> in ordered table (ordering is based on the orderBy parameter)
#' @param orderBy can be "pValue","GeneRatio","Counts"
#' @param col_low defining color
#' @param col_high defining color
#' @return a ggplot2 object, add + facet_grid(~Cluster) to display the clustered result. remove + facet_grid(~othergroup) for unclustered results
#' @export
#' @author Zhiming Ye
#' @examples
setMethod("PathwayStatsPlot",signature(x="BioThemeFinder.ORA"),function(x,clusterType="NetworkResult",showStart=1,showEnd=20,orderBy,col_low="#1e6091",col_high="#99d98c",...){
  require(forcats)
  if(showEnd>nrow(x@Results)){
    if(nrow(x@Results)>30){
      showEnd<-30
    }
    else{
      showEnd<-nrow(x@Results)
    }
  }
  clusterType <- match.arg(clusterType, c("MatrixResult", "NetworkResult"))
  orderBy<-match.arg(orderBy,c("pValue","GeneRatio","Counts"))
  if(sum(c("GSCluster","NetworkCluster")%in%colnames(x@Results))==0){
    FigDF<-data.frame(Terms=paste0(x@Results$Database,":",x@Results$Description),pValue=x@Results$p.adjust,Counts=x@Results$Count,GeneRatio=x@Results$GeneRatio,Cluster=NA)
  }
  if("GSCluster"%in%colnames(x@Results)&clusterType=="MatrixResult"){
    FigDF<-data.frame(Terms=paste0(x@Results$Database,":",x@Results$Description),pValue=x@Results$p.adjust,Counts=x@Results$Count,GeneRatio=x@Results$GeneRatio,Cluster=x@Results$GSCluster)%>%na.omit()
  }
  if("NetworkCluster"%in%colnames(x@Results)&clusterType=="NetworkResult"){
    FigDF<-data.frame(Terms=paste0(x@Results$Database,":",x@Results$Description),pValue=x@Results$p.adjust,Counts=x@Results$Count,GeneRatio=x@Results$GeneRatio,Cluster=x@Results$NetworkCluster)%>%na.omit()
  }
  if(nrow(FigDF)<2){
    stop("Enrichment result is too few!\n")
  }
  FigDF<-FigDF%>%dplyr::arrange(desc(!!sym(orderBy)))
  if(showEnd<showStart){
    stop("End is smaller than start!\n")
  }
  keep<-c(1:nrow(FigDF))[c(1:nrow(FigDF))%in%c(showStart:showEnd)]
  if(length(keep)==0){
    stop("NO overlap!\n")
  }
  p<-ggplot(FigDF[keep,],aes(x=parse_ratio_(GeneRatio), y=fct_reorder(Terms,!!sym(orderBy)), size=Counts, color=pValue))+geom_point()+scale_color_continuous(low=col_low, high=col_high, name = "adjustedP",guide=guide_colorbar(reverse=F))+scale_size(range=c(2, 8))+xlab("Gene Ratio")+ylab("Gene Set")+theme_bw()
  print(p)
  return(p)
})


#' @rdname PathwayStatsPlot
#' @param x BioThemeFinder object
#' @param clusterType can be one of "MatrixResult", "NetworkResult"
#' @param showStart showing from <showStart> to <showEnd> in ordered table (ordering is based on the orderBy parameter)
#' @param showEnd showing from <showStart> to <showEnd> in ordered table (ordering is based on the orderBy parameter)
#' @param orderBy can be "pValue","GeneRatio","Counts"
#' @param col_low defining color
#' @param col_high defining color
#' @export
#' @author Zhiming Ye
#' @examples
setMethod("PathwayStatsPlot",signature(x="BioThemeFinder.ORA_FC"),function(x,clusterType="NetworkResult",showStart=1,showEnd=20,orderBy,col_low="#1e6091",col_high="#99d98c",...){
  require(forcats)
  if(showEnd>nrow(x@Results)){
    if(nrow(x@Results)>30){
      showEnd<-30
    }
    else{
      showEnd<-nrow(x@Results)
    }
  }
  clusterType <- match.arg(clusterType, c("MatrixResult", "NetworkResult"))
  orderBy<-match.arg(orderBy,c("pValue","GeneRatio","Counts"))
  if(sum(c("GSCluster","NetworkCluster")%in%colnames(x@Results))==0){
    FigDF<-data.frame(Terms=paste0(x@Results$Database,":",x@Results$Description),pValue=x@Results$p.adjust,Counts=x@Results$Count,GeneRatio=x@Results$GeneRatio,Cluster=NA,Regulation=x@Results$RegType)
  }
  if("GSCluster"%in%colnames(x@Results)&clusterType=="MatrixResult"){
    FigDF<-data.frame(Terms=paste0(x@Results$Database,":",x@Results$Description),pValue=x@Results$p.adjust,Counts=x@Results$Count,GeneRatio=x@Results$GeneRatio,Cluster=x@Results$GSCluster,Regulation=x@Results$RegType)%>%na.omit()
  }
  if("NetworkCluster"%in%colnames(x@Results)&clusterType=="NetworkResult"){
    FigDF<-data.frame(Terms=paste0(x@Results$Database,":",x@Results$Description),pValue=x@Results$p.adjust,Counts=x@Results$Count,GeneRatio=x@Results$GeneRatio,Cluster=x@Results$NetworkCluster,Regulation=x@Results$RegType)%>%na.omit()
  }
  if(nrow(FigDF)<2){
    stop("Enrichment result is too few!\n")
  }
  FigDF<-FigDF%>%dplyr::arrange(desc(!!sym(orderBy)))
  if(showEnd<showStart){
    stop("End is smaller than start!\n")
  }
  keep<-c(1:nrow(FigDF))[c(1:nrow(FigDF))%in%c(showStart:showEnd)]
  if(length(keep)==0){
    stop("NO overlap!\n")
  }
  p<-ggplot(FigDF[keep,],aes(x=parse_ratio_(GeneRatio), y=fct_reorder(Terms,!!sym(orderBy)), size=Counts, color=pValue,shape=Regulation))+geom_point()+scale_color_continuous(low=col_low, high=col_high, name = "adjustedP",guide=guide_colorbar(reverse=F))+scale_size(range=c(2, 8))+xlab("Gene Ratio")+ylab("Gene Set")+theme_bw()
  print(p)
  return(p)
})



#' @rdname PathwayStatsPlot
#' @param x BioThemeFinder object
#' @param clusterType can be one of "MatrixResult", "NetworkResult"
#' @param showStart showing from <showStart> to <showEnd> in ordered table (ordering is based on the NES)
#' @param showEnd showing from <showStart> to <showEnd> in ordered table (ordering is based on the NES)
#' @param col_low defining color
#' @param col_high defining color
#' @export
#' @author Zhiming Ye
#' @examples
setMethod("PathwayStatsPlot",signature(x="BioThemeFinder.GSEA"),function(x,clusterType="NetworkResult",showStart=1,showEnd=20,col_low="#1e6091",col_high="#99d98c",...){
  require(forcats)
  if(showEnd>nrow(x@Results)){
    if(nrow(x@Results)>30){
      showEnd<-30
    }
    else{
      showEnd<-nrow(x@Results)
    }
  }
  clusterType <- match.arg(clusterType, c("MatrixResult", "NetworkResult"))
  orderBy<-"NES"
  if(sum(c("GSCluster","NetworkCluster")%in%colnames(x@Results))==0){
    FigDF<-data.frame(Terms=paste0(x@Results$Database,":",x@Results$Description),pValue=x@Results$p.adjust,NES=x@Results$NES,Cluster=NA)
  }
  if("GSCluster"%in%colnames(x@Results)&clusterType=="MatrixResult"){
    FigDF<-data.frame(Terms=paste0(x@Results$Database,":",x@Results$Description),pValue=x@Results$p.adjust,NES=x@Results$NES,Cluster=x@Results$GSCluster)%>%na.omit()
  }
  if("NetworkCluster"%in%colnames(x@Results)&clusterType=="NetworkResult"){
    FigDF<-data.frame(Terms=paste0(x@Results$Database,":",x@Results$Description),pValue=x@Results$p.adjust,NES=x@Results$NES,Cluster=x@Results$NetworkCluster)%>%na.omit()
  }
  if(nrow(FigDF)<2){
    stop("Enrichment result is too few!\n")
  }
  FigDF<-FigDF%>%dplyr::arrange(desc(!!sym(orderBy)))
  if(showEnd<showStart){
    stop("End is smaller than start!\n")
  }
  keep<-c(1:nrow(FigDF))[c(1:nrow(FigDF))%in%c(showStart:showEnd)]
  if(length(keep)==0){
    stop("NO overlap!\n")
  }
  p<-ggplot(FigDF[keep,],aes(x=NES, y=fct_reorder(Terms,NES), size=NES, color=pValue))+geom_point()+scale_color_continuous(low=col_low, high=col_high, name = "adjustedP",guide=guide_colorbar(reverse=F))+scale_size(range=c(2, 8))+xlab("Gene Ratio")+ylab("Gene Set")+theme_bw()
  print(p)
  return(p)
})



#parse_ratio is Cited from enrichplot package
#Copyright Guangchuang Yu Lab @ SMU
parse_ratio_<-function (ratio)
{
  gsize <- as.numeric(sub("/\\d+$", "", as.character(ratio)))
  gcsize <- as.numeric(sub("^\\d+/", "", as.character(ratio)))
  return(gsize/gcsize)
}
