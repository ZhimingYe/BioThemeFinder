DupClustering.Fuzzy<-function(x,k,...){
  require(fclust)
  res.fcm <- fclust::Fclust(x@DupMatrix, k=k)
  kk<-res.fcm[["clus"]][,-2]%>%as.data.frame()
  colnames(kk)[1]<-"Cluster"
  return(kk)
}
DupClustering.NMF<-function(x,k,method="brunet",...){
  require(NMF)
  if(nrow(x@DupMatrix)>=1000){
    message("Too many pathways! NMF can be very slow!\n")
  }
  if(nrow(x@DupMatrix)>2500){
    stop("Please reduce the number of result before NMF.\n")
  }
  NMFresult<-NMF::nmf(x@DupMatrix,rank = 10,method = method,seed=123)
  NMFgroup <-NMF::predict(NMFresult)
  res0<-as.data.frame(NMFgroup)
  colnames(res0)[1]<-"Cluster"
  return(res0)
}
DupClustering.Hclust<-function(x,k,method1="euclidean",method2="ward.D2",...){
  hc<-hclust(dist(x@DupMatrix,method = method1),method =method2)
  clusters<-cutree(hc,k=k)
  res0<-as.data.frame(clusters)
  colnames(res0)[1]<-"Cluster"
  return(res0)
}
DupClustering.PAM<-function(x,k,method1="euclidean",...){
  PAMresult<-pam(x=x@DupMatrix,k=k,diss=F,stand = F,metric=method1)
  res0<-as.data.frame(PAMresult[["clustering"]])
  colnames(res0)[1]<-"Cluster"
  return(res0)
}
#' @rdname MatrixClustering
#' @title Cluster repetition rate matrix based on matrix-based methods
#'
#' @param x a BioThemeFinder object
#' @param k number of clusters. You can draw PathwayHeatmap to visualize preliminarily and determine the number of clusters
#' @param method can be "nmf","hc","pam","fuzzy"
#' @param dist_method default:euclidean
#' @param cluster_method default:NMF~brunet, hc~ward.D2
#' @param self_def whether re-define the dist_method or cluster_method or not
#' @return a BioThemeFinder object
#' @export
#' @author Zhiming Ye
#'
#' @examples
MatrixClustering<-function(x,k,method="nmf",dist_method,cluster_method,self_def=F){
  method<-match.arg(method,c("nmf","hc","pam","fuzzy"))
  if("GSCluster"%in%colnames(x@Results)){
    x@Results<-x@Results%>%dplyr::select(-GSCluster)
    message("Old Cluster has been removed.\n")
  }
  if(method=="nmf"){
    if(!self_def){
      cluster_method<-"brunet"
    }
    clsResult<-DupClustering.NMF(x=x,k=k,method=cluster_method)
  }
  if(method=="hc"){
    if(!self_def){
      cluster_method<-"ward.D2"
      dist_method<-"euclidean"
    }
    clsResult<-DupClustering.Hclust(x=x,k=k,method1=dist_method,method2=cluster_method)
  }
  if(method=="pam"){
    if(!self_def){
      dist_method<-"euclidean"
    }
    clsResult<-DupClustering.PAM(x=x,k=k,method1=dist_method)
  }
  if(method=="fuzzy"){
    clsResult<-DupClustering.Fuzzy(x=x,k=k)
  }
  # return(clsResult)
  clsResult<-clsResult%>%dplyr::mutate(Method=method)
  if(sum((colnames(x@DupMatrix)==rownames(x@DupMatrix))==0)!=0){
    warning("Unknown ERROR, colnames can't map rownames! condiction=1\n")
  }
  if(nrow(clsResult)!=nrow(x@Results)){
    warning("Unknown ERROR, colnames can't map rownames! condiction=2\n")
  }
  clsResult<-mapTBs.(clsResult,x@Results$Description)
  x@Results<-x@Results%>%dplyr::mutate(GSCluster=clsResult$Cluster)
  return(x)
}

