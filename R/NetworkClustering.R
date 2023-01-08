#' @rdname NetworkClustering
#' @title Cluster repetition rate matrix based on network-based greedy optimization of modularity method
#' @param x a BioThemeFinder object
#' @param EdgeCutoff cutoff of repetition rate
#' @param ...
#'
#' @return BioThemeFinder object
#' @export
#' @author Zhiming Ye
#' @examples
NetworkClustering<-function(x,EdgeCutoff=0.5,...){
  if("NetworkCluster"%in%colnames(x@Results)){
    x@Results<-x@Results%>%dplyr::select(-NetworkCluster)
    DFempty<-data.frame(from=c("1","2","3"),to=c("3","2","1"))
    x@Network<-igraph::graph.data.frame(DFempty, directed = FALSE)
    message("Old Cluster has been removed.\n")
  }
  OriginDF<-x@DupMatrix
  colnames(OriginDF)<-paste0(x@Results$Database,":",x@Results$Description)
  rownames(OriginDF)<-paste0(x@Results$Database,":",x@Results$Description)
  ind <- which(upper.tri(OriginDF, diag = TRUE), arr.ind = TRUE)
  NwkDF<-cbind(ind, OriginDF[ind])
  colnames(NwkDF)<-c("from","to","DupRate")
  NwkDF<-as.data.frame(NwkDF)
  NwkDF$from<-rownames(OriginDF)[NwkDF$from]
  NwkDF$to<-colnames(OriginDF)[NwkDF$to]
  NwkDF<-NwkDF%>%dplyr::filter(DupRate>EdgeCutoff)
  NwkDF<-NwkDF%>%dplyr::filter(from!=to)
  x@Network<-igraph::graph.data.frame(NwkDF, directed = FALSE)
  x@Communities<-igraph::cluster_fast_greedy(x@Network)#greedy optimization of modularity
  FavorsUpPathwaysdf<-x@Results[x@Results$RegType=="Favor_UpReg",]
  FavorsUpPathways<-paste0(FavorsUpPathwaysdf$Database,":",FavorsUpPathwaysdf$Description)
  FavorsDnPathwaysdf<-x@Results[x@Results$RegType=="Favor_DnReg",]
  FavorsDnPathways<-paste0(FavorsDnPathwaysdf$Database,":",FavorsDnPathwaysdf$Description)
  V(x@Network)$color<-case_when(V(x@Network)$name%in%FavorsUpPathways~"#e5989b",V(x@Network)$name%in%FavorsDnPathways~"#023e8a",T~"#dad7cd")
  V(x@Network)$ggcolor<-case_when(V(x@Network)$name%in%FavorsUpPathways~"UP",V(x@Network)$name%in%FavorsDnPathways~"DOWN",T~"UNKNOWN")
  ResDF<-data.frame(Description2=x@Communities$names,NetworkCluster=x@Communities$membership)
  suppressMessages(x@Results<-x@Results%>%dplyr::mutate(Description2=paste0(x@Results$Database,":",x@Results$Description))%>%left_join(ResDF)%>%dplyr::select(-Description2))
  #NetworkCls.Argument
  x@IsNetworkClustered<-T
  return(x)
}

#' @rdname PlotNetwork
#' @title plot network
#' @param x a BioThemeFinder object
#' @param method igraph or ggplot2
#' @param Label
#' @param ...
#'
#' @return
#' @export
#' @author Zhiming Ye
#' @examples
PlotNetwork<-function(x,method="igraph",Label=T,...){
  if(x@IsNetworkClustered==T){
    message("For better viewing, please zoom the network, or output it as a PDF and zoom in to find clusters")
    if(method=="igraph"&Label==T){
      message("if can NOT display TEXT, please stop and try [ggplot2] method.\n")
      plot(x@Communities,x@Network,edge.color="#457b9d",edge.width=5,e=TRUE,vertex.size=6,vertex.label.color="#1d3557",vertex.label.cex=0.8,vertex.label.font=2,col=V(x@Network)$color,...)
    }
    if(method=="igraph"&Label==F){
      message("if can not display texts, please stop and try [ggplot2] method.\n")
      plot(x@Communities,x@Network,edge.color="#457b9d",edge.width=5,e=TRUE,vertex.size=6,vertex.label=NA,vertex.label.color="#1d3557",vertex.label.cex=0.8,vertex.label.font=2,col=V(x@Network)$color,...)
    }
    if(method=="ggplot2"&Label==T){
      xr<-x
      V(xr@Network)$ggcommunity<-xr@Communities$membership
      p<-ggraph(xr@Network,layout = "fr")+
        geom_edge_arc(strength = 0.2, width = 0.5, alpha = 0.8)+
        geom_node_point(aes(color = factor(ggcommunity),shape=factor(ggcolor,levels = c("UNKNOWN","UP","DOWN"))),size=10) +
        geom_node_text(aes(label = name), repel = F,alpha = 0.8) +theme_void()+labs(color="Communitiy Cluster",shape="UP/DOWN regulated")
      return(p)
    }
    if(method=="ggplot2"&Label==F){
      xr<-x
      V(xr@Network)$ggcommunity<-xr@Communities$membership
      p<-ggraph(xr@Network,layout = "fr")+
        geom_edge_arc(strength = 0.2, width = 0.5, alpha = 0.8)+
        geom_node_point(aes(color = factor(ggcommunity),shape=factor(ggcolor,levels = c("UNKNOWN","UP","DOWN"))),size=10)+theme_void()+labs(color="Communitiy Cluster",shape="UP/DOWN regulated")
      return(p)
    }
  }
  else{
    warning("Can't find network clustered result! Please do cluster first.\n")
  }
}
