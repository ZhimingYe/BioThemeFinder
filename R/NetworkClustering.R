ToNetworkDF<-function(x,EdgeCutoff=0.5,...){
  OriginDF<-x@DupMatrix
  #name加入DB
  ind <- which(upper.tri(OriginDF, diag = TRUE), arr.ind = TRUE)
  NwkDF<-cbind(ind, OriginDF[ind])
  colnames(NwkDF)<-c("from","to","Prob")
  NwkDF<-as.data.frame(NwkDF)
  NwkDF$from<-rownames(x@DupMatrix)[NwkDF$from]
  NwkDF$to<-colnames(x@DupMatrix)[NwkDF$to]
  NwkDF<-NwkDF%>%dplyr::filter(Prob>EdgeCutoff)
  NwkDF<-NwkDF%>%dplyr::filter(from!=to)
  g <- graph.data.frame(NwkDF, directed = FALSE)
  return(g)
}


NwkOrigin<-NwkOrigin[,-3]
g <- graph.data.frame(NwkOrigin, directed = FALSE)
