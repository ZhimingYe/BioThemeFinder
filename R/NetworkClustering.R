ToNetworkDF<-function(x,EdgeCutoff=0.5,nGeneCutOff=3,...){

  ind <- which(upper.tri(df, diag = TRUE), arr.ind = TRUE)
  NwkDF<-cbind(ind, df[ind])
  colnames(NwkDF)<-c("from","to","Prob")
  NwkDF$from<-rownames(x@dup)[NwkDF$from]
  NwkDF$from<-colnames(x@dup)[NwkDF$to]
  NwkDF<-NwkDF%>%dplyr::filter(Prob>EdgeCutoff)
  return(NwkDF)
}
a<-c(1,2,3,4,5)
b<-c(0,5,6,7,8)
d<-c(0,0,4,6,75)
e<-c(0,0,0,43,53)
f<-c(0,0,0,0,13)
df<-rbind(a,b,d,e,f)
ind <- which(upper.tri(df, diag = TRUE), arr.ind = TRUE)
cbind(ind, df[ind])

