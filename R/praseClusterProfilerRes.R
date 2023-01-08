prase.clusterProfiler.result<-function(x,GOBP=NULL,GOCC=NULL,GOMF=NULL,GOALL=NULL,KEGG=NULL,MKEGG=NULL,Reactome=NULL,Self=NULL,nGeneCutOff=3){
  if(class(x)[1]%in%c("BioThemeFinder.ORA_FC","BioThemeFinder.ORA")){
    Type<-"ORA"
  }
  else{
    Type<-"GSEA"
  }
  Collection<-c()
  Result00<-data.frame()
  for(i in c("GOBP","GOCC","GOMF","GOALL","KEGG","MKEGG","Reactome","Self")){
    if(!is.null(get(i))){
      Collection<-c(Collection,i)
      cpResult<-(get(i))%>%as.data.frame()%>%dplyr::mutate(Database=i)
      if(Type=="ORA"){
        if(!"geneID"%in%colnames(cpResult)){
          stop(paste0(i," result is NOT ORA result!\n"))
        }
      }
      if(Type=="GSEA"){
        if(!"core_enrichment"%in%colnames(cpResult)){
          stop(paste0(i," result is NOT GSEA result!\n"))
        }
      }
      Result00<-rbind(Result00,cpResult)
    }
  }
  if(nrow(Result00)>0){
    x@Results<-Result00
    NumFiltered<-sum((getNumOfGENEs(x))<=nGeneCutOff)
    x@Results<-x@Results%>%dplyr::filter((getNumOfGENEs(x))>nGeneCutOff)
    cat("Included result from clusterProfiler : ",Collection,"\nFiltered ",NumFiltered," results with less than ",nGeneCutOff," genes enriched.\n")
    x@IsAnalysed<-T
    x@IsClustered<-F
    x@IsNetworkClustered<-F
    x@dbName<-Collection
    message("Finished!\n")
    return(x)
  }
  else{
    stop("Can't find any pathways. Please re-check.\n")
  }
}


