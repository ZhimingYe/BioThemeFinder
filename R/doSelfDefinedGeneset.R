
setGeneric("doSelfDefinedGeneset",function(x,PValCutOff,Term2GENE,...) standardGeneric("doSelfDefinedGeneset"))

setMethod("doSelfDefinedGeneset",signature(x="BioThemeFinder.ORA"),function(x,PValCutOff,Term2GENE,QValCutOff,...){
  if(x@Species=="human"){
    require(org.Hs.eg.db)
    OrgDB<-org.Hs.eg.db
  }
  if(x@Species=="mouse"){
    require(org.Mm.eg.db)
    OrgDB<-org.Mm.eg.db
  }
  if(sum(unlist(Term2GENE[,2])%in%x@WholeGenes)<20){
    stop("Please re-check your TERM2GENE table! \n")
  }
  if(x@Species%in%c("human","mouse")){
    SDres<-clusterProfiler::enricher(x@WholeGenes,TERM2GENE=Term2GENE,pvalueCutoff = PValCutOff,qvalueCutoff = QValCutOff)%>%DOSE::setReadable(OrgDb = OrgDB,keyType = "ENTREZID")%>%as.data.frame()%>%dplyr::mutate(Database="Self")
  }
  else{
    SDres<-clusterProfiler::enricher(x@WholeGenes,TERM2GENE=Term2GENE,pvalueCutoff = PValCutOff,qvalueCutoff = QValCutOff)%>%as.data.frame()%>%dplyr::mutate(Database="Self")
  }
  RESdf<-SDres
  return(RESdf)
})


setMethod("doSelfDefinedGeneset",signature(x="BioThemeFinder.ORA_FC"),function(x,PValCutOff,Term2GENE,QValCutOff,...){
  if(x@Species=="human"){
    require(org.Hs.eg.db)
    OrgDB<-org.Hs.eg.db
  }
  if(x@Species=="mouse"){
    require(org.Mm.eg.db)
    OrgDB<-org.Mm.eg.db
  }
  if(sum(unlist(Term2GENE[,2])%in%x@WholeGenes)<20){
    stop("Please re-check your TERM2GENE table! \n")
  }
  if(x@Species%in%c("human","mouse")){
    SDres<-clusterProfiler::enricher(x@WholeGenes,TERM2GENE=Term2GENE,pvalueCutoff = PValCutOff,qvalueCutoff = QValCutOff)%>%DOSE::setReadable(OrgDb = OrgDB,keyType = "ENTREZID")%>%as.data.frame()%>%dplyr::mutate(Database="Self")
  }
  else{
    SDres<-clusterProfiler::enricher(x@WholeGenes,TERM2GENE=Term2GENE,pvalueCutoff = PValCutOff,qvalueCutoff = QValCutOff)%>%as.data.frame()%>%dplyr::mutate(Database="Self")
  }

  ResultsDF<-SDres
  GeneRegType<-rep("UnKnown",nrow(ResultsDF))
  geneIDcol<-which(colnames(ResultsDF)=="geneID")
  if(x@Species%in%c("human","mouse")){
    for(i in 1:nrow(ResultsDF)){
      GeneRegType[i]<-DetermineDirection(ResultsDF[i,geneIDcol],NumOfDiff,x@UpRegGenes,x@DownRegGenes,OrgDB,T)
    }
  }
  else{
    for(i in 1:nrow(ResultsDF)){
      GeneRegType[i]<-DetermineDirection(ResultsDF[i,geneIDcol],NumOfDiff,x@UpRegGenes,x@DownRegGenes,needConvert=F)
    }
  }
  ResultsDF<-ResultsDF%>%dplyr::mutate(RegType=GeneRegType)
  RESdf<-ResultsDF
  return(RESdf)
})

setMethod("doSelfDefinedGeneset",signature(x="BioThemeFinder.GSEA"),function(x,PValCutOff,Term2GENE,...){
  if(x@Species=="human"){
    require(org.Hs.eg.db)
    OrgDB<-org.Hs.eg.db
  }
  if(x@Species=="mouse"){
    require(org.Mm.eg.db)
    OrgDB<-org.Mm.eg.db
  }
  if(sum(unlist(Term2GENE[,2])%in%names(x@RankedGenes))<20){
    stop("Please re-check your TERM2GENE table! \n")
  }
  if(x@Species%in%c("human","mouse")){
    SDres<-clusterProfiler::GSEA(x@RankedGenes,TERM2GENE=Term2GENE,pvalueCutoff = PValCutOff)%>%DOSE::setReadable(OrgDb = OrgDB,keyType = "ENTREZID")%>%as.data.frame()%>%dplyr::mutate(Database="Self")
  }
  else{
    SDres<-clusterProfiler::GSEA(x@RankedGenes,TERM2GENE=Term2GENE,pvalueCutoff = PValCutOff)%>%as.data.frame()%>%dplyr::mutate(Database="Self")
  }
  ResultsDF<-SDres
  ResultsDF<-ResultsDF%>%dplyr::mutate(RegType=case_when(NES>1~"Favor_UpReg",NES<(-1)~"Favor_DnReg",T~"UnKnown"))
  RESdf<-ResultsDF
  return(RESdf)
})
