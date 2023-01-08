
setGeneric("doKEGGanalysis",function(x,PValCutOff,MKEGG,...) standardGeneric("doKEGGanalysis"))

setMethod("doKEGGanalysis",signature(x="BioThemeFinder.ORA"),function(x,PValCutOff,MKEGG,QValCutOff,...){
  if(x@Specics=="human"){
    orgid<-"hsa"
    require(org.Hs.eg.db)
    OrgDB<-org.Hs.eg.db
  }
  if(x@Specics=="mouse"){
    orgid<-"mmu"
    require(org.Mm.eg.db)
    OrgDB<-org.Mm.eg.db
  }
  KEGGres<-clusterProfiler::enrichKEGG(x@WholeGenes,organism=orgid,pvalueCutoff = PValCutOff,qvalueCutoff = QValCutOff)%>%DOSE::setReadable(OrgDb = OrgDB,keyType = "ENTREZID")%>%as.data.frame()%>%dplyr::mutate(Database="KEGG")
  KEGGres2<-data.frame()
  if(MKEGG){
    KEGGres2<-clusterProfiler::enrichMKEGG(x@WholeGenes,organism=orgid,pvalueCutoff = PValCutOff,qvalueCutoff = QValCutOff)%>%DOSE::setReadable(OrgDb = OrgDB,keyType = "ENTREZID")%>%as.data.frame()%>%dplyr::mutate(Database="KEGG")
  }
  RESdf<-rbind(KEGGres,KEGGres2)
  return(RESdf)
})


setMethod("doKEGGanalysis",signature(x="BioThemeFinder.ORA_FC"),function(x,PValCutOff,MKEGG,QValCutOff,...){
  NumOfDiff<-x@CutOff_Reg
  if(x@Specics=="human"){
    orgid<-"hsa"
    require(org.Hs.eg.db)
    OrgDB<-org.Hs.eg.db
  }
  if(x@Specics=="mouse"){
    orgid<-"mmu"
    require(org.Mm.eg.db)
    OrgDB<-org.Mm.eg.db
  }
  KEGGres<-clusterProfiler::enrichKEGG(x@WholeGenes,organism=orgid,pvalueCutoff = PValCutOff,qvalueCutoff = QValCutOff)%>%DOSE::setReadable(OrgDb = OrgDB,keyType = "ENTREZID")%>%as.data.frame()%>%dplyr::mutate(Database="KEGG")
  KEGGres2<-data.frame()
  if(MKEGG){
    KEGGres2<-clusterProfiler::enrichMKEGG(x@WholeGenes,organism=orgid,pvalueCutoff = PValCutOff,qvalueCutoff = QValCutOff)%>%DOSE::setReadable(OrgDb = OrgDB,keyType = "ENTREZID")%>%as.data.frame()%>%dplyr::mutate(Database="KEGG")
  }
  ResultsDF<-rbind(KEGGres,KEGGres2)
  GeneRegType<-rep("UnKnown",nrow(ResultsDF))
  geneIDcol<-which(colnames(ResultsDF)=="geneID")
  for(i in 1:nrow(ResultsDF)){
    GeneRegType[i]<-DetermineDirection(ResultsDF[i,geneIDcol],NumOfDiff,x@UpRegGenes,x@DownRegGenes,OrgDB,T)
  }
  ResultsDF<-ResultsDF%>%dplyr::mutate(RegType=GeneRegType)
  RESdf<-ResultsDF
  return(RESdf)
})

setMethod("doKEGGanalysis",signature(x="BioThemeFinder.GSEA"),function(x,PValCutOff,MKEGG,...){
  if(x@Specics=="human"){
    orgid<-"hsa"
    require(org.Hs.eg.db)
    OrgDB<-org.Hs.eg.db
  }
  if(x@Specics=="mouse"){
    orgid<-"mmu"
    require(org.Mm.eg.db)
    OrgDB<-org.Mm.eg.db
  }
  KEGGres<-clusterProfiler::gseKEGG(x@RankedGenes,organism=orgid,pvalueCutoff = PValCutOff)%>%DOSE::setReadable(OrgDb = OrgDB,keyType = "ENTREZID")%>%as.data.frame()%>%dplyr::mutate(Database="KEGG")
  KEGGres2<-data.frame()
  if(MKEGG){
    KEGGres2<-clusterProfiler::gseMKEGG(x@RankedGenes,organism=orgid,pvalueCutoff = PValCutOff)%>%DOSE::setReadable(OrgDb = OrgDB,keyType = "ENTREZID")%>%as.data.frame()%>%dplyr::mutate(Database="KEGG")
  }
  ResultsDF<-rbind(KEGGres,KEGGres2)
  ResultsDF<-ResultsDF%>%dplyr::mutate(RegType=ifelse(.$NES>0,"Favor_UpReg","Favor_DnReg"))
  RESdf<-ResultsDF
  return(RESdf)
})
