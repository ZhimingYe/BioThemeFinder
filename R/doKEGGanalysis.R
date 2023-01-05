
setGeneric("doKEGGanalysis",function(x,PValCutOff,...) standardGeneric("doKEGGanalysis"))

setMethod("doKEGGanalysis",signature(x="BioThemeFinder.ORA"),function(x,PValCutOff,QValCutOff,...){
  if(x@Specics=="human"){
    orgid<-"hsa"
  }
  if(x@Specics=="mouse"){
    orgid<-"mmu"
  }
  KEGGres<-clusterProfiler::enrichKEGG(x@WholeGenes,organism=orgid,pvalueCutoff = PValCutOff,qvalueCutoff = QValCutOff)%>%as.data.frame()%>%dplyr::mutate(Database="KEGG")
  RESdf<-KEGGres
  return(RESdf)
})


setMethod("doKEGGanalysis",signature(x="BioThemeFinder.ORA_FC"),function(x,PValCutOff,QValCutOff,...){
  NumOfDiff<-x@CutOff_Reg
  if(x@Specics=="human"){
    orgid<-"hsa"
  }
  if(x@Specics=="mouse"){
    orgid<-"mmu"
  }
  KEGGres<-clusterProfiler::enrichKEGG(x@WholeGenes,organism=orgid,pvalueCutoff = PValCutOff,qvalueCutoff = QValCutOff)%>%as.data.frame()%>%dplyr::mutate(Database="KEGG")
  ResultsDF<-KEGGres
  GeneRegType<-rep("UnKnown",nrow(ResultsDF))
  geneIDcol<-which(colnames(ResultsDF)=="geneID")
  for(i in 1:nrow(ResultsDF)){
    GeneRegType[i]<-DetermineDirection(ResultsDF[i,geneIDcol],NumOfDiff,x@UpRegGenes,x@DownRegGenes)
  }
  ResultsDF<-ResultsDF%>%dplyr::mutate(RegType=GeneRegType)
  RESdf<-ResultsDF
  return(RESdf)
})

setMethod("doKEGGanalysis",signature(x="BioThemeFinder.GSEA"),function(x,PValCutOff,...){
  if(x@Specics=="human"){
    orgid<-"hsa"
  }
  if(x@Specics=="mouse"){
    orgid<-"mmu"
  }
  KEGGres<-clusterProfiler::gseKEGG(x@RankedGenes,organism=orgid,pvalueCutoff = PValCutOff)%>%as.data.frame()%>%dplyr::mutate(Database="KEGG")
  ResultsDF<-KEGGres
  ResultsDF<-ResultsDF%>%dplyr::mutate(RegType=ifelse(.$NES>0,"Favor_UpReg","Favor_DnReg"))
  RESdf<-ResultsDF
  return(RESdf)
})
