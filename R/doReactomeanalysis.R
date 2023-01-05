
setGeneric("doReactomeanalysis",function(x,PValCutOff,...) standardGeneric("doReactomeanalysis"))

setMethod("doReactomeanalysis",signature(x="BioThemeFinder.ORA"),function(x,PValCutOff,QValCutOff,...){
  if(x@Specics=="human"){
    orgid<-"human"
  }
  if(x@Specics=="mouse"){
    orgid<-"mouse"
  }
  RAres<-ReactomePA::enrichPathway(x@WholeGenes,organism=orgid,pvalueCutoff = PValCutOff,qvalueCutoff = QValCutOff)%>%as.data.frame()%>%dplyr::mutate(Database="Reactome")
  RESdf<-RAres
  return(RESdf)
})


setMethod("doReactomeanalysis",signature(x="BioThemeFinder.ORA_FC"),function(x,PValCutOff,QValCutOff,...){
  NumOfDiff<-x@CutOff_Reg
  if(x@Specics=="human"){
    orgid<-"human"
  }
  if(x@Specics=="mouse"){
    orgid<-"mouse"
  }
  RAres<-ReactomePA::enrichPathway(x@WholeGenes,organism=orgid,pvalueCutoff = PValCutOff,qvalueCutoff = QValCutOff)%>%as.data.frame()%>%dplyr::mutate(Database="Reactome")
  ResultsDF<-RAres
  GeneRegType<-rep("UnKnown",nrow(ResultsDF))
  geneIDcol<-which(colnames(ResultsDF)=="geneID")
  for(i in 1:nrow(ResultsDF)){
    GeneRegType[i]<-DetermineDirection(ResultsDF[i,geneIDcol],NumOfDiff,x@UpRegGenes,x@DownRegGenes)
  }
  ResultsDF<-ResultsDF%>%dplyr::mutate(RegType=GeneRegType)
  RESdf<-ResultsDF
  return(RESdf)
})

setMethod("doReactomeanalysis",signature(x="BioThemeFinder.GSEA"),function(x,PValCutOff,...){
  if(x@Specics=="human"){
    orgid<-"human"
  }
  if(x@Specics=="mouse"){
    orgid<-"mouse"
  }
  RAres<-ReactomePA::gsePathway(x@RankedGenes,organism=orgid,pvalueCutoff = PValCutOff)%>%as.data.frame()%>%dplyr::mutate(Database="Reactome")
  ResultsDF<-RAres
  ResultsDF<-ResultsDF%>%dplyr::mutate(RegType=ifelse(.$NES>0,"Favor_UpReg","Favor_DnReg"))
  RESdf<-ResultsDF
  return(RESdf)
})