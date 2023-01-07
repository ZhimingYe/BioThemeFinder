
setGeneric("doSelfDefinedGeneset",function(x,PValCutOff,Term2GENE,...) standardGeneric("doSelfDefinedGeneset"))

setMethod("doSelfDefinedGeneset",signature(x="BioThemeFinder.ORA"),function(x,PValCutOff,Term2GENE,QValCutOff,...){
  if(sum(unlist(Term2GENE[,2])%in%x@WholeGenes)<20){
    stop("Please re-check your TERM2GENE table! \n")
  }
  SDres<-clusterProfiler::enricher(x@WholeGenes,TERM2GENE=Term2GENE,pvalueCutoff = PValCutOff,qvalueCutoff = QValCutOff)%>%as.data.frame()%>%dplyr::mutate(Database="Self")
  RESdf<-SDres
  return(RESdf)
})


setMethod("doSelfDefinedGeneset",signature(x="BioThemeFinder.ORA_FC"),function(x,PValCutOff,Term2GENE,QValCutOff,...){
  if(sum(unlist(Term2GENE[,2])%in%x@WholeGenes)<20){
    stop("Please re-check your TERM2GENE table! \n")
  }
  SDres<-clusterProfiler::enricher(x@WholeGenes,TERM2GENE=Term2GENE,pvalueCutoff = PValCutOff,qvalueCutoff = QValCutOff)%>%as.data.frame()%>%dplyr::mutate(Database="Self")
  ResultsDF<-SDres
  GeneRegType<-rep("UnKnown",nrow(ResultsDF))
  geneIDcol<-which(colnames(ResultsDF)=="geneID")
  for(i in 1:nrow(ResultsDF)){
    GeneRegType[i]<-DetermineDirection(ResultsDF[i,geneIDcol],NumOfDiff,x@UpRegGenes,x@DownRegGenes)
  }
  ResultsDF<-ResultsDF%>%dplyr::mutate(RegType=GeneRegType)
  RESdf<-ResultsDF
  return(RESdf)
})

setMethod("doSelfDefinedGeneset",signature(x="BioThemeFinder.GSEA"),function(x,PValCutOff,Term2GENE,...){
  if(sum(unlist(Term2GENE[,2])%in%names(x@RankedGenes))<20){
    stop("Please re-check your TERM2GENE table! \n")
  }
  SDres<-clusterProfiler::enricher(x@RankedGenes,TERM2GENE=Term2GENE,pvalueCutoff = PValCutOff)%>%as.data.frame()%>%dplyr::mutate(Database="Self")
  ResultsDF<-SDres
  ResultsDF<-ResultsDF%>%dplyr::mutate(RegType=ifelse(.$NES>0,"Favor_UpReg","Favor_DnReg"))
  RESdf<-ResultsDF
  return(RESdf)
})