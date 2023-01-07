
setGeneric("doGOanalysis",function(x,PValCutOff,...) standardGeneric("doGOanalysis"))

setMethod("doGOanalysis",signature(x="BioThemeFinder.ORA"),function(x,PValCutOff,QValCutOff,simplifycutoff,...){
  if(x@Specics=="human"){
    require(org.Hs.eg.db)
    OrgDB = org.Hs.eg.db
  }
  if(x@Specics=="mouse"){
    require(org.Mm.eg.db)
    OrgDB = org.Mm.eg.db
  }
  GOBP<-clusterProfiler::enrichGO(x@WholeGenes,OrgDb = OrgDB,keyType = "ENTREZID",ont = "BP",pvalueCutoff = PValCutOff,qvalueCutoff = QValCutOff)%>%clusterProfiler::simplify(cutoff=simplifycutoff)%>%DOSE::setReadable(OrgDb = OrgDB,keyType = "ENTREZID")%>%as.data.frame()%>%dplyr::mutate(Database="GOBP")
  GOCC<-clusterProfiler::enrichGO(x@WholeGenes,OrgDb = OrgDB,keyType = "ENTREZID",ont = "CC",pvalueCutoff = PValCutOff,qvalueCutoff = QValCutOff)%>%clusterProfiler::simplify(cutoff=simplifycutoff)%>%DOSE::setReadable(OrgDb = OrgDB,keyType = "ENTREZID")%>%as.data.frame()%>%dplyr::mutate(Database="GOCC")
  GOMF<-clusterProfiler::enrichGO(x@WholeGenes,OrgDb = OrgDB,keyType = "ENTREZID",ont = "MF",pvalueCutoff = PValCutOff,qvalueCutoff = QValCutOff)%>%clusterProfiler::simplify(cutoff=simplifycutoff)%>%DOSE::setReadable(OrgDb = OrgDB,keyType = "ENTREZID")%>%as.data.frame()%>%dplyr::mutate(Database="GOMF")
  RESdf<-rbind(GOBP,GOCC,GOMF)
  return(RESdf)
})


setMethod("doGOanalysis",signature(x="BioThemeFinder.ORA_FC"),function(x,PValCutOff,QValCutOff,simplifycutoff,...){
  NumOfDiff<-x@CutOff_Reg
  if(x@Specics=="human"){
    require(org.Hs.eg.db)
    OrgDB = org.Hs.eg.db
  }
  if(x@Specics=="mouse"){
    require(org.Mm.eg.db)
    OrgDB = org.Mm.eg.db
  }
  GOBP<-clusterProfiler::enrichGO(x@WholeGenes,OrgDb = OrgDB,keyType = "ENTREZID",ont = "BP",pvalueCutoff = PValCutOff,qvalueCutoff = QValCutOff)%>%clusterProfiler::simplify(cutoff=simplifycutoff)%>%DOSE::setReadable(OrgDb = OrgDB,keyType = "ENTREZID")%>%as.data.frame()%>%dplyr::mutate(Database="GOBP")
  GOCC<-clusterProfiler::enrichGO(x@WholeGenes,OrgDb = OrgDB,keyType = "ENTREZID",ont = "CC",pvalueCutoff = PValCutOff,qvalueCutoff = QValCutOff)%>%clusterProfiler::simplify(cutoff=simplifycutoff)%>%DOSE::setReadable(OrgDb = OrgDB,keyType = "ENTREZID")%>%as.data.frame()%>%dplyr::mutate(Database="GOCC")
  GOMF<-clusterProfiler::enrichGO(x@WholeGenes,OrgDb = OrgDB,keyType = "ENTREZID",ont = "MF",pvalueCutoff = PValCutOff,qvalueCutoff = QValCutOff)%>%clusterProfiler::simplify(cutoff=simplifycutoff)%>%DOSE::setReadable(OrgDb = OrgDB,keyType = "ENTREZID")%>%as.data.frame()%>%dplyr::mutate(Database="GOMF")
  ResultsDF<-rbind(GOBP,GOCC,GOMF)
  GeneRegType<-rep("UnKnown",nrow(ResultsDF))
  geneIDcol<-which(colnames(ResultsDF)=="geneID")
  for(i in 1:nrow(ResultsDF)){
    GeneRegType[i]<-DetermineDirection(ResultsDF[i,geneIDcol],NumOfDiff,x@UpRegGenes,x@DownRegGenes)
  }
  ResultsDF<-ResultsDF%>%dplyr::mutate(RegType=GeneRegType)
  RESdf<-ResultsDF
  return(RESdf)
})

setMethod("doGOanalysis",signature(x="BioThemeFinder.GSEA"),function(x,PValCutOff,...){
  if(x@Specics=="human"){
    require(org.Hs.eg.db)
    OrgDB = org.Hs.eg.db
  }
  if(x@Specics=="mouse"){
    require(org.Mm.eg.db)
    OrgDB = org.Mm.eg.db
  }
  GOBP<-clusterProfiler::gseGO(x@RankedGenes,OrgDb = OrgDB,keyType = "ENTREZID",ont = "BP",pvalueCutoff = PValCutOff)%>%DOSE::setReadable(OrgDb = OrgDB,keyType = "ENTREZID")%>%as.data.frame()%>%dplyr::mutate(Database="GOBP")
  GOCC<-clusterProfiler::gseGO(x@RankedGenes,OrgDb = OrgDB,keyType = "ENTREZID",ont = "CC",pvalueCutoff = PValCutOff)%>%DOSE::setReadable(OrgDb = OrgDB,keyType = "ENTREZID")%>%as.data.frame()%>%dplyr::mutate(Database="GOCC")
  GOMF<-clusterProfiler::gseGO(x@RankedGenes,OrgDb = OrgDB,keyType = "ENTREZID",ont = "MF",pvalueCutoff = PValCutOff)%>%DOSE::setReadable(OrgDb = OrgDB,keyType = "ENTREZID")%>%as.data.frame()%>%dplyr::mutate(Database="GOMF")
  ResultsDF<-rbind(GOBP,GOCC,GOMF)
  ResultsDF<-ResultsDF%>%dplyr::mutate(RegType=ifelse(.$NES>0,"Favor_UpReg","Favor_DnReg"))
  RESdf<-ResultsDF
  return(RESdf)
})
