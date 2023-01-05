
setClass("BioThemeFinder.ORA",slots=list(WholeGenes="character",
                                         Specics="character",
                                         Results="data.frame",
                                         dbName="character",
                                         IsAnalysed="logical",
                                         IsClustered="logical",
                                         DupMatrix="data.frame"))

setClass("BioThemeFinder.ORA_FC",slots=list(UpRegGenes="character",
                                            DownRegGenes="character",
                                            WholeGenes="character",
                                            Specics="character",
                                            Results="data.frame",
                                            dbName="character",
                                            IsAnalysed="logical",
                                            IsClustered="logical",
                                            DupMatrix="data.frame",
                                            CutOff_Reg="numeric"))

setClass("BioThemeFinder.GSEA",slots=list(RankedGenes="numeric",
                                          Specics="character",
                                          Results="data.frame",
                                          dbName="character",
                                          IsAnalysed="logical",
                                          IsClustered="logical",
                                          DupMatrix="data.frame"))

Create.newBioThemeFinder.ORA<-function (Gene,FromType = "SYMBOL", Specics="human"){
  if(Specics=="human"){
    OrgDB = org.Hs.eg.db
  }
  if(Specics=="mouse"){
    OrgDB = org.Mm.eg.db
  }
  if(!Specics%in%c("human","mouse")){
    warning("unsupported specics.\n")
  }
  Genetable <- data.frame(Gene=Gene)
  ENTREZIDtable <- clusterProfiler::bitr(Genetable$Gene, fromType = FromType,toType = "ENTREZID", OrgDb = OrgDB)
  GSresult<-new("BioThemeFinder.ORA")
  GSresult@WholeGenes<-ENTREZIDtable$ENTREZID
  GSresult@Specics<-Specics
  GSresult@IsAnalysed<-F
  GSresult@IsClustered<-F
  return(GSresult)
}

Create.newBioThemeFinder.ORAwithFC<-function (Gene, log2FC,Pvalue,FCcutoff=1,PvalueCutOff=0.05,FromType = "SYMBOL", Specics="human",CntsOfDiffGene=2){
  if(Specics=="human"){
    OrgDB = org.Hs.eg.db
  }
  if(Specics=="mouse"){
    OrgDB = org.Mm.eg.db
  }
  if(!Specics%in%c("human","mouse")){
    warning("unsupported specics.\n")
  }
  Genetable <- data.frame(Gene=Gene,log2FC=log2FC,PVal=Pvalue)%>%dplyr::filter(Pvalue<PvalueCutOff)%>%dplyr::filter(abs(log2FC)>FCcutoff)
  ENTREZIDtable <- clusterProfiler::bitr(Genetable$Gene, fromType = FromType,toType = "ENTREZID", OrgDb = OrgDB)
  colnames(ENTREZIDtable)[1] <- "Gene"
  Genetable <- Genetable %>% left_join(ENTREZIDtable) %>% arrange(desc(log2FC))%>%na.omit()
  GSElist <- as.numeric(Genetable$log2FC)
  names(GSElist) <- Genetable$ENTREZID
  GSElist = sort(GSElist, decreasing = TRUE)
  GSresult<-new("BioThemeFinder.ORA_FC")
  GSresult@UpRegGenes<-names(GSElist[GSElist>0])
  GSresult@DownRegGenes<-names(GSElist[GSElist<0])
  GSresult@WholeGenes<-names(GSElist)
  GSresult@Specics<-Specics
  GSresult@IsAnalysed<-F
  GSresult@IsClustered<-F
  GSresult@CutOff_Reg<-CntsOfDiffGene
  return(GSresult)
}

Create.newBioThemeFinder.GSEA<-function (Gene, log2FC,FromType = "SYMBOL", Specics="human"){
  if(Specics=="human"){
    OrgDB = org.Hs.eg.db
  }
  if(Specics=="mouse"){
    OrgDB = org.Mm.eg.db
  }
  if(!Specics%in%c("human","mouse")){
    warning("unsupported specics.\n")
  }
  Genetable <- data.frame(Gene=Gene,log2FC=log2FC)
  ENTREZIDtable <- clusterProfiler::bitr(Genetable$Gene, fromType = FromType,toType = "ENTREZID", OrgDb = OrgDB)
  colnames(ENTREZIDtable)[1] <- "Gene"
  Genetable <- Genetable %>% left_join(ENTREZIDtable) %>% arrange(desc(log2FC))%>%na.omit()
  GSElist <- as.numeric(Genetable$log2FC)
  names(GSElist) <- Genetable$ENTREZID
  GSElist = sort(GSElist, decreasing = TRUE)
  GSresult<-new("BioThemeFinder.GSEA")
  GSresult@RankedGenes<-GSElist
  GSresult@Specics<-Specics
  GSresult@IsAnalysed<-F
  GSresult@IsClustered<-F
  return(GSresult)
}
