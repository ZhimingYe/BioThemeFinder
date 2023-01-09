setOldClass("igraph")
setOldClass("communities")
setClass("BioThemeFinder.ORA",slots=list(WholeGenes="character",
                                         Species="character",
                                         Results="data.frame",
                                         dbName="character",
                                         IsAnalysed="logical",
                                         IsClustered="logical",
                                         IsNetworkClustered="logical",
                                         DupMatrix="data.frame",
                                         Network="igraph",
                                         Communities="communities",
                                         SelectedResultNames="character"))

setClass("BioThemeFinder.ORA_FC",slots=list(UpRegGenes="character",
                                            DownRegGenes="character",
                                            WholeGenes="character",
                                            Species="character",
                                            Results="data.frame",
                                            dbName="character",
                                            IsAnalysed="logical",
                                            IsClustered="logical",
                                            IsNetworkClustered="logical",
                                            DupMatrix="data.frame",
                                            Network="igraph",
                                            Communities="communities",
                                            CutOff_Reg="numeric",
                                            SelectedResultNames="character"))

setClass("BioThemeFinder.GSEA",slots=list(RankedGenes="numeric",
                                          Species="character",
                                          Results="data.frame",
                                          dbName="character",
                                          IsAnalysed="logical",
                                          IsClustered="logical",
                                          IsNetworkClustered="logical",
                                          DupMatrix="data.frame",
                                          Network="igraph",
                                          Communities="communities",
                                          SelectedResultNames="character"))
#' @rdname Create.newBioThemeFinder.ORA
#' @title When only gene names are available, create BioThemeFinder object using ORA.
#'
#' @param Gene A vector containing genes
#' @param FromType Type of genes, can be the columns supported by org.DBs like SYMBOL for gene symbol or ENSEMBL for Ensembl ID
#' @param Species one of the human, mouse of non-gene
#'
#' @return A BioThemeFinder.ORA object
#' @export
#' @author Zhiming Ye
#'
#' @examples
Create.newBioThemeFinder.ORA<-function (Gene,FromType = "SYMBOL", Species="human"){
  require(clusterProfiler)
  require(ReactomePA)
  Species <- match.arg(Species, c("human", "mouse","non-gene"))
  if(Species=="human"){
    require(org.Hs.eg.db)
    OrgDB = org.Hs.eg.db
  }
  if(Species=="mouse"){
    require(org.Mm.eg.db)
    OrgDB = org.Mm.eg.db
  }
  Genetable <- data.frame(Gene=Gene)
  if(Species!="non-gene"){
    ENTREZIDtable <- clusterProfiler::bitr(Genetable$Gene, fromType = FromType,toType = "ENTREZID", OrgDb = OrgDB)
  }
  else{
    ENTREZIDtable<-data.frame(ENTREZID=Genetable$Gene)
  }
  GSresult<-new("BioThemeFinder.ORA")
  GSresult@WholeGenes<-ENTREZIDtable$ENTREZID
  GSresult@Species<-Species
  GSresult@IsAnalysed<-F
  GSresult@IsClustered<-F
  return(GSresult)
}
#' @rdname Create.newBioThemeFinder.ORAwithFC
#' @title When gene names and vector measuring the chenge of specific genes, like logFC, are available, create BioThemeFinder object with FC using ORA.
#'
#' @param Gene A vector containing genes
#' @param log2FC A vector containing genes
#' @param Pvalue A vector containing pvalue, the order of the above three must be the same
#' @param FCcutoff the cut off of abs(FC) value, default as 1
#' @param PvalueCutOff the cut off of P value, default as 0.05
#' @param FromType Type of genes, can be the columns supported by org.DBs like SYMBOL for gene symbol or ENSEMBL for Ensembl ID
#' @param Species one of the human, mouse of non-gene
#' @param CntsOfDiffGene Threshold for judging the proportion of downregulated genes on the same pathway. For example, when set this value equals to 2, in a particular gene set, there are at least two up-regulated genes more than down regulated genes to judge that there is an up trend for that gene set.
#'
#' @return A BioThemeFinder.ORA_FC object
#' @export
#' @author Zhiming Ye
#'
#' @examples
Create.newBioThemeFinder.ORAwithFC<-function (Gene, log2FC,Pvalue,FCcutoff=1,PvalueCutOff=0.05,FromType = "SYMBOL", Species="human",CntsOfDiffGene=2){
  require(clusterProfiler)
  require(ReactomePA)
  Species <- match.arg(Species, c("human", "mouse","non-gene"))
  if(Species=="human"){
    require(org.Hs.eg.db)
    OrgDB = org.Hs.eg.db
  }
  if(Species=="mouse"){
    require(org.Mm.eg.db)
    OrgDB = org.Mm.eg.db
  }
  if(!is.null(PvalueCutOff)){
    Genetable <- data.frame(Gene=Gene,log2FC=log2FC,PVal=Pvalue)%>%dplyr::filter(Pvalue<PvalueCutOff)%>%dplyr::filter(abs(log2FC)>FCcutoff)
  }
  else{
    Genetable <- data.frame(Gene=Gene,log2FC=log2FC)%>%dplyr::filter(abs(log2FC)>FCcutoff)
  }
  if(Species!="non-gene"){
    ENTREZIDtable <- clusterProfiler::bitr(Genetable$Gene, fromType = FromType,toType = "ENTREZID", OrgDb = OrgDB)
  }
  else{
    ENTREZIDtable<-data.frame(Gene=Genetable$Gene,ENTREZID=Genetable$Gene)
  }
  colnames(ENTREZIDtable)[1] <- "Gene"
  Genetable <- Genetable %>% left_join(ENTREZIDtable) %>% arrange(desc(log2FC))%>%na.omit()
  GSElist <- as.numeric(Genetable$log2FC)
  names(GSElist) <- Genetable$ENTREZID
  GSElist = sort(GSElist, decreasing = TRUE)
  GSresult<-new("BioThemeFinder.ORA_FC")
  GSresult@UpRegGenes<-names(GSElist[GSElist>0])
  GSresult@DownRegGenes<-names(GSElist[GSElist<0])
  GSresult@WholeGenes<-names(GSElist)
  GSresult@Species<-Species
  GSresult@IsAnalysed<-F
  GSresult@IsClustered<-F
  GSresult@CutOff_Reg<-CntsOfDiffGene
  return(GSresult)
}
#' @rdname Create.newBioThemeFinder.GSEA
#' @title When gene names and vector measuring the chenge of specific genes, like logFC, are available, create BioThemeFinder object with FC using ORA.
#'
#' @param Gene A vector containing genes
#' @param log2FC A vector containing genes, the order of the above two must be the same
#' @param FromType Type of genes, can be the columns supported by org.DBs like SYMBOL for gene symbol or ENSEMBL for Ensembl ID
#' @param Species one of the human, mouse of non-gene
#'
#' @return A BioThemeFinder.GSEA object
#' @export
#' @author Zhiming Ye
#'
#' @examples
Create.newBioThemeFinder.GSEA<-function (Gene, log2FC,FromType = "SYMBOL", Species="human"){
  require(clusterProfiler)
  require(ReactomePA)
  Species <- match.arg(Species, c("human", "mouse","non-gene"))
  if(Species=="human"){
    require(org.Hs.eg.db)
    OrgDB = org.Hs.eg.db
  }
  if(Species=="mouse"){
    require(org.Mm.eg.db)
    OrgDB = org.Mm.eg.db
  }
  Genetable <- data.frame(Gene=Gene,log2FC=log2FC)
  if(Species!="non-gene"){
    ENTREZIDtable <- clusterProfiler::bitr(Genetable$Gene, fromType = FromType,toType = "ENTREZID", OrgDb = OrgDB)
  }
  else{
    ENTREZIDtable<-data.frame(Gene=Genetable$Gene,ENTREZID=Genetable$Gene)
  }
  colnames(ENTREZIDtable)[1] <- "Gene"
  Genetable <- Genetable %>% left_join(ENTREZIDtable) %>% arrange(desc(log2FC))%>%na.omit()
  GSElist <- as.numeric(Genetable$log2FC)
  names(GSElist) <- Genetable$ENTREZID
  GSElist = sort(GSElist, decreasing = TRUE)
  GSresult<-new("BioThemeFinder.GSEA")
  GSresult@RankedGenes<-GSElist
  GSresult@Species<-Species
  GSresult@IsAnalysed<-F
  GSresult@IsClustered<-F
  return(GSresult)
}
