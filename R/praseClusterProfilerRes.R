#' @rdname prase.clusterProfiler.result
#' @title prase result from clusterProfiler
#'
#' @param x A BioThemeFinder object
#' @param GOBP clusterProfiler object
#' @param GOCC clusterProfiler object
#' @param GOMF clusterProfiler object
#' @param GOALL clusterProfiler object
#' @param KEGG clusterProfiler object
#' @param MKEGG clusterProfiler object
#' @param Reactome clusterProfiler object
#' @param Self self-defined gene collection
#' @param nGeneCutOff Enriched pathways with a gene counts less than that value were filtered
#' @param IsSymbol True~SYMBOL (have done "setreadable" function), not~ENTREZID. NOTE:It is recommended to perform setreadable operation before importing BioThemeFinder
#' @param OrgDB the org.DB using to convert ID
#'
#' @return A BioThemeFinder object
#' @export
#' @author Zhiming Ye
#'
#' @examples
prase.clusterProfiler.result<-function(x,GOBP=NULL,GOCC=NULL,GOMF=NULL,GOALL=NULL,KEGG=NULL,MKEGG=NULL,Reactome=NULL,Self=NULL,nGeneCutOff=3,IsSymbol=F,OrgDB=NULL){
  if(class(x)[1]%in%c("BioThemeFinder.ORA_FC","BioThemeFinder.ORA")){
    Type<-"ORA"
  }
  else{
    Type<-"GSEA"
  }
  if(x@Species=="non-gene"){
    message("Please only input in Self because of non-gene input.\n")
  }
  Collection<-c()
  Result00<-data.frame()
  for(i in c("GOBP","GOCC","GOMF","GOALL","KEGG","MKEGG","Reactome","Self")){
    if(!is.null(get(i))){
      Collection<-c(Collection,i)
      if(i=="MKEGG"){
        dbName0<-"KEGG"
      }
      else{
        dbName0<-i
      }
      cpResult<-(get(i))%>%as.data.frame()%>%dplyr::mutate(Database=dbName0)
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
      if((x@Species!="non-gene")|(x@Species=="non-gene"&i=="Self")){
        Result00<-rbind(Result00,cpResult)
      }
    }
  }
  if(nrow(Result00)>0){
    if(sum(duplicated(Result00$Description))>0){
      cat("Duplicated Terms:        will only keep the first\n")
      print(Result00$Description[duplicated(Result00$Description)])
      Result00<-Result00%>%dplyr::filter(!duplicated(Result00$Description))
    }
    if(nrow(Result00)<=2){
      stop("Not enough enrichment result, please consider a lower P or Q value!\nOr wrong gene set. \n")
    }
    if(class(x)[1]%in%c("BioThemeFinder.ORA_FC")){
      geneIDcol<-which(colnames(Result00)=="geneID")
      for(i in 1:nrow(Result00)){
        GeneRegType[i]<-DetermineDirection(Result00[i,geneIDcol],NumOfDiff,x@UpRegGenes,x@DownRegGenes,OrgDB,needConvert=IsSymbol)
      }
      Result00<-Result00%>%dplyr::mutate(RegType=GeneRegType)
    }
    if(class(x)[1]%in%c("BioThemeFinder.GSEA")){
      Result00<-Result00%>%dplyr::mutate(RegType=ifelse(.$NES>0,"Favor_UpReg","Favor_DnReg"))
    }
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


