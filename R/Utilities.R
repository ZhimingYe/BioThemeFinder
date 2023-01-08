DetermineDirection<-function(x,NumOfDiff,UPset,DNset){
  GeneList<-unlist(strsplit(x,"/"))
  UPsetNum<-sum(GeneList%in%UPset)
  DNsetNum<-sum(GeneList%in%DNset)
  Type0<-"UnKnown"
  if((UPsetNum>DNsetNum)&(abs(UPsetNum-DNsetNum))>=NumOfDiff){
    Type0<-"Favor_UpReg"
    return(Type0)
  }
  if((UPsetNum<DNsetNum)&(abs(UPsetNum-DNsetNum))>=NumOfDiff){
    Type0<-"Favor_DnReg"
    return(Type0)
  }
  return(Type0)
}

MultiDBanalysis0<-function(x,PVal = 0.05,QVal = 0.1,DBlist= c("GO","KEGG"),nGeneCutOff=3,simplify_cutoff=0.7,useMKEGG=T,Term2GENE){
  message("***BioThemeFinder***\nEnrichment analysis will be finished by clusterProfiler ",packageVersion('clusterProfiler')," and ReactomePA ",packageVersion('ReactomePA')," .\nPlease cite article The Innovation. 2021, 2(3):100141. doi: 10.1016/j.xinn.2021.100141 when using them. \n")
  if(x@Specics=="non-gene"){
    DBlist<-"SelfDefinedGS"
    message("No-Genes! please use self-defined collection to enrich\n")
  }
  Res0<-data.frame()
  if(sum(c("GO","KEGG","Reactome")%in%DBlist)<1){
    stop("Can't find any database!\n")
  }
  if("GO"%in%DBlist){
    cat("Applying GO analysis...\n")
    Res0<-rbind(Res0,doGOanalysis(x=x,PValCutOff = PVal,QValCutOff=QVal,simplifycutoff=simplify_cutoff))
  }
  if("KEGG"%in%DBlist){
    cat("Applying KEGG analysis...\nKEGG database is only available online for non-commercial use. \nSo there might be an error when unable to reach to the KEGG server. \n")
    Res0<-rbind(Res0,doKEGGanalysis(x=x,PValCutOff = PVal,QValCutOff=QVal,MKEGG=useMKEGG))
  }
  if("Reactome"%in%DBlist){
    cat("Applying Reactome analysis...\n")
    Res0<-rbind(Res0,doReactomeanalysis(x=x,PValCutOff = PVal,QValCutOff=QVal))
  }
  if("SelfDefinedGS"%in%DBlist){
    if(!is.null(Term2GENE)){
      cat("Applying analysis in self-defined gene set...\n")
      Res0<-rbind(Res0,doSelfDefinedGeneset(x=x,PValCutOff = PVal,QValCutOff=QVal,Term2GENE=Term2GENE))
    }
    else{
      stop("Please provide self-defined gene set in Term2GENE argument. \n")
    }
  }
  if(sum(duplicated(Res0$Description))>0){
    cat("Duplicated Terms:        will only keep the first\n")
    print(Res0$Description[duplicated(Res0$Description)])
    Res0<-Res0%>%dplyr::filter(!duplicated(Res0$Description))
  }
  if(nrow(Res0)<=2){
    stop("Not enough enrichment result, please consider a lower P or Q value!\nOr wrong gene set. \n")
  }
  x@Results<-Res0%>%as.data.frame()
  NumFiltered<-sum((getNumOfGENEs(x))<=nGeneCutOff)
  x@Results<-x@Results%>%dplyr::filter((getNumOfGENEs(x))>nGeneCutOff)
  cat("Filtered ",NumFiltered," results with less than ",nGeneCutOff," genes enriched. \n PCutoff:",PVal,"  QCutoff (in ORA)",QVal," GO ORA simplify cutoff:",simplify_cutoff,". \n")
  x@IsAnalysed<-T
  x@IsClustered<-F
  x@IsNetworkClustered<-F
  x@dbName<-DBlist[DBlist%in%c("GO","KEGG","Reactome","SelfDefinedGS")]
  message("Finished!\n")
  return(x)
}

CalcDUP<-function(x,y){
  GeneListA<-unlist(strsplit(x,"/"))
  GeneListB<-unlist(strsplit(y,"/"))
  if(length(GeneListA)>=length(GeneListB)){
    largegs<-GeneListA
    smallgs<-GeneListB
  }
  else{
    largegs<-GeneListB
    smallgs<-GeneListA
  }
  dupNum<-sum(smallgs%in%largegs)
  dupRate<-dupNum/length(largegs)
  return(dupRate)
}


mapTBs.<-function (Mat, FilterList, mode = "V2")
{
  if (mode == "V2") {
    FilterList <- FilterList[FilterList %in% rownames(Mat)]
  }
  Mat <- as.data.frame(Mat) %>% dplyr::filter(rownames(Mat) %in% FilterList)
  rn1 <- rownames(Mat)
  rn2 <- FilterList
  Mat <- Mat[rn1[match(rn2, rn1)], ]
  return(Mat)
}



setGeneric("getNumOfGENEs",function(x,...) standardGeneric("getNumOfGENEs"))


setMethod("getNumOfGENEs",signature(x="BioThemeFinder.ORA"),function(x,...){
  res<-unlist(lapply(strsplit(x@Results$geneID,"/"),length))
  return(res)
})


setMethod("getNumOfGENEs",signature(x="BioThemeFinder.ORA_FC"),function(x,...){
  res<-unlist(lapply(strsplit(x@Results$geneID,"/"),length))
  return(res)
})

setMethod("getNumOfGENEs",signature(x="BioThemeFinder.GSEA"),function(x,...){
  res<-unlist(lapply(strsplit(x@Results$core_enrichment,"/"),length))
  return(res)
})




RemoveTerms<-function(x,Item){
  if(nrow(x@Results)<1){
    stop("Please first do analysis before remove.\n")
  }
  if(length(Item)==1){
    x@Results<-x@Results%>%dplyr::filter(!grepl(Item,.$Description))
  }
  if(length(Item)>1){
    x@Results<-x@Results%>%dplyr::filter(!Description%in%Item)#add ! ,fix bug
  }
  x@DupMatrix<-NULL
  x@IsClustered<-F
  return(x)
}# need to rewrite. if clustered, generate a new frame. remove selectterms.

# SelectTerms<-function(x,Item){
#   cat("TIPS:\nThis function only used after calculated duplicated matrix, to select a few pathways for following drawing.\nIf you want to modifiy the object before caluclating duplicate matrix, please use RemoveTerms function instead.\n")
#   x@SelectedResultNames<-Item[Item%in%((x@Results)$Description)]
#   return(x)
# }

ExtractGenes<-function(x,Cluster){
  clusterType <- match.arg(clusterType, c("MatrixResult", "NetworkResult"))
  if(Cluster=="MatrixResult"){
    ClsID<-"GSCluster"
    cat("Cluster list:\n")
    print(table(x@Results$GSCluster))
  }
  if(Cluster=="NetworkResult"){
    ClsID<-"NetworkCluster"
    cat("Cluster list:\n")
    print(table(x@Results$NetworkCluster))
  }
  if(sum(c("GSCluster","NetworkCluster")%in%colnames(x@Results))!=0){
    DF0<-x@Results%>%dplyr::filter(!!sym(ClsID)==Cluster)
    if("geneID"%in%colnames(DF0)){
      genelst<-unlist(strsplit(DF0$geneID,"/"))
    }
    if("core_enrichment"%in%colnames(DF0)){
      genelst<-unlist(strsplit(DF0$core_enrichment,"/"))
    }
    doAnn<-F
    if(x@Specics=="human"){
      require(org.Hs.eg.db)
      OrgDB = org.Hs.eg.db
      doAnn<-T
    }
    if(x@Specics=="mouse"){
      require(org.Mm.eg.db)
      OrgDB = org.Mm.eg.db
      doAnn<-T
    }
    if(doAnn==T){
      Names0<-mapIds(x = OrgDB, keys = genelst,column = "GENENAME", keytype = "SYMBOL", multiVals = "first")
      REsDF<-data.frame(Gene=genelst,Name=Names0,Cluster=Cluster)
    }
    else{
      REsDF<-data.frame(Gene=genelst,Cluster=Cluster)
    }
    #casewhen
    return(REsDF)
  }
  else{
    stop("unclustered.\n")
  }
}

write.resultDF<-function(x,file){
  write.csv(x@Results,file=file)
}

