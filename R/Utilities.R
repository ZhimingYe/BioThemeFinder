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

MultiDBanalysis0<-function(x,PVal = 0.05,QVal = 0.1,DBlist= c("GO","KEGG"),simplify_cutoff=0.5){
  message("***BioThemeFinder***\nEnrichment analysis will be finished by clusterProfiler ",packageVersion('clusterProfiler')," and ReactomePA ",packageVersion('ReactomePA')," .\nPlease cite article The Innovation. 2021, 2(3):100141. doi: 10.1016/j.xinn.2021.100141 when using them. \n")
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
    Res0<-rbind(Res0,doKEGGanalysis(x=x,PValCutOff = PVal,QValCutOff=QVal))
  }
  if("Reactome"%in%DBlist){
    cat("Applying Reactome analysis...\n")
    Res0<-rbind(Res0,doReactomeanalysis(x=x,PValCutOff = PVal,QValCutOff=QVal))
  }
  x@Results<-Res0%>%as.data.frame()
  x@IsAnalysed<-T
  x@dbName<-DBlist
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
  # warning("in Version 5.6 update, new V2 mode is enabled. you can choose the V1 mode to change to the origin mode.\n")
  if (mode == "V2") {
    FilterList <- FilterList[FilterList %in% rownames(Mat)]
  }
  Mat <- as.data.frame(Mat) %>% dplyr::filter(rownames(Mat) %in% FilterList)
  rn1 <- rownames(Mat)
  rn2 <- FilterList
  Mat <- Mat[rn1[match(rn2, rn1)], ]
  return(Mat)
}