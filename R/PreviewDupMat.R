setGeneric("PathwayHeatmap",function(x,using_cluster=T,...) standardGeneric("PathwayHeatmap"))


setMethod("PathwayHeatmap",signature(x="BioThemeFinder.ORA"),function(x,using_cluster=T,clusterType="NetworkResult",...){
  do_cluster<-F
  if(("GSCluster"%in%colnames(x@Results)&using_cluster)|("NetworkCluster"%in%colnames(x@Results)&using_cluster)){
    message("Drawing using existed matrix cluster.\nif want to remove cluster, Please use removeCluster function.\n")
    clusterType <- match.arg(clusterType, c("MatrixResult", "NetworkResult"))
    do_cluster<-T
    if(clusterType=="MatrixResult"){
      x@Results<-x@Results%>%dplyr::mutate(geneSetCluster=GSCluster)
    }
    else{
      x@Results<-x@Results%>%dplyr::mutate(geneSetCluster=case_when((!is.na(NetworkCluster))~as.character(NetworkCluster),is.na(NetworkCluster)~"TooLowRate"))
    }
  }
  if(x@IsClustered!=T){
    stop("unClustered! Please run GenerateDupMat before preview it.\n")
  }
  dupmat0<-x@DupMatrix%>%as.matrix()
  rowANN<-rowAnnotation(Database=x@Results$Database,Change=x@Results$RegType,show_legend=F,col=list(Database=c(`GOBP`="#264653",
                                                                                                                 `GOCC`="#2a9d8f",
                                                                                                                 `GOMF`="#e9c46a",
                                                                                                                 `KEGG`="#f4a261",
                                                                                                                 `Reactome`="#e76f51",
                                                                                                               `Self`="#e5989b")))
  colANN<-columnAnnotation(Database=x@Results$Database,Change=x@Results$RegType,col=list(Database=c(`GOBP`="#264653",
                                                                                                      `GOCC`="#2a9d8f",
                                                                                                      `GOMF`="#e9c46a",
                                                                                                      `KEGG`="#f4a261",
                                                                                                      `Reactome`="#e76f51",
                                                                                                    `Self`="#e5989b")))
  if(!do_cluster){
    Heatmap(dupmat0,name = "Duplicate Rate",col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),row_names_gp = gpar(fontsize = 9),column_names_gp = gpar(fontsize = 9),row_names_max_width = max_text_width(rownames(dupmat0),gp = gpar(fontsize = 9)),column_names_max_height = max_text_width(rownames(dupmat0),gp = gpar(fontsize = 3)),right_annotation = rowANN,top_annotation = colANN)
  }
  else{
    warning(paste0("using ",clusterType," to display result.\n"))
    Heatmap(dupmat0,name = "Duplicate Rate",col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),row_names_gp = gpar(fontsize = 9),column_names_gp = gpar(fontsize = 9),row_names_max_width = max_text_width(rownames(dupmat0),gp = gpar(fontsize = 9)),column_names_max_height = max_text_width(rownames(dupmat0),gp = gpar(fontsize = 3)),right_annotation = rowANN,top_annotation = colANN,column_split = x@Results$geneSetCluster,row_split = x@Results$geneSetCluster)
  }
})


setMethod("PathwayHeatmap",signature(x="BioThemeFinder.ORA_FC"),function(x,using_cluster=T,clusterType="NetworkResult",...){
  do_cluster<-F
  if(("GSCluster"%in%colnames(x@Results)&using_cluster)|("NetworkCluster"%in%colnames(x@Results)&using_cluster)){
    message("Drawing using existed matrix cluster.\nif want to remove cluster, Please use removeCluster function.\n")
    clusterType <- match.arg(clusterType, c("MatrixResult", "NetworkResult"))
    do_cluster<-T
    if(clusterType=="MatrixResult"){
      x@Results<-x@Results%>%dplyr::mutate(geneSetCluster=GSCluster)
    }
    else{
      x@Results<-x@Results%>%dplyr::mutate(geneSetCluster=case_when((!is.na(NetworkCluster))~as.character(NetworkCluster),is.na(NetworkCluster)~"TooLowRate"))
    }
  }
  if(x@IsClustered!=T){
    stop("unClustered! Please run GenerateDupMat before preview it.\n")
  }
  dupmat0<-x@DupMatrix%>%as.matrix()
  rowANN<-rowAnnotation(Database=x@Results$Database,Change=x@Results$RegType,show_legend=F,col=list(Database=c(`GOBP`="#264653",
                                                                                                                 `GOCC`="#2a9d8f",
                                                                                                                 `GOMF`="#e9c46a",
                                                                                                                 `KEGG`="#f4a261",
                                                                                                                 `Reactome`="#e76f51",
                                                                                                               `Self`="#e5989b"),
                                                                                                      Change=c(`Favor_UpReg`="#e63946",
                                                                                                               `Favor_DnReg`="#1d3557",
                                                                                                               `UnKnown`="#fefae0")))
  colANN<-columnAnnotation(Database=x@Results$Database,Change=x@Results$RegType,col=list(Database=c(`GOBP`="#264653",
                                                                                                      `GOCC`="#2a9d8f",
                                                                                                      `GOMF`="#e9c46a",
                                                                                                      `KEGG`="#f4a261",
                                                                                                      `Reactome`="#e76f51",
                                                                                                    `Self`="#e5989b"),
                                                                                           Change=c(`Favor_UpReg`="#e63946",
                                                                                                    `Favor_DnReg`="#1d3557",
                                                                                                    `UnKnown`="#fefae0")))
  if(!do_cluster){
    Heatmap(dupmat0,name = "Duplicate Rate",col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),row_names_gp = gpar(fontsize = 9),column_names_gp = gpar(fontsize = 9),row_names_max_width = max_text_width(rownames(dupmat0),gp = gpar(fontsize = 9)),column_names_max_height = max_text_width(rownames(dupmat0),gp = gpar(fontsize = 3)),right_annotation = rowANN,top_annotation = colANN)
  }
  else{
    warning(paste0("using ",clusterType," to display result.\n"))
    Heatmap(dupmat0,name = "Duplicate Rate",col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),row_names_gp = gpar(fontsize = 9),column_names_gp = gpar(fontsize = 9),row_names_max_width = max_text_width(rownames(dupmat0),gp = gpar(fontsize = 9)),column_names_max_height = max_text_width(rownames(dupmat0),gp = gpar(fontsize = 3)),right_annotation = rowANN,top_annotation = colANN,column_split = x@Results$geneSetCluster,row_split = x@Results$geneSetCluster)
  }
})

setMethod("PathwayHeatmap",signature(x="BioThemeFinder.GSEA"),function(x,using_cluster=T,clusterType="NetworkResult",...){
  do_cluster<-F
  if(("GSCluster"%in%colnames(x@Results)&using_cluster)|("NetworkCluster"%in%colnames(x@Results)&using_cluster)){
    message("Drawing using existed matrix cluster.\nif want to remove cluster, Please use removeCluster function.\n")
    clusterType <- match.arg(clusterType, c("MatrixResult", "NetworkResult"))
    do_cluster<-T
    if(clusterType=="MatrixResult"){
      x@Results<-x@Results%>%dplyr::mutate(geneSetCluster=GSCluster)
    }
    else{
      x@Results<-x@Results%>%dplyr::mutate(geneSetCluster=case_when((!is.na(NetworkCluster))~as.character(NetworkCluster),is.na(NetworkCluster)~"TooLowRate"))
    }
  }
  if(x@IsClustered!=T){
    stop("unClustered! Please run GenerateDupMat before preview it.\n")
  }
  dupmat0<-x@DupMatrix%>%as.matrix()
  rowANN<-rowAnnotation(Database=x@Results$Database,NES=x@Results$NES,show_legend=F,col=list(Database=c(`GOBP`="#264653",
                                                                                                                 `GOCC`="#2a9d8f",
                                                                                                                 `GOMF`="#e9c46a",
                                                                                                                 `KEGG`="#f4a261",
                                                                                                                 `Reactome`="#e76f51",
                                                                                                        `Self`="#e5989b"),
                                                                                                      Change=c(`Favor_UpReg`="#e63946",
                                                                                                               `Favor_DnReg`="#1d3557",
                                                                                                               `UnKnown`="#fefae0")))
  colANN<-columnAnnotation(Database=x@Results$Database,NES=x@Results$NES,col=list(Database=c(`GOBP`="#264653",
                                                                                                      `GOCC`="#2a9d8f",
                                                                                                      `GOMF`="#e9c46a",
                                                                                                      `KEGG`="#f4a261",
                                                                                                      `Reactome`="#e76f51",
                                                                                             `Self`="#e5989b"),
                                                                                            `NES` = colorRamp2(c(min(x@Results$NES), max(x@Results$NES)), c("#0e6ba8","#e56b6f"))))
  if(!do_cluster){
    Heatmap(dupmat0,name = "Duplicate Rate",col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),row_names_gp = gpar(fontsize = 9),column_names_gp = gpar(fontsize = 9),row_names_max_width = max_text_width(rownames(dupmat0),gp = gpar(fontsize = 9)),column_names_max_height = max_text_width(rownames(dupmat0),gp = gpar(fontsize = 3)),right_annotation = rowANN,top_annotation = colANN)
  }
  else{
    warning(paste0("using ",clusterType," to display result.\n"))
    Heatmap(dupmat0,name = "Duplicate Rate",col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),row_names_gp = gpar(fontsize = 9),column_names_gp = gpar(fontsize = 9),row_names_max_width = max_text_width(rownames(dupmat0),gp = gpar(fontsize = 9)),column_names_max_height = max_text_width(rownames(dupmat0),gp = gpar(fontsize = 3)),right_annotation = rowANN,top_annotation = colANN,column_split = x@Results$geneSetCluster,row_split = x@Results$geneSetCluster)
  }

})

