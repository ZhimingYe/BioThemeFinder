setGeneric("PreviewDupMat",function(x,...) standardGeneric("PreviewDupMat"))


setMethod("PreviewDupMat",signature(x="BioThemeFinder.ORA"),function(x,...){
  if(x@IsClustered!=T){
    stop("unClustered! Please run GenerateDupMat before preview it.\n")
  }
  dupmat0<-x@DupMatrix%>%as.matrix()
  rowANN<-rowAnnotation(Database=x@Results$Database,Change=x@Results$RegType,show_legend=F,col=list(Database=c(`GOBP`="#264653",
                                                                                                                 `GOCC`="#2a9d8f",
                                                                                                                 `GOMF`="#e9c46a",
                                                                                                                 `KEGG`="#f4a261",
                                                                                                                 `Reactome`="#e76f51")))
  colANN<-columnAnnotation(Database=x@Results$Database,Change=x@Results$RegType,col=list(Database=c(`GOBP`="#264653",
                                                                                                      `GOCC`="#2a9d8f",
                                                                                                      `GOMF`="#e9c46a",
                                                                                                      `KEGG`="#f4a261",
                                                                                                      `Reactome`="#e76f51")))
  Heatmap(dupmat0,name = "Duplicate Rate",col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),row_names_gp = gpar(fontsize = 9),column_names_gp = gpar(fontsize = 9),row_names_max_width = max_text_width(rownames(dupmat0),gp = gpar(fontsize = 3)),column_names_max_height = max_text_width(rownames(dupmat0),gp = gpar(fontsize = 3)),right_annotation = rowANN,top_annotation = colANN)
})


setMethod("PreviewDupMat",signature(x="BioThemeFinder.ORA_FC"),function(x,...){
  if(x@IsClustered!=T){
    stop("unClustered! Please run GenerateDupMat before preview it.\n")
  }
  dupmat0<-x@DupMatrix%>%as.matrix()
  rowANN<-rowAnnotation(Database=x@Results$Database,Change=x@Results$RegType,show_legend=F,col=list(Database=c(`GOBP`="#264653",
                                                                                                                 `GOCC`="#2a9d8f",
                                                                                                                 `GOMF`="#e9c46a",
                                                                                                                 `KEGG`="#f4a261",
                                                                                                                 `Reactome`="#e76f51"),
                                                                                                      Change=c(`Favor_UpReg`="#e63946",
                                                                                                               `Favor_DnReg`="#1d3557",
                                                                                                               `UnKnown`="#fefae0")))
  colANN<-columnAnnotation(Database=x@Results$Database,Change=x@Results$RegType,col=list(Database=c(`GOBP`="#264653",
                                                                                                      `GOCC`="#2a9d8f",
                                                                                                      `GOMF`="#e9c46a",
                                                                                                      `KEGG`="#f4a261",
                                                                                                      `Reactome`="#e76f51"),
                                                                                           Change=c(`Favor_UpReg`="#e63946",
                                                                                                    `Favor_DnReg`="#1d3557",
                                                                                                    `UnKnown`="#fefae0")))
  Heatmap(dupmat0,name = "Duplicate Rate",col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),row_names_gp = gpar(fontsize = 9),column_names_gp = gpar(fontsize = 9),row_names_max_width = max_text_width(rownames(dupmat0),gp = gpar(fontsize = 3)),column_names_max_height = max_text_width(rownames(dupmat0),gp = gpar(fontsize = 3)),right_annotation = rowANN,top_annotation = colANN)
})

setMethod("PreviewDupMat",signature(x="BioThemeFinder.GSEA"),function(x,...){
  if(x@IsClustered!=T){
    stop("unClustered! Please run GenerateDupMat before preview it.\n")
  }
  dupmat0<-x@DupMatrix%>%as.matrix()
  rowANN<-rowAnnotation(Database=x@Results$Database,NES=x@Results$NES,show_legend=F,col=list(Database=c(`GOBP`="#264653",
                                                                                                                 `GOCC`="#2a9d8f",
                                                                                                                 `GOMF`="#e9c46a",
                                                                                                                 `KEGG`="#f4a261",
                                                                                                                 `Reactome`="#e76f51"),
                                                                                                      Change=c(`Favor_UpReg`="#e63946",
                                                                                                               `Favor_DnReg`="#1d3557",
                                                                                                               `UnKnown`="#fefae0")))
  colANN<-columnAnnotation(Database=x@Results$Database,NES=x@Results$NES,col=list(Database=c(`GOBP`="#264653",
                                                                                                      `GOCC`="#2a9d8f",
                                                                                                      `GOMF`="#e9c46a",
                                                                                                      `KEGG`="#f4a261",
                                                                                                      `Reactome`="#e76f51"),
                                                                                            `NES` = colorRamp2(c(min(x@Results$NES), max(x@Results$NES)), c("#0e6ba8","#e56b6f"))))
  Heatmap(dupmat0,name = "Duplicate Rate",col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),row_names_gp = gpar(fontsize = 9),column_names_gp = gpar(fontsize = 9),row_names_max_width = max_text_width(rownames(dupmat0),gp = gpar(fontsize = 3)),column_names_max_height = max_text_width(rownames(dupmat0),gp = gpar(fontsize = 3)),right_annotation = rowANN,top_annotation = colANN)
})

