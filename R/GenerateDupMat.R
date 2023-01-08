setGeneric("GenerateDupMat",function(x,...) standardGeneric("GenerateDupMat"))


#' @rdname GenerateDupMat
#' @title generate repetition rate matrix
#' @param x BioThemeFinder.ORA.
#'
#' @return BioThemeFinder object
#' @export
#' @author Zhiming Ye
#' @examples
setMethod("GenerateDupMat",signature(x="BioThemeFinder.ORA"),function(x,...){
  dupmat0<-matrix(nrow=nrow(x@Results),ncol=nrow(x@Results))
  geneIDcol<-which(colnames(x@Results)=="geneID")
  for(i in 1:nrow(x@Results)){
    for(j in 1:nrow(x@Results)){
      dupmat0[i,j]<-CalcDUP(unlist(x@Results[i,geneIDcol]),unlist(x@Results[j,geneIDcol]))
    }
  }
  colnames(dupmat0)<-x@Results$Description
  rownames(dupmat0)<-x@Results$Description
  x@DupMatrix<-dupmat0%>%as.data.frame()
  x@IsClustered<-T
  return(x)
})


#' @rdname GenerateDupMat
#' @title generate repetition rate matrix
#' @param x BioThemeFinder.ORA_FC.
#'
#' @return BioThemeFinder object
#' @export
#' @author Zhiming Ye
#' @examples
setMethod("GenerateDupMat",signature(x="BioThemeFinder.ORA_FC"),function(x,...){
  dupmat0<-matrix(nrow=nrow(x@Results),ncol=nrow(x@Results))
  geneIDcol<-which(colnames(x@Results)=="geneID")
  for(i in 1:nrow(x@Results)){
    for(j in 1:nrow(x@Results)){
      dupmat0[i,j]<-CalcDUP(unlist(x@Results[i,geneIDcol]),unlist(x@Results[j,geneIDcol]))
    }
  }
  colnames(dupmat0)<-x@Results$Description
  rownames(dupmat0)<-x@Results$Description
  x@DupMatrix<-dupmat0%>%as.data.frame()
  x@IsClustered<-T
  return(x)
})

#' @rdname GenerateDupMat
#' @title generate repetition rate matrix
#' @param x BioThemeFinder.GSEA.
#'
#' @return BioThemeFinder object
#' @export
#' @author Zhiming Ye
#' @examples
setMethod("GenerateDupMat",signature(x="BioThemeFinder.GSEA"),function(x,...){
  dupmat0<-matrix(nrow=nrow(x@Results),ncol=nrow(x@Results))
  geneIDcol<-which(colnames(x@Results)=="core_enrichment")
  for(i in 1:nrow(x@Results)){
    for(j in 1:nrow(x@Results)){
      dupmat0[i,j]<-CalcDUP(unlist(x@Results[i,geneIDcol]),unlist(x@Results[j,geneIDcol]))
    }
  }
  colnames(dupmat0)<-x@Results$Description
  rownames(dupmat0)<-x@Results$Description
  x@DupMatrix<-dupmat0%>%as.data.frame()
  x@IsClustered<-T
  return(x)
})

