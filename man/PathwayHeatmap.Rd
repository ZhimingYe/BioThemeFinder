% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PreviewDupMat.R
\name{PathwayHeatmap}
\alias{PathwayHeatmap}
\alias{PathwayHeatmap,BioThemeFinder.ORA-method}
\alias{PathwayHeatmap,BioThemeFinder.ORA_FC-method}
\alias{PathwayHeatmap,BioThemeFinder.GSEA-method}
\title{Draw heatmap of repetition rate matrix}
\usage{
PathwayHeatmap(x, using_cluster = T, ...)

\S4method{PathwayHeatmap}{BioThemeFinder.ORA}(x, using_cluster = T, clusterType = "NetworkResult", ...)

\S4method{PathwayHeatmap}{BioThemeFinder.ORA_FC}(x, using_cluster = T, clusterType = "NetworkResult", ...)

\S4method{PathwayHeatmap}{BioThemeFinder.GSEA}(x, using_cluster = T, clusterType = "NetworkResult", ...)
}
\arguments{
\item{x}{BioThemeFinder object}

\item{using_cluster}{Draw figure using cluster result or not}

\item{clusterType}{can be MatrixResult or NetworkResult}
}
\description{
Draw heatmap of repetition rate matrix
}
\author{
Zhiming Ye
}
