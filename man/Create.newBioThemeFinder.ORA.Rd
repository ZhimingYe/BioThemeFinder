% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ObjectFrame.R
\name{Create.newBioThemeFinder.ORA}
\alias{Create.newBioThemeFinder.ORA}
\title{When only gene names are available, create BioThemeFinder object using ORA.}
\usage{
Create.newBioThemeFinder.ORA(Gene, FromType = "SYMBOL", Species = "human")
}
\arguments{
\item{Gene}{A vector containing genes}

\item{FromType}{Type of genes, can be the columns supported by org.DBs like SYMBOL for gene symbol or ENSEMBL for Ensembl ID}

\item{Species}{one of the human, mouse of non-gene}
}
\value{
A BioThemeFinder.ORA object
}
\description{
When only gene names are available, create BioThemeFinder object using ORA.
}
\author{
Zhiming Ye
}
