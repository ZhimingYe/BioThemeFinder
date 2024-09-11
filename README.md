# BioThemeFinder

This package is an expansion pack of the famous package "clusterProfiler". BioThemeFinder presents an enrichment analysis method for multi-database union (like GO, KEGG, Reactome and self-defined gene sets, etc. Also, redundant enrichment results can be turned into more explanatory biological stories using methods such as NMF or network-based community clustering.

## Description

![](figs/BioThemeFinderhomepage.png)

## Installation

    library(devtools)
    install_github("ZhimingYe/BioThemeFinder")
    
## Demo

### Loading R Packages

```r
library(BioThemeFinder)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
```

Here, the Gene names and differential expression values are vectors. The differential expression values can come from the log2FC of DEG analysis or other metrics like ROC power (e.g., from Seurat). `LUAD_DEG` is differential gene expression data for TCGA lung adenocarcinoma and is provided as sample data with the package.

### Creating a BioThemeFinder Object

1. `Create.newBioThemeFinder.ORA` – Only Gene names, using over-representation analysis (ORA).
2. `Create.newBioThemeFinder.ORAwithFC` – Contains Gene names and differential expression values, using over-representation analysis (ORA).
3. `Create.newBioThemeFinder.GSEA` – Contains Gene names and differential expression values, using GSEA analysis.

We will use the second scenario as the demo.

`Species` specifies the organism, supporting human and mouse, and `FromType` indicates the type of gene ID provided, such as `ENSEMBL`, `ENTREZID`, or `SYMBOL`.

You can use `?` before each function to view detailed documentation.

```r
data("LUAD_DEG")
BTFdemo <- Create.newBioThemeFinder.ORAwithFC(
  Gene = rownames(LUAD_DEG),
  log2FC = LUAD_DEG$log2FoldChange,
  Pvalue = LUAD_DEG$padj,
  FCcutoff = 1,
  PvalueCutOff = 0.05,
  Species = "human",
  FromType = "ENSEMBL"
)
```

### Performing Multi-Database Enrichment Analysis

`DBlist` can include GO, KEGG, Reactome, or SelfDefinedGS.

When using a custom database (`SelfDefinedGS`), you need to pass a custom gene set using `Term2GENE`.

Check detailed documentation using `?MultiDBanalysis`.

Additionally, you can use `parse.clusterProfiler.result` to directly import results from clusterProfiler.

```r
BTFdemo <- MultiDBanalysis(
  BTFdemo,
  PVal = 0.05,
  QVal = 0.05,
  DBlist = c("GO", "KEGG", "Reactome"),
  nGeneCutOff = 5,
  simplify_cutoff = 0.7,
  useMKEGG = TRUE
)
```

Preview enrichment results:

```r
PathwayStatsPlot(BTFdemo, orderBy = "pValue")
```

### Calculating Gene Redundancy Between Terms

You can pass the `EdgeCutoff` parameter to decide the threshold for redundancy between pathways in enrichment analysis to be included in the network plot. The default is 0.5.

```r
BTFdemo <- GenerateDupMat(BTFdemo)
```

You can preview the result using `PathwayHeatmap(BTFdemo)`. Note that the image might be too large for R's preview. It is recommended to export it as a PDF of approximately 50x50 inches and view it in a PDF viewer.

```r
PathwayHeatmap(BTFdemo)
```

Image interpretation:

### Network-Based Clustering

```r
BTFdemo <- NetworkClustering(BTFdemo)
```

### Network Clustering Visualization

```r
PlotNetwork(BTFdemo, method = "igraph", Label = TRUE)
PlotNetwork(BTFdemo, method = "igraph", Label = FALSE)
PlotNetwork(BTFdemo, method = "ggplot2", Label = TRUE)
```

### Matrix-Based Clustering

The available clustering tools are "nmf", "hc", "pam", and "fuzzy". You can specify the tool via the `method` parameter. `nmf` refers to non-negative matrix factorization, `hc` to hierarchical clustering, `pam` to Partitioning Around Medoids, and `fuzzy` to fuzzy clustering. For detailed help, see `?MatrixClustering`.

The value of `k` can be determined roughly by inspecting the `PathwayHeatmap`.

```r
BTFdemo <- MatrixClustering(BTFdemo, k = 20, method = "nmf")
```

### Visualization Based on Heatmap

`clusterType` can be one of "MatrixResult" or "NetworkResult".

```r
PathwayHeatmap(BTFdemo, using_cluster = TRUE, clusterType = "NetworkResult")
```

### Faceted Dot Plot

For detailed help, see `?PathwayStatsPlot`. Add `+ facet_grid(~Cluster)` for faceting.

You can observe that each module's function is quite clear and makes good use of multiple databases.

```r
PathwayStatsPlot(BTFdemo, orderBy = "pValue", clusterType = "NetworkResult") + facet_grid(~Cluster)
PathwayStatsPlot(BTFdemo, orderBy = "pValue", clusterType = "NetworkResult") + facet_grid(~Cluster)
```

### Extracting Results

```r
Genes <- ExtractGenes(BTFdemo, clusterType = "NetworkResult")
Result <- resultDF(BTFdemo)
```
