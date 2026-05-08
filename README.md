# final-project-pvarelad

# Huntington's Disease Prefrontal Cortex — RNA-seq Explorer

An interactive R Shiny application for exploring RNA-seq data from post-mortem prefrontal cortex samples of Huntington's Disease (HD) patients compared to neurologically healthy controls (Original data and study: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810).

---

## About the Data

The dataset contains 69 post-mortem prefrontal cortex samples (HD patients and controls) sourced from the NCBI GEO repository. The analysis pipeline consisted of:

- **Normalization:** FPKM normalized counts from original published study (GSE64810)
- **Differential Expression:** DESeq2 (default parameters), HD vs. Control; significance defined as FDR < 0.05 and |log2FC| > 1
- **Gene Set Enrichment:** fgsea against the MSigDB Hallmark gene set collection using DESeq2 Wald statistics as the ranking metric

---

## Features

The app is organized into four tabs:

### Samples
- Column-by-column summary (mean ± SD for numeric, distinct values for categorical)
- Full sortable metadata table
- Histograms of Age, RIN, and PMI colored by grouping variable
- Sample count bar chart

### Counts
- Variance percentile slider to filter genes in real time
- Filter summary table (genes passing/failing, % breakdown)
- Diagnostic scatter plots (median count vs. variance; median count vs. zero-count samples)
- Clustered heatmap of top N highest-variance genes with optional log2 transformation
- PCA plot with selectable principal components, colored by condition

### Differential Expression
- Adjustable padj and |log2FC| thresholds
- Summary of up/down/NS gene counts
- Full sortable and searchable results table
- Volcano plot with swappable axes
- MA plot

### GSEA
- **Top Pathways:** Bar chart of top N pathways by adjusted p-value; click a bar to see its full table entry
- **Table:** Filterable by padj and NES direction; downloadable as CSV
- **Scatter Plot:** NES vs. −log10(padj) with a reactive significance threshold line

---

## Input File Formats

| Tab | File | Required Columns |
|-----|------|-----------------|
| Samples | `metadata.csv` | `sample_id`, `geo_accession`, `condition` |
| Counts | `normalized_counts.csv` | Gene ID column + one column per sample (named by GEO accession) |
| DE | `deseq2_results.csv` | `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj` |
| GSEA | `gsea_results.csv` | `pathway`, `NES`, `padj`, `leadingEdge` (semicolon-separated) |

All files must be `.csv` or `.tsv`. 

---

## Installation

```r
install.packages(c("shiny", "tidyverse", "DT", "bslib", "pheatmap"))
```

```bash
git clone https://github.com/BF591-R/final-project-pvarelad.git
cd your-repo-name
```

```r
shiny::runApp("app.R")
```

---

## Repository Structure

```
.
├── FileProcessing.Rmd    # File read in and data pre-processing
├── app.R                 # Main Shiny application (UI + server)
├── README.md             # This file
└── data/                 # Example input files
```

---

## Known Issues

- Counts column names must match the `geo_accession` column in metadata, not `sample_id`
- GSEA `leadingEdge` must be semicolon-separated for correct parsing
- Upload the **unfiltered** normalized counts for best PCA and heatmap results

---

## Built With

[R Shiny](https://shiny.posit.co/) · [tidyverse](https://www.tidyverse.org/) · [DT](https://rstudio.github.io/DT/) · [bslib](https://rstudio.github.io/bslib/) · [pheatmap](https://cran.r-project.org/package=pheatmap) · [DESeq2](https://bioconductor.org/packages/DESeq2) · [fgsea](https://bioconductor.org/packages/fgsea)
