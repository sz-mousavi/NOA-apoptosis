# Single-cell analysis of human spermatogenesis

R Markdown notebook for the human spermatogenesis scRNA-seq workflow:

1. Build Seurat objects from sparse count matrices  
2. Merge related samples (SPG, SPM, SPC, SPT)  
3. Quality control (MAD-based filters and mitochondrial RNA)  
4. Normalization and integration  
5. PCA, UMAP, and clustering  
6. Marker visualization and cell-cycle scoring  
7. Differential expression  
8. Monocle 3 pseudotime  

Notebook file: [`SPM.Rmd`](SPM.Rmd)

## Repository contents

| File | Description |
|------|-------------|
| `SPM.Rmd` | Analysis notebook |
| `README.md` | This file |
| `.gitignore` | Excludes large data and Seurat objects |

Raw counts and `.h5seurat` checkpoints are **not** included. They are too large for GitHub.

## Expected data layout

Place local files like this before running the notebook:

```text
.
├── SPM.Rmd
├── README.md
└── data/
    ├── Mitogene.csv
    ├── counts/
    │   ├── Sparse_SRR6459187.tsv.gz
    │   ├── Sparse_SRR6459188.tsv.gz
    │   ├── Sparse_SRR6459189.tsv.gz
    │   ├── Sparse_SRR6459190.tsv.gz
    │   ├── Sparse_SRR6459191.tsv.gz
    │   ├── Sparse_SRR6459192.tsv.gz
    │   ├── Sparse_SRR6459193.tsv.gz
    │   ├── Sparse_SRR6459194.tsv.gz
    │   ├── Sparse_SRR7553900.tsv.gz
    │   └── Sparse_SRR7553901.tsv.gz
    └── state_objects/    # optional saved objects
```

Input accessions used in the notebook:

- `SRR6459187`–`SRR6459194`
- `SRR7553900`, `SRR7553901`

## How to run

This notebook is computationally heavy and is not executed by GitHub.

```r
# install once
install.packages(c(
  "Seurat", "dplyr", "ggplot2", "patchwork",
  "clusterProfiler", "enrichplot", "ggVennDiagram", "rmarkdown"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("EnsDb.Hsapiens.v86", "org.Hs.eg.db", "monocle3"))

# SeuratDisk / SeuratWrappers as needed for your Seurat version
```

Then:

```r
rmarkdown::render("SPM.Rmd", output_format = "github_document")
```

or open `SPM.Rmd` in RStudio and knit it.

## Interactive steps

These chunks are set to `eval=FALSE` because they need a local session:

- `CellSelector()` for manual outlier removal
- Monocle 3 `order_cells()` for choosing the trajectory root
- optional embedding-outlier cleanup in spermatogonia

Run those chunks interactively after the upstream objects exist.

## Outputs written locally

If the export chunks are enabled, results go to:

- `results/DEG.xlsx`
- `CData_cells_per_cluster.xlsx`
- `SeuratProject.h5Seurat`
- `GCData.h5Seurat`

Those files should stay out of git.
