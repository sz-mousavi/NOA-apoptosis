# NOA-apoptosis

Transcriptomic analysis of **non-obstructive azoospermia (NOA)**, with emphasis on spermatogenesis and apoptosis-related genes. The repository contains three independent analysis modules, each an R Markdown notebook.

| Folder | Analysis | Public dataset |
|--------|----------|----------------|
| `NOA_bulk_RNA-seq/` | Bulk RNA-seq, differential expression (DESeq2) | [GSE190752](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190752) |
| `Microarray_NOA/` | Two-color Agilent microarray (limma) | [GSE145467](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145467) |
| `Spermatogenesis/` | Single-cell RNA-seq and trajectory (Seurat, Monocle3) | [GSE109037](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109037) |

Knitted GitHub-flavored reports (`.md`) of the DESeq2 and limma notebooks are included for viewing. The Spermatogenesis notebook is computationally heavy and is not knitted here; see `Spermatogenesis/README_spermatogenesis.md`.

## Installation

R 4.5.0 was used for the DESeq2 and limma notebooks.

```r
install.packages(c("BiocManager", "ggplot2", "pheatmap", "rmarkdown"))
BiocManager::install(c("DESeq2", "limma"))
```

Package versions from the original sessions:

| Notebook | Packages |
|----------|----------|
| `DESeq2_GSE190752_analysis.Rmd` | DESeq2 1.50.2, ggplot2 4.0.3, pheatmap 1.0.13 |
| `Agilent_GSE145467_limma.Rmd` | limma 3.66.0, ggplot2 4.0.3 |
| `SPM.Rmd` | Seurat, SeuratDisk, SeuratWrappers, monocle3, clusterProfiler, enrichplot, org.Hs.eg.db, EnsDb.Hsapiens.v86, ggVennDiagram, dplyr, patchwork (versions: run `sessionInfo()` in your own session; installation notes are in `Spermatogenesis/README_spermatogenesis.md`) |

## Data

Raw data are **not** stored in this repository. Download them from GEO and place them under `data/` next to the corresponding notebook. All notebooks use paths relative to their own folder; the input location can be changed through the `data_dir` variable at the top of each notebook. Output tables are written to `results/`.

### Bulk RNA-seq: GSE190752

Download `GSE190752_gene_count_matrix.csv.gz` from the GEO supplementary files and prepare two tab-separated files:

```text
NOA_bulk_RNA-seq/
├── DESeq2_GSE190752_analysis.Rmd
└── data/
    └── GSE190752/
        ├── counts.txt   # count matrix: genes as rows, samples as columns (Control1-3, NOA1-3)
        └── Pheno.txt    # sample table: samples as rows, with a "status" column (Normal / NOA)
```

`counts.txt` is the GEO count matrix saved as a tab-separated file with the gene identifier column used as row names.

### Microarray: GSE145467

Download `GSE145467_RAW.tar` from the GEO supplementary files and extract it into `GSE145467_RAW/`. The targets file lists file names ending in `.txt`, so decompress the `.txt.gz` files.

```text
Microarray_NOA/
├── Agilent_GSE145467_limma.Rmd
└── data/
    └── GSE145467/
        ├── pheno.txt        # limma targets file: FileName, Cy3, Cy5 columns
        ├── spot.txt         # spot-type file for control probes (limma readSpotTypes format)
        └── GSE145467_RAW/   # raw Agilent files listed in pheno.txt
```

### Single-cell: Spermatogenesis

See `Spermatogenesis/README_spermatogenesis.md` for the expected file layout (`data/counts/Sparse_SRR*.tsv.gz`, `data/Mitogene.csv`).

## How to reproduce

The three notebooks are independent and can be run in any order. Run each from its own folder (in RStudio, open the `.Rmd` and knit, or use `rmarkdown::render`):

1. `NOA_bulk_RNA-seq/DESeq2_GSE190752_analysis.Rmd`: differential expression, PCA, heatmap, volcano plot.
2. `Microarray_NOA/Agilent_GSE145467_limma.Rmd`: normalization, differential expression, PCA.
3. `Spermatogenesis/SPM.Rmd`: Seurat integration, markers, Monocle3 pseudotime (needs interactive steps, see its README).

```r
rmarkdown::render("NOA_bulk_RNA-seq/DESeq2_GSE190752_analysis.Rmd")
```

## License

MIT, see [LICENSE](LICENSE).
