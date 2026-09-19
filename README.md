# NOA-apoptosis

Transcriptomic analysis of **non-obstructive azoospermia (NOA)**, with emphasis on **spermatogenesis** process.

This repository contains three analysis modules.

---

## Repository contents

| Folder | Analysis |
|--------|----------|
| `NOA_bulk_RNA-seq/` | Bulk RNA-seq (DESeq2) |
| `Microarray_NOA/` | Dual-channel NOA microarray (limma) |
| `Spermatogenesis/` | Single-cell RNA-seq and trajectory analysis (Seurat, Monocle3) |

---

## 1. Bulk RNA-seq (`NOA_bulk_RNA-seq/`)

**Inputs**

- raw count matrix  
- phenotype table with a `status` column (e.g. `Normal` vs `NOA`)

Public dataset used here: [GSE190752](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190752).

Run the scripts in that folder in R after installing Bioconductor packages (`DESeq2`, …).

---

## 2. Microarray (`Microarray_NOA/`)

Limma-based differential expression for a dual-channel NOA microarray dataset.

---

## 3. Spermatogenesis (`Spermatogenesis/`)

Single-cell RNA-seq processing and trajectory analysis with **Seurat** and **Monocle3**.

---

## How to use

1. Clone the repository.  
2. Open the folder that matches the analysis you want.  
3. Follow the comments inside the scripts in that folder.

