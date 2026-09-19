# NOA-apoptosis

Transcriptomic analysis of **non-obstructive azoospermia (NOA)**, with emphasis on **spermatogenesis** process.

This repository contains three analysis modules.

---

## Repository contents

| Folder | Analysis | Public dataset |
|--------|----------|----------------|
| `NOA_bulk_RNA-seq/` | Bulk RNA-seq (DESeq2) | [GSE190752](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190752) |
| `Microarray_NOA/` | Dual-channel NOA microarray (limma) | [GSE145467](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145467) |
| `Spermatogenesis/` | scRNA-seq and trajectory (Seurat, Monocle3) | [GSE109037](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109037) |

---

## 1. Bulk RNA-seq (`NOA_bulk_RNA-seq/`)

**Inputs**

- raw count matrix  
- phenotype table with a `status` column (e.g. `Normal` vs `NOA`)

Dataset: [GSE190752](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190752)  
Altered mRNA profile in testis of patients with secondary idiopathic NOA.

Run the scripts in that folder in R after installing Bioconductor packages (`DESeq2`, …).

---

## 2. Microarray (`Microarray_NOA/`)

Limma-based differential expression for a dual-channel NOA microarray dataset.

Dataset: [GSE145467](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145467)  
Transcriptome changes in patients with severely impaired spermatogenesis.

---

## 3. Spermatogenesis (`Spermatogenesis/`)

Single-cell RNA-seq processing and trajectory analysis with **Seurat** and **Monocle3**.

Dataset: [GSE109037](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109037)  
10x Genomics / Drop-seq of adult human spermatogonia, spermatocytes, spermatids, and steady-state spermatogenic cells.

---

## How to use

1. Clone the repository.  
2. Open the folder that matches the analysis you want.  
3. Follow the comments inside the scripts in that folder.
