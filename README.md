# Capriotti_et_al_2026_Rprolixus_CNS_RNASeq
This repository contains the code, experiment description, count tables, and differential expression analysis for the RNA-Seq study investigating the transcriptional modulation in the nervous system of Rhodnius prolixus fifth-instar nymphs following a blood meal.

# Transcriptomics of the nervous system from Rhodnius prolixus, and its modulation after a blood meal in an immature stage

R code used to identify differentially expressed genes (DEGs) between fed and unfed insects in the manuscript:

> [Full citation of the manuscript: authors, title, journal, year, DOI]

[![DOI](https://zenodo.org/badge/DOI/[ZENODO-DOI].svg)](https://doi.org/[ZENODO-DOI])

## Overview

We analyze gene-level counts generated with STAR (16,138 annotated genes) using DESeq2. Genes with zero counts across all samples are excluded automatically by DESeq2 (14,992 genes tested). No additional gene-filtering step is applied.

- **Contrast:** Fed vs Unfed (reference level = Unfed). Positive log2 fold changes indicate higher expression in Fed.
- **Model:** negative binomial GLM, `~ treatment`, Wald test against |log2FC| > 1 (`lfcThreshold = 1`).
- **Effect-size shrinkage:** apeglm (`lfcShrink`).
- **DEG criterion:** apeglm s-value < 0.05, evaluated against the |log2FC| > 1 threshold.
- **Other settings:** `independentFiltering = TRUE`, `cooksCutoff = TRUE`. Independent filtering only affects adjusted p-values and does not influence s-value-based calls. Genes flagged as outliers by Cook's distance have their s-value set to `NA`.
- **Exploratory analyses:** rlog transformation (`blind = TRUE`), sample-distance heatmap and PCA (500 most variable genes).

## Repository contents

```
.
├── Read trimming code with Trimmomatic and mapping/counting code with STAR   # analysis script
├── Differential_expression_analysis.R   # analysis script
├── input/
│   ├── Raw_count_table.txt                # raw counts (see format below)
│   └── Experiment_description.txt             # sample information
├── output/                             # created when the script is run
│   ├── Deseq2_statistics.txt                # Differential expression analysis output
│   ├── GSR_Biological_Process.txt                # GO-enrichment output for Biological process
│   ├── GSR_Molecular_Function.txt                # GO-enrichment output for Molecular Function
│   ├── GSR_Cellular_Component.txt                # GO-enrichment output for Cellular Component
└── README.md
```

## Requirements

- R RStudio 2026.07.1+147 "Pacific Dogwood" Release (49299327da3b03a79e7aa615c2388bcd05a1261a, 2026-07-15) for Ubuntu Jammy
Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) rstudio/2026.07.1+147 Chrome/146.0.7680.216 Electron/41.9.0 Safari/537.36, Quarto 1.9.38
- R packages: `DESeq2` [1.42.1], `apeglm` [version 1.24], `pheatmap` [version 1.0.13], `RColorBrewer` [version 1.1-13]

Exact versions used are recorded in `output/sessionInfo.txt` after running the script.

Installation:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "apeglm"))
install.packages(c("pheatmap", "RColorBrewer"))
```

## Input format

Both files are tab-delimited with a header row.

**`input/Raw_count_table.txt`**: first column = gene IDs; remaining columns = raw (unnormalized) integer counts, one column per sample.

| gene_id | sample1 | sample2 | ... |
|---------|---------|---------|-----|
| gene_A  | 120     | 98      | ... |

**`input/Experiment_description.txt`**: one row per sample, in the same order as the count columns. Must contain a column named `treatment` with the values `Fed` and `Unfed`.

## Usage

From the repository folder:

```bash
Rscript Differential_expression_analysis.R
```

or open the script in RStudio and run it. Paths, thresholds, and the decimal separator are set in the "Parameters" section at the top of the script.

## Outputs (written to `output/`)

| File | Description |
|------|-------------|
| `DESeq2_results_all_genes.txt` | All tested genes: baseMean, shrunken log2FoldChange, lfcSE, s-value, Wald p-value and adjusted p-value (both against the |log2FC| > 1 threshold) |
| `DEGs_Fed_vs_Unfed.txt` | DEGs (s-value < 0.05), with regulation direction |
| `MAplot_Fed_vs_Unfed.tif` | MA plot |
| `rlog_values.txt` | rlog-transformed expression values |
| `Sample_distance_heatmap.tif` | Sample-to-sample distance heatmap |
| `PCA_rlog.tif` | PCA plot |
| `sessionInfo.txt` | Software versions used |

## Data availability

- Raw sequencing data: [repository and accession number, e.g., NCBI SRA/GEO: XXXX]
- Count table used in this analysis: [included in `input/` / deposited at XXXX]

## Citation

If you use this code, please cite the manuscript above and the archived version of this repository:

> [Natalia Capriotti; Lucila Traverso; Jose Manuel Latorre Estivalis; Ivana Sierra; Juan P. Ianowski; Sheila Ons]. [Capriotti_et_al_2026_Rprolixus_CNS_RNASeq]. Zenodo. [Year]. https://doi.org/[ZENODO-DOI]

DESeq2: Love MI, Huber W, Anders S (2014). Genome Biology 15:550.
apeglm: Zhu A, Ibrahim JG, Love MI (2019). Bioinformatics 35:2084–2092.

## License

[e.g., MIT License; see `LICENSE`]

## Contact

Jose Manuel Latorre Estivalis, IBBEA - CONICET UBA, jmlatorre@conicet.gov.ar
