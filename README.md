# NiV_bioinformatics_analysis

This repository contains the analysis code and processed data tables used in our study on Nipah virus (NiV)-induced host transcriptional responses and host-directed drug repurposing.

## Overview

Nipah virus infection induces dynamic host transcriptional remodeling. In this study, we integrated time-point differential expression analysis, temporal gene analysis, core temporal gene identification, host-response signature construction, CMap-based gene-level reversal analysis, GSEA, and network proximity-based drug prioritization to identify candidate host-directed therapeutic agents against NiV infection.

## Public datasets

The main NiV transcriptomic dataset analyzed in this study was obtained from the NCBI Gene Expression Omnibus (GEO) under accession number GSE166707.

Additional public transcriptomic datasets used for drug-related GSEA analyses were also obtained from GEO:

- `GSE46263`: metformin-related transcriptomic dataset.
- `GSE15483`: fenofibrate-related transcriptomic dataset.

This repository does not include raw sequencing data. Users should download the original datasets from GEO and modify file paths according to their local computing environment before running the scripts.

## Repository structure

- `code/`: R scripts used for the main bioinformatics analyses.
- `data/`: key processed data tables and input files used in this study.

## Code files

- `01_timepoint_DEG_and_enrichment_analysis.R`: time-point differential expression analysis and enrichment analysis.
- `02_temporal_and_cluster_enrichment_analysis.R`: temporal analysis, cluster enrichment analysis, enrichment visualization, and combined plotting.
- `03_core_temporal_genes_and_intersection_plot.R`: identification of core temporal genes and generation of the intersection plot.
- `04_network_proximity_drug_prioritization.R`: network proximity-based drug prioritization analysis.
- `05_cmap_gene_level_analysis.R`: CMap gene-level reversal analysis for NiV host-response signatures.
- `06_gsea_analysis.R`: GSEA analysis for candidate drug-related transcriptional effects.

## Data files

The `data/` folder contains processed input and result tables used in the analyses, including:

- core temporal genes;
- time-point DEG union genes;
- drug information and drug-target information;
- human PPI and NiV-related PPI files;
- temporal cluster gene mapping;
- gene-disease association results;
- core temporal gene annotations with expression direction and virus-host factor labels.

Large raw datasets and original database files are not included in this repository. Processed tables are provided to support reproducibility of the main analyses described in the manuscript.

## Notes

File paths in the R scripts should be modified according to the user's local computing environment before running the analyses.
