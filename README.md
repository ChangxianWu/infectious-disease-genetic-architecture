# Analysis code for: Shared host genetic architecture of susceptibility to bacterial and viral infections

Code accompanying:

> **Shared host genetic architecture of susceptibility to bacterial and viral infections: a landscape of
> genetic correlation, pleiotropy and colocalisation**
> Changxian Wu, Haoyi Ma, Hong Su, Xiaoman Zheng, Bo Zhang, Hui Zhao
> College of Agriculture and Animal Husbandry Science and Technology,
> Ningxia Vocational and Technical College for Minorities, Wuzhong 751100, Ningxia, China

This repository contains the **analysis code only**. Every intermediate table, every figure source table and
the frozen keyed-value file from which each reported number is read are provided as **Data Set S1**
(archived with the submission, SHA256 manifest included). Keeping data out of the repository is deliberate:
the archive is the citable record of the numbers, and the repository is the executable record of how they
were produced.

## What is here

| Path | Contents |
|---|---|
| `code/` | Analysis scripts in execution order: data acquisition, harmonisation, LD score regression (300 pairs), locus atlas, clumping, colocalisation (coloc.abf, coloc.susie), gene sets, MHC sensitivity, canonical-value freeze, figure and table builders, sensitivity layers |
| `code/ldsc_py3/` | Local Python 3 port of the published LD score regression estimators used for every genetic correlation in the paper |

## Reproduction order

1. `code/00_download_assets.sh` (LD score regression reference files, FinnGen release 12 endpoints) and
   `code/01_download_covidhgi.sh` (COVID-19 Host Genetics Initiative freeze 7, European ancestry).
   Both are sharded, resumable, and verify the reported content length before renaming a partial file.
2. `code/10_munge_finngen.py <PHENOCODE> <CASES> <CONTROLS>` and `code/10b_munge_hgi.py <C2|A2>`.
3. `code/12_matrix_runner.py` — pairwise genetic correlation, 300 pairs.
4. `14_clump_leads.py`, `15_locus_atlas.py`, `16_annotate_leads.py`, `17_coloc_prep.py`,
   `17b_coloc_N.py`, `18_coloc.R`, `22_hgi_coloc_prep.py`, `20_prep_susie_ld.py`,
   `19d_coloc_susie_final.R`, `19_gene_sets.py`, `24_pathway_enrichment.py`,
   `23_mhc_sensitivity.py`, `21_figures.py`, `27_tables.py`, `32_supp_xlsx.py`.
5. `25_canonical_numbers.py` freezes the reported values into one keyed table; the manuscript reads
   every number from that table and nothing is typed by hand.

Scripts whose names start with 39-54 are the additional evidence layers (inborn errors of immunity,
MHC subregions, cross-dataset triangulation, polygenic-score transfer, autoimmune anchors, ontology
panels, sensitivity analyses).

## Software

Python 3.13 (pandas, numpy, scipy, matplotlib); R 4.6.1 and R 4.4.2 (coloc, susieR); PLINK 1.9 for
clumping and LD matrices. The manuscript text was assembled from a single configuration file with Node.js
(docx, jszip); that assembly code is not needed to reproduce any number and is not included here.

## Estimator provenance

Genetic correlations use a local Python 3 port of the published LD score regression estimators. The port is
a syntax-level translation only (the estimator algebra is unchanged). Before any atlas analysis it was
validated against three controls:

- positive: type 2 diabetes x coronary heart disease, rg 0.366 (SE 0.031), p 4.2e-33, inside the
  literature band;
- negative: sign randomisation of one trait, rg 0.008 (SE 0.009), p 0.41, near zero;
- guard: full row permutation of one trait, where the published out-of-bounds rule triggers as expected.

## Data

All inputs are public summary statistics: FinnGen release 12, COVID-19 Host Genetics Initiative freeze 7
(European ancestry), the 1000 Genomes Phase 3 European reference panel, HGNC gene groups, Reactome
pathways and the PanelApp primary-immunodeficiency panel. Access dates are recorded in the manuscript's
Data availability section. No individual-level data are used or redistributed here.

## Citation

If you use this code, please cite the associated article:

> Wu C, Ma H, Su H, Zheng X, Zhang B, Zhao H. Shared host genetic architecture of susceptibility to
> bacterial and viral infections: a landscape of genetic correlation, pleiotropy and colocalisation.
> (Manuscript submitted; the journal citation and DOI will be added here once assigned.)

The frozen keyed-value table from which every number in the article is read is in Data Set S1, not in this
repository, so that the code and the numeric record stay separate.

## Licence

No licence file is included. The corresponding author chooses the licence at upload; until one is chosen,
default copyright applies and reuse requires permission. This is a deliberate choice, not an omission.
