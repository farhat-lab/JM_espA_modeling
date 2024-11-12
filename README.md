# _espA_ Mathematical Modeling
## Introduction
This repository contains the codes of Jiazheng Miao's (MBI Class of 2025) Capstone Project, 
a mathematical modeling on the expression level of _espA_ in _Mycobacterium tuberculosis_.

## Pipeline
Execute the scripts according to the following order:
1. `sbatch_download.sh`: Download data listed in `dna_accession.txt` and `rna_accession.txt`
2. `sbatch_dnaseq_pe.sh`: Process paired-ended DNA-seq data
3. `sbatch_dnaseq_se.sh`: Process single-ended DNA-seq data
4. `sbatch_rnaseq.sh`: Process RNA-seq data
5. `organize_data.py`: Aggregate `VCF` files to a `CSV` file, convert RNA read counts to LogFKPM, and screen for RD8/RD236a deletions
6. `merge_replicates.R`: Combine identified variants and average the expression level across the technique replicates
7. `pca.R`: Decompose the variant matrix (_espA_ regulatory region excluded) into PCs
8. `modeling.Rmd`: Perform the modeling