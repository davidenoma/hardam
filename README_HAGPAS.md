# HAGPAS Project Documentation

## 1. Training with `har-weight-train.R`

### Gene Annotation File
- The input gene annotation file (`gene_annot_file`) must be downloaded from GENCODE: [GENCODE Release 26](https://www.gencodegenes.org/human/release_26.html).
- Format: GTF, build GRCh38.

### Example Usage
```
Rscript har-weight-train.R {genotype_prefix} hars.bed gencode.v26.annotation.gtf {tissue-expression}.bed {chromosome}
```
- `{genotype_prefix}`: Prefix for PLINK genotype files (e.g., `gtex_geno`)
- `hars.bed`: BED file with HAR locations
- `gencode.v26.annotation.gtf`: GENCODE gene annotation file
- `{tissue-expression}.bed`: Tissue-specific gene expression file
- `{chromosome}`: Chromosome number (e.g., `1`)

### R Version and Dependencies
- R version: 4.4.2
- Required R packages:
  - bigsnpr
  - glmnet
  - grpreg
  - dplyr
  - data.table
  - stringr

## 2. Running HAGPAS

### Example Command
```
bash exec_hagpas.sh scz "SCZ_sumstats_hg38_chr22.tsv.gz" ID A1 A2 BETA SE
```
- `scz`: Disease/trait label
- `SCZ_sumstats_hg38_chr22.tsv.gz`: GWAS summary statistics (chromosome 22)
- `ID`, `A1`, `A2`, `BETA`, `SE`: Column names in the GWAS file

### What Happens
- The script runs S-PrediXcan using the MetaXcan framework on the specified GWAS summary statistics, using the provided S-PrediXcan database and covariance files for the tissue.
- The S-PrediXcan database files should be called HAR trained weight databases.
- Covariance files are calculated from the GTEx genotypes.
- Output is written to the `TWAS_results/` directory.

## 3. Environment Setup

### Python Environment
- The default environment is `/opt/anaconda3/envs/oldnumpy`.
- For reproducibility, create a dedicated environment for HAGPAS with Python 3.9:

```
conda create -n hagpas_env python=3.9
conda activate hagpas_env
pip install numpy scipy pandas pyarrow h5py statsmodels cyvcf2 bgen_reader pyliftover
```
- Or use the provided `conda_env.yaml` as a template and update Python version to 3.9.

### R Environment
- Use R 4.4.2 and install the required packages:
```
install.packages(c("bigsnpr", "glmnet", "grpreg", "dplyr", "data.table", "stringr"))
```

## 4. Test Data
- HAR trained S-PrediXcan styled database weight files: `database/`
- Covariance files (from GTEx genotypes): `covariances/`
- Example summary statistics: `sumstats/SCZ_sumstats_hg38_chr22.tsv.gz`

## 5. References
- [MetaXcan Wiki](https://github.com/hakyimlab/MetaXcan/wiki)
- [GENCODE Human Release 26](https://www.gencodegenes.org/human/release_26.html)

## 6. License and Credits

This project uses the MetaXcan and S-PrediXcan software, which are licensed under the MIT License (see: https://github.com/hakyimlab/MetaXcan). 

> Copyright (c) 2016 Hakymlab
> 
> Permission is hereby granted, free of charge, to any person obtaining a copy
> of this software and associated documentation files (the "Software"), to deal
> in the Software without restriction, including without limitation the rights
> to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
> copies of the Software, and to permit persons to whom the Software is
> furnished to do so, subject to the following conditions:
> 
> The above copyright notice and this permission notice shall be included in all
> copies or substantial portions of the Software.

*Last updated: January 20, 2026*
