# HARDAM Project Documentation

HARDAM trains HAR-informed expression prediction models and runs S-PrediXcan/MetaXcan-style TWAS against GWAS summary statistics. The current tested local workflow uses a chromosome 1 GTEx genotype subset named `test_geno` and the PFC HAR+CIS model.

Run all commands from the repository root. Paths below are repo-relative.

## 1. Tested Local Files

The cleaned test workflow keeps these files:

- Genotype prefix: `geno_data/test_geno`
  - `geno_data/test_geno.bed`
  - `geno_data/test_geno.bim` - 200 chromosome 1 SNPs
  - `geno_data/test_geno.fam` - 300 individuals
- Expression file: `geno_data/Brain_Frontal_Cortex_BA9.v8.normalized_expression.chr1_test20_v2.bed`
- HAR BED file: `hars.bed`
- GENCODE annotation: `annotations/gencode.v26.annotation.gtf`
- Test GWAS summary statistics: `sumstats/test_with_cis_sumstats.tsv.gz`
- Final HAR+CIS model DB: `database/Brain_Frontal_Cortex_BA9_with_cis.db`
- Final HAR+CIS covariance: `covariances/snp_weights_Brain_Frontal_Cortex_BA9_with_cis_cov.tsv.gz`
- Final test TWAS output: `TWAS_results/test_with_cis_PFC_with_cis.csv`

No HAR expansion CSV is required for this tested workflow.

## 2. Train HAR+CIS Weights

Run `har-weight-train.R` in `har_cis` mode using the `test_geno` PLINK prefix:

```bash
Rscript har-weight-train.R \
  geno_data/test_geno \
  hars.bed \
  annotations/gencode.v26.annotation.gtf \
  geno_data/Brain_Frontal_Cortex_BA9.v8.normalized_expression.chr1_test20_v2.bed \
  1 \
  har_cis
```

The training script accepts these modes:

- `har_only`
- `har_cis`
- `cis_only`

For the tested with_cis run, the important outputs are:

- `output_results_seed/model_results_Brain_Frontal_Cortex_BA9.chr1_test20_v2_with_cis.txt`
- `output_results_seed/snp_weights_Brain_Frontal_Cortex_BA9.chr1_test20_v2_with_cis.txt`
- `output_results_seed/predicted_expression_Brain_Frontal_Cortex_BA9.chr1_test20_v2_with_cis.txt`

## 3. Create Model Database

Create the MetaXcan-compatible SQLite DB from the with_cis weights and model results:

```bash
python3 software_deps/create_db.py \
  --weights output_results_seed/snp_weights_Brain_Frontal_Cortex_BA9.chr1_test20_v2_with_cis.txt \
  --results output_results_seed/model_results_Brain_Frontal_Cortex_BA9.chr1_test20_v2_with_cis.txt \
  --expr geno_data/Brain_Frontal_Cortex_BA9.v8.normalized_expression.chr1_test20_v2.bed \
  --out database/Brain_Frontal_Cortex_BA9_with_cis.db
```

Expected test DB contents:

- `weights`: 161 rows
- `extra`: 10 genes
- `sample_info`: 175 expression samples

The model weights are trained on PLINK `--recode A` dosages. For this genotype, those dosages count BIM `allele1`, so `create_db.py` writes `eff_allele=allele1` and `ref_allele=allele2`. This allele convention is required for correct GWAS beta/z-score alignment in S-PrediXcan.

## 4. Create SNP Covariance

Create the long-form SNP covariance file from `test_geno` and the with_cis model weights:

```bash
python3 software_deps/create_covariance.py \
  --plink geno_data/test_geno \
  --weights output_results_seed/snp_weights_Brain_Frontal_Cortex_BA9.chr1_test20_v2_with_cis.txt \
  --expr geno_data/Brain_Frontal_Cortex_BA9.v8.normalized_expression.chr1_test20_v2.bed \
  --out covariances/snp_weights_Brain_Frontal_Cortex_BA9_with_cis_cov.tsv.gz
```

Expected test covariance size:

- `covariances/snp_weights_Brain_Frontal_Cortex_BA9_with_cis_cov.tsv.gz`: 1,888 lines including header

## 5. Run HARDAM

`exec_HARDAM.sh` reads GWAS files from `sumstats/`, model DBs from `database/`, covariance files from `covariances/`, and writes outputs to `TWAS_results/`.

Use `HARDAM_DB_GLOB` to select the with_cis DB explicitly:

```bash
PYTHON_BIN=python3 HARDAM_DB_GLOB=database/Brain_Frontal_Cortex_BA9_with_cis.db \
  bash exec_HARDAM.sh \
  test_with_cis \
  'test_with_cis_sumstats.tsv.gz' \
  ID A1 A2 BETA SE \
  --additional_output \
  --overwrite
```

This writes:

```text
TWAS_results/test_with_cis_PFC_with_cis.csv
```

The tested output has 11 lines: 1 header plus 10 gene results.

For another GWAS file, place the file in `sumstats/` and keep the same DB selector:

```bash
PYTHON_BIN=python3 HARDAM_DB_GLOB=database/Brain_Frontal_Cortex_BA9_with_cis.db \
  bash exec_HARDAM.sh \
  <trait_name> \
  '<gwas_file_pattern.tsv.gz>' \
  <snp_column> <effect_allele_column> <non_effect_allele_column> <beta_column> <se_column> \
  --additional_output \
  --overwrite
```

If `--gwas_h2` and `--gwas_N` are not passed through the wrapper, MetaXcan reports that p-values and z-scores are uncalibrated for inflation. That warning is expected for the local test command above.

## 6. Environment Setup

### Python

The wrapper defaults to `python3.9`, but the tested local command overrides this with `PYTHON_BIN=python3`.

Recommended setup:

```bash
conda create -n HARDAM_env python=3.9
conda activate HARDAM_env
pip install numpy scipy pandas pyarrow h5py statsmodels cyvcf2 bgen_reader pyliftover
```

The vendored MetaXcan code has been patched for current pandas/numpy compatibility in the local repo.

### R

Use R 4.4.2 and install:

```r
install.packages(c("bigsnpr", "glmnet", "grpreg", "dplyr", "data.table", "stringr"))
```

`plink` must also be available on `PATH`.

## 7. Quick Verification

Check the retained test genotype and final TWAS output:

```bash
wc -l geno_data/test_geno.bim geno_data/test_geno.fam TWAS_results/test_with_cis_PFC_with_cis.csv
```

Expected counts:

```text
200 geno_data/test_geno.bim
300 geno_data/test_geno.fam
11 TWAS_results/test_with_cis_PFC_with_cis.csv
```

Check the with_cis DB:

```bash
sqlite3 database/Brain_Frontal_Cortex_BA9_with_cis.db \
  'select "weights", count(*) from weights union all select "extra", count(*) from extra union all select "sample_info", count(*) from sample_info;'
```

Expected:

```text
weights|161
extra|10
sample_info|1
```

## 8. References

- [MetaXcan Wiki](https://github.com/hakyimlab/MetaXcan/wiki)
- [GENCODE Human Release 26](https://www.gencodegenes.org/human/release_26.html)

## 9. License and Credits

This project uses the MetaXcan and S-PrediXcan software, which are licensed under the MIT License. See: https://github.com/hakyimlab/MetaXcan

Copyright (c) 2016 Hakymlab

Last updated: May 24, 2026
