# HARDAM — HARs-Directed Association Mapping

HARDAM trains expression prediction models informed by human accelerated regions (HARs) and runs transcriptome-wide association studies (TWAS) using S-PrediXcan and GWAS summary statistics. This guide describes a chromosome 1 test workflow for frontal cortex, BA9 (PFC), using the HAR+CIS model.

<img width="800" height="600" alt="hardam-method-figure-1" src="https://github.com/user-attachments/assets/948644bc-c453-4a3b-a7d0-517fce683643" />


Author-provided test data allow users to exercise the workflow without accessing raw genotypes. The test example demonstrates model training and association testing

Run all commands from the repository root. Paths below are repo-relative.

## 1. Test Inputs and Generated Files

Before running, place the author-provided test inputs at the following paths:

- Genotype prefix: `geno_data/test_geno`
  - `geno_data/test_geno.bed`
  - `geno_data/test_geno.bim` - 200 chromosome 1 SNPs
  - `geno_data/test_geno.fam` - 300 individuals
- Expression file: `geno_data/Brain_Frontal_Cortex_BA9.v8.normalized_expression.chr1_test20_v2.bed`
- HAR BED file: `hars.bed`
- GENCODE annotation: `annotations/gencode.v26.annotation.gtf.gz`
- Test GWAS summary statistics: `sumstats/test_with_cis_sumstats.tsv.gz`

The workflow generates:

- HAR+CIS model DB: `database/Brain_Frontal_Cortex_BA9_with_cis.db`
- HAR+CIS covariance: `covariances/snp_weights_Brain_Frontal_Cortex_BA9_with_cis_cov.tsv.gz`
- Test TWAS output: `TWAS_results/test_with_cis_PFC_with_cis.csv`

HAR expansion is computed during training; no separate expansion CSV is required. Complete [Environment Setup](#6-environment-setup) before running the commands below.

## 2. Train HAR+CIS Weights

Run `har-weight-train.R` in `har_cis` mode using the `test_geno` PLINK prefix:

```bash
Rscript har-weight-train.R \
  geno_data/test_geno \
  hars.bed \
  annotations/gencode.v26.annotation.gtf.gz \
  geno_data/Brain_Frontal_Cortex_BA9.v8.normalized_expression.chr1_test20_v2.bed \
  1 \
  har_cis
```

The arguments are positional:

| Argument | Description |
| --- | --- |
| `plink_prefix` | Common path prefix for the PLINK `.bed`, `.bim`, and `.fam` files. |
| `har_bed_file` | Tab-separated HAR intervals with four columns: chromosome, start, end, and HAR ID. |
| `gene_annot_file` | GENCODE GTF annotation; the supplied `.gtf.gz` can be read directly. |
| `gene_expr_file` | Tab-separated expression BED with chromosome, start, end, and versioned gene ID in the first four columns, followed by sample columns matching genotype IIDs. |
| `chrom` | Chromosome to analyze, such as `1` or `chr1`. |
| `mode` | Optional: `har_only` (default), `har_cis`, or `cis_only`. |

### Training outputs

For the HAR+CIS example, training writes:

- `output_results_seed/model_results_Brain_Frontal_Cortex_BA9.chr1_test20_v2_with_cis.txt`
- `output_results_seed/snp_weights_Brain_Frontal_Cortex_BA9.chr1_test20_v2_with_cis.txt`
- `output_results_seed/predicted_expression_Brain_Frontal_Cortex_BA9.chr1_test20_v2_with_cis.txt`

| File | Contents |
| --- | --- |
| `model_results_*.txt` | Retained genes, gene names, out-of-fold squared correlation (`R2`), selected lambda (`BestLambda`), and nonzero SNP coefficient count (`NumSNPs`). |
| `snp_weights_*.txt` | Coefficients fitted on all training samples at the selected lambda, with SNP alleles, gene IDs, and mode-specific annotations. Intercept rows can occur here; the database converter excludes them. |
| `predicted_expression_*.txt` | Observed expression and held-out elastic-net predictions for each retained gene and sample. |


## 3. Create Model Database

Create the MetaXcan-compatible SQLite DB from the with_cis weights and model results:

```bash
python3 software_deps/create_db.py \
  --weights output_results_seed/snp_weights_Brain_Frontal_Cortex_BA9.chr1_test20_v2_with_cis.txt \
  --results output_results_seed/model_results_Brain_Frontal_Cortex_BA9.chr1_test20_v2_with_cis.txt \
  --expr geno_data/Brain_Frontal_Cortex_BA9.v8.normalized_expression.chr1_test20_v2.bed \
  --out database/Brain_Frontal_Cortex_BA9_with_cis.db
```

The database contains:

- `weights`: nonzero SNP coefficients and their effect/reference alleles.
- `extra`: retained gene metadata and prediction performance.
- `sample_info`: one row containing `n_samples`. With `--expr`, this value is the number of expression sample columns (175 in the test expression file); the converter does not independently count the genotype–expression intersection.

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

The output has columns `GENE`, `RSID1`, `RSID2`, and `VALUE`, with upper-triangular SNP covariance entries, including the diagonal, for each gene. Its row count depends on the newly retained model SNPs.

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

The output has one header row followed by gene association results. The number of results depends on retained models and SNP overlap with the GWAS.

The wrapper currently recognizes the PFC database names `Brain_Frontal_Cortex_BA9.db` and `Brain_Frontal_Cortex_BA9_with_cis.db`. Other tissue or mode database names require extending the wrapper's tissue mapping or invoking S-PrediXcan directly.

For another GWAS file, place the file in `sumstats/` and keep the same DB selector. Replace the uppercase placeholders below with your trait label, file pattern, and column names:

```bash
PYTHON_BIN=python3 HARDAM_DB_GLOB=database/Brain_Frontal_Cortex_BA9_with_cis.db \
  bash exec_HARDAM.sh \
  TRAIT_NAME \
  'GWAS_FILE_PATTERN.tsv.gz' \
  SNP_COLUMN EFFECT_ALLELE_COLUMN NON_EFFECT_ALLELE_COLUMN BETA_COLUMN SE_COLUMN \
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

### PLINK

HARDAM calls the PLINK 1.9 executable as `plink`. The commands below download the 7 September 2026 build from the [official PLINK downloads](https://www.cog-genomics.org/plink/1.9/), extract the executable, add it to `PATH` immediately and in future shell sessions, and print its version. They require `wget` and `unzip` and do not require administrator privileges.

**Linux x86-64 with Bash** — paste this as one command:

```bash
mkdir -p "$HOME/.local/bin/plink-1.9" && \
  wget -O "$HOME/.local/bin/plink-1.9/plink.zip" https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20260907.zip && \
  unzip -o "$HOME/.local/bin/plink-1.9/plink.zip" plink -d "$HOME/.local/bin/plink-1.9" && \
  chmod +x "$HOME/.local/bin/plink-1.9/plink" && \
  export PATH="$HOME/.local/bin/plink-1.9:$PATH" && \
  plink --version && \
  printf '\nexport PATH="$HOME/.local/bin/plink-1.9:$PATH"\n' >> "$HOME/.bashrc"
```

**macOS with Zsh** — use the macOS archive and shell configuration instead:

```zsh
mkdir -p "$HOME/.local/bin/plink-1.9" && \
  wget -O "$HOME/.local/bin/plink-1.9/plink.zip" https://s3.amazonaws.com/plink1-assets/plink_mac_20260907.zip && \
  unzip -o "$HOME/.local/bin/plink-1.9/plink.zip" plink -d "$HOME/.local/bin/plink-1.9" && \
  chmod +x "$HOME/.local/bin/plink-1.9/plink" && \
  export PATH="$HOME/.local/bin/plink-1.9:$PATH" && \
  plink --version && \
  printf '\nexport PATH="$HOME/.local/bin/plink-1.9:$PATH"\n' >> "$HOME/.zshrc"
```

Run HARDAM from the same terminal after installation. On cluster jobs, also include `export PATH="$HOME/.local/bin/plink-1.9:$PATH"` in the job script, since batch shells may not load your interactive shell configuration. For other platforms or architecture compatibility, consult the official download page.

The previously checked local PLINK version was 1.90b7.2; the newer download above has not been validated with a complete HARDAM run. Record the version printed by `plink --version` with your analysis. The database inspection commands use the `sqlite3` command-line tool.


## 7. References

- [MetaXcan Wiki](https://github.com/hakyimlab/MetaXcan/wiki)
- [GENCODE Human Release 26](https://www.gencodegenes.org/human/release_26.html)

## 8. Citation

If you use HARDAM in your research, please cite:

> Enoma, D. O., Wang, D., Weeraman, J., Gordon, P. M. K., de Koning, A. P. J., Cao, B., Long, Q., & Cao, C. (2026). *Human accelerated regions-directed association mapping discovers trans-regulatory associations in brain disorders*. Manuscript.

This citation describes the manuscript; journal, volume, pages, and DOI will be added when available.

```bibtex
@unpublished{enoma2026hardam,
  title = {Human accelerated regions-directed association mapping discovers trans-regulatory associations in brain disorders},
  author = {Enoma, David O. and Wang, Dinghao and Weeraman, Janith and Gordon, Paul M. K. and de Koning, A. P. Jason and Cao, Bo and Long, Quan and Cao, Chen},
  year = {2026},
  note = {Manuscript}
}
```

Machine-readable citation metadata are provided in [CITATION.cff](CITATION.cff). For reproducibility, also record the software release or Git commit used in your analysis:

```bash
git rev-parse HEAD
```

Describe any local modifications used for the analysis. Please also cite the upstream methods and data resources used in your workflow, including S-PrediXcan/MetaXcan, GTEx, GENCODE, and the relevant GWAS studies.

## 9. License and Credits

HARDAM's original source code and documentation are distributed under the [MIT License](LICENSE), copyright © 2026 David O. Enoma and HARDAM contributors.


Vendored MetaXcan/S-PrediXcan code retains its upstream MIT license and attribution: copyright © 2015 hakyimlab, with software mostly written by heroico. The [included MetaXcan license](software_deps/MetaXcan/LICENSE) reproduces the [upstream notice](https://github.com/hakyimlab/MetaXcan/blob/master/LICENSE). Other dependencies retain their respective licenses.


Last updated: September 8, 2026
