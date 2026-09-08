#!/usr/bin/env bash
set -euo pipefail

# -------------------------------------------------------------------
# exec_HARDAM.sh

#
# Usage:
#   ./exec_HARDAM.sh <DISEASE> <GWAS_PATTERN> <SNP_COL> <EFF_ALLELE> \
#                        <NON_EFF_ALLELE> <BETA_COL> <SE_COL> [--additional_output]
# -------------------------------------------------------------------

if [[ $# -lt 7 ]]; then
  echo "Usage: $0 <DISEASE> <GWAS_PATTERN> <SNP_COL> <EFF_ALLELE> <NON_EFF_ALLELE> <BETA_COL> <SE_COL> [--additional_output]"
  exit 1
fi

DISEASE=$1
GWAS_PATTERN=$2
SNP_COL=$3
EFF_ALLELE=$4
NON_EFF=$5
BETA_COL=$6
SE_COL=$7
shift 7
EXTRA_ARGS="$*"

GWAS_FOLDER="sumstats"
DB_FOLDER="database"
COV_FOLDER="covariances"
OUT_FOLDER="TWAS_results"
DB_GLOB=${HARDAM_DB_GLOB:-"$DB_FOLDER/*.db"}
PYTHON_BIN=${PYTHON_BIN:-python3.9}
mkdir -p "$OUT_FOLDER"

for DB_PATH in $DB_GLOB; do
  tissue=$(basename "$DB_PATH" .db)
  case "$tissue" in
    Brain_Frontal_Cortex_BA9)
      COV_FILE="snp_weights_Brain_Frontal_Cortex_BA9_cov.tsv.gz"
      SUFFIX="PFC"
      MODEL_LABEL="HAR only" ;;
    Brain_Frontal_Cortex_BA9_with_cis)
      COV_FILE="snp_weights_Brain_Frontal_Cortex_BA9_with_cis_cov.tsv.gz"
      SUFFIX="PFC_with_cis"
      MODEL_LABEL="HAR+CIS" ;;
    *)
      echo "Skipping unknown tissue file: $DB_PATH" >&2
      continue ;;
  esac

  OUT_FILE="$OUT_FOLDER/${DISEASE}_${SUFFIX}.csv"

  echo ">>> Running S-PrediXcan for $DISEASE on $tissue ($MODEL_LABEL)"
  "$PYTHON_BIN" software_deps/MetaXcan/software/SPrediXcan.py \
    --model_db_path            "$DB_PATH" \
    --covariance               "$COV_FOLDER/$COV_FILE" \
    --gwas_folder              "$GWAS_FOLDER" \
    --gwas_file_pattern        "$GWAS_PATTERN" \
    --snp_column               "$SNP_COL" \
    --effect_allele_column     "$EFF_ALLELE" \
    --non_effect_allele_column "$NON_EFF" \
    --beta_column              "$BETA_COL" \
    --se_column                "$SE_COL" \
    --output_file              "$OUT_FILE" \
    $EXTRA_ARGS

  echo
done

echo "All done - outputs in $OUT_FOLDER"
