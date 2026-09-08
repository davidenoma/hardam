#!/usr/bin/env python3
import argparse
import gzip
import os
import re
import subprocess
import sys
import tempfile
from datetime import datetime

import numpy as np
import pandas as pd


def now():
    return datetime.now().strftime("%H:%M:%S")


def die(message):
    sys.stderr.write(message + "\n")
    sys.exit(1)


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def read_keep_ids(expr_bed):
    with open_text(expr_bed) as handle:
        header = handle.readline().rstrip("\n").split("\t")
    if len(header) < 5:
        die(f"Expression BED has fewer than 5 columns: {expr_bed}")
    return header[4:]


def write_keep_file(sample_ids, tmp_dir):
    path = os.path.join(tmp_dir, "keep_samples.txt")
    with open(path, "w") as handle:
        for iid in sample_ids:
            handle.write(f"{iid}\t{iid}\n")
    return path


def load_bim(plink_prefix):
    bim = pd.read_csv(
        plink_prefix + ".bim",
        sep=r"\s+",
        header=None,
        usecols=[0, 1, 3],
        names=["chr", "snp", "pos"],
        dtype={"chr": str, "snp": str, "pos": int},
    )
    bim["chr"] = bim["chr"].str.replace(r"^chr", "", regex=True)
    return bim.drop_duplicates(subset=["snp"], keep="first")


def load_weight_snps(path):
    weights = pd.read_csv(path, sep="\t", dtype=str, usecols=["Gene", "SNP", "Weight"])
    weights = weights[weights["SNP"].notna() & (weights["SNP"] != "(Intercept)")].copy()
    weights["Weight"] = pd.to_numeric(weights["Weight"], errors="coerce")
    weights = weights[weights["Weight"].notna() & (weights["Weight"] != 0)]
    if weights.empty:
        die("No nonzero SNP weights found")
    return weights[["Gene", "SNP"]].drop_duplicates()


def compute_covariance(raw_file):
    df = pd.read_csv(raw_file, sep=r"\s+")
    for col in ("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"):
        if col in df.columns:
            df = df.drop(columns=col)
    if df.empty:
        return [], None

    snps = [re.sub(r"_[ACGTacgt]+$", "", col) for col in df.columns]
    x = df.to_numpy(dtype=float)
    x -= np.nanmean(x, axis=0)
    if x.shape[1] == 1:
        cov = np.array([[float(np.nanvar(x[:, 0], ddof=1))]])
    else:
        cov = np.cov(x, rowvar=False, bias=False)
    return snps, cov


def write_flat_covariance(gene, snps, cov, handle):
    seen = set()
    for i, snp_i in enumerate(snps):
        for j in range(i, len(snps)):
            snp_j = snps[j]
            key = (snp_i, snp_j)
            if key in seen:
                continue
            handle.write(f"{gene}\t{snp_i}\t{snp_j}\t{cov[i, j]:.8g}\n")
            seen.add(key)


def run_plink_recode(plink_bin, plink_prefix, keep_file, snp_file, out_prefix):
    return subprocess.run(
        [
            plink_bin,
            "--bfile",
            plink_prefix,
            "--keep",
            keep_file,
            "--extract",
            snp_file,
            "--recode",
            "A",
            "--out",
            out_prefix,
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )


def main():
    parser = argparse.ArgumentParser(
        description="Create MetaXcan long-form SNP covariance from a PLINK prefix and model weights."
    )
    parser.add_argument("--plink", required=True, help="PLINK prefix")
    parser.add_argument("--weights", "--wgt", dest="weights", required=True, help="Weights TSV from har-weight-train.R")
    parser.add_argument("--expr", required=True, help="Expression BED/BED.GZ used for training")
    parser.add_argument("--out", required=True, help="Output covariance TSV/TSV.GZ")
    parser.add_argument("--plink-bin", default="plink", help="PLINK executable")
    parser.add_argument("--no-gzip", action="store_true", help="Write plain text even if --out ends in .gz")
    parser.add_argument("--gtf", help="Accepted for backward compatibility; not needed for model-SNP covariance")
    parser.add_argument("--window", type=int, default=500000, help="Accepted for backward compatibility")
    args = parser.parse_args()

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    model_snps = load_weight_snps(args.weights)
    bim = load_bim(args.plink)
    model_snps = model_snps[model_snps["SNP"].isin(set(bim["snp"]))]
    if model_snps.empty:
        die("No weighted SNPs were present in the PLINK BIM")

    opener = gzip.open if (args.out.endswith(".gz") and not args.no_gzip) else open

    with tempfile.TemporaryDirectory(prefix="hardam_cov_") as tmp_dir:
        keep_file = write_keep_file(read_keep_ids(args.expr), tmp_dir)
        with opener(args.out, "wt") as out_handle:
            out_handle.write("GENE\tRSID1\tRSID2\tVALUE\n")
            for gene, gene_df in model_snps.groupby("Gene", sort=True):
                snps = gene_df["SNP"].drop_duplicates().tolist()
                snp_file = os.path.join(tmp_dir, f"{re.sub(r'[^A-Za-z0-9_.-]', '_', gene)}.snps")
                with open(snp_file, "w") as handle:
                    handle.write("\n".join(snps) + "\n")

                out_prefix = os.path.join(tmp_dir, f"plink_{re.sub(r'[^A-Za-z0-9_.-]', '_', gene)}")
                result = run_plink_recode(args.plink_bin, args.plink, keep_file, snp_file, out_prefix)
                raw_file = out_prefix + ".raw"
                if result.returncode != 0 or not os.path.exists(raw_file):
                    print(f"[{now()}] {gene}: PLINK failed", file=sys.stderr)
                    continue

                raw_snps, cov = compute_covariance(raw_file)
                if cov is None or len(raw_snps) == 0:
                    print(f"[{now()}] {gene}: no covariance written", file=sys.stderr)
                    continue
                write_flat_covariance(gene, raw_snps, cov, out_handle)
                print(f"[{now()}] {gene}: {len(raw_snps)} SNPs")

    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
