#!/usr/bin/env python3
import argparse
import gzip
import math
import os
import sqlite3
import sys

import pandas as pd

try:
    from scipy.stats import t as student_t
except Exception:  # pragma: no cover - optional runtime dependency
    student_t = None


def die(message):
    sys.stderr.write(message + "\n")
    sys.exit(1)


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def count_expression_samples(expr_bed):
    with open_text(expr_bed) as handle:
        header = handle.readline().rstrip("\n").split("\t")
    if len(header) < 5:
        die(f"Expression BED has fewer than 5 columns: {expr_bed}")
    return len(header) - 4


def bh_adjust(pvalues):
    indexed = [(i, p) for i, p in enumerate(pvalues) if p is not None and not math.isnan(p)]
    adjusted = [None] * len(pvalues)
    if not indexed:
        return adjusted

    indexed.sort(key=lambda x: x[1])
    m = len(indexed)
    running = 1.0
    for rank_from_end, (idx, pvalue) in enumerate(reversed(indexed), start=1):
        rank = m - rank_from_end + 1
        running = min(running, pvalue * m / rank)
        adjusted[idx] = min(running, 1.0)
    return adjusted


def r2_pvalue(r2, n_samples):
    if n_samples is None or n_samples <= 2 or r2 is None or math.isnan(r2):
        return None
    r2 = max(0.0, min(float(r2), 1.0))
    if r2 >= 1.0:
        return 0.0
    if student_t is None:
        return None
    r = math.sqrt(r2)
    t_value = r * math.sqrt((n_samples - 2) / max(1e-300, 1 - r2))
    return float(2 * student_t.sf(abs(t_value), df=n_samples - 2))


def load_weights(path):
    weights = pd.read_csv(path, sep="\t", dtype=str)
    required = {"SNP", "Weight", "Gene", "allele1", "allele2"}
    missing = required - set(weights.columns)
    if missing:
        die(f"Weights file is missing required columns: {', '.join(sorted(missing))}")

    weights = weights[weights["SNP"].notna() & (weights["SNP"] != "(Intercept)")].copy()
    weights["Weight"] = pd.to_numeric(weights["Weight"], errors="coerce")
    weights = weights[weights["Weight"].notna() & (weights["Weight"] != 0)]
    weights = weights[weights["Gene"].notna()]
    weights = weights.drop_duplicates(subset=["SNP", "Gene"], keep="first")

    # har-weight-train.R fits weights on PLINK --recode A dosages, which count BIM allele1.
    db_weights = pd.DataFrame(
        {
            "rsid": weights["SNP"],
            "gene": weights["Gene"],
            "weight": weights["Weight"],
            "ref_allele": weights["allele2"],
            "eff_allele": weights["allele1"],
        }
    )
    return db_weights


def load_extra(path, n_samples):
    results = pd.read_csv(path, sep="\t")
    required = {"Gene", "GeneName", "R2", "NumSNPs"}
    missing = required - set(results.columns)
    if missing:
        die(f"Results file is missing required columns: {', '.join(sorted(missing))}")

    results["R2"] = pd.to_numeric(results["R2"], errors="coerce")
    results["NumSNPs"] = pd.to_numeric(results["NumSNPs"], errors="coerce").fillna(0).astype(int)
    results = results.drop_duplicates(subset=["Gene"], keep="first").copy()
    pvalues = [r2_pvalue(r2, n_samples) for r2 in results["R2"]]
    qvalues = bh_adjust(pvalues)

    return pd.DataFrame(
        {
            "gene": results["Gene"],
            "genename": results["GeneName"],
            "pred.perf.R2": results["R2"],
            "n.snps.in.model": results["NumSNPs"],
            "pred.perf.pval": pvalues,
            "pred.perf.qval": qvalues,
        }
    )


def write_db(out_path, weights, extra, n_samples):
    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    if os.path.exists(out_path):
        os.remove(out_path)

    with sqlite3.connect(out_path) as conn:
        extra.to_sql("extra", conn, index=False)
        weights.to_sql("weights", conn, index=False)
        pd.DataFrame({"n_samples": [n_samples if n_samples is not None else 0]}).to_sql(
            "sample_info", conn, index=False
        )
        conn.execute("CREATE INDEX extra_gene ON extra(gene)")
        conn.execute("CREATE INDEX weights_rsid ON weights(rsid)")
        conn.execute("CREATE INDEX weights_gene ON weights(gene)")
        conn.execute("CREATE INDEX weights_rsid_gene ON weights(rsid, gene)")


def main():
    parser = argparse.ArgumentParser(description="Create a MetaXcan-compatible SQLite model DB.")
    parser.add_argument("--weights", "--wgt", dest="weights", required=True, help="SNP weights TSV from har-weight-train.R")
    parser.add_argument("--results", required=True, help="Model results TSV from har-weight-train.R")
    parser.add_argument("--out", required=True, help="Output SQLite DB path")
    parser.add_argument("--expr", help="Expression BED/BED.GZ used for training; used to count samples")
    parser.add_argument("--n-samples", type=int, help="Training sample count if --expr is not provided")
    args = parser.parse_args()

    n_samples = args.n_samples
    if args.expr:
        n_samples = count_expression_samples(args.expr)

    weights = load_weights(args.weights)
    extra = load_extra(args.results, n_samples)
    write_db(args.out, weights, extra, n_samples)

    print(f"wrote {args.out}")
    print(f"weights: {len(weights)} rows")
    print(f"genes: {len(extra)} rows")
    print(f"n_samples: {n_samples if n_samples is not None else 0}")


if __name__ == "__main__":
    main()
