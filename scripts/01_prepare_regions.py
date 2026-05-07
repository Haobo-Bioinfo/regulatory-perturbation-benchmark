#!/usr/bin/env python3
"""Prepare a balanced small set of K562 regulatory regions from cCRE class BED files."""
from pathlib import Path
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
RAW = BASE / "data/raw/ccre"
OUT = BASE / "data/processed"
OUT.mkdir(parents=True, exist_ok=True)

N_PER_CLASS = 150
WINDOW_BP = 1000
RANDOM_STATE = 1


def read_bed(path):
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    if df.shape[1] < 3:
        raise ValueError(f"BED file has fewer than three columns: {path}")
    return df


def main():
    pls_file = RAW / "pls_hg38.bed"
    pels_file = RAW / "pels_hg38.bed"
    if not pls_file.exists() or not pels_file.exists():
        raise FileNotFoundError("Expected data/raw/ccre/pls_hg38.bed and pels_hg38.bed")

    standard_chrs = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

    pls = read_bed(pls_file)
    pels = read_bed(pels_file)
    pls = pls[pls[0].isin(standard_chrs)].copy()
    pels = pels[pels[0].isin(standard_chrs)].copy()

    pls = pls.sample(n=min(N_PER_CLASS, len(pls)), random_state=RANDOM_STATE).copy()
    pels = pels.sample(n=min(N_PER_CLASS, len(pels)), random_state=RANDOM_STATE).copy()
    pls["region_type"] = "PLS"
    pels["region_type"] = "pELS"
    regions = pd.concat([pls, pels], axis=0).reset_index(drop=True)

    out = pd.DataFrame()
    out["region_id"] = [f"region_{i:05d}" for i in range(len(regions))]
    out["chr"] = regions[0]
    start_raw = regions[1].astype(int)
    end_raw = regions[2].astype(int)
    out["center"] = ((start_raw + end_raw) // 2).astype(int)
    half = WINDOW_BP // 2
    out["start"] = (out["center"] - half).clip(lower=0)
    out["end"] = out["start"] + WINDOW_BP
    out["region_type"] = regions["region_type"].values
    out = out[["region_id", "chr", "start", "end", "center", "region_type"]]

    path = OUT / "regions_k562_mvp.tsv"
    out.to_csv(path, sep="\t", index=False)
    print(out.head())
    print(f"Saved {path} n={len(out)}")


if __name__ == "__main__":
    main()
