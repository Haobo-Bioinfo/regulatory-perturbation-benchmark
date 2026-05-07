#!/usr/bin/env python3
"""Fetch hg38 DNA sequences for prepared regions using the UCSC REST API."""
from pathlib import Path
import time
import requests
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "data/processed"
IN_FILE = OUT / "regions_k562_mvp.tsv"
OUT_FILE = OUT / "regions_k562_mvp_with_sequence.tsv"


def fetch_ucsc_sequence(chrom, start, end, genome="hg38", max_retries=3):
    url = (
        "https://api.genome.ucsc.edu/getData/sequence?"
        f"genome={genome};chrom={chrom};start={int(start)};end={int(end)}"
    )
    last_err = None
    for attempt in range(max_retries):
        try:
            r = requests.get(url, timeout=30)
            r.raise_for_status()
            data = r.json()
            seq = data.get("dna", "").upper()
            if not seq:
                raise RuntimeError(f"No DNA returned for {chrom}:{start}-{end}")
            return seq
        except Exception as err:
            last_err = err
            time.sleep(1 + attempt)
    raise RuntimeError(f"Failed to fetch {chrom}:{start}-{end}: {last_err}")


def main():
    df = pd.read_csv(IN_FILE, sep="\t")
    seqs = []
    for i, row in df.iterrows():
        if i % 25 == 0:
            print(f"Fetching {i}/{len(df)}", flush=True)
        seqs.append(fetch_ucsc_sequence(row["chr"], row["start"], row["end"]))
        time.sleep(0.03)

    df["sequence_ref"] = seqs
    df.to_csv(OUT_FILE, sep="\t", index=False)
    print(df["sequence_ref"].str.len().describe())
    print(f"Saved {OUT_FILE}")


if __name__ == "__main__":
    main()
