#!/usr/bin/env python3
"""Download reference input files for the perturbation benchmark pilot."""
from pathlib import Path
import requests

BASE = Path(__file__).resolve().parents[1]
CCRE_DIR = BASE / "data/raw/ccre"
PWM_DIR = BASE / "data/raw/jaspar_pwm"
CCRE_DIR.mkdir(parents=True, exist_ok=True)
PWM_DIR.mkdir(parents=True, exist_ok=True)

CCRE_URLS = {
    "pls_hg38.bed": "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.PLS.bed",
    "pels_hg38.bed": "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.pELS.bed",
}

JASPAR_IDS = {
    "GATA1_MA0035.5.pfm": "MA0035.5",
    "TAL1_MA0091.1.pfm": "MA0091.1",
    "SPI1_MA0080.1.pfm": "MA0080.1",
}


def download_text(url, out_path):
    out_path = Path(out_path)
    if out_path.exists() and out_path.stat().st_size > 0:
        print(f"Exists, skipping: {out_path}")
        return
    print(f"Downloading {url}")
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    out_path.write_text(r.text, encoding="utf-8")
    print(f"Saved: {out_path}")


def main():
    for name, url in CCRE_URLS.items():
        download_text(url, CCRE_DIR / name)

    for filename, matrix_id in JASPAR_IDS.items():
        url = f"https://jaspar.elixir.no/api/v1/matrix/{matrix_id}.pfm"
        download_text(url, PWM_DIR / filename)

    print("Done.")


if __name__ == "__main__":
    main()
