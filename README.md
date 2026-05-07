# Regulatory perturbation benchmark pilot

This repository contains a small pilot benchmark for comparing in silico motif perturbation strategies in regulatory sequence interpretation.

## Question

Many regulatory sequence models use in silico motif mutation or replacement to estimate regulatory effects. However, the perturbation design itself can introduce non-motif artifacts such as GC shifts, local k-mer distortion, excessive edit distance, or creation of new off-target motifs. This pilot asks:

> Before using model predictions, do different perturbation strategies differ in their sequence-level realism and artifact profiles?

## Current pilot design

- **System**: human hg38 ENCODE SCREEN candidate cis-regulatory elements (cCREs).
- **Regions**: balanced sample of promoter-like cCREs (PLS) and proximal enhancer-like cCREs (pELS).
- **Note**: the current pilot is not yet cell-type-specific. K562- or hematopoietic-context validation will require additional filtering or validation using K562 ATAC/ChIP/CRISPRi evidence.
- **TF motifs**: GATA1, TAL1, SPI1 from JASPAR PWMs.
- **Perturbation strategies**:
  1. random motif replacement;
  2. GC-matched replacement;
  3. minimal PWM-score-reducing edit;
  4. motif-preserving flank control.
- **QC metrics**:
  - GC shift;
  - relative PWM score reduction;
  - edit distance inside target motif;
  - local dinucleotide Jensen-Shannon divergence;
  - off-target motif creation rate.

This is a perturbation-realism audit, not yet a model-prediction benchmark. The next step is to feed the same reference and mutant sequences into pretrained regulatory sequence models such as Enformer or Sei.

## Repository structure

```text
regulatory-perturbation-benchmark/
├── scripts/
│   ├── 00_download_reference_inputs.py
│   ├── 01_prepare_regions.py
│   ├── 02_fetch_sequences_ucsc.py
│   ├── 03_scan_and_perturb_lite.py
│   ├── 04_plot_qc.py
│   └── motif_perturb_utils.py
├── data/
│   ├── raw/
│   └── processed/
└── results/
    ├── figures/
    └── tables/
```

Raw data and generated results are not tracked by git.

## Quick start

```bash
pip install -r requirements.txt

# Download SCREEN cCRE BED files and JASPAR PFMs
python scripts/00_download_reference_inputs.py

# Prepare a balanced hg38 cCRE region panel
python scripts/01_prepare_regions.py

# Fetch hg38 DNA sequences from UCSC
python scripts/02_fetch_sequences_ucsc.py

# Scan motifs and generate perturbations
python scripts/03_scan_and_perturb_lite.py

# Plot QC figures
python scripts/04_plot_qc.py
```

## Notes

- The current script keeps at most 80 high-scoring motif hits for each TF × region-type group. This keeps the pilot balanced and avoids short motifs dominating the analysis.
- The minimal PWM-edit strategy is intentionally simple. In the current pilot, it can introduce GC shifts, which motivates a GC-aware minimal-edit variant in future work.
- The motif-preserving flank control is a negative control: it changes nearby sequence while keeping the target motif intact.
- Earlier local file names may contain `k562_mvp` because this pilot was originally motivated by K562/hematopoietic regulatory biology. In the public repository, the current region set should be interpreted as a lightweight hg38 cCRE panel rather than a K562-specific regulatory element set.

## Main outputs

- `data/processed/motif_hits_mvp.tsv`
- `data/processed/perturbation_table_mvp_q1.tsv`
- `results/tables/qc_summary_by_strategy.tsv`
- `results/figures/fig*_*.png`
