#!/usr/bin/env python3
"""Plot perturbation QC metrics."""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parents[1]
IN_FILE = BASE / "data/processed/perturbation_table_mvp_q1.tsv"
FIG_DIR = BASE / "results/figures"
TAB_DIR = BASE / "results/tables"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TAB_DIR.mkdir(parents=True, exist_ok=True)

ORDER = [
    "gc_matched_replacement",
    "minimal_pwm_edit",
    "motif_preserving_flank_control",
    "random_replacement",
]


def boxplot_by_strategy(df, col, ylabel, title, filename):
    data = [df.loc[df["strategy"] == s, col].values for s in ORDER]
    plt.figure(figsize=(8, 4.8))
    plt.boxplot(data, labels=ORDER, showfliers=True)
    plt.xticks(rotation=30, ha="right")
    plt.ylabel(ylabel)
    plt.xlabel("Perturbation strategy")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(FIG_DIR / filename, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    df = pd.read_csv(IN_FILE, sep="\t")

    summary = (
        df.groupby("strategy")
        .agg(
            n=("perturb_id", "size"),
            mean_delta_gc=("delta_gc", "mean"),
            mean_delta_pwm_rel=("delta_pwm_rel", "mean"),
            mean_edit_distance=("edit_distance", "mean"),
            mean_local_dinuc_jsd=("local_dinuc_jsd", "mean"),
            off_target_rate=("off_target_motif_created", "mean"),
        )
        .reset_index()
    )
    summary.to_csv(TAB_DIR / "qc_summary_by_strategy.tsv", sep="\t", index=False)
    print(summary)

    boxplot_by_strategy(df, "delta_gc", "Delta GC", "Delta GC by perturbation strategy", "fig1_delta_gc_by_strategy.png")
    boxplot_by_strategy(df, "delta_pwm_rel", "Delta relative PWM score", "Delta relative PWM score by perturbation strategy", "fig2_delta_pwm_rel_by_strategy.png")
    boxplot_by_strategy(df, "edit_distance", "Edit distance in target motif", "Edit distance in target motif by perturbation strategy", "fig3_edit_distance_by_strategy.png")
    boxplot_by_strategy(df, "local_dinuc_jsd", "Local dinucleotide JSD", "Local dinucleotide JSD by perturbation strategy", "fig4_local_dinuc_jsd_by_strategy.png")

    off = summary.set_index("strategy").loc[ORDER].reset_index()
    plt.figure(figsize=(8, 4.5))
    plt.bar(off["strategy"], off["off_target_rate"])
    plt.xticks(rotation=30, ha="right")
    plt.ylabel("Off-target motif creation rate")
    plt.xlabel("Perturbation strategy")
    plt.title("Off-target motif creation by strategy")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig5_off_target_creation_rate.png", dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Figures saved to {FIG_DIR}")


if __name__ == "__main__":
    main()
