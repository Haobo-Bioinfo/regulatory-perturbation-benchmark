#!/usr/bin/env python3
"""Scan TF motifs and generate balanced perturbation table for the MVP benchmark."""
from pathlib import Path
import random
import pandas as pd
from motif_perturb_utils import (
    parse_jaspar_pfm, pfm_to_pwm, pwm_length, max_pwm_score,
    scan_pwm_both_strands, revcomp, random_replacement, gc_matched_replacement,
    minimal_pwm_edit, motif_preserving_flank_control, mutate_region_sequence,
    score_seq_pwm, rel_pwm_score, gc_content, edit_distance_same_length,
    local_kmer_jsd, check_off_target_motif_creation,
)

BASE = Path(__file__).resolve().parents[1]
PWM_DIR = BASE / "data/raw/jaspar_pwm"
OUT = BASE / "data/processed"
OUT.mkdir(parents=True, exist_ok=True)

MAX_HITS_PER_TF_REGION = 80
MIN_REL_SCORE = 0.80
RANDOM_SEED = 7

MOTIF_FILES = {
    "GATA1": PWM_DIR / "GATA1_MA0035.5.pfm",
    "TAL1": PWM_DIR / "TAL1_MA0091.1.pfm",
    "SPI1": PWM_DIR / "SPI1_MA0080.1.pfm",
}


def main():
    rng = random.Random(RANDOM_SEED)
    pwms = {tf: pfm_to_pwm(parse_jaspar_pfm(path)) for tf, path in MOTIF_FILES.items()}
    print("PWM summary:")
    print({tf: (pwm_length(pwm), max_pwm_score(pwm)) for tf, pwm in pwms.items()})

    regions = pd.read_csv(OUT / "regions_k562_mvp_with_sequence.tsv", sep="\t")
    print(f"Loaded regions: {regions.shape}")

    # Scan motifs
    hits = []
    for idx, row in regions.iterrows():
        if idx % 50 == 0:
            print(f"Scanning region {idx}/{len(regions)}", flush=True)
        seq = row.sequence_ref.upper()
        for tf, pwm in pwms.items():
            for j, (s, e, strand, mseq, score, rel) in enumerate(scan_pwm_both_strands(seq, pwm, MIN_REL_SCORE)):
                hits.append({
                    "hit_id": f"{row.region_id}_{tf}_{j}",
                    "region_id": row.region_id,
                    "tf_name": tf,
                    "chr": row.chr,
                    "motif_start_local": s,
                    "motif_end_local": e,
                    "motif_start_genome": int(row.start) + s,
                    "motif_end_genome": int(row.start) + e,
                    "strand": strand,
                    "motif_seq_ref": mseq,
                    "pwm_score_ref": score,
                    "pwm_score_rel": rel,
                    "region_type": row.region_type,
                })

    motif_hits_all = pd.DataFrame(hits)
    if motif_hits_all.empty:
        raise RuntimeError("No motif hits found. Consider lowering MIN_REL_SCORE.")

    print("All motif hits:", motif_hits_all.shape)
    print(motif_hits_all.groupby(["tf_name", "region_type"]).size())

    # Keep balanced high-scoring motif hits.
    motif_hits = (
        motif_hits_all.sort_values("pwm_score_rel", ascending=False)
        .groupby(["tf_name", "region_type"], group_keys=False)
        .head(MAX_HITS_PER_TF_REGION)
        .reset_index(drop=True)
    )
    motif_hits["hit_id"] = [f"hit_{i:05d}_{row.tf_name}_{row.region_type}" for i, row in motif_hits.iterrows()]
    motif_hits.to_csv(OUT / "motif_hits_mvp.tsv", sep="\t", index=False)
    print("MVP motif hits:", motif_hits.shape)
    print(motif_hits.groupby(["tf_name", "region_type"]).size())

    seqmap = dict(zip(regions.region_id, regions.sequence_ref))
    rows = []

    for idx, hit in motif_hits.iterrows():
        if idx % 50 == 0:
            print(f"Generating perturbations {idx}/{len(motif_hits)}", flush=True)
        tf = hit.tf_name
        pwm = pwms[tf]
        region_seq = seqmap[hit.region_id].upper()
        s = int(hit.motif_start_local)
        e = int(hit.motif_end_local)
        motif_ref = hit.motif_seq_ref.upper()  # oriented relative to PWM

        muts = {
            "random_replacement": random_replacement(motif_ref, rng),
            "gc_matched_replacement": gc_matched_replacement(motif_ref, rng),
        }
        minmut, _, _ = minimal_pwm_edit(motif_ref, pwm, target_rel_score=0.35, max_edits=2)
        muts["minimal_pwm_edit"] = minmut

        for strategy, mut_oriented in muts.items():
            mut_genomic = mut_oriented if hit.strand == "+" else revcomp(mut_oriented)
            seq_mut = mutate_region_sequence(region_seq, s, e, mut_genomic)
            off, off_tfs = check_off_target_motif_creation(region_seq, seq_mut, s, e, pwms, 0.85, 50)

            rows.append({
                "perturb_id": f"{hit.hit_id}__{strategy}",
                "hit_id": hit.hit_id,
                "region_id": hit.region_id,
                "tf_name": tf,
                "region_type": hit.region_type,
                "strategy": strategy,
                "chr": hit.chr,
                "strand": hit.strand,
                "motif_start_local": s,
                "motif_end_local": e,
                "motif_start_genome": hit.motif_start_genome,
                "motif_end_genome": hit.motif_end_genome,
                "motif_seq_ref": motif_ref,
                "motif_seq_mut": mut_oriented,
                "sequence_ref": region_seq,
                "sequence_mut": seq_mut,
                "pwm_score_ref": score_seq_pwm(motif_ref, pwm),
                "pwm_score_mut": score_seq_pwm(mut_oriented, pwm),
                "pwm_rel_ref": rel_pwm_score(motif_ref, pwm),
                "pwm_rel_mut": rel_pwm_score(mut_oriented, pwm),
                "delta_pwm_score": score_seq_pwm(motif_ref, pwm) - score_seq_pwm(mut_oriented, pwm),
                "delta_pwm_rel": rel_pwm_score(motif_ref, pwm) - rel_pwm_score(mut_oriented, pwm),
                "gc_ref": gc_content(motif_ref),
                "gc_mut": gc_content(mut_oriented),
                "delta_gc": gc_content(mut_oriented) - gc_content(motif_ref),
                "edit_distance": edit_distance_same_length(motif_ref, mut_oriented),
                "local_dinuc_jsd": local_kmer_jsd(region_seq, seq_mut, s, e, 30, 2),
                "off_target_motif_created": off,
                "off_target_motif_tfs": off_tfs,
            })

        # Motif-preserving flank control
        flank_seq, _ = motif_preserving_flank_control(region_seq, s, e, flank=10, n_edits=2, rng=rng)
        off, off_tfs = check_off_target_motif_creation(region_seq, flank_seq, s, e, pwms, 0.85, 50)
        rows.append({
            "perturb_id": f"{hit.hit_id}__motif_preserving_flank_control",
            "hit_id": hit.hit_id,
            "region_id": hit.region_id,
            "tf_name": tf,
            "region_type": hit.region_type,
            "strategy": "motif_preserving_flank_control",
            "chr": hit.chr,
            "strand": hit.strand,
            "motif_start_local": s,
            "motif_end_local": e,
            "motif_start_genome": hit.motif_start_genome,
            "motif_end_genome": hit.motif_end_genome,
            "motif_seq_ref": motif_ref,
            "motif_seq_mut": motif_ref,
            "sequence_ref": region_seq,
            "sequence_mut": flank_seq,
            "pwm_score_ref": score_seq_pwm(motif_ref, pwm),
            "pwm_score_mut": score_seq_pwm(motif_ref, pwm),
            "pwm_rel_ref": rel_pwm_score(motif_ref, pwm),
            "pwm_rel_mut": rel_pwm_score(motif_ref, pwm),
            "delta_pwm_score": 0.0,
            "delta_pwm_rel": 0.0,
            "gc_ref": gc_content(motif_ref),
            "gc_mut": gc_content(motif_ref),
            "delta_gc": 0.0,
            "edit_distance": 0,
            "local_dinuc_jsd": local_kmer_jsd(region_seq, flank_seq, s, e, 30, 2),
            "off_target_motif_created": off,
            "off_target_motif_tfs": off_tfs,
        })

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "perturbation_table_mvp_q1.tsv", sep="\t", index=False)
    print("Saved perturbation table:", df.shape)
    print(df.groupby(["strategy", "tf_name"]).size())


if __name__ == "__main__":
    main()
