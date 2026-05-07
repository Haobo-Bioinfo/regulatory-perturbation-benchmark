"""Utility functions for motif scanning and perturbation design.

The functions in this file are intentionally lightweight and dependency-minimal.
They are designed for a small pilot benchmark rather than genome-scale motif
analysis.
"""

from __future__ import annotations

import math
import random as random_module
from itertools import combinations, product
from collections import Counter
import numpy as np

BASES = "ACGT"
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1].upper()


def parse_jaspar_pfm(path):
    """Parse a JASPAR PFM file.

    Supports both labeled rows:
        A [ 1 2 3 ]
        C [ ... ]
    and unlabeled four-row format returned by the JASPAR API, where rows are
    assumed to be A, C, G, T.
    """
    path = str(path)
    raw_lines = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(">"):  # skip header
                continue
            raw_lines.append(line)

    if not raw_lines:
        raise ValueError(f"No matrix rows found in {path}")

    counts = {}
    # Labeled format
    if raw_lines[0][0].upper() in BASES and ("[" in raw_lines[0] or raw_lines[0].split()[0].upper() in BASES):
        for line in raw_lines:
            base = line[0].upper()
            if base not in BASES:
                continue
            nums = (
                line.replace("[", " ")
                .replace("]", " ")
                .replace(base, " ", 1)
                .split()
            )
            counts[base] = [float(x) for x in nums]
    else:
        # Unlabeled four-row format, assumed A/C/G/T.
        if len(raw_lines) != 4:
            raise ValueError(f"Expected four unlabeled PFM rows in {path}, got {len(raw_lines)}")
        for base, line in zip(BASES, raw_lines):
            counts[base] = [float(x) for x in line.split()]

    if set(counts) != set(BASES):
        raise ValueError(f"PFM must contain A/C/G/T rows. Got {counts.keys()} in {path}")

    lengths = {len(v) for v in counts.values()}
    if len(lengths) != 1:
        raise ValueError(f"PFM rows have inconsistent lengths in {path}: {lengths}")

    return counts


def pfm_to_pwm(pfm, pseudocount=0.8, bg=None):
    if bg is None:
        bg = {b: 0.25 for b in BASES}
    L = len(pfm["A"])
    pwm = {b: [] for b in BASES}
    for i in range(L):
        col_sum = sum(pfm[b][i] for b in BASES) + pseudocount * 4
        for b in BASES:
            p = (pfm[b][i] + pseudocount) / col_sum
            pwm[b].append(math.log2(p / bg[b]))
    return pwm


def pwm_length(pwm) -> int:
    return len(pwm["A"])


def score_seq_pwm(seq: str, pwm) -> float:
    seq = seq.upper()
    if len(seq) != pwm_length(pwm):
        raise ValueError("Sequence length and PWM length do not match")
    score = 0.0
    for i, b in enumerate(seq):
        if b not in BASES:
            return -1e9
        score += pwm[b][i]
    return float(score)


def max_pwm_score(pwm) -> float:
    L = pwm_length(pwm)
    return float(sum(max(pwm[b][i] for b in BASES) for i in range(L)))


def min_pwm_score(pwm) -> float:
    L = pwm_length(pwm)
    return float(sum(min(pwm[b][i] for b in BASES) for i in range(L)))


def rel_pwm_score(seq: str, pwm) -> float:
    """Relative PWM score scaled to [0, 1] using min and max PWM scores."""
    s = score_seq_pwm(seq, pwm)
    smin = min_pwm_score(pwm)
    smax = max_pwm_score(pwm)
    if smax == smin:
        return 0.0
    return float((s - smin) / (smax - smin))


def scan_pwm_both_strands(seq: str, pwm, min_rel_score=0.80):
    """Scan a sequence on both strands.

    Returns tuples: (start, end, strand, motif_seq_oriented, raw_score, rel_score)
    The returned motif sequence is oriented relative to the PWM. For negative
    strand hits, it is the reverse-complement of the genomic sequence window.
    """
    seq = seq.upper()
    L = pwm_length(pwm)
    hits = []
    for i in range(0, len(seq) - L + 1):
        genomic_sub = seq[i : i + L]
        if any(b not in BASES for b in genomic_sub):
            continue
        for strand, oriented in [('+', genomic_sub), ('-', revcomp(genomic_sub))]:
            rel = rel_pwm_score(oriented, pwm)
            if rel >= min_rel_score:
                hits.append((i, i + L, strand, oriented, score_seq_pwm(oriented, pwm), rel))
    return hits


def gc_content(seq: str) -> float:
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / max(1, len(seq))


def edit_distance_same_length(a: str, b: str) -> int:
    if len(a) != len(b):
        raise ValueError("Edit distance helper expects equal-length strings")
    return sum(x != y for x, y in zip(a.upper(), b.upper()))


def random_replacement(seq: str, rng=random_module) -> str:
    return "".join(rng.choice(BASES) for _ in seq)


def gc_matched_replacement(seq: str, rng=random_module) -> str:
    seq = seq.upper()
    n_gc = seq.count("G") + seq.count("C")
    n_at = len(seq) - n_gc
    letters = [rng.choice("GC") for _ in range(n_gc)] + [rng.choice("AT") for _ in range(n_at)]
    rng.shuffle(letters)
    return "".join(letters)


def minimal_pwm_edit(seq: str, pwm, target_rel_score=0.35, max_edits=2):
    """Find a small edit that reduces relative PWM score.

    Returns (mutant_sequence, number_of_edits, edited_positions).
    This version is intentionally simple and does not constrain GC content.
    """
    seq = seq.upper()
    best_seq = seq
    best_score = rel_pwm_score(seq, pwm)
    best_positions = []

    for k in range(1, max_edits + 1):
        for positions in combinations(range(len(seq)), k):
            choices = [[b for b in BASES if b != seq[pos]] for pos in positions]
            for repls in product(*choices):
                arr = list(seq)
                for pos, b in zip(positions, repls):
                    arr[pos] = b
                mut = "".join(arr)
                rel = rel_pwm_score(mut, pwm)
                if rel < best_score:
                    best_seq = mut
                    best_score = rel
                    best_positions = list(positions)
                if rel <= target_rel_score:
                    return mut, k, list(positions)
    return best_seq, len(best_positions), best_positions


def mutate_region_sequence(region_seq: str, motif_start: int, motif_end: int, motif_seq_mut_genomic: str) -> str:
    return region_seq[:motif_start] + motif_seq_mut_genomic + region_seq[motif_end:]


def motif_preserving_flank_control(region_seq: str, motif_start: int, motif_end: int, flank=10, n_edits=2, rng=random_module):
    """Mutate flanking sequence while leaving the target motif unchanged."""
    region_seq = region_seq.upper()
    candidate_positions = list(range(max(0, motif_start - flank), motif_start)) + list(
        range(motif_end, min(len(region_seq), motif_end + flank))
    )
    candidate_positions = [p for p in candidate_positions if region_seq[p] in BASES]
    if not candidate_positions:
        return region_seq, []

    n_edits = min(n_edits, len(candidate_positions))
    positions = rng.sample(candidate_positions, n_edits)
    arr = list(region_seq)
    for p in positions:
        choices = [b for b in BASES if b != arr[p]]
        arr[p] = rng.choice(choices)
    return "".join(arr), positions


def kmer_counts(seq: str, k=2):
    seq = seq.upper()
    counts = Counter()
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if all(b in BASES for b in kmer):
            counts[kmer] += 1
    return counts


def js_divergence_from_counts(c1, c2):
    alphabet = sorted(set(c1) | set(c2))
    if not alphabet:
        return 0.0
    p = np.array([c1.get(a, 0) for a in alphabet], dtype=float)
    q = np.array([c2.get(a, 0) for a in alphabet], dtype=float)
    p = p / p.sum() if p.sum() > 0 else np.ones_like(p) / len(p)
    q = q / q.sum() if q.sum() > 0 else np.ones_like(q) / len(q)
    m = 0.5 * (p + q)

    def kl(a, b):
        mask = a > 0
        return float(np.sum(a[mask] * np.log2(a[mask] / b[mask])))

    return 0.5 * kl(p, m) + 0.5 * kl(q, m)


def local_kmer_jsd(seq_ref: str, seq_mut: str, motif_start: int, motif_end: int, flank=30, k=2) -> float:
    left = max(0, motif_start - flank)
    right = min(len(seq_ref), motif_end + flank)
    return js_divergence_from_counts(kmer_counts(seq_ref[left:right], k), kmer_counts(seq_mut[left:right], k))


def check_off_target_motif_creation(seq_ref: str, seq_mut: str, motif_start: int, motif_end: int, pwms, min_rel_score=0.85, flank=50):
    """Flag new high-scoring motif hits in a local window after perturbation."""
    left = max(0, motif_start - flank)
    right = min(len(seq_ref), motif_end + flank)
    ref_window = seq_ref[left:right].upper()
    mut_window = seq_mut[left:right].upper()

    created_tfs = []
    for tf, pwm in pwms.items():
        ref_hits = scan_pwm_both_strands(ref_window, pwm, min_rel_score)
        mut_hits = scan_pwm_both_strands(mut_window, pwm, min_rel_score)
        ref_set = set((h[0], h[1], h[2], h[3]) for h in ref_hits)
        mut_set = set((h[0], h[1], h[2], h[3]) for h in mut_hits)
        if len(mut_set - ref_set) > 0:
            created_tfs.append(tf)
    created_tfs = sorted(set(created_tfs))
    return len(created_tfs) > 0, ",".join(created_tfs)
