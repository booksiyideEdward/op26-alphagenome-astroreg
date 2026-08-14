#!/usr/bin/env python3
"""Final OP26 QA: structure-aware shuffled-AlphaGenome stacking null.

This diagnostic holds the frozen cohort, EGrf probabilities, labels, fold
assignments, calibration model and metrics fixed. It changes only the mapping
between AlphaGenome scores and enhancer-gene relations.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from run_targeted_audit import (
    EGRF,
    SCORES,
    dual_group_splits,
    fit_logistic,
    grouped_stack,
    metrics,
)


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
PERMUTATION_ROWS = RESULTS / "permutation_null_iterations.csv"
PERMUTATION_SUMMARY = RESULTS / "permutation_null_summary.json"
SEED = 20260813
PERMUTATIONS = 1_000


def reconcile_inputs() -> pd.DataFrame:
    scores = pd.read_csv(SCORES)
    author = pd.read_csv(EGRF)
    required = {"Pair", "EnsID", "Enh", "HitPermissive_NegZ", "EGrf"}
    if not required.issubset(author.columns):
        raise ValueError("Author EGrf table lacks required columns")
    if len(scores) != 2_307 or scores["relation_id"].nunique() != 2_307:
        raise ValueError("Frozen AstroREG cohort cardinality changed")
    merged = scores.merge(
        author[["Pair", "EnsID", "Enh", "HitPermissive_NegZ", "EGrf"]],
        left_on="relation_id",
        right_on="Pair",
        how="left",
        validate="one_to_one",
        suffixes=("", "_author"),
    )
    if merged["EGrf"].isna().any() or set(scores["relation_id"]) != set(author["Pair"]):
        raise ValueError("EGrf rows do not exactly cover the frozen cohort")
    if (merged["gene_id"] != merged["EnsID"]).any() or (merged["Enh"] != merged["Enh_author"]).any():
        raise ValueError("EGrf gene or enhancer identifiers do not reconcile")
    if (merged["crispri_positive"].astype(bool) != merged["HitPermissive_NegZ"].astype(bool)).any():
        raise ValueError("EGrf labels do not reconcile")
    return merged.reset_index(drop=True)


def relation_order(frame: pd.DataFrame, indices: np.ndarray) -> np.ndarray:
    """Stable within-enhancer order used only to transfer whole score vectors."""
    local = frame.loc[indices, ["gene_id", "relation_id"]].copy()
    return local.sort_values(["gene_id", "relation_id"], kind="mergesort").index.to_numpy()


def permute_score_blocks(frame: pd.DataFrame, rng: np.random.Generator) -> np.ndarray:
    """Permute complete score vectors among equally sized enhancer blocks.

    Enhancers with the same number of tested genes exchange their entire score
    vectors through a random cyclic derangement. This preserves the exact score
    distribution, enhancer multiplicity and within-donor block structure. Five
    rare multiplicities occur for only one enhancer each; for those blocks, a
    non-zero cyclic rotation breaks the within-block relation mapping while
    conservatively retaining the enhancer mean.
    """
    observed = frame["alphagenome_deletion_strength"].to_numpy(float)
    shuffled = np.full(len(frame), np.nan)
    blocks = [z.index.to_numpy() for _, z in frame.groupby("Enh", sort=True)]
    by_size: dict[int, list[np.ndarray]] = {}
    for block in blocks:
        by_size.setdefault(len(block), []).append(block)

    for size, size_blocks in sorted(by_size.items()):
        n_blocks = len(size_blocks)
        if n_blocks > 1:
            cycle = rng.permutation(n_blocks)
            shift = int(rng.integers(1, n_blocks))
            for position, recipient_number in enumerate(cycle):
                donor_number = cycle[(position + shift) % n_blocks]
                recipient = relation_order(frame, size_blocks[recipient_number])
                donor = relation_order(frame, size_blocks[donor_number])
                shuffled[recipient] = observed[donor]
        else:
            recipient = relation_order(frame, size_blocks[0])
            if size < 2:
                raise RuntimeError("A unique singleton block cannot be permuted")
            shift = int(rng.integers(1, size))
            shuffled[recipient] = np.roll(observed[recipient], shift)

    if not np.isfinite(shuffled).all():
        raise RuntimeError("Structure-aware permutation left missing scores")
    if not np.array_equal(np.sort(shuffled), np.sort(observed)):
        raise RuntimeError("Structure-aware permutation changed the score multiset")
    return shuffled


def plus_predictions(
    y: np.ndarray,
    genes: np.ndarray,
    enhancers: np.ndarray,
    egrf: np.ndarray,
    ag: np.ndarray,
) -> np.ndarray:
    """Run the same EGrf-plus-AlphaGenome calibration with fixed folds."""
    clipped = np.clip(egrf, 1e-6, 1 - 1e-6)
    egrf_logit = np.log(clipped / (1 - clipped))
    probability = np.full(len(y), np.nan)
    for train, test in dual_group_splits(y, genes, enhancers):
        x_train = np.c_[egrf_logit[train], ag[train]]
        x_test = np.c_[egrf_logit[test], ag[test]]
        probability[test] = fit_logistic(x_train, y[train], x_test)
    if not np.isfinite(probability).all():
        raise RuntimeError("Null stack left missing predictions")
    return probability


def summarize(null: np.ndarray, observed: float) -> dict[str, float | int]:
    return {
        "observed": float(observed),
        "null_mean": float(np.mean(null)),
        "null_median": float(np.median(null)),
        "null_p2_5": float(np.quantile(null, 0.025)),
        "null_p97_5": float(np.quantile(null, 0.975)),
        "empirical_upper_tail_p": float((1 + np.sum(null >= observed)) / (len(null) + 1)),
        "observed_percentile": float(100 * np.mean(null < observed)),
        "observed_ascending_rank": int(1 + np.sum(null < observed)),
        "permutations": int(len(null)),
    }


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    frame = reconcile_inputs()
    y = frame["crispri_positive"].astype(int).to_numpy()
    genes = frame["gene_id"].to_numpy()
    enhancers = frame["Enh"].to_numpy()
    egrf = frame["EGrf"].to_numpy(float)
    ag = frame["alphagenome_deletion_strength"].to_numpy(float)

    base_probability, observed_plus_probability = grouped_stack(y, genes, enhancers, egrf, ag)
    base_ap, base_auc, base_loss = metrics(y, base_probability)
    plus_ap, plus_auc, plus_loss = metrics(y, observed_plus_probability)
    observed = np.array([plus_ap - base_ap, plus_auc - base_auc, base_loss - plus_loss])

    existing = pd.read_csv(RESULTS / "model_oof_predictions.csv")
    expected_base = existing["dual_group_calibrated_egrf_probability"].to_numpy(float)
    expected_plus = existing["dual_group_egrf_plus_ag_probability"].to_numpy(float)
    if not np.allclose(base_probability, expected_base, rtol=0, atol=1e-12):
        raise RuntimeError("Observed EGrf-only stack does not reproduce the targeted audit")
    if not np.allclose(observed_plus_probability, expected_plus, rtol=0, atol=1e-12):
        raise RuntimeError("Observed EGrf-plus-AlphaGenome stack does not reproduce the targeted audit")

    rng = np.random.default_rng(SEED)
    rows: list[dict[str, float | int]] = []
    for permutation in range(1, PERMUTATIONS + 1):
        shuffled_ag = permute_score_blocks(frame, rng)
        null_probability = plus_predictions(y, genes, enhancers, egrf, shuffled_ag)
        null_ap, null_auc, null_loss = metrics(y, null_probability)
        rows.append({
            "permutation": permutation,
            "ap_increment": null_ap - base_ap,
            "auc_increment": null_auc - base_auc,
            "log_loss_improvement": base_loss - null_loss,
            "null_plus_ap": null_ap,
            "null_plus_auc": null_auc,
            "null_plus_log_loss": null_loss,
        })

    null_frame = pd.DataFrame(rows)
    null_frame.to_csv(PERMUTATION_ROWS, index=False)
    summary = {
        "design": {
            "label": "structure-aware shuffled-AlphaGenome stacking null",
            "cohort": "frozen 2,307-relation AstroREG cohort",
            "permutation_unit": "whole AlphaGenome score vectors exchanged among enhancer blocks with identical relation multiplicity",
            "rare_block_rule": "within-block cyclic rotation for the five multiplicities represented by one enhancer",
            "fixed_components": "labels, EGrf probabilities, gene/enhancer folds, logistic calibration, metrics",
            "seed": SEED,
            "permutations": PERMUTATIONS,
        },
        "observed_models": {
            "egrf_only_ap": base_ap,
            "egrf_plus_alphagenome_ap": plus_ap,
            "egrf_only_auc": base_auc,
            "egrf_plus_alphagenome_auc": plus_auc,
            "egrf_only_log_loss": base_loss,
            "egrf_plus_alphagenome_log_loss": plus_loss,
        },
        "ap_increment": summarize(null_frame["ap_increment"].to_numpy(), observed[0]),
        "auc_increment": summarize(null_frame["auc_increment"].to_numpy(), observed[1]),
        "log_loss_improvement": summarize(null_frame["log_loss_improvement"].to_numpy(), observed[2]),
        "interpretation_boundary": "This diagnostic tests an arbitrary-second-feature explanation under the same stacking pipeline. It does not resolve the non-nested provenance of the fixed EGrf probabilities or prove absence of leakage.",
    }
    PERMUTATION_SUMMARY.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("OP26_PERMUTATION_NULL_PASS")


if __name__ == "__main__":
    main()
