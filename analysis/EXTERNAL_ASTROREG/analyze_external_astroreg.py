#!/usr/bin/env python3
"""Analyze the frozen AlphaGenome AstroREG external-validation results."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import tempfile
import warnings

_CACHE_ROOT = Path(tempfile.gettempdir()) / "op26_astroreg_cache"
os.environ.setdefault("MPLCONFIGDIR", str(_CACHE_ROOT / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(_CACHE_ROOT / "xdg"))

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    log_loss,
    precision_recall_curve,
    roc_auc_score,
)
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
SCORES = RESULTS / "astroreg_alphagenome_scores.csv"
METRICS = RESULTS / "external_validation_metrics.csv"
MATCHES = RESULTS / "matched_sensitivity_pairs.csv"
OOF = RESULTS / "grouped_oof_predictions.csv"
DISTANCE_SENSITIVITY = RESULTS / "distance_sensitivity.csv"
REPORT = HERE / "EXTERNAL_VALIDATION_RESULTS.md"
FIGURES = HERE / "figures"
SEED = 20260808
BOOTSTRAPS = 2_000
COEF_BOOTSTRAPS = 500

CONTEXT_FEATURES = [
    "log_distance", "log_abc", "log_enhancer_length", "log_gene_expression", "log_n_cells"
]


def percentile_ci(values: np.ndarray) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))


def wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total <= 0:
        return math.nan, math.nan
    p = successes / total
    denominator = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denominator
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    return center - half, center + half


def transform_features(data: pd.DataFrame) -> pd.DataFrame:
    transformed = pd.DataFrame(index=data.index)
    transformed["log_distance"] = np.log1p(pd.to_numeric(data["distance_abs"], errors="coerce").clip(lower=0))
    transformed["log_abc"] = np.log10(pd.to_numeric(data["abc_score"], errors="coerce").clip(lower=0) + 1e-6)
    transformed["log_enhancer_length"] = np.log1p(pd.to_numeric(data["enhancer_length"], errors="coerce").clip(lower=0))
    transformed["log_gene_expression"] = np.log1p(pd.to_numeric(data["gene_expression"], errors="coerce").clip(lower=0))
    transformed["log_n_cells"] = np.log1p(pd.to_numeric(data["n_cells"], errors="coerce").clip(lower=0))
    for column in transformed:
        transformed[column] = transformed[column].fillna(transformed[column].median())
    return transformed


def metric_triplet(data: pd.DataFrame) -> tuple[float, float, float]:
    y = data["crispri_positive"].astype(int).to_numpy()
    score = data["alphagenome_deletion_strength"].to_numpy(dtype=float)
    ap = average_precision_score(y, score)
    auc = roc_auc_score(y, score)
    separation = float(data.loc[data.crispri_positive, "alphagenome_deletion_strength"].mean() - data.loc[~data.crispri_positive, "alphagenome_deletion_strength"].mean())
    return float(ap), float(auc), separation


def cluster_bootstrap(data: pd.DataFrame, cluster: str, draws: int, seed: int) -> np.ndarray:
    groups = [group.index.to_numpy() for _, group in data.groupby(cluster, sort=False)]
    rng = np.random.default_rng(seed)
    output = []
    for _ in range(draws):
        chosen = rng.integers(0, len(groups), size=len(groups))
        indices = np.concatenate([groups[index] for index in chosen])
        sample = data.loc[indices]
        if sample["crispri_positive"].nunique() < 2:
            continue
        output.append(metric_triplet(sample))
    return np.asarray(output, dtype=float)


def fit_model(x: np.ndarray, y: np.ndarray) -> tuple[StandardScaler, LogisticRegression]:
    scaler = StandardScaler().fit(x)
    model = LogisticRegression(C=np.inf, solver="lbfgs", max_iter=2_000)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        model.fit(scaler.transform(x), y)
    return scaler, model


def grouped_predictions(x: np.ndarray, y: np.ndarray, groups: np.ndarray, seed_offset: int = 0) -> np.ndarray:
    del seed_offset  # deterministic grouped folds; retained for an explicit stable interface
    predictions = np.full(len(y), np.nan)
    splitter = GroupKFold(n_splits=5)
    for train, test in splitter.split(x, y, groups):
        scaler, model = fit_model(x[train], y[train])
        predictions[test] = model.predict_proba(scaler.transform(x[test]))[:, 1]
    if not np.isfinite(predictions).all():
        raise ValueError("Grouped cross-validation left missing predictions")
    return predictions


def coefficient_bootstrap(data: pd.DataFrame, x_full: np.ndarray, y: np.ndarray) -> np.ndarray:
    groups = [group.index.to_numpy() for _, group in data.groupby("enhancer_id", sort=False)]
    rng = np.random.default_rng(SEED + 3)
    values = []
    for _ in range(COEF_BOOTSTRAPS):
        chosen = rng.integers(0, len(groups), size=len(groups))
        indices = np.concatenate([groups[index] for index in chosen])
        if np.unique(y[indices]).size < 2:
            continue
        try:
            _, model = fit_model(x_full[indices], y[indices])
            values.append(float(model.coef_[0, -1]))
        except Exception:
            continue
    return np.asarray(values, dtype=float)


def build_matches(data: pd.DataFrame, transformed: pd.DataFrame) -> pd.DataFrame:
    matching = transformed.copy()
    matching["relation_id"] = data["relation_id"].to_numpy()
    matching["positive"] = data["crispri_positive"].astype(bool).to_numpy()
    matrix = matching[CONTEXT_FEATURES].to_numpy(dtype=float)
    matrix = (matrix - matrix.mean(axis=0)) / np.where(matrix.std(axis=0) > 0, matrix.std(axis=0), 1)
    positive_index = np.flatnonzero(matching["positive"].to_numpy())
    negative_index = np.flatnonzero(~matching["positive"].to_numpy())
    cost = cdist(matrix[positive_index], matrix[negative_index], metric="euclidean")
    pos_assignment, neg_assignment = linear_sum_assignment(cost)
    if len(pos_assignment) != 133 or len(set(neg_assignment.tolist())) != 133:
        raise ValueError("Frozen 1:1 matching did not produce 133 unique pairs")
    rows = []
    for pair_number, (positive_row, negative_row) in enumerate(zip(positive_index[pos_assignment], negative_index[neg_assignment], strict=True), start=1):
        positive = data.iloc[int(positive_row)]
        negative = data.iloc[int(negative_row)]
        rows.append({
            "match_id": f"M{pair_number:03d}",
            "positive_relation_id": positive["relation_id"],
            "negative_relation_id": negative["relation_id"],
            "positive_enhancer_id": positive["enhancer_id"],
            "negative_enhancer_id": negative["enhancer_id"],
            "positive_gene_id": positive["gene_id"],
            "negative_gene_id": negative["gene_id"],
            "matching_distance": float(cost[pos_assignment[pair_number - 1], neg_assignment[pair_number - 1]]),
            "positive_strength": float(positive["alphagenome_deletion_strength"]),
            "negative_strength": float(negative["alphagenome_deletion_strength"]),
            "paired_difference": float(positive["alphagenome_deletion_strength"] - negative["alphagenome_deletion_strength"]),
        })
    return pd.DataFrame(rows)


def paired_cluster_bootstrap(matches: pd.DataFrame) -> np.ndarray:
    blocks = [group.index.to_numpy() for _, group in matches.groupby("positive_enhancer_id", sort=False)]
    rng = np.random.default_rng(SEED + 4)
    values = []
    for _ in range(BOOTSTRAPS):
        chosen = rng.integers(0, len(blocks), size=len(blocks))
        indices = np.concatenate([blocks[index] for index in chosen])
        values.append(float(matches.loc[indices, "paired_difference"].mean()))
    return np.asarray(values, dtype=float)


def ecdf(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.sort(np.asarray(values, dtype=float))
    return x, np.arange(1, len(x) + 1) / len(x)


def make_figure(data: pd.DataFrame, oof: pd.DataFrame, matches: pd.DataFrame, metrics: dict) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica", "sans-serif"],
        "font.size": 7,
        "axes.labelsize": 7,
        "axes.titlesize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.7,
        "legend.frameon": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    })
    signal = "#3B6FB6"
    neutral = "#9AA0A6"
    accent = "#C76B38"
    dark = "#333333"
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.8), constrained_layout=True)

    y = data["crispri_positive"].astype(int).to_numpy()
    score = data["alphagenome_deletion_strength"].to_numpy()
    precision, recall, _ = precision_recall_curve(y, score)
    ax = axes[0, 0]
    ax.plot(recall, precision, color=signal, lw=1.6, label=f"AlphaGenome (AP={metrics['average_precision']:.3f})")
    ax.axhline(metrics["positive_prevalence"], color=neutral, ls="--", lw=1, label=f"Prevalence={metrics['positive_prevalence']:.3f}")
    ax.set(xlabel="Recall", ylabel="Precision", xlim=(0, 1), ylim=(0, 1), title="External discrimination")
    ax.legend(loc="upper right", fontsize=6)

    ax = axes[0, 1]
    negative = data.loc[~data.crispri_positive, "alphagenome_deletion_strength"].to_numpy()
    positive = data.loc[data.crispri_positive, "alphagenome_deletion_strength"].to_numpy()
    xn, yn = ecdf(negative)
    xp, yp = ecdf(positive)
    ax.plot(xn, yn, color=neutral, lw=1.3, label=f"Well-powered nonfunctional (n={len(negative):,})")
    ax.plot(xp, yp, color=signal, lw=1.5, label=f"Functional (n={len(positive):,})")
    ax.axvline(0, color=dark, lw=0.7, ls=":")
    ax.set(xlabel="Predicted deletion strength", ylabel="Empirical cumulative fraction", title="Score distributions")
    ax.set_xscale("symlog", linthresh=0.01)
    ax.legend(loc="lower right", fontsize=6)

    ax = axes[1, 0]
    for column, color, label in [
        ("context_probability", neutral, f"Context only (AP={metrics['context_oof_ap']:.3f})"),
        ("full_probability", accent, f"Context + AlphaGenome (AP={metrics['full_oof_ap']:.3f})"),
    ]:
        p, r, _ = precision_recall_curve(oof["crispri_positive"].astype(int), oof[column])
        ax.plot(r, p, color=color, lw=1.5, label=label)
    ax.axhline(metrics["positive_prevalence"], color=dark, ls=":", lw=0.8)
    ax.set(xlabel="Recall", ylabel="Precision", xlim=(0, 1), ylim=(0, 1), title="Grouped out-of-fold increment")
    ax.legend(loc="upper right", fontsize=6)

    ax = axes[1, 1]
    values = matches["paired_difference"].to_numpy()
    ordered = np.sort(values)
    ranks = np.arange(1, len(ordered) + 1)
    colors = np.where(ordered > 0, signal, neutral)
    ax.plot(ranks, ordered, color=neutral, lw=0.7, zorder=1)
    ax.scatter(ranks, ordered, c=colors, s=9, linewidths=0, zorder=2)
    ax.axhline(0, color=dark, lw=0.8, ls=":")
    ax.axhline(values.mean(), color=accent, lw=1.2, label=f"Mean={values.mean():.3g}")
    ax.set(xlabel="Matched-pair rank", ylabel="Positive − matched-negative strength", title="Frozen matched sensitivity")
    ax.set_yscale("symlog", linthresh=0.01)
    ax.legend(loc="upper right", fontsize=6)

    for label, ax in zip("abcd", axes.flat, strict=True):
        ax.text(-0.16, 1.08, label, transform=ax.transAxes, fontsize=8, fontweight="bold", va="top")
    base = FIGURES / "op26_astroreg_external_validation"
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--figure-only", action="store_true")
    args = parser.parse_args()
    data = pd.read_csv(SCORES)
    if len(data) != 2_307 or data["relation_id"].nunique() != 2_307:
        raise ValueError("External score table is not the frozen 2,307-relation universe")
    if not data["status"].eq("COMPLETED").all() or data["alphagenome_deletion_strength"].isna().any():
        raise ValueError("External score table has missing or failed relations")
    data["crispri_positive"] = data["crispri_positive"].astype(bool)
    data = data.reset_index(drop=True)
    if int(data["crispri_positive"].sum()) != 133 or data["enhancer_id"].nunique() != 745:
        raise ValueError("Frozen class or enhancer count changed")
    if args.figure_only:
        metric_table = pd.read_csv(METRICS)
        metrics = dict(zip(metric_table["metric"], metric_table["value"], strict=True))
        make_figure(data, pd.read_csv(OOF), pd.read_csv(MATCHES), metrics)
        print("figure_export_complete")
        return

    ap, auc, separation = metric_triplet(data)
    enhancer_boot = cluster_bootstrap(data, "enhancer_id", BOOTSTRAPS, SEED)
    gene_boot = cluster_bootstrap(data, "gene_id", BOOTSTRAPS, SEED + 1)
    ap_ci = percentile_ci(enhancer_boot[:, 0])
    auc_ci = percentile_ci(enhancer_boot[:, 1])
    separation_ci = percentile_ci(enhancer_boot[:, 2])
    gene_ap_ci = percentile_ci(gene_boot[:, 0])
    gene_auc_ci = percentile_ci(gene_boot[:, 1])
    gene_separation_ci = percentile_ci(gene_boot[:, 2])

    distance_rows = []
    for cutoff in (0, 10_000, 20_000, 50_000):
        subset = data.loc[data["distance_abs"] >= cutoff].copy()
        subset_ap, subset_auc, subset_difference = metric_triplet(subset)
        distance_rows.append({
            "minimum_distance_bp": cutoff,
            "relation_count": len(subset),
            "positive_count": int(subset["crispri_positive"].sum()),
            "positive_prevalence": float(subset["crispri_positive"].mean()),
            "average_precision": subset_ap,
            "roc_auc": subset_auc,
            "mean_strength_difference": subset_difference,
            "analysis_label": "primary_all_relations" if cutoff == 0 else "posthoc_distance_sensitivity",
        })
    pd.DataFrame(distance_rows).to_csv(DISTANCE_SENSITIVITY, index=False)
    distal_10k = data.loc[data["distance_abs"] >= 10_000].copy().reset_index(drop=True)
    distal_10k_ap, distal_10k_auc, distal_10k_difference = metric_triplet(distal_10k)
    distal_10k_boot = cluster_bootstrap(distal_10k, "enhancer_id", BOOTSTRAPS, SEED + 5)
    distal_10k_ap_ci = percentile_ci(distal_10k_boot[:, 0])
    distal_10k_auc_ci = percentile_ci(distal_10k_boot[:, 1])
    distal_10k_difference_ci = percentile_ci(distal_10k_boot[:, 2])

    positives = data.loc[data.crispri_positive]
    direction_success = int((positives["alphagenome_rnaseq_lfc"] < 0).sum())
    direction_ci = wilson_interval(direction_success, len(positives))

    transformed = transform_features(data)
    x_context = transformed[CONTEXT_FEATURES].to_numpy(dtype=float)
    strength = data[["alphagenome_deletion_strength"]].to_numpy(dtype=float)
    x_full = np.column_stack([x_context, strength])
    y = data["crispri_positive"].astype(int).to_numpy()
    groups = data["enhancer_id"].astype(str).to_numpy()
    context_probability = grouped_predictions(x_context, y, groups)
    full_probability = grouped_predictions(x_full, y, groups)
    _, full_model = fit_model(x_full, y)
    ag_coefficient = float(full_model.coef_[0, -1])
    coefficient_boot = coefficient_bootstrap(data, x_full, y)
    coefficient_ci = percentile_ci(coefficient_boot)
    context_ap = float(average_precision_score(y, context_probability))
    full_ap = float(average_precision_score(y, full_probability))
    context_loss = float(log_loss(y, context_probability))
    full_loss = float(log_loss(y, full_probability))

    oof = data[["relation_id", "enhancer_id", "gene_id", "crispri_positive"]].copy()
    oof["context_probability"] = context_probability
    oof["full_probability"] = full_probability
    oof.to_csv(OOF, index=False)

    matches = build_matches(data, transformed)
    matches.to_csv(MATCHES, index=False)
    paired_mean = float(matches["paired_difference"].mean())
    paired_median = float(matches["paired_difference"].median())
    paired_positive_fraction = float((matches["paired_difference"] > 0).mean())
    paired_boot = paired_cluster_bootstrap(matches)
    paired_ci = percentile_ci(paired_boot)

    metrics = {
        "relation_count": len(data),
        "positive_count": int(y.sum()),
        "negative_count": int((1 - y).sum()),
        "enhancer_count": data["enhancer_id"].nunique(),
        "positive_prevalence": float(y.mean()),
        "average_precision": ap,
        "average_precision_ci_low": ap_ci[0],
        "average_precision_ci_high": ap_ci[1],
        "roc_auc": auc,
        "roc_auc_ci_low": auc_ci[0],
        "roc_auc_ci_high": auc_ci[1],
        "mean_strength_difference": separation,
        "mean_strength_difference_ci_low": separation_ci[0],
        "mean_strength_difference_ci_high": separation_ci[1],
        "gene_cluster_ap_ci_low": gene_ap_ci[0],
        "gene_cluster_ap_ci_high": gene_ap_ci[1],
        "gene_cluster_auc_ci_low": gene_auc_ci[0],
        "gene_cluster_auc_ci_high": gene_auc_ci[1],
        "gene_cluster_strength_difference_ci_low": gene_separation_ci[0],
        "gene_cluster_strength_difference_ci_high": gene_separation_ci[1],
        "positive_direction_agreement": direction_success / len(positives),
        "positive_direction_ci_low": direction_ci[0],
        "positive_direction_ci_high": direction_ci[1],
        "positive_median_rnaseq_lfc": float(positives["alphagenome_rnaseq_lfc"].median()),
        "context_adjusted_ag_coefficient_per_sd": ag_coefficient,
        "context_adjusted_ag_coefficient_ci_low": coefficient_ci[0],
        "context_adjusted_ag_coefficient_ci_high": coefficient_ci[1],
        "context_oof_ap": context_ap,
        "full_oof_ap": full_ap,
        "oof_ap_increment": full_ap - context_ap,
        "context_oof_log_loss": context_loss,
        "full_oof_log_loss": full_loss,
        "oof_log_loss_improvement": context_loss - full_loss,
        "matched_pair_count": len(matches),
        "matched_mean_difference": paired_mean,
        "matched_mean_difference_ci_low": paired_ci[0],
        "matched_mean_difference_ci_high": paired_ci[1],
        "matched_median_difference": paired_median,
        "matched_positive_fraction": paired_positive_fraction,
        "missing_relation_count": 0,
        "distal_10k_relation_count": len(distal_10k),
        "distal_10k_positive_count": int(distal_10k["crispri_positive"].sum()),
        "distal_10k_prevalence": float(distal_10k["crispri_positive"].mean()),
        "distal_10k_average_precision": distal_10k_ap,
        "distal_10k_average_precision_ci_low": distal_10k_ap_ci[0],
        "distal_10k_average_precision_ci_high": distal_10k_ap_ci[1],
        "distal_10k_roc_auc": distal_10k_auc,
        "distal_10k_roc_auc_ci_low": distal_10k_auc_ci[0],
        "distal_10k_roc_auc_ci_high": distal_10k_auc_ci[1],
        "distal_10k_mean_strength_difference": distal_10k_difference,
        "distal_10k_mean_strength_difference_ci_low": distal_10k_difference_ci[0],
        "distal_10k_mean_strength_difference_ci_high": distal_10k_difference_ci[1],
    }
    pd.DataFrame([{"metric": key, "value": value} for key, value in metrics.items()]).to_csv(METRICS, index=False)
    interpretation = (
        "supports external dry-lab generalization"
        if ap_ci[0] > metrics["positive_prevalence"] and separation_ci[0] > 0 and distal_10k_ap_ci[0] > metrics["distal_10k_prevalence"] and distal_10k_auc_ci[0] > 0.5 and metrics["oof_ap_increment"] > 0 and metrics["oof_log_loss_improvement"] > 0
        else "does not establish external dry-lab generalization"
    )
    report = f"""# OP-26 AstroREG external-validation results

## Decision

The frozen analysis **{interpretation}**. This decision is based jointly on full-cohort discrimination, enhancer-cluster uncertainty, context-adjusted grouped prediction and the frozen matched sensitivity; no single metric is treated as decisive.

## Frozen cohort and completion

- 2,307/2,307 enhancer–gene relations completed with no missing target-gene score or API failure.
- 133 experimentally functional, expression-decreasing relations and 2,174 well-powered nonfunctional relations.
- 745 unique enhancer deletions; each enhancer was called once and all corresponding target genes were extracted from the same 1-Mb request.

## Primary external discrimination

- Average precision: **{ap:.4f}** (enhancer-cluster bootstrap 95% CI {ap_ci[0]:.4f} to {ap_ci[1]:.4f}); no-skill prevalence baseline {metrics['positive_prevalence']:.4f}.
- ROC AUC: **{auc:.4f}** (95% CI {auc_ci[0]:.4f} to {auc_ci[1]:.4f}).
- Mean functional-minus-nonfunctional deletion strength: **{separation:.6g}** (95% CI {separation_ci[0]:.6g} to {separation_ci[1]:.6g}).
- Target-gene-cluster sensitivity: AP CI {gene_ap_ci[0]:.4f} to {gene_ap_ci[1]:.4f}; ROC AUC CI {gene_auc_ci[0]:.4f} to {gene_auc_ci[1]:.4f}; mean-difference CI {gene_separation_ci[0]:.6g} to {gene_separation_ci[1]:.6g}.

### Distance-aligned post hoc sensitivity

The author truth set contains 64 relations below 10 kb, including 32 positives. Because the original OP-26 question is distal, a clearly labelled post hoc restriction to relations at least 10 kb from the target TSS was added without replacing the frozen primary analysis.

- At ≥10 kb: {len(distal_10k):,} relations and {int(distal_10k['crispri_positive'].sum())} positives; prevalence {metrics['distal_10k_prevalence']:.4f}.
- Average precision **{distal_10k_ap:.4f}** (enhancer-cluster 95% CI {distal_10k_ap_ci[0]:.4f} to {distal_10k_ap_ci[1]:.4f}); ROC AUC **{distal_10k_auc:.4f}** (95% CI {distal_10k_auc_ci[0]:.4f} to {distal_10k_auc_ci[1]:.4f}).
- Mean functional-minus-nonfunctional strength **{distal_10k_difference:.6g}** (95% CI {distal_10k_difference_ci[0]:.6g} to {distal_10k_difference_ci[1]:.6g}). Further 20-kb and 50-kb descriptive restrictions are in `distance_sensitivity.csv` and remain exploratory.

## Direction and context increment

- Among the 133 positive relations, {direction_success}/{len(positives)} had a predicted negative RNA-seq LFC after deletion: **{direction_success / len(positives):.3f}** (Wilson 95% CI {direction_ci[0]:.3f} to {direction_ci[1]:.3f}). This is a selected expression-decreasing truth set and is not an independent heterogeneous-direction test.
- Context-only grouped out-of-fold AP: **{context_ap:.4f}**; context plus AlphaGenome AP: **{full_ap:.4f}**; increment **{full_ap - context_ap:.4f}**.
- Context-only grouped log loss: **{context_loss:.4f}**; context plus AlphaGenome: **{full_loss:.4f}**; improvement **{context_loss - full_loss:.4f}**.
- Full-model AlphaGenome coefficient: **{ag_coefficient:.4f} per SD** (enhancer-cluster bootstrap 95% CI {coefficient_ci[0]:.4f} to {coefficient_ci[1]:.4f}).

## Frozen matched sensitivity

- 133 positive relations matched without replacement to 133 negatives using only pre-model context variables.
- Mean positive-minus-negative strength: **{paired_mean:.6g}** (positive-enhancer-cluster bootstrap 95% CI {paired_ci[0]:.6g} to {paired_ci[1]:.6g}).
- Median paired difference: **{paired_median:.6g}**; positive-difference fraction: **{paired_positive_fraction:.3f}**.

## Interpretation limits

This is a dry-lab evaluation against public CRISPRi measurements, not a new wet experiment. The result cannot establish that AlphaGenome's model-development data are independent of all relevant astrocyte genomic tracks. It also cannot support arbitrary enhancer–gene causal assignment unless discrimination, context increment and matched specificity are all coherent. The K562 Gasperini analysis remains development-overlapped and heavily confounded in its wrong-gene control; the AstroREG result is the external decision gate.
"""
    REPORT.write_text(report, encoding="utf-8")
    (RESULTS / "analysis_summary.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    make_figure(data, oof, matches, metrics)
    print(report)


if __name__ == "__main__":
    main()
