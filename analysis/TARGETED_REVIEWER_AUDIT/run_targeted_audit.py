#!/usr/bin/env python3
"""Bounded OP26 reviewer audit: EGrf, gene confounding, nonlinear expression."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, log_loss, roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import SplineTransformer, StandardScaler


HERE = Path(__file__).resolve().parent
OP26 = HERE.parent
SCORES = OP26 / "EXTERNAL_ASTROREG/results/astroreg_alphagenome_scores.csv"
FROZEN_OOF = OP26 / "EXTERNAL_ASTROREG/results/grouped_oof_predictions.csv"
EGRF = HERE / "inputs/Astrocyte_trainingData_pluspredictions.csv"
RESULTS = HERE / "results"
SEED = 20260811
BOOTSTRAPS = 5_000
PERMUTATIONS = 20_000


def ci(x: np.ndarray) -> tuple[float, float]:
    x = np.asarray(x, dtype=float)
    return float(np.quantile(x, 0.025)), float(np.quantile(x, 0.975))


def metrics(y: np.ndarray, p: np.ndarray) -> tuple[float, float, float]:
    p = np.asarray(p, dtype=float)
    return (
        float(average_precision_score(y, p)),
        float(roc_auc_score(y, p)),
        float(log_loss(y, np.clip(p, 1e-8, 1 - 1e-8))),
    )


def fit_logistic(x_train: np.ndarray, y_train: np.ndarray, x_test: np.ndarray) -> np.ndarray:
    scaler = StandardScaler().fit(x_train)
    model = LogisticRegression(C=1e6, solver="lbfgs", max_iter=3_000)
    model.fit(scaler.transform(x_train), y_train)
    return model.predict_proba(scaler.transform(x_test))[:, 1]


def dual_group_splits(y: np.ndarray, genes: np.ndarray, enhancers: np.ndarray):
    """Hold out both a gene fold and an enhancer fold, matching the EGrf design."""
    gene_fold = np.full(len(y), -1, dtype=int)
    enhancer_fold = np.full(len(y), -1, dtype=int)
    for fold, (_, test) in enumerate(GroupKFold(n_splits=5).split(y, y, genes)):
        gene_fold[test] = fold
    for fold, (_, test) in enumerate(GroupKFold(n_splits=5).split(y, y, enhancers)):
        enhancer_fold[test] = fold
    if gene_fold.min() < 0 or enhancer_fold.min() < 0:
        raise RuntimeError("Dual-group fold assignment failed")
    for gf in range(5):
        for ef in range(5):
            test = np.flatnonzero((gene_fold == gf) & (enhancer_fold == ef))
            if not len(test):
                continue
            train = np.flatnonzero((gene_fold != gf) & (enhancer_fold != ef))
            if np.unique(y[train]).size < 2:
                raise RuntimeError("Dual-group training fold has only one class")
            yield train, test


def grouped_stack(y: np.ndarray, genes: np.ndarray, enhancers: np.ndarray, egrf: np.ndarray, ag: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    clipped = np.clip(egrf, 1e-6, 1 - 1e-6)
    egrf_logit = np.log(clipped / (1 - clipped))
    p_egrf = np.full(len(y), np.nan)
    p_plus = np.full(len(y), np.nan)
    for train, test in dual_group_splits(y, genes, enhancers):
        p_egrf[test] = fit_logistic(egrf_logit[train, None], y[train], egrf_logit[test, None])
        x_train = np.c_[egrf_logit[train], ag[train]]
        x_test = np.c_[egrf_logit[test], ag[test]]
        p_plus[test] = fit_logistic(x_train, y[train], x_test)
    if not np.isfinite(p_egrf).all() or not np.isfinite(p_plus).all():
        raise RuntimeError("Grouped EGrf stack left missing predictions")
    return p_egrf, p_plus


def nonlinear_expression_oof(data: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    log_distance = np.log1p(data["distance_abs"].clip(lower=0).to_numpy(float))
    log_abc = np.log10(data["abc_score"].clip(lower=0).to_numpy(float) + 1e-6)
    log_length = np.log1p(data["enhancer_length"].clip(lower=0).to_numpy(float))
    log_cells = np.log1p(data["n_cells"].clip(lower=0).to_numpy(float))
    log_expr = np.log1p(data["gene_expression"].clip(lower=0).to_numpy(float))[:, None]
    ag = data["alphagenome_deletion_strength"].to_numpy(float)
    fixed = np.c_[log_distance, log_abc, log_length, log_cells]
    y = data["crispri_positive"].astype(int).to_numpy()
    genes = data["gene_id"].to_numpy()
    enhancers = data["Enh"].to_numpy()
    base = np.full(len(y), np.nan)
    plus = np.full(len(y), np.nan)
    for train, test in dual_group_splits(y, genes, enhancers):
        spline = SplineTransformer(n_knots=4, degree=3, include_bias=False, knots="quantile")
        spline_train = spline.fit_transform(log_expr[train])
        spline_test = spline.transform(log_expr[test])
        x_train = np.c_[fixed[train], spline_train]
        x_test = np.c_[fixed[test], spline_test]
        base[test] = fit_logistic(x_train, y[train], x_test)
        plus[test] = fit_logistic(np.c_[x_train, ag[train]], y[train], np.c_[x_test, ag[test]])
    return base, plus


def cluster_bootstrap_deltas(frame: pd.DataFrame, base: str, plus: str) -> pd.DataFrame:
    blocks = [z.index.to_numpy() for _, z in frame.groupby("enhancer_id", sort=False)]
    rng = np.random.default_rng(SEED + sum(map(ord, base + plus)))
    rows = []
    y_all = frame["crispri_positive"].astype(int).to_numpy()
    for _ in range(BOOTSTRAPS):
        chosen = rng.integers(0, len(blocks), size=len(blocks))
        idx = np.concatenate([blocks[j] for j in chosen])
        y = y_all[idx]
        if np.unique(y).size < 2:
            continue
        mb = metrics(y, frame[base].to_numpy()[idx])
        mp = metrics(y, frame[plus].to_numpy()[idx])
        rows.append((mp[0] - mb[0], mp[1] - mb[1], mb[2] - mp[2]))
    values = np.asarray(rows)
    point_base = metrics(y_all, frame[base].to_numpy())
    point_plus = metrics(y_all, frame[plus].to_numpy())
    names = ["ap_increment", "auc_increment", "log_loss_improvement"]
    points = [point_plus[0] - point_base[0], point_plus[1] - point_base[1], point_base[2] - point_plus[2]]
    out = []
    for j, name in enumerate(names):
        lo, hi = ci(values[:, j])
        out.append({"comparison": f"{plus}_minus_{base}", "metric": name, "estimate": points[j], "ci_low": lo, "ci_high": hi, "cluster": "enhancer", "resamples": len(values)})
    return pd.DataFrame(out)


def within_gene_analysis(data: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    per_gene = []
    gene_arrays = []
    for gene, z in data.groupby("gene_id", sort=True):
        y = z["crispri_positive"].astype(bool).to_numpy()
        if not y.any() or y.all():
            continue
        score = z["alphagenome_deletion_strength"].to_numpy(float)
        pos, neg = score[y], score[~y]
        wins = float((pos[:, None] > neg).sum() + 0.5 * (pos[:, None] == neg).sum())
        pairs = int(len(pos) * len(neg))
        per_gene.append({
            "gene_id": gene,
            "relations": len(z),
            "positives": len(pos),
            "negatives": len(neg),
            "mean_contrast": float(pos.mean() - neg.mean()),
            "conditional_auc": wins / pairs,
            "comparison_pairs": pairs,
            "wins_with_half_ties": wins,
        })
        gene_arrays.append((score, int(y.sum())))
    per_gene = pd.DataFrame(per_gene)
    observed_contrast = float(per_gene["mean_contrast"].mean())
    observed_equal_auc = float(per_gene["conditional_auc"].mean())
    observed_pair_auc = float(per_gene["wins_with_half_ties"].sum() / per_gene["comparison_pairs"].sum())

    rng = np.random.default_rng(SEED + 91)
    boot = np.empty((BOOTSTRAPS, 3))
    for i in range(BOOTSTRAPS):
        take = rng.integers(0, len(per_gene), size=len(per_gene))
        b = per_gene.iloc[take]
        boot[i] = [b.mean_contrast.mean(), b.conditional_auc.mean(), b.wins_with_half_ties.sum() / b.comparison_pairs.sum()]

    null = np.empty((PERMUTATIONS, 3))
    for i in range(PERMUTATIONS):
        contrasts, aucs, total_wins, total_pairs = [], [], 0.0, 0
        for score, npos in gene_arrays:
            pos_index = rng.choice(len(score), size=npos, replace=False)
            mask = np.zeros(len(score), dtype=bool)
            mask[pos_index] = True
            pos, neg = score[mask], score[~mask]
            wins = float((pos[:, None] > neg).sum() + 0.5 * (pos[:, None] == neg).sum())
            pairs = len(pos) * len(neg)
            contrasts.append(pos.mean() - neg.mean())
            aucs.append(wins / pairs)
            total_wins += wins
            total_pairs += pairs
        null[i] = [np.mean(contrasts), np.mean(aucs), total_wins / total_pairs]

    estimates = [observed_contrast, observed_equal_auc, observed_pair_auc]
    metric_names = ["equal_gene_mean_contrast", "equal_gene_conditional_auc", "pair_weighted_conditional_auc"]
    summary = []
    for j, name in enumerate(metric_names):
        lo, hi = ci(boot[:, j])
        p = float((1 + np.sum(null[:, j] >= estimates[j])) / (PERMUTATIONS + 1))
        summary.append({"metric": name, "estimate": estimates[j], "ci_low": lo, "ci_high": hi, "permutation_p_one_sided": p, "genes": len(per_gene), "bootstrap_resamples": BOOTSTRAPS, "permutations": PERMUTATIONS})
    return per_gene, pd.DataFrame(summary)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(SCORES)
    frozen = pd.read_csv(FROZEN_OOF)
    author = pd.read_csv(EGRF)
    required = {"Pair", "Enh", "EnsID", "HitPermissive_NegZ", "EGrf"}
    if not required.issubset(author.columns):
        raise ValueError("Author EGrf table lacks required columns")
    if len(data) != 2307 or data["relation_id"].nunique() != 2307:
        raise ValueError("Frozen OP26 cohort cardinality changed")
    merged = data.merge(author[["Pair", "EnsID", "Enh", "HitPermissive_NegZ", "EGrf"]], left_on="relation_id", right_on="Pair", how="left", validate="one_to_one", suffixes=("", "_author"))
    if merged["EGrf"].isna().any() or set(data["relation_id"]) != set(author["Pair"]):
        raise ValueError("EGrf rows do not exactly cover the frozen cohort")
    if (merged["gene_id"] != merged["EnsID"]).any() or (merged["Enh"] != merged["Enh_author"]).any():
        raise ValueError("EGrf gene/enhancer identifiers do not reconcile")
    if (merged["crispri_positive"].astype(bool) != merged["HitPermissive_NegZ"].astype(bool)).any():
        raise ValueError("EGrf labels do not reconcile")

    merged = merged.merge(frozen, on=["relation_id", "gene_id"], how="left", validate="one_to_one", suffixes=("", "_frozen"))
    y = merged["crispri_positive"].astype(int).to_numpy()
    genes = merged["gene_id"].to_numpy()
    enhancers = merged["Enh"].to_numpy()
    ag = merged["alphagenome_deletion_strength"].to_numpy(float)
    egrf = merged["EGrf"].to_numpy(float)
    egrf_cal, egrf_plus_ag = grouped_stack(y, genes, enhancers, egrf, ag)
    spline_context, spline_plus_ag = nonlinear_expression_oof(merged)

    oof = pd.DataFrame({
        "relation_id": merged["relation_id"],
        "enhancer_id": merged["Enh"],
        "gene_id": merged["gene_id"],
        "crispri_positive": y,
        "context_probability": merged["context_probability"],
        "context_plus_ag_probability": merged["full_probability"],
        "author_egrf_oof_probability": egrf,
        "dual_group_calibrated_egrf_probability": egrf_cal,
        "dual_group_egrf_plus_ag_probability": egrf_plus_ag,
        "spline_expression_context_probability": spline_context,
        "spline_expression_context_plus_ag_probability": spline_plus_ag,
    })
    oof.to_csv(RESULTS / "model_oof_predictions.csv", index=False)

    model_rows = []
    for name in oof.columns[4:]:
        ap, auc, loss = metrics(y, oof[name].to_numpy())
        model_rows.append({"model": name, "average_precision": ap, "roc_auc": auc, "log_loss": loss, "relations": len(oof), "positives": int(y.sum()), "evaluation_unit": "relation; grouped/cross-fitted as documented"})
    model_metrics = pd.DataFrame(model_rows)
    model_metrics.to_csv(RESULTS / "model_comparison_metrics.csv", index=False)

    deltas = pd.concat([
        cluster_bootstrap_deltas(oof, "context_probability", "context_plus_ag_probability"),
        cluster_bootstrap_deltas(oof, "dual_group_calibrated_egrf_probability", "dual_group_egrf_plus_ag_probability"),
        cluster_bootstrap_deltas(oof, "spline_expression_context_probability", "spline_expression_context_plus_ag_probability"),
    ], ignore_index=True)
    deltas.to_csv(RESULTS / "paired_cluster_increment_metrics.csv", index=False)

    per_gene, gene_summary = within_gene_analysis(merged)
    per_gene.to_csv(RESULTS / "within_gene_per_gene.csv", index=False)
    gene_summary.to_csv(RESULTS / "within_gene_summary.csv", index=False)

    relation_corr = spearmanr(ag, merged["gene_expression"].to_numpy(float))
    by_gene = merged.groupby("gene_id").agg(mean_ag=("alphagenome_deletion_strength", "mean"), expression=("gene_expression", "median"), relations=("relation_id", "size"))
    gene_corr = spearmanr(by_gene["mean_ag"], by_gene["expression"])
    label_corr = spearmanr(y, merged["gene_expression"].to_numpy(float))
    expression_ap = float(average_precision_score(y, merged["gene_expression"]))
    expression_auc = float(roc_auc_score(y, merged["gene_expression"]))
    confounding = pd.DataFrame([
        {"metric": "relation_spearman_ag_vs_measured_expression", "estimate": relation_corr.statistic, "p_value_two_sided": relation_corr.pvalue, "n": len(merged)},
        {"metric": "gene_equal_spearman_mean_ag_vs_measured_expression", "estimate": gene_corr.statistic, "p_value_two_sided": gene_corr.pvalue, "n": len(by_gene)},
        {"metric": "relation_spearman_functional_label_vs_measured_expression", "estimate": label_corr.statistic, "p_value_two_sided": label_corr.pvalue, "n": len(merged)},
        {"metric": "expression_only_average_precision", "estimate": expression_ap, "p_value_two_sided": np.nan, "n": len(merged)},
        {"metric": "expression_only_roc_auc", "estimate": expression_auc, "p_value_two_sided": np.nan, "n": len(merged)},
    ])
    confounding.to_csv(RESULTS / "gene_expression_confounding_metrics.csv", index=False)

    summary = {
        "relations": len(merged),
        "positives": int(y.sum()),
        "enhancers": int(merged["Enh"].nunique()),
        "genes": int(merged["gene_id"].nunique()),
        "mixed_label_genes": int(len(per_gene)),
        "mixed_label_gene_relations": int(merged["gene_id"].isin(per_gene["gene_id"]).sum()),
        "mixed_label_gene_positives": int(merged.loc[merged["gene_id"].isin(per_gene["gene_id"]), "crispri_positive"].sum()),
        "model_metrics": model_metrics.set_index("model")[["average_precision", "roc_auc", "log_loss"]].to_dict(orient="index"),
        "increments": deltas.to_dict(orient="records"),
        "within_gene": gene_summary.to_dict(orient="records"),
        "expression_confounding": confounding.astype(object).where(pd.notna(confounding), None).to_dict(orient="records"),
        "author_egrf_fairness": {
            "row_coverage": "2307/2307 exact Pair identifiers",
            "author_prediction_method": "gene-by-enhancer held-out cross-validation",
            "caveat": "author RF hyperparameters were tuned before its gene-by-enhancer held-out prediction step; our calibration/stacking also holds out both gene and enhancer folds and treats published cross-fitted scores as fixed inputs",
        },
    }
    (RESULTS / "analysis_summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False), encoding="utf-8")
    print("OP26_TARGETED_AUDIT_ANALYSIS_PASS")


if __name__ == "__main__":
    main()
