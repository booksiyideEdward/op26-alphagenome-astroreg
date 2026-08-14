#!/usr/bin/env python3
"""Post-hoc same-request AstroREG target-relation specificity analysis."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SEED = 20260808
B = 20_000
HERE = Path(__file__).resolve().parent
SCORES = HERE / "EXTERNAL_ASTROREG/results/astroreg_alphagenome_scores.csv"


def main() -> None:
    out = HERE / "results"; figs = HERE / "figures"
    out.mkdir(exist_ok=True); figs.mkdir(exist_ok=True)
    data = pd.read_csv(SCORES)
    mixed_ids = data.groupby("Enh").crispri_positive.nunique().loc[lambda x: x == 2].index
    mixed = data[data.Enh.isin(mixed_ids)].copy()
    if (len(mixed_ids), len(mixed), int(mixed.crispri_positive.sum())) != (106, 471, 115):
        raise ValueError("Frozen mixed-label enhancer cardinality mismatch")

    rows = []
    for enh, g in mixed.groupby("Enh", sort=True):
        pos = g.loc[g.crispri_positive, "alphagenome_deletion_strength"].to_numpy()
        neg = g.loc[~g.crispri_positive, "alphagenome_deletion_strength"].to_numpy()
        comparisons = (pos[:, None] > neg).sum() + 0.5 * (pos[:, None] == neg).sum()
        rows.append({"Enh": enh, "n_relations": len(g), "n_positive": len(pos), "n_negative": len(neg),
                     "positive_minus_negative_mean": pos.mean() - neg.mean(),
                     "conditional_auc": comparisons / (len(pos) * len(neg)),
                     "positive_negative_comparisons": len(pos) * len(neg)})
    per = pd.DataFrame(rows)
    per.to_csv(out / "astroreg_within_enhancer_per_enhancer.csv", index=False)
    mixed.to_csv(out / "astroreg_mixed_enhancer_rows.csv", index=False)

    rng = np.random.default_rng(SEED)
    idx = rng.integers(0, len(per), size=(B, len(per)))
    diff = per.positive_minus_negative_mean.to_numpy()
    auc = per.conditional_auc.to_numpy()
    weights = per.positive_negative_comparisons.to_numpy()
    boot_diff = diff[idx].mean(axis=1)
    boot_auc_equal = auc[idx].mean(axis=1)
    boot_auc_weighted = (auc[idx] * weights[idx]).sum(axis=1) / weights[idx].sum(axis=1)

    observed_diff = diff.mean()
    observed_auc_equal = auc.mean()
    observed_auc_weighted = np.average(auc, weights=weights)

    perm_diff = np.zeros(B)
    perm_auc_equal = np.zeros(B)
    perm_auc_weighted_num = np.zeros(B)
    total_pair_weight = weights.sum()
    for _, g in mixed.groupby("Enh", sort=True):
        scores = g.alphagenome_deletion_strength.to_numpy()
        n, k = len(scores), int(g.crispri_positive.sum())
        selected = np.argpartition(rng.random((B, n)), kth=k-1, axis=1)[:, :k]
        pos_sum = scores[selected].sum(axis=1)
        contrast = pos_sum / k - (scores.sum() - pos_sum) / (n - k)
        ranks = pd.Series(scores).rank(method="average").to_numpy()
        rank_sum = ranks[selected].sum(axis=1)
        pair_weight = k * (n - k)
        c_auc = (rank_sum - k * (k + 1) / 2) / pair_weight
        perm_diff += contrast / len(per)
        perm_auc_equal += c_auc / len(per)
        perm_auc_weighted_num += c_auc * pair_weight
    perm_auc_weighted = perm_auc_weighted_num / total_pair_weight

    metrics = pd.DataFrame([
        {"metric": "equal_enhancer_mean_positive_minus_negative_strength", "estimate": observed_diff,
         "bootstrap_ci_low": np.quantile(boot_diff, .025), "bootstrap_ci_high": np.quantile(boot_diff, .975),
         "permutation_null": 0.0, "one_sided_permutation_p": (1 + np.sum(perm_diff >= observed_diff)) / (B + 1)},
        {"metric": "equal_enhancer_conditional_auc", "estimate": observed_auc_equal,
         "bootstrap_ci_low": np.quantile(boot_auc_equal, .025), "bootstrap_ci_high": np.quantile(boot_auc_equal, .975),
         "permutation_null": 0.5, "one_sided_permutation_p": (1 + np.sum(perm_auc_equal >= observed_auc_equal)) / (B + 1)},
        {"metric": "pair_weighted_conditional_auc", "estimate": observed_auc_weighted,
         "bootstrap_ci_low": np.quantile(boot_auc_weighted, .025), "bootstrap_ci_high": np.quantile(boot_auc_weighted, .975),
         "permutation_null": 0.5, "one_sided_permutation_p": (1 + np.sum(perm_auc_weighted >= observed_auc_weighted)) / (B + 1)},
    ])
    metrics["n_mixed_enhancers"] = len(per); metrics["n_relations"] = len(mixed)
    metrics["n_positive"] = int(mixed.crispri_positive.sum()); metrics["n_negative"] = int((~mixed.crispri_positive).sum())
    metrics["bootstrap_resamples"] = B; metrics["permutations"] = B; metrics["seed"] = SEED
    metrics.to_csv(out / "astroreg_within_enhancer_specificity.csv", index=False)

    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                         "font.size": 7.5, "svg.fonttype": "none", "pdf.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.65), constrained_layout=True)
    ax = axes[0]
    y = np.arange(1, len(per) + 1)
    order = np.argsort(diff)
    ax.scatter(diff[order], y, s=12, alpha=.55, color="#7A7A7A", linewidths=0)
    lo, hi = np.quantile(boot_diff, [.025, .975])
    ax.errorbar(observed_diff, len(per) + 8, xerr=[[observed_diff-lo], [hi-observed_diff]], fmt="o", color="#2878B5", capsize=2)
    ax.axvline(0, color="#7A7A7A", ls="--", lw=.8); ax.set_ylim(0, len(per) + 14)
    ax.set_ylabel("Enhancers ranked by contrast")
    ax.set_xlabel("Within-enhancer positive − negative strength")
    ax.set_title("A  Same-request contrasts (106 enhancers)", loc="left", fontweight="bold")
    ax = axes[1]
    vals = metrics.iloc[1:][::-1]
    ax.errorbar(vals.estimate, [0, 1], xerr=[vals.estimate-vals.bootstrap_ci_low, vals.bootstrap_ci_high-vals.estimate],
                fmt="o", color="#C73E3A", capsize=2)
    ax.axvline(.5, color="#7A7A7A", ls="--", lw=.8)
    ax.set_yticks([0, 1], ["Pair-weighted", "Equal-enhancer"]); ax.set_xlim(.35, .83)
    ax.set_xlabel("Conditional AUC within enhancer")
    ax.set_title("B  Relation ranking", loc="left", fontweight="bold")
    fig.savefig(figs / "op26_astroreg_within_enhancer_specificity.svg", bbox_inches="tight")
    fig.savefig(figs / "op26_astroreg_within_enhancer_specificity.pdf", bbox_inches="tight")
    fig.savefig(figs / "op26_astroreg_within_enhancer_specificity.png", bbox_inches="tight", dpi=300)
    fig.savefig(figs / "op26_astroreg_within_enhancer_specificity.tiff", bbox_inches="tight", dpi=600)
    plt.close(fig)
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
