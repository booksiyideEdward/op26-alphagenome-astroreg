#!/usr/bin/env python3
"""Retrieve/reconcile the non-redistributed AstroREG EGrf source table.

This script joins the official author table to the minimized OP26 AlphaGenome
score table. It does not call AlphaGenome and stops on any cohort mismatch.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from urllib.request import urlretrieve

import pandas as pd


REPO = Path(__file__).resolve().parents[1]
EGRF_URL = (
    "https://raw.githubusercontent.com/Voineagulab/astrocyte_crispri/main/"
    "6_PredictiveModels/3.Predictions/RF_Results/Astrocytes/"
    "Astrocyte_trainingData_pluspredictions.csv"
)
DEFAULT_EGRF = REPO / "analysis/TARGETED_REVIEWER_AUDIT/inputs/Astrocyte_trainingData_pluspredictions.csv"
DEFAULT_MINIMAL = REPO / "data/derived/astroreg_alphagenome_relation_scores.csv"
DEFAULT_SCORES = REPO / "analysis/EXTERNAL_ASTROREG/results/astroreg_alphagenome_scores.csv"


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--egrf", type=Path, default=DEFAULT_EGRF)
    parser.add_argument("--minimal-scores", type=Path, default=DEFAULT_MINIMAL)
    parser.add_argument("--out", type=Path, default=DEFAULT_SCORES)
    parser.add_argument("--download-egrf", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = arguments()
    if args.download_egrf:
        args.egrf.parent.mkdir(parents=True, exist_ok=True)
        urlretrieve(EGRF_URL, args.egrf)
    if not args.egrf.exists():
        raise FileNotFoundError(
            f"Missing official EGrf source table: {args.egrf}. "
            "Use --download-egrf or provide --egrf PATH."
        )

    author = pd.read_csv(args.egrf)
    scores = pd.read_csv(args.minimal_scores)
    required_author = {
        "Pair", "Enh", "EnsID", "Gene", "Enh.Pos", "Enh.chr", "Enh.start",
        "Enh.end", "Enh.size", "Gene.TSS", "Distance", "Gene.Exp", "nCells",
        "ABCScore_ENCODE", "HitPermissive_NegZ", "EGrf",
    }
    required_scores = {
        "relation_id", "Enh", "enhancer_id", "gene_id", "model_ontology",
        "model_track_name", "model_track_source", "alphagenome_rnaseq_lfc",
        "alphagenome_deletion_strength", "status",
    }
    if not required_author.issubset(author.columns):
        raise ValueError("Official EGrf table lacks required columns")
    if not required_scores.issubset(scores.columns):
        raise ValueError("OP26 minimized score table lacks required columns")
    if len(author) != 2307 or author["Pair"].nunique() != 2307:
        raise ValueError("Official EGrf table no longer matches the frozen 2,307-relation cohort")
    if len(scores) != 2307 or scores["relation_id"].nunique() != 2307:
        raise ValueError("OP26 score table no longer matches the frozen 2,307-relation cohort")

    joined = scores.merge(
        author,
        left_on="relation_id",
        right_on="Pair",
        how="inner",
        validate="one_to_one",
        suffixes=("", "_author"),
    )
    if len(joined) != 2307:
        raise ValueError("Author and OP26 relation identifiers do not reconcile")
    if (joined["Enh_author"] != joined["Enh"]).any():
        raise ValueError("Enhancer identifiers do not reconcile")
    if (joined["EnsID"] != joined["gene_id"]).any():
        raise ValueError("Gene identifiers do not reconcile")

    output = pd.DataFrame({
        "relation_id": joined["relation_id"],
        "Enh": joined["Enh"],
        "enhancer_region_hg38": joined["Enh.Pos"],
        "chromosome": joined["Enh.chr"],
        "enhancer_start_hg38": joined["Enh.start"],
        "enhancer_end_hg38": joined["Enh.end"],
        "enhancer_length": joined["Enh.size"],
        "gene_id": joined["gene_id"],
        "gene_name": joined["Gene"],
        "gene_tss_hg38": joined["Gene.TSS"],
        "distance_abs": joined["Distance"].abs(),
        "crispri_positive": joined["HitPermissive_NegZ"].astype(bool),
        "abc_score": joined["ABCScore_ENCODE"],
        "gene_expression": joined["Gene.Exp"],
        "n_cells": joined["nCells"],
        "enhancer_id": joined["enhancer_id"],
        "model_ontology": joined["model_ontology"],
        "model_track_name": joined["model_track_name"],
        "model_track_source": joined["model_track_source"],
        "status": joined["status"],
        "alphagenome_rnaseq_lfc": joined["alphagenome_rnaseq_lfc"],
        "alphagenome_deletion_strength": joined["alphagenome_deletion_strength"],
    })
    if int(output["crispri_positive"].sum()) != 133:
        raise ValueError("Frozen positive count changed")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.out, index=False)
    print(f"Wrote {len(output)} reconciled relations to {args.out}")


if __name__ == "__main__":
    main()
