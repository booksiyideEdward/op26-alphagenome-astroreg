# OP26 targeted reviewer-response audit — Phase 1 decision

**Date:** 2026-08-11
**Scope:** one expert comparator (EGrf), gene-level confounding, timeline wording, then a value-of-information decision on at most one additional external dataset. No model zoo or pipeline rewrite.

## Decisions before new analysis

### EGrf — MUST DO, using a bounded comparison

The objection is scientifically valid but the third-party framing is too simple. EGrf is a dataset-specific supervised random forest trained on the AstroREG CRISPRi labels, whereas AlphaGenome is applied without fitting to those labels. The meaningful question is not whether AlphaGenome must beat EGrf; it is whether AlphaGenome provides held-out information beyond EGrf.

The author repository makes a bounded comparison feasible:

- `Astrocyte_trainingData_pluspredictions.csv` contains exactly the same 2,307 `Pair` identifiers as the frozen OP26 cohort, with no missing EGrf score and no label mismatch.
- The author code defines EGrf using H3K4me3, H3K27ac, nearest-gene status, distance, ATAC signal, the number of tested TSSs per enhancer and gene stability.
- The supplied EGrf predictions are not training-set fitted scores. `crossFoldGeneLevel()` holds out both target-gene and enhancer groups before predicting each test intersection.
- Caveat: the author code tunes RF hyperparameters with ordinary full-dataset caret cross-validation before generating the gene-by-enhancer held-out predictions. The published scores are therefore substantially cross-fitted but not a fully nested, untouched evaluation.

Required comparison:

1. retain the author EGrf score as the faithful expert-comparator readout;
2. report its raw AP, ROC AUC and log loss on the identical frozen cohort;
3. use the same five enhancer-grouped outer folds to calibrate EGrf alone and EGrf plus AlphaGenome;
4. quantify paired enhancer-bootstrap changes in AP and log loss;
5. state the asymmetric training conditions and the hyperparameter-selection caveat.

### Gene-level confounding — MUST DO

The same-enhancer analysis conditions out enhancer/request effects, but a mirror analysis is needed because some target genes may be intrinsically easier for the model. The frozen cohort contains enough repeated structure for a bounded test: genes with at least one functional and one well-powered nonfunctional enhancer can support a within-gene contrast, conditional AUC and label permutation.

Required checks:

- relation-weighted and gene-equal Spearman association between AlphaGenome deletion strength and author-measured gene expression;
- expression-only discrimination/enrichment of functional labels;
- equal-gene functional-minus-nonfunctional contrast and conditional AUC;
- within-gene permutation and gene bootstrap;
- a cubic-spline sensitivity for measured expression in the context model.

The frozen AlphaGenome score table does not retain an absolute reference RNA-seq or promoter-CAGE output. Author-measured gene expression is the available baseline-expression covariate; no new AlphaGenome endpoint will be selected post hoc.

### Timeline wording — MUST CORRECT, not strengthen

The third-party claim that AstroREG CRISPRi labels were not public at AlphaGenome release is contradicted by primary-source dates:

- AstroREG manuscript received: 2024-09-24;
- the authors' public GitHub history places the processed AstroREG training table, EGrf predictions and RF code in commit `0cda01f` on 2024-02-22, with an update on 2024-10-16;
- GEO GSE236057 status: public on 2025-01-01;
- AlphaGenome public API announcement and bioRxiv v1: 2025-06-25;
- AstroREG accepted: 2025-10-23 and published online: 2025-12-18;
- AlphaGenome Nature publication: 2026-01-28.

Therefore the manuscript must not claim that AstroREG labels postdated AlphaGenome's public release or were necessarily unseen. Exact processed relation labels and author EGrf predictions were already present in a public repository before AlphaGenome's June 2025 release. The defensible statement is only that AstroREG is external to OP26 and is not identified as a named AlphaGenome benchmark in the public materials reviewed. Public availability does not establish actual use in model training or evaluation, but it removes any defensible post-release-label claim.

## Phase-1 stop decision

Proceed with EGrf and gene-level analyses. Do not search for or run a second validation dataset until those results are known. A second dataset will be considered only by the prespecified value-of-information questions and the one-dataset stopping rule.
