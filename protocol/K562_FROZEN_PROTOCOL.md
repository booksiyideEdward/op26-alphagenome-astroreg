# OP-26 frozen full-study protocol

**Freeze date:** 2026-08-07, before any new 161-pair AlphaGenome output was generated or inspected.

## Scientific question

In K562 endogenous CRISPRi enhancer–gene pairs, does AlphaGenome's distal perturbation response contain target-gene-specific, experimental-direction-concordant information beyond ABC and simple sequence/context variables?

## Frozen primary set and sealed subset

- Primary set: the existing 161 positive–matched-tested-negative pairs in `FEASIBILITY_2026-08/OP26/matched_pairs_hg38.csv`, unchanged.
- Independent unit: target gene/regulatory-locus cluster, not rows, TSS bins, tracks, or model calls.
- Sealed confirmation subset: the existing 33 pairs from 24 genes labelled `confirmation` in that manifest.
- The previous 20-pair smoke test used development genes only. It must have zero gene overlap with the sealed subset.
- The repository's original row split is descriptive only because target genes overlap across train and test.

## Frozen model and endpoint

- Model: AlphaGenome API, Python client `alphagenome==0.7.0`; local client repository commit `95cdbfce7981411453e5e094519bcf0605720199`.
- Sequence length: 1,048,576 bp.
- Cell ontology: K562, `EFO:0002067`.
- Output: the strand-matched `hCAGE EFO:0002067` track.
- Target window: CAGE sum within +/-500 bp of the experimentally linked target-gene TSS.
- Perturbation: VCF-style deletion with the left anchor retained, identical to the completed smoke test.
- Response: `delta_log1p_cage = log1p(alternate) - log1p(reference)`; deletion reliance/strength is `-delta_log1p_cage`.
- Each pair uses the same fixed 1-Mb interval construction as the smoke test: the smallest span containing the positive element, matched negative element, and linked TSS is centered in a 1-Mb request.

## Frozen analyses

1. Primary: positive minus matched-negative deletion strength, averaged within gene and then across genes. Report a gene-cluster bootstrap 95% CI, gene-level sign test, and gene-level sign-flip permutation result.
2. Sealed confirmation: report the same estimate on the 33 pairs/24 genes separately; do not pool it into a claim of independent generalization because exact model-training-label overlap is unresolved.
3. CRISPR direction: for positive elements, compare the sign of AlphaGenome `delta_log1p_cage` with the sign of experimental `gene_expr_change`; use genes as the inference units.
4. Target specificity: compare the same positive-element deletion at the linked target with a geometry-selected, distance-matched non-target TSS. To keep the control symmetric, the linked-control and wrong-gene-control calls each use a 1-Mb interval constructed from that element and the corresponding TSS only. Wrong genes are selected before new model output and cannot be selected using predicted response. These control calls are separate from, and cannot replace, the frozen primary positive-versus-negative calls.
5. ABC adjustment: matched conditional logistic analysis using within-pair differences. Compare a baseline using ABC plus frozen matching covariates with a full model that additionally includes AlphaGenome deletion strength. Report cluster-bootstrap uncertainty and development-to-sealed held-out log-likelihood increment.
6. GC sensitivity: construct a separate outcome-blind within-gene match using distance, ABC, region length, and sequence GC. It is a sensitivity set and never replaces the frozen primary 161 pairs.
7. Overlap/contamination: report source-label/model-development overlap as known, excluded, or unresolved. Do not call the sealed subset independent model generalization unless overlap is excluded.

## Interpretation boundary

Positive-versus-negative separation alone supports at most selective sensitivity to functional CRE sequence. A claim of correct distal use additionally requires coherent CRISPR direction, linked-target superiority over the same-element wrong-gene control, incremental information beyond ABC/context, and no direction reversal in the sealed subset.

No second DNA model, attribution expansion, genome-wide expansion, sequence optimization, or post-result hypothesis change is permitted in this execution.
