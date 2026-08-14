# OP-26 external AstroREG dry-lab validation: frozen protocol

Freeze date: 2026-08-08, before any AlphaGenome outputs for this cohort were generated or inspected.

## Purpose

Test whether the narrow K562 result generalizes to a different cell context and a different CRISPRi truth set. This is a dry-lab external validation using the public author data from Green et al., *Nature Neuroscience* (2026), DOI `10.1038/s41593-025-02154-3`. It does not add or imply new wet-lab work.

## Frozen cohort

- Source: the authors' `TrainingDataframe_Astrocytes.csv` and `CRISPRi_Power Simulation - Power per EGP.csv` from `Voineagulab/astrocyte_crispri`.
- Genome build: hg38, as specified by the author analysis repository.
- Positive relation: `HitPermissive_NegZ == TRUE`, matching the authors' definition of a significant CRISPRi hit with negative expression change.
- Negative relation: `HitPermissive_NegZ == FALSE` and `WellPowered015 == TRUE`, matching the authors' well-powered nonfunctional definition.
- Expected frozen universe: 2,307 unique enhancer-gene relations, comprising 133 positives and 2,174 negatives across 745 unique enhancers. Any mismatch stops execution.
- No records will be added, removed, or relabelled after AlphaGenome results are viewed except for a documented technical failure such as an unresolvable gene identifier.

## Frozen model request

- AlphaGenome client: the same local v0.7.0 client used in Step 6.
- Model output: `GeneMaskLFCScorer(OutputType.RNA_SEQ)`.
- Biological context: primary astrocyte total RNA-seq, ontology `CL:0000127`.
- Perturbation: VCF-style deletion of the full tested enhancer interval with its immediately preceding reference base retained as the anchor.
- One API request per unique enhancer. Each 1,048,576-bp request interval must contain the enhancer and every eligible target-gene TSS linked to that enhancer. The returned gene-level score is then extracted separately for each eligible target ENSG identifier.
- Primary score: deletion strength for the target gene, defined as the negative of AlphaGenome's ALT-versus-REF RNA-seq gene-mask log-fold-change score. Larger values therefore mean a larger predicted loss of target-gene expression after enhancer deletion.
- API credentials are loaded from the existing local `.env`; the key is never printed or stored in an output.

## Frozen analyses

1. Discrimination across all 2,307 relations: area under the precision-recall curve (primary) and ROC area (secondary), with the observed positive prevalence shown as the no-skill PR baseline.
2. Direction on the 133 positive relations: proportion with predicted ALT-versus-REF RNA-seq score below zero; also report the median score and its distribution.
3. Effect separation: difference in mean deletion strength between positives and negatives. Use an enhancer-cluster bootstrap as primary uncertainty and a target-gene-cluster bootstrap as sensitivity.
4. Context-adjusted association: logistic regression of CRISPRi label on standardized log distance, log ABC score, enhancer length, author-measured gene expression, number of assayed cells, and AlphaGenome deletion strength. Report the AlphaGenome coefficient and grouped cross-validated change in log loss and average precision when it is added to the context-only model.
5. Matched sensitivity: one negative per positive, selected without replacement before AlphaGenome outputs are used by minimum standardized distance on log genomic distance, log ABC score, enhancer length, gene expression and assayed-cell count. Matching must not use any AlphaGenome value.
6. Missingness: report target genes absent from the AlphaGenome score output and technical API failures explicitly. Do not treat them as zero.

Bootstrap and grouped cross-validation resampling are seeded with `20260808`. The enhancer is the primary dependence unit because one deletion can generate several enhancer-gene rows; target-gene clustering is reported as a sensitivity because genes can be linked to several enhancers.

## Interpretation boundary

A positive result would show external, cell-context-relevant computational association with experimentally supported enhancer-gene effects. It would not prove causal enhancer-gene assignment for an arbitrary locus, would not establish independence from AlphaGenome's track-level training data, and would not replace experimental validation. A null result would be treated as a failed external generalization, not hidden by returning to the K562-only analysis.

The publication novelty claim is limited to the matched-tested-negative, target-relation specificity and context-increment questions. It must explicitly acknowledge the 2026 public Koo Lab benchmark of AlphaGenome on Fulco and Gasperini CRISPRi screens and AlphaGenome's own ENCODE-rE2G enhancer-gene evaluation.
