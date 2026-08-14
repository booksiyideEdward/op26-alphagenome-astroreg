# Data sources and reconstruction

## AstroREG source data

The decisive cohort comes from Green, Sutton, Pérez-Burillo et al., *Nature Neuroscience* 29, 703–716 (2026), DOI `10.1038/s41593-025-02154-3`, with public study accessions including GEO `GSE236057`.

The exact EGrf table used in OP26 is maintained by the authors at:

`https://github.com/Voineagulab/astrocyte_crispri/blob/main/6_PredictiveModels/3.Predictions/RF_Results/Astrocytes/Astrocyte_trainingData_pluspredictions.csv`

Exact source path:

`6_PredictiveModels/3.Predictions/RF_Results/Astrocytes/Astrocyte_trainingData_pluspredictions.csv`

The AstroREG repository was public at release audit but did not declare a repository-wide license. The table is not duplicated here. Run `python scripts/prepare_public_inputs.py --download-egrf` to retrieve it directly from the official repository and reconstruct the local analysis input.

The join key is the author's `Pair` field to OP26 `relation_id`. The reconstruction script additionally requires exact agreement of `EnsID` with `gene_id`, `Enh` with the enhancer identifier and `HitPermissive_NegZ` with the frozen binary label. It stops if the 2,307 rows do not reconcile one-to-one.

## OP26-derived AstroREG scores

`derived/astroreg_alphagenome_relation_scores.csv` contains only OP26 relation identifiers, model-track provenance, relation-level AlphaGenome RNA-seq log-fold-change scores, sign-reversed deletion strengths and completion status. It excludes the original author covariates and raw AlphaGenome tracks. The scores were produced for non-commercial research use with the AlphaGenome primary-astrocyte total-RNA-seq output and one request per enhancer, as documented in `protocol/ASTROREG_FROZEN_PROTOCOL.md`.

## K562 sources

The supporting K562 analysis traces to Gasperini et al., GEO `GSE120861`, DOI `10.1016/j.cell.2018.11.029`, and a DNALONGBENCH-derived matched table. Public benchmark-overlap audits used the Murphy–Koo archive, DOI `10.5281/zenodo.21383205`, and the ENCODE-rE2G resources accompanying DOI `10.1038/s41586-026-10781-4`.

The K562 third-party source tables are not redistributed. `derived/k562_*_scores.csv` provides only minimized OP26 model outputs needed to inspect the supporting contrasts. Machine-readable summary results are under `results/k562/`.

## Credentials

No AlphaGenome API key, GitHub token, institutional-access credential or browser session is included. Reproducing already published statistical analyses uses the committed numerical scores and does not require an AlphaGenome API call.
