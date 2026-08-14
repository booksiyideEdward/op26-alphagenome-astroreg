# OP26 source and timeline audit

## Conclusion

The third-party assertion that AstroREG supplied labels created only after AlphaGenome's release is **factually incorrect**. The exact processed AstroREG training table, the 2,307-row author EGrf prediction table and the relevant random-forest code were present in the public `Voineagulab/astrocyte_crispri` repository by **22 February 2024**. GEO GSE236057 became public on **1 January 2025**. AlphaGenome was publicly announced and made available by API on **25 June 2025**.

This chronology does not show that AstroREG labels were used by AlphaGenome. It does show that OP26 cannot call AstroREG a post-release label set or infer temporal independence.

## Evidence table

| Event | Date | Primary evidence | Implication |
|---|---|---|---|
| Public repository commit adds processed AstroREG training data, EGrf predictions and model scripts | 2024-02-22 | Git commit `0cda01f3b07c1eff4802a1381b66f4cde2daf3a0`; files include `TrainingDataframe_Astrocytes.csv`, `Astrocyte_trainingData_pluspredictions.csv`, `Functions.R` and `RunRFmodels.R` | Exact labels and comparator outputs were public before AlphaGenome release |
| EGrf prediction file updated | 2024-10-16 | Git commit `05d7fa3cd82b286e7a041baed69d6644bb916894` | The final public comparator also predates AlphaGenome release |
| AstroREG manuscript received | 2024-09-24 | Nature Neuroscience article history | Manuscript submission also predates AlphaGenome public release |
| GEO GSE236057 submitted | 2023-06-28 | NCBI GEO accession record | Raw/processed-data submission began well before release |
| GEO GSE236057 public | 2025-01-01 | NCBI GEO accession record | Public experimental data predate AlphaGenome release |
| AlphaGenome manuscript received | 2025-05-16 | Nature article history | Public AstroREG repository/GEO availability predates even this submission date |
| AlphaGenome public announcement and API | 2025-06-25 | Google DeepMind announcement; bioRxiv v1 identifier `10.1101/2025.06.25.661532` | Reference date for public model availability |
| AstroREG accepted | 2025-10-23 | Nature Neuroscience article history | Acceptance occurred after AlphaGenome release, but label availability did not |
| AstroREG version of record | 2025-12-18 | Nature Neuroscience article history | Publication date is not the label-availability date |

## EGrf source-method audit

The author `crossFoldGeneLevel()` function:

- randomly divides unique genes into ten folds and unique enhancers into ten folds;
- predicts each gene-fold × enhancer-fold intersection;
- trains after excluding every row sharing the held-out gene fold or enhancer fold;
- uses `ranger` probability predictions.

The author workflow first tunes RF hyperparameters with `caret`, then reuses the selected settings inside the gene×enhancer held-out loop. Therefore the released predictions are cross-fitted at relation level but not the product of fully nested hyperparameter selection.

## Sources used

- Public repository: `https://github.com/Voineagulab/astrocyte_crispri`
- Decisive commit: `https://github.com/Voineagulab/astrocyte_crispri/commit/0cda01f3b07c1eff4802a1381b66f4cde2daf3a0`
- GEO: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE236057`
- AstroREG article: `https://www.nature.com/articles/s41593-025-02154-3`
- AlphaGenome public announcement: `https://deepmind.google/blog/alphagenome-ai-for-better-understanding-the-genome/`
- AlphaGenome article: `https://www.nature.com/articles/s41586-025-10014-0`

The repository was cloned only to a temporary local directory for read-only history inspection. The OP26 package retains the exact 2,307-row author prediction table used for the comparator, not a redundant repository copy.
