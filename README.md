# OP26: AlphaGenome enhancer–gene relation evaluation

This repository accompanies the OP26 preprint on AlphaGenome deletion responses and enhancer–gene relation prediction in primary human astrocytes. It contains the manuscript, analysis records, focused analysis code, OP26-derived score tables, machine-readable results and final figures.

The main result uses 2,307 public AstroREG enhancer–gene relations (133 functional and 2,174 well-powered nonfunctional relations). AlphaGenome deletion strength is evaluated alone, beyond a frozen activity/contact/context model, and as a complementary feature to the authors' supervised EGrf score. K562 analyses are retained as supporting evidence because exact joins establish overlap with released Gasperini and ENCODE-rE2G resources.

## Repository map

- `manuscript/`: release-candidate manuscript.
- `protocol/`: project-internal freeze records and the later targeted audit plan.
- `data/derived/`: minimized OP26-derived AlphaGenome relation-level scores. Third-party source tables are not redistributed.
- `analysis/`: focused AstroREG, EGrf, within-enhancer and permutation-null code plus machine-readable outputs.
- `results/k562/`: supporting K562 summaries and derived results.
- `figures/`: the three final manuscript figures in PDF.
- `docs/`: chronology, policy and reproducibility notes.

## Reproduce the AstroREG analyses

Python 3.11 was used for the release calculations. Install the recorded environment:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Retrieve the exact public EGrf input from the authors' repository and build the local analysis table:

```bash
python scripts/prepare_public_inputs.py --download-egrf
```

Then run the focused analyses:

```bash
python analysis/EXTERNAL_ASTROREG/analyze_external_astroreg.py
python analysis/analyze_astroreg_within_enhancer.py
python analysis/TARGETED_REVIEWER_AUDIT/run_targeted_audit.py
python analysis/TARGETED_REVIEWER_AUDIT/run_permutation_null.py
```

The last command runs the fixed 1,000-permutation structure-aware shuffled-AlphaGenome stacking diagnostic and does not call the AlphaGenome API.

## Analysis chronology

The K562 and AstroREG protocol files were maintained inside the project before the corresponding result inspection. The EGrf and gene-confounding plan was recorded later as a targeted post-hoc audit before those module outputs were inspected. These records were maintained contemporaneously within the project but were not externally preregistered or third-party timestamped. The first public GitHub commit is a release archive, not retrospective evidence of preregistration.

## Data and rights boundary

AstroREG source tables are public but the source repository does not declare a repository-wide license. They are therefore retrieved directly from the authors rather than copied here. AlphaGenome client code is Apache-2.0; use of the AlphaGenome API remains subject to Google DeepMind's non-commercial terms. This repository contains relation-level OP26 numerical outputs, not model weights, raw prediction tracks, source repositories or credentials. See `data/README.md` and `docs/REPRODUCIBILITY_AND_RIGHTS.md`.

No project license is asserted in this release candidate. License selection remains a human release action if the authors intend downstream reuse beyond GitHub's default copyright terms.
