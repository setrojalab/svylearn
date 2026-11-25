# svylearn 0.3.0

- Add `svy_attrition()` specifications and `svy_attrition_weights()` for sequential
  retention weighting from complete baseline panel rosters, with automatically
  lagged histories, optional wave-specific stabilization and clipping audits.
- `svy_fit(attrition = ...)` refits retention and preprocessing in every inner and
  outer training partition for all four outcome families. It retains full-roster
  folds and repeat-specific evaluation weights. Full-workflow replication also
  refits the retention models from baseline-roster replicate weights.
- Clarify sequential-MAR assumptions, positivity, pooled-wave stabilization targets,
  estimated-weight uncertainty and the distinction from MNAR sensitivity models.

# svylearn 0.2.0

* Froze 12 core estimators and 17 reporting rows in `FROZEN_METRICS.md`, preserving
  the supplied specifications and exact programmatic keys.
* Added `svy_eval()` / `evaluate()` and `svy_glossary()` with explicit metric
  availability, conditional uncertainty and plain-language interpretations.
* Added population utility/count/ratio intervals, raking residual audits,
  transport balance, split-conformal construction and independent test evaluation,
  sharp bounded-weight mean limits and separately labeled conditional CI tipping.
* Added joint-covariance domain calibration contrasts, conservative maximum-gap
  envelopes, SHAP intervals and trajectory-preserving within-PSU importance.
* Added training-only weighted thin-SVD PCA, explicit time-cutoff trajectory
  features and transitive validation dependence groups that preserve inference IDs.
* Added separate finite-population and conformal simulation evidence. Historical
  0.1.0 results retain their original attribution; negative findings are retained.
* Confirmed ECLS-K:2011 as the real longitudinal companion to NHANES and added
  official-source download and integrity tooling.

# svylearn 0.1.0

* Preserved the exact software title from the Figshare reference in the supplied
  method article; the method article keeps its own title in citations.
* Added gaussian, poisson and multinomial nested workflows, family-specific
  scores, calibration, conditional loss inference and full replicate refits.
  Added independent engine checks, outer-outcome isolation tests and a separate
  response-family simulation study.

* Added a dependency-light R implementation of whole-PSU nested validation,
  weighted preprocessing, learner selection, convex stacking and calibration.
* Added explicit cross-sectional and repeated-subject design metadata.
* Delegated multistage, FPC and replicate-weight inference to survey.
* Added conditional fixed-score inference, optional complete replicate refits,
  subgroup support checks, within-PSU permutation summaries and SHAP aggregation.
* Added reproducible simulation, source-integrity and R/Python kernel checks.
* Added the originating method citation and software manuscript sources.

## Reproducibility verification

- Every nested training pipeline has an isolated, restored random seed, including
  stochastic learner calls and outcome-dependent fallback paths.
- Portable XGBoost model bytes support prediction after loading a saved fit in a
  fresh R process while keeping engine namespaces lazy at package attachment.
- Nonconverged candidates are recorded and excluded; invalid candidates cannot
  contaminate selection or stacking with missing predictions.
