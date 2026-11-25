# Capability evidence through the 0.3.0 local implementation

This ledger supplements the frozen requirements; it does not declare the entire
editorial wish list implemented. The software title and authors remain unchanged.

## Twelve-estimator evaluation contract

The stable keys are exposed by `svy_eval()` / `evaluate()` and documented by
`svy_glossary()` (17 reporting rows). All supplied source texts are preserved
verbatim with SHA-256 hashes in `validation/metric_contract_audit.json`.

| Core | Implemented interface | Evidence and precise scope |
|---|---|---|
| 1 | `svy_eval(..., iid_predictions=)` | Observed IID-minus-design OOF loss contrast; caller establishes IID prediction provenance |
| 2–3 | `svy_eval()` | Frozen unweighted SRS variance reference and predictive effective n; exact PSU-sum oracle and survey comparison |
| 4 | `svy_calibration_error()`, `svy_eval()` | Weighted spline ICI/E90; supplied replicates refit the smoother; no training refits |
| 5 | `svy_rake()`, `svy_rake_audit()` | Separate demographic constraints and known event-total discrepancy; direct margin checks |
| 6 | `svy_transport()` | Before/after SMD for supplied transport weights and target-row prediction divergence; automatic IOSW/AIPW fitting remains in the broader specification |
| 7 | `svy_conformal()`, `svy_conformal_eval()` | Regression and classification split construction and separate empirical test evaluation; 600,000 test predictions across three assumption settings; cross-fitted and censored-outcome constructions remain separate work |
| 8 | `svy_weight_bounds()`, `svy_sensitivity()` | Sharp bounded-weight mean extrema, tested by exhaustive box vertices; explicitly labeled mean-minimizing CI path. Global optimization of the CI lower bound is not claimed |
| 9–10 | `svy_decision_infer()` | Conditional utility, counts and ratio intervals; both supplied follow-up and population-screening denominators; 8,000 finite-population samples |
| 11 | `svy_domain_calibration()` | Joint covariance, supported-domain coefficients, pairwise contrasts and Bonferroni maximum-gap envelopes; descriptive calibration audit |
| 12 | `svy_shap()`, `svy_importance_infer()` | Supplied SHAP aggregation and OOF within-PSU permutation with whole compatible panel trajectories; conditional design intervals and separate permutation Monte Carlo SE |

## Additional frozen extensions advanced here

- A07: weighted thin-SVD PCA is trained inside each recipe; covariance reconstruction
  and P greater than N checks pass. Embeddings can enter as externally extracted
  numeric predictors; image/text encoders themselves are not included.
- A08: trajectory summaries require an explicit cutoff and reject duplicate
  subject-wave records. Tests intervene on future observations.
- A09: optional dependence-group columns are joined transitively with sampling
  clusters for validation. Inference identifiers remain unchanged.
- D01: ECLS-K:2011 is downloaded from NCES, verified and used in two completed
  workflows: 8,487 fifth-grade observations and 25,453 grade 3–5 child-wave records.
  The official longitudinal weight, 80 JK2 replicates and 80 agency-specified
  degrees of freedom are recorded. Public original PSU IDs are suppressed;
  validation uses released variance units.

## Broader scope still requiring implementation and separate validation

The 0.3.0 extension adds `svy_attrition()` and `svy_attrition_weights()` for
estimated sequential monotone retention. Passing the specification into
`svy_fit()` refits the nuisance models inside every nested partition, for all
four response families. Baseline-roster replicate weights trigger retention
refitting as part of complete-workflow replication. Tests intervene on held-out
outcomes and response histories and compare direct weighted logistic fits.
This is an observed-history MAR model, distinct from MNAR identification or
automatic population transport-weight estimation.

Ordinal/Dirichlet, Huber/quantile/distributional, survival/competing-risk and
negative-binomial/zero-inflated model pipelines are not relabeled as supported
by these metric additions. Automatic IOSW/AIPW, DR/R causal learning, E-values,
Rosenbaum and missing-outcome tipping models, federated estimation and full
Python model parity remain in the frozen extension specification. Public GitHub,
CRAN, arXiv and JSS publication is separate from a checked local release.

The source records retain the observed undercoverage of sparse Wald utility
intervals and failure of arbitrary conformal weights under shift. These findings
are part of the evidence, not omitted exceptions.
