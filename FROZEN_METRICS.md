# Frozen evaluation contract: 24 November 2026

The user's 12-estimator mapping and both supplied specifications are preserved
verbatim in `specifications/metric-suite-2026-11-24.txt` and
`specifications/metric-equations-2026-11-24.txt`. The package and software-paper
title remains **svylearn: Design-consistent statistical learning, calibration,
and inference for complex surveys**. ECLS-K:2011 is the confirmed longitudinal
companion to NHANES. The broader contract is `FROZEN_EXTENSION_SPEC.md`.

| Core | Stable return keys | Operational definition |
|---|---|---|
| 1 | `optimism_bias` | Unweighted loss of supplied IID out-of-fold predictions minus weighted loss of supplied design-aware out-of-fold predictions |
| 2 | `deff_loss` | Design variance of mean loss divided by unweighted sample loss variance / n (SRS with replacement) |
| 3 | `n_eff_pred` | n / `deff_loss`; undefined for a zero SRS reference variance |
| 4 | `ici_svy`, `e90_svy` | Weighted mean and inverse-ECDF 90th percentile of absolute spline-smoothed calibration error |
| 5 | `rake_residual` | (Weighted predicted event total - externally known event total) / known event total |
| 6 | `transport_smd`, `transport_div` | Covariate SMD before/after supplied transport weighting, with a fixed pooled reference SD; weighted squared difference between two prediction vectors on the target |
| 7 | `conf_coverage`, `conf_width` | Weighted empirical test coverage and interval width or classification set cardinality |
| 8 | `sensitivity_gamma` | Explicitly labeled bounded-weight sensitivity tipping functional; report the target and whether the crossing is bracketed |
| 9 | `net_benefit_svy` | (Weighted TP - threshold odds times weighted FP) / sum(weights), also reported per 100,000 |
| 10 | `pop_tp`, `pop_fp`, `pop_tn`, `pop_fn`, `nns_svy` | Weighted confusion totals; supplied NNS formula (TP+FP)/TP retained as positive follow-ups per true positive; add `nns_population` = sum(weights)/TP |
| 11 | `domain_cal_gap` | Maximum absolute domain joint-logistic intercept and maximum absolute slope deviation from one; not a pairwise fairness test |
| 12 | `shap_svy`, `vip_psu` | Weighted mean absolute supplied SHAP contributions and mean within-PSU out-of-fold permutation loss increase |

Foundational keys include `loss_svy`, `ipa_svy`, `cal_intercept`, `cal_slope`,
`singleton_audit`, `ci_lower`, and `ci_upper`. `cal_intercept` is calibration in
the large (offset model with slope fixed to one); the free intercept of the
joint intercept/slope model is additionally identified. For continuous outcomes,
calibration uses the corresponding linear model. Binary-only metrics must be
marked inapplicable to other outcomes, not fabricated. Missing optional inputs
must have explicit availability reasons, not numerical zero substitutes.

## Scientific clarifications that form part of this contract

- Normalized weighted means are Hajek ratio estimators; weighted totals are
  Horvitz-Thompson estimators only when the weights have the appropriate inclusion
  interpretation. Rescaled weights cannot identify population headcounts without
  a supplied population size. Longitudinal person-wave totals count records,
  not necessarily unique people.
- The IID-versus-survey contrast combines weighting, partitioning, fitting and
  calibration changes. It is not proof of true finite-population optimism or
  attribution to a single mechanism. Its sign may be either positive or negative.
- The supplied unweighted SRS variance reference is frozen explicitly; it can
  differ from `survey::deff(..., deff="replace")`, whose reference is different.
  Effective sample size is an estimated variance equivalence, not a power theorem.
- ICI/E90 depend on the chosen smoother. Replication must refit the smoother.
  Calibration intercept in the large and the free joint intercept are distinct.
- Raking covariates does not ensure event-total calibration. A known event total
  must be provided separately. `transport_div` follows the supplied displayed
  prediction-versus-prediction formula; prediction-versus-outcome Brier loss is
  a separate quantity. Small SMD does not prove transportability.
- Weighted split conformal retains the test-point atom at infinity. Exact
  marginal theory requires exchangeability or appropriate known density ratios
  and a separately fitted prediction rule. Arbitrary design weights, pooled
  cross-validation residuals and dependent panel rows do not automatically meet
  those conditions. Empirical coverage is evaluated on a separate test set.
- Minimizing a bounded-weight mean is a linear-fractional problem. Minimizing its
  recomputed Wald lower bound is a different, generally nonconvex problem.
  A mean-minimizing tilt with a conditional interval must not be advertised as a
  globally optimal confidence-bound tipping point. Both targets must be labeled.
- Taylor/replicate intervals are estimated, usually asymptotic, conditional on
  supplied predictions unless the entire fitting procedure is repeated. Pointwise
  intervals across thresholds are not simultaneous confidence bands. Maxima,
  quantiles and adaptive smoothers need separately documented inference.
- Aggregate SHAP depends on the supplied model, background and explanation
  algorithm. Within-PSU permutation is associational. Panel trajectory permutation
  requires compatible whole-subject blocks, not arbitrary wave shuffling.
- Novelty, exact finite-population coverage, regulatory endorsement and journal
  acceptance are not established by these specifications or by implementation.

Primary references: the `survey` package documentation
(https://r-survey.r-forge.r-project.org/pkgdown/docs/reference/surveysummary.html),
Austin and Steyerberg's ICI paper (https://pmc.ncbi.nlm.nih.gov/articles/PMC6771733/),
and Tibshirani et al.'s weighted conformal paper
(https://proceedings.neurips.cc/paper_files/paper/2019/hash/8fb21ee7a2207526da55a679f0332de2-Abstract.html).
