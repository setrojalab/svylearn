# Frozen extension specification: 24 November 2026

The user requested that the supplied editorial brief be frozen and implemented.
The verbatim source is `specifications/editorial-brief-2026-11-24.txt`.
It records the user's requested scope; it is not independent evidence of journal
acceptance, regulatory endorsement or mathematical guarantees.

The title remains exactly **svylearn: Design-consistent statistical learning,
calibration, and inference for complex surveys**. Both cross-sectional and
longitudinal designs are required. The published method retains its own title
in citations. Existing released evidence remains attributed to the implementation
that generated it; extensions need their own validation.

## Fixed scientific interpretation

1. Scores must identify their target. Logarithmic and Brier scores assess
   probabilities, RPS assesses ordered category distributions, and quantile and
   Huber losses assess specified point functionals. Partial likelihood is not a
   complete survival-distribution scoring rule.
2. Conditional design-based intervals, repeated-sampling empirical coverage,
   conformal marginal coverage and conditional coverage are different objects.
   Universal nontrivial distribution-free conditional coverage is not attainable
   without restrictions: https://arxiv.org/abs/1903.04684.
3. Survey inclusion weights are not automatically the density ratios required by
   weighted conformal theory. Any coverage guarantee must name the sampling or
   exchangeability assumptions and account for the test point's probability mass.
4. Causal and transport estimators must state consistency, overlap, exchangeability,
   nuisance-model and target-population assumptions. Sensitivity analyses vary
   specified assumptions; they do not establish absence of hidden bias.
5. Aggregate sharing is distinct from differential privacy. No regulatory
   certification, legal compliance, clinical benefit, universal optimality or
   adoption guarantee follows from installing a package.

## Required interfaces and acceptance evidence

| ID | Frozen capability | Required evidence before declaring complete |
|---|---|---|
| A01 | Shared outcome/scoring/calibration protocol; preserve design/split/tune/calibrate/inference workflow | Explicit response and prediction schema; unsupported inputs fail clearly; stable existing API and binary results |
| A02 | Binary/rare-event log loss, Brier/IPA and Platt/isotonic/spline calibration | Existing independent score checks and sparse-fold evidence remain passing |
| A03 | Multiclass and ordinal cross-entropy/RPS; temperature and Dirichlet scaling | Probability/simplex checks, known ordinal RPS oracle, independent optimizer checks and outer-outcome isolation |
| A04 | Continuous/skewed Gaussian, Huber and quantile losses; distributional PIT diagnostics | Hand calculations, quantile ordering and finite-support checks; distinguish mean from distribution calibration |
| A05 | Right-censored survival and competing risks; hazards, IPCW Brier, cause-specific/Fine-Gray cumulative incidence evaluation | Known survival/CIF curves, censoring left-limit convention, competing-event accounting and leakage-safe nuisance estimation |
| A06 | Poisson, negative-binomial and zero-inflated counts | Independent likelihood/deviance calculations, zero-only boundary and dispersion checks |
| A07 | Weighted training-only PCA/SVD, mixed predictors and pre-extracted embeddings/omics | Frozen transformations, P much greater than N test, weighted covariance reconstruction and validation-set intervention |
| A08 | Longitudinal trajectory features | Subject/wave alignment, explicit time eligibility, duplicate rejection and no use of future records when a cutoff is specified |
| A09 | Survey, site/subject, spatial-block and graph dependency boundaries | Entire declared dependency groups stay in one fold; inference retains the original sampling design; impossible partitions are reported |
| A10 | External weighted learner protocol | Explicit sample-weight and prediction contracts; custom and installed-adapter checks; no engine loaded at attachment |
| B01 | IOSW and augmented transport/generalizability estimation | Known-density-ratio oracle, target-weighted balance/overlap diagnostics and independent-source/target variance calculation |
| B02 | Weighted split and cross-fitted conformal prediction | Test-point mass, infinite-boundary behavior, marginal coverage simulations and explicit conditions; no universal conditional-coverage claim |
| B03 | Orthogonal DR/R learners, ATE/CATE and policy evaluation | Cross-fitted nuisance provenance, positivity checks, randomized-treatment oracle, nuisance misspecification simulations and honest policy-value target |
| B04 | E-values, Rosenbaum bounds and MNAR tipping curves | Published formula checks; explicit supported estimand/scale; observed-data invariance and missing-data sensitivity grid |
| B05 | Federated summaries, influence-function aggregation and meta-analysis | No row pooling, agreement with centralized sufficient-statistic calculations, independent-shard assumptions and explicit privacy properties |
| C01 | Adaptive evaluate/inference reports | Outcome-appropriate accuracy, calibration, uncertainty, utility and subgroup sections; no binary metrics silently assigned to unrelated outcomes |
| C02 | Calibration/transportability diagnostics | Mean/slope or distribution-specific calibration, ICI/E90 where defined, weighted target/source SMD |
| C03 | Uncertainty and DEFF | Agreement with survey linearization/replication; conditional nuisance treatment disclosed; no exact universal finite-sample claim |
| C04 | Decision utility, screening/follow-up counts and budget inputs | Hand-calculated weighted TP/FP/FN and utility checks; explicit NNS denominator and user-supplied cost/benefit assumptions |
| C05 | Equity and explainability | Subgroup support checks, descriptive parity/odds gaps, covariance-aware contrasts where intervals are reported, SHAP intervals |
| D01 | Real longitudinal application alongside NHANES | Authorized public-use source, exact data hashes, documented weights/design, repeated-person workflow and separate empirical results |
| D02 | Release and manuscript | Reproducible tests/simulations, current capability ledger, small lazy-loading package, preserved title and method citation, truthful release status |

## Baseline and completion accounting

The baseline is local commit `4550f69`, with four response families, 262 passing
test expectations, 94.60% measured line coverage and the previously documented
simulation studies. This specification does not relabel unimplemented extensions
as completed. Acceptance evidence will be recorded by capability ID as work proceeds.

For the longitudinal application, the user confirmed ECLS-K:2011 as the open
public-use panel companion (nine rounds; weights and design documentation).
The 12-estimator evaluation contract is additionally frozen in `FROZEN_METRICS.md`.
Source descriptions:
https://nces.ed.gov/ecls/dataproducts.asp and
https://hrsdata.isr.umich.edu/data-products/rand.
