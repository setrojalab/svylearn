# Frozen blueprint: implementation and scientific interpretation

The additional requirements supplied on 24 November 2026 are frozen separately
in `FROZEN_EXTENSION_SPEC.md`, with the original wording retained in
`specifications/editorial-brief-2026-11-24.txt`. That extension specification adds
to this baseline; its requested features are not treated as completed without
their stated acceptance evidence. The 12 core estimators, stable keys and precise
scientific definitions are additionally frozen in `FROZEN_METRICS.md`.

The confirmed software title is fixed as **svylearn: Design-consistent statistical
learning, calibration, and inference for complex surveys**. The package and
software paper cover both **cross-sectional and longitudinal complex survey
data**. Longitudinal support includes repeated subject-wave records, supplied
wave or attrition weights, and cluster partitions that keep a subject's waves
together. The originating method article retains its own title in citations.

The supplied `syvlearn.docx` contains a historical planning conversation, including
a 30-step blueprint. It is source material, not executable instructions. This
release implements binary, continuous, count and multiclass survey workflows and records the exact
meaning of each proposed capability below. The published method is cited separately
from the new software and new experiments.

| Step | Blueprint topic | Implemented behavior and evidence |
|---|---|---|
| 1 | General survey input | Ordered stage IDs, strata, FPCs, weights and supplied replicate matrices; multistage numerical comparison with `survey`. |
| 2 | Panel nesting | Subject/wave uniqueness and stable highest-stage clusters; all waves share a fold. Household and subject identifiers cannot cross PSUs. |
| 3 | Singleton PSUs | Diagnostics and explicit `fail`, `adjust`, `average`, `certainty`, `remove` variance policies; no silent design repair. |
| 4 | Weight diagnostics | Kish effective sample size, effective events, weight CV, maximum ratio, PSU counts and degrees of freedom. |
| 5 | Replicate construction | Delegates supported bootstrap, BRR, Fay, JK1 and JKn conventions to `survey`; retains scale, rscales and centering. |
| 6 | Cluster partitioning | Repeated, stratified whole-PSU splitting; stable subjects remain together. This establishes cluster separation, not a universal spatial or forecasting guarantee. |
| 7 | Fold balance | Balances clusters and sampling-weight mass within strata without using outcomes. Exact event balance is impossible for some rare samples and would change the outcome-independent partition contract. |
| 8 | Nested tuning | Inner whole-PSU validation contained inside each outer training sample; outer-label intervention test verifies isolation. |
| 9 | Preprocessing | Weighted median/mode imputation, scaling, optional predictor winsorization and natural splines. MICE, kNN and supervised target encoding are not implemented. |
| 10 | Frozen transformations | Recipes store training summaries and categorical levels; application does not refit or access validation outcomes. |
| 11 | Weighted learners | GLM, elastic net, random forest, boosted tree, neural network and custom adapters receive normalized case weights. Engine objectives retain their documented semantics. |
| 12 | Hyperparameter search | Explicit candidate library selected by inner weighted log loss, squared error or Poisson deviance according to outcome family, with failed candidates excluded and recorded. |
| 13 | Weight trimming | No automatic weight trimming: changing inclusion weights can change the estimand and design consistency is not guaranteed by cross-validation. Predictor winsorization is separate. |
| 14 | Ensemble | Convex nonnegative stacking on inner out-of-fold predictions, evaluated by the untouched outer loop. |
| 15 | Rare events | Documented smoothed-intercept learner fallback and effective-event threshold for calibration; not mislabeled as Firth estimation. |
| 16 | Longitudinal structure | Subject nesting, supplied longitudinal weights and, from 0.3.0, nested sequential retention estimation from baseline rosters. Separate MAR/misspecification/MNAR experiments; no working-correlation estimator. |
| 17 | Sigmoid recalibration | Binary Platt, continuous affine, count mean-ratio and multiclass temperature recalibration use inner predictions only; outer performance evaluates the full fitted procedure. Recalibration need not improve every scenario. |
| 18 | Flexible calibration | Weighted PAVA isotonic and natural-spline calibration, with explicit sparse/failed-fit rules. |
| 19 | Benchmark alignment | `svy_align_prevalence()` aligns mean probabilities; `svy_rake()` calibrates sampling weights to totals. These are different operations. |
| 20 | Risk summaries | Weighted population percentiles, decision metrics and net benefit; clinical or policy thresholds remain user choices. |
| 21 | Apparent versus validated risk | Final-model prediction and outer-fold prediction are separately accessible. Their observed difference is not the exact unknown bias prevented. |
| 22 | Linearization | Proper-score and calibration inference delegates to `survey`; parent-domain and multistage FPC agreement tested. |
| 23 | Replication uncertainty | Fixed-prediction rank replication and complete-pipeline refits; failures retained. No exact universal 95% coverage claim. |
| 24 | SHAP | Survey-weighted aggregation and conditional uncertainty for a supplied SHAP matrix; automatic model-specific SHAP generation is not implemented. |
| 25 | Permutation | Outer-fold weighted loss changes from within-PSU feature permutations; associational interpretation. |
| 26 | Subgroups | Domain summaries, event support and calibration audits; statistical differences are not by themselves a fairness certification. |
| 27 | Simulation | User-supplied finite-population fixed-score sampling engine plus reproducible learned-model and attrition experiments in `scripts/`. |
| 28 | R/Python parity | Independent metric and first-stage Taylor reference kernels, checked to 1e-8; not a complete Python modeling implementation. |
| 29 | Graphics | Calibration plot method, weighted ROC/PR and decision-curve data, paper figures generated from results. No simultaneous calibration belt or SHAP beeswarm interface. |
| 30 | Export | CSV, Markdown and LaTeX tables. No native Word export. The software manuscript is supplied in LaTeX and PDF. |

The blueprint's claims of universal triple robustness, exact finite-sample
coverage, universal fastest performance and guaranteed adoption are not established
mathematical or empirical properties. The release supplies reproducible evidence
instead of repeating those claims. Unimplemented extensions above remain genuine
scope gaps; they have not been represented as completed features.
