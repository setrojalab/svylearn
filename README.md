# svylearn: Design-consistent statistical learning, calibration, and inference for complex surveys

**Learning, calibration and inference across cross-sectional and longitudinal complex survey datasets.**

The title is preserved from the Figshare software citation in the published method
paper: DOI `10.6084/m9.figshare.33437179`. NHANES is a motivating application;
the survey interfaces accept the analyst's own data, variable names and design.

Authors: Sunday A. Adetunji (Oregon State University) and Oluwatimileyin A. Adetunji (Bowen University). Paper correspondence: Sunday A. Adetunji <adetunjs@oregonstate.edu>. Package maintainer: Sunday A. Adetunji <sundayadetunjisa@gmail.com>. Coauthor contact: <eniola16.to@gmail.com>. Project home: <https://setrojalabs.github.io/>.

The package implements reusable elements of the method of Adetunji et al. (2026), *Design-consistent statistical learning for rare survey outcomes: a cardiovascular-kidney-metabolic multimorbidity plasmode benchmark in women*, Current Biomedicine 5(1), 7–23, DOI: [10.29244/currbiomed.5.1.7](https://doi.org/10.29244/currbiomed.5.1.7). The new software adds nested calibration, explicit rare-event fallbacks, general survey and panel specifications, and separate software validation experiments.

## Install from GitHub

```r
install.packages("remotes")
remotes::install_github("setrojalabs/svylearn", dependencies = NA)
library(svylearn)
```

Optional engines can be installed when needed, for example
`install.packages(c("survey", "glmnet"))`. Both cross-sectional and longitudinal
workflows use the same package. Report issues at
<https://github.com/setrojalabs/svylearn/issues>.

## Install the local release

From an R session whose working directory is this project:

For a source checkout, first run `R CMD build .` in a terminal. The prepared
release archive already contains the resulting source tarball.

```r
install.packages("svylearn_0.3.1.tar.gz", repos = NULL, type = "source")
install.packages(c("survey", "glmnet"))
library(svylearn)
```

The core attaches without loading optional modeling or survey packages. Nothing downloads, installs, starts a worker, or accesses the network when the package loads. The package has no compiled code and only imports packages shipped with R. CRAN installation by package name becomes available only after CRAN acceptance.

## Frozen evaluation suite

The 12 core estimators and their stable return keys are defined in
[FROZEN_METRICS.md](FROZEN_METRICS.md). `svy_glossary()` returns the 17 reporting
rows and their interpretations. After fitting, run `svy_eval(fit, thresholds = .1)`
to evaluate consensus out-of-fold predictions, or supply a design, outcomes and
independent prediction vector. `print(audit, interpret = TRUE)` includes definitions.

Transport samples, known event totals, conformal test sets, sensitivity targets,
SHAP contributions and importance audits require their own inputs. Missing inputs
are reported explicitly in `audit$availability`. Binary decision metrics are not
assigned to continuous or multiclass outcomes. The supplied `nns_svy` formula is
positive follow-ups per true positive; `nns_population` reports total screens per
true positive when the whole target population is screened.

`svy_conformal()` retains the test-point mass at infinity. Its exchangeable and
covariate-shift interpretations have explicit assumptions; arbitrary design
weights do not establish finite-sample coverage. `svy_weight_bounds()` computes
sharp box-weight mean limits. `svy_sensitivity(target = "tilted_lower_ci")` follows
a mean-minimizing tilt with conditional intervals and does not claim global
optimization of the confidence bound. Full definitions and examples are in the
installed help files.

## Cross-sectional workflow

```r
d <- svy_example(1000)
design <- svy_design(d, weights = "weight", ids = "psu", strata = "stratum")
learners <- list(
  svy_learner("glmnet", list(alpha = 0.5, lambda = 0.005)),
  svy_learner("glmnet", list(alpha = 0.5, lambda = 0.02)))
fit <- svy_fit(design, "y", c("x1", "x2"), learners,
               v = 5, inner_v = 3, repeats = 5, seed = 17)
summary(fit)
predict(fit, d[1:5, ])
plot(fit)
svy_export(fit, "performance.csv")
```

Each outer fold is held out of preprocessing, candidate selection, stacking and calibration. The calibrator is trained on inner out-of-fold scores from the outer training sample. Its base learner is then refitted on that training sample. Rare training folds use a documented, smoothed intercept fallback and report it in `fit$audit`.

`predict(fit, newdata)` uses the final model. Evaluate internal performance using `fit$metrics` or `fit$consensus`; evaluating final-model predictions on its training sample estimates apparent performance.

## Panel workflow

```r
panel <- svy_example(1000, waves = 3)
design <- svy_design(panel, "weight", "psu", "stratum",
                     subject = "subject", wave = "wave")
fit <- svy_fit(design, "y", c("x1", "x2", "wave"),
               learners = learners, v = 5, inner_v = 3)
```

All waves of a subject stay in one PSU and one fold. Supply valid wave-specific or longitudinal weights for your intended estimand. The illustrated target is classification of person-wave outcomes in unseen clusters; it is not a forward-time forecasting design. Cross-PSU movers must be assigned a defensible stable cluster before fitting.

## Inference and diagnostics

```r
fixed <- svy_infer(design, panel$y, fit$consensus)
fixed$scores
fixed$calibration
replicates <- svy_replicates(design, type = "bootstrap", replicates = 500)
svy_infer(design, panel$y, fit$consensus, replicate_design = replicates)$ranks
svy_decision(panel$y, fit$consensus, panel$weight)
svy_curves(panel$y, fit$consensus, panel$weight)
svy_calibration_error(panel$y, fit$consensus, panel$weight)
svy_importance(fit, repeats = 5)
```

`svy_infer()` estimates survey uncertainty **conditional on the supplied predictions**. It does not account for fitting, tuning or calibration uncertainty. `svy_refit_replicates()` reruns the complete nested statistic under supplied replicate weights and retains failed replicates. Replication validity depends on the sampling design and learning procedure; neither interface guarantees exact finite-sample coverage.

Multistage IDs, finite population corrections, supplied replicate weights and raking delegate to `survey`. Supported constructed schemes include bootstrap, BRR, Fay, JK1 and JKn as appropriate for the design. Singleton corrections are explicit and restored after each call. Domains should use a parent design and `domain=`.

## Optional learners

`svy_learner()` supports `glm`, `glmnet`, `ranger`, `xgboost`, `nnet`, and custom weighted fit/predict callbacks. Supply candidate hyperparameters as separate learners. `stack = TRUE` estimates convex ensemble weights from inner out-of-fold predictions. Third-party engines supply their own optimized numerical code.

## Continuous, count and multiclass responses

The same design, recipe, nested validation, stacking, permutation and replicate-refit
interfaces accept four explicitly declared outcome families. Designs can describe
independent records, stratified or multistage samples, clusters and repeated subjects.
Numeric and categorical predictors may be mixed; training-fold recipes handle missing
predictors and previously unseen predictor categories.

| `family` | Outcome | Prediction | Selection loss | Default calibration |
|---|---|---|---|---|
| `binomial` | Numeric 0/1 | Probability vector | Log loss | Platt |
| `gaussian` | Finite numeric | Mean-response vector | Squared error | Weighted affine |
| `poisson` | Nonnegative integer | Nonnegative mean vector | Poisson deviance | Weighted mean ratio |
| `multinomial` | Factor with declared class levels | Named class-probability matrix | Multiclass log loss | Temperature scaling |

```r
study <- svy_example(500)
set.seed(42)
study$amount <- 2 + study$x1 + rnorm(nrow(study))
study$count <- rpois(nrow(study), exp(0.3 + 0.3 * study$x1))
study_design <- svy_design(study, "weight", "psu", "stratum")
amount_fit <- svy_fit(study_design, "amount", c("x1", "x2"),
                      family = "gaussian", v = 3, inner_v = 2)
count_fit <- svy_fit(study_design, "count", c("x1", "x2"),
                     family = "poisson", v = 3, inner_v = 2)
svy_infer(study_design, study$amount, amount_fit$consensus,
          family = "gaussian")$scores

flowers <- datasets::iris
flowers$weight <- 1
flower_design <- svy_design(flowers, "weight")
class_fit <- svy_fit(flower_design, "Species", names(flowers)[1:4],
                     learners = list(svy_learner("glm", list(lambda = 0.01))),
                     family = "multinomial", v = 3, inner_v = 2)
predict(class_fit, flowers[1:3, ])
```

Use `calibration = "none"` to retain raw mean/probability predictions. Multiclass
factor levels specify the class schema before resampling; prediction columns always
retain those levels, including classes absent in a training fold. Sparse class
training uses the recorded smoothed-intercept fallback. Multiclass predictions have
dimensions observations by classes by repeats; scalar predictions are observations
by repeats. `fit$consensus` averages repeats and `summary(fit)` reports each target.

The `glm`, `glmnet`, `ranger` and `xgboost` adapters cover all four families. `nnet`
covers binary, continuous and multiclass responses. Count random forests estimate
the conditional mean by regression; candidate comparison still uses Poisson deviance.
The dependency-free multinomial GLM uses baseline-category softmax with a small
slope ridge penalty (`lambda = 1e-6` by default). Custom multiclass callbacks must
return class-named probability matrices. Missing outcomes are rejected. Ordinal
factors are treated as nominal classes; censoring and survival targets require a
separately specified model and score, and are not silently treated as ordinary counts.

## Reproduce the software paper

Run the scripts from the project root. Source NHANES files and the supplied published-result directory remain separate from the installable package.

For a clean checkout, run `Rscript scripts/setup_replication.R` once, then
`Rscript scripts/replicate.R`. The latter explicitly downloads and verifies
the public NHANES inputs, performs the numerical checks and simulations, and
generates the paper tables. Set the `PYTHON` environment variable if Python
is not available under that name.

```text
Rscript scripts/validate.R
python scripts/audit_sources.py
Rscript scripts/simulate.R
Rscript scripts/simulate_comparative.R 500 8
Rscript scripts/simulate_responses.R 500 8
Rscript scripts/nhanes_workflow.R
Rscript scripts/benchmark_large.R 1000000 500
```

The scripts expect the development installation under `work/library`; create it with `R CMD INSTALL --library=work/library .`. Required analysis dependencies and the exact software environment are listed in `validation/R_package_versions.csv`. The NHANES reconstruction additionally uses `haven`, `dplyr`, `tibble`, `forcats`, `digest`, and `preventr`. Replication scripts never overwrite the supplied original results.

The `paper` directory contains JSS-style LaTeX sources. The `validation` directory separates new computations, reproduction of supplied results, and retained development failure records. The Python directory contains independent reference kernels for score and first-stage Taylor-variance parity; it does not provide the complete R modeling API.

## Citation and interpretation

Run `citation("svylearn")` for both the method and software citations. SHAP summaries and within-PSU permutation summaries describe model behavior; they are not causal effects. Decision-curve output requires a defensible decision context and prospective validation before clinical deployment. A software release does not establish universal design consistency, triple robustness, exact coverage, or guaranteed adoption.
## Sequential retention weighting

For a full scheduled panel roster with baseline survey weights, construct an
unfitted specification and pass it into nested learning:

```r
retention <- svy_attrition("observed", history = c("y", "x"), analysis = "target")
design <- svy_design(panel, "w0", "psu", "stratum", subject = "id", wave = "wave")
fit <- svy_fit(design, "y", "x", family = "gaussian", attrition = retention)
summary(fit)$attrition
```

Keep a row for every scheduled subject-wave, including dropout rows. `observed`
is monotone 0/1 starting at one; `target` is a prespecified eligibility indicator.
History columns are lagged automatically. Every inner and outer training fold
refits its retention models and preprocessing. `svy_attrition_weights()` also
provides a standalone helper with an explicit whole-group training mask.

Use baseline weights for estimated retention adjustment; do not re-adjust agency
longitudinal weights for the same dropout. Identification requires sequential MAR,
positivity and adequate recorded history/model specification. Optional stabilization
preserves within-wave normalized means but changes pooled wave mass. Full-pipeline
replication requires baseline-roster replicate weights and refits retention too;
fixed-score intervals condition on the estimated weights.

