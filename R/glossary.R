#' Frozen evaluation names, keys and interpretations
#' @return A data frame with the 17 reporting rows covering 12 core estimators
#'   and foundational summaries. Definitions distinguish empirical diagnostics,
#'   design-based uncertainty and assumptions needed for stronger guarantees.
#' @export
svy_glossary <- function() {
  data.frame(
    pillar = c(rep("Scoring", 3), rep("Survey geometry", 3), rep("Calibration and transport", 4),
      rep("Uncertainty", 3), rep("Decision utility", 2), rep("Equity and explanation", 2)),
    core = c(NA, NA, 1, 2, 3, NA, NA, 4, 5, 6, NA, 7, 8, 9, 10, 11, 12),
    name = c("Design-Weighted Loss", "Scaled Brier Score (IPA)", "Design-Induced Optimism Contrast",
      "Loss Design Effect", "Predictive Effective Sample Size", "Singleton-PSU Variance Audit",
      "Calibration-in-the-Large and Slope", "Design-Weighted Integrated Calibration Index and E90",
      "Post-Raking Marginal Residual", "Target-Population Transportability Overlap",
      "Design-Based Confidence Intervals", "Design-Weighted Conformal Coverage and Set Width",
      "Informative Sampling Sensitivity Bound", "Design-Weighted Net Benefit",
      "Projected Population Headcounts and Screening Ratios", "Subgroup Calibration Disparity",
      "Design-Weighted SHAP and Within-PSU Importance"),
    key = c("loss_svy", "ipa_svy", "optimism_bias", "deff_loss", "n_eff_pred", "singleton_audit",
      "cal_intercept, cal_slope", "ici_svy, e90_svy", "rake_residual", "transport_smd, transport_div",
      "ci_lower, ci_upper", "conf_coverage, conf_width", "sensitivity_gamma", "net_benefit_svy",
      "pop_tp, pop_fp, pop_tn, pop_fn, nns_svy, nns_population", "domain_cal_gap", "shap_svy, vip_psu"),
    interpretation = c(
      "Weighted outcome-appropriate evaluation loss; lower is better.",
      "One minus Brier loss divided by prevalence*(1-prevalence); undefined for a single class.",
      "IID unweighted OOF loss minus design-weighted OOF loss; an observed contrast, not known true bias.",
      "Conditional design variance divided by unweighted sample loss variance/n (SRS with replacement).",
      "Nominal n divided by the loss design effect; an estimated variance equivalence, possibly above n.",
      "Detected single-PSU strata and the explicitly chosen variance policy.",
      "Offset intercept and free calibration slope; the joint-model intercept is reported separately.",
      "Weighted mean and 90th percentile of absolute errors from a fitted calibration spline.",
      "Relative predicted-event-total discrepancy against a separately known positive event total.",
      "Before/after weighted covariate SMD and target-row prediction divergence; not proof of transportability.",
      "Estimated conditional Taylor or replicate intervals; pointwise unless explicitly labeled otherwise.",
      "Empirical weighted test-set coverage and width/cardinality; construction assumptions are separate.",
      "Tipping under a stated bounded-weight model; sharp point bounds are distinct from a tilted CI path.",
      "True-positive proportion minus threshold odds times false-positive proportion, also per 100,000.",
      "Population-scale confusion totals; nns_svy=(TP+FP)/TP, nns_population=population/TP.",
      "Maximum absolute joint intercept and slope departure from one across supported domains.",
      "Weighted supplied absolute SHAP values and within-PSU OOF permutation loss increases; associational."),
    required_inputs = c("outcomes, predictions, design", "binary outcomes and probabilities",
      "separate IID out-of-fold predictions", "loss vector and valid variance design", "loss design effect",
      "PSUs, strata and singleton policy", "nonconstant outcomes and predictions", "binary outcomes and probabilities",
      "external event total and optional raked design", "source and target covariates and weights",
      "estimable design and metric-specific inference", "independent test outcomes and prediction sets",
      "row-level utility, sensitivity target and gamma search range", "binary outcomes, probabilities and thresholds",
      "population-scale weights or known population size", "domain labels and adequate support",
      "external SHAP matrix and/or fitted OOF pipeline"), stringsAsFactors = FALSE)
}
