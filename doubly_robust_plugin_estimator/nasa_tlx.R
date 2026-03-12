################################################################################
# Doubly Robust Plug-in Estimator: NASA-TLX Subscales -> DCS
#
# Stabilized IPW + G-formula (linear outcome model)
#   Treatment model : ordinal logistic regression (polr)
#   Outcome model   : IPW-weighted linear regression (lm)
#   Inference       : nonparametric bootstrap (1,000 iterations)
#
# Exposure : NASA-TLX subscales (6 items; each scored 0-10)
# Outcome  : DCS Overall (range 0-100)
#
# Assumes `analysis_df` exists with columns named per README.md.
################################################################################

library(tidyverse)
library(MASS)
library(ggplot2)

################################################################################
# Covariates and Exposures
################################################################################

# Minimal sufficient adjustment set for NASA-TLX models
# (derived from DAG via backdoor criterion; see Appendix B, Table B1)
covariates <- c(
  "age",
  "gender_male",

  # Education Level (ref: High School)
  "middle_school",
  "vocational_school",
  "university",
  "graduate_school",
  "education_other",
  "education_no_response",

  # Occupation (ref: Self-Employed / Student / Other)
  "company_employee",
  "government_employee",
  "part_time_worker",
  "unemployed",

  # Time Since Diagnosis (ref: 4-5 years)
  "diagnosis_lt1y",
  "diagnosis_1to2y",
  "diagnosis_2to3y",
  "diagnosis_3to4y",

  # Cancer Type (not mutually exclusive)
  "colorectal",
  "lung",
  "breast",
  "prostate",

  # Comorbidities
  "brain_disorder",
  "mental_disorders",
  "n_comorbidities",

  # Decision-Making Style (factor: shared / patient / physician / pattern)
  "dm_style_actual",

  # Decision Content (not mutually exclusive)
  "whether_to_undergo_treatment",
  "choosing_among_options",
  "additional_treatment_or_tests",
  "determining_timing",
  "changing_plan",
  "other_decision_content",

  # Number of Decision Aids Used
  "n_decision_aids",

  # Treatment Selection (not mutually exclusive)
  "surgery",
  "pharmacotherapy",
  "radiation_therapy",
  "palliative_care",
  "clinical_trial",
  "alternative_therapy",
  "no_treatment",
  "other_treatment",

  # HVICS Subscales (confounders in TLX models)
  "horizontal_individualism",
  "vertical_individualism",
  "horizontal_collectivism",
  "vertical_collectivism"
)

exposures <- c(
  "mental_demand",
  "physical_demand",
  "temporal_demand",
  "own_performance",
  "effort",
  "frustration"
)

for (e in exposures) {
  analysis_df[[e]] <- factor(
    analysis_df[[e]],
    levels  = sort(unique(analysis_df[[e]])),
    ordered = TRUE
  )
}

################################################################################
# Stabilized Inverse Probability Weighting
#
# GPS estimated via ordinal logistic regression.
# Polynomial terms in the treatment model:
#   age^3, vertical_individualism^3, vertical_collectivism^3
# Weight: SW_i = P(T=t) / P(T=t|X)
# Truncation: 99.5th percentile, hard cap at 30.
################################################################################

calculate_weights <- function(df, exposure, covariates,
                              truncate_pctl = 0.995,
                              max_weight    = 30) {

  poly_terms <- c(
    "I(age^2)", "I(age^3)",
    "I(vertical_individualism^2)", "I(vertical_individualism^3)",
    "I(vertical_collectivism^2)",  "I(vertical_collectivism^3)"
  )

  fml       <- as.formula(paste(
    exposure, "~", paste(c(covariates, poly_terms), collapse = " + ")
  ))
  gps_model <- polr(fml, data = df, Hess = TRUE)
  gps_probs <- predict(gps_model, type = "probs")

  levs      <- levels(df[[exposure]])
  marginals <- prop.table(table(df[[exposure]]))

  sw <- vapply(seq_len(nrow(df)), function(i) {
    t_i <- as.character(df[[exposure]][i])
    marginals[t_i] / gps_probs[i, which(levs == t_i)]
  }, numeric(1))

  sw <- pmin(sw, quantile(sw, truncate_pctl, na.rm = TRUE))
  sw <- pmin(sw, max_weight)
  sw
}

################################################################################
# Bootstrap Doubly Robust Estimation
#
# Each bootstrap iteration:
#   1. Resample with replacement
#   2. Estimate stabilized IPW (polr)
#   3. Fit IPW-weighted linear outcome model
#   4. Compute G-formula counterfactual predictions
#
# Polynomial terms in the outcome model (from specification checks):
#   horizontal_individualism^3, vertical_individualism^2,
#   vertical_collectivism^2, n_decision_aids^2
################################################################################

bootstrap_dr <- function(df, exposures, covariates, n_iter = 1000) {

  outcome_poly <- list(
    horizontal_individualism = 3,
    vertical_individualism   = 2,
    vertical_collectivism    = 2,
    n_decision_aids          = 2
  )

  results <- list()

  for (exposure in exposures) {
    cat("Exposure:", exposure, "\n")
    store  <- vector("list", n_iter)
    valid  <- 0
    attempt <- 0

    while (valid < n_iter && attempt < n_iter * 10) {
      attempt <- attempt + 1
      boot_df <- df[sample(nrow(df), replace = TRUE), ]

      if (any(sapply(Filter(is.factor, boot_df), nlevels) < 2)) next

      w <- tryCatch(
        calculate_weights(boot_df, exposure, covariates),
        error = function(e) NULL
      )
      if (is.null(w)) next

      boot_df$t_num <- as.numeric(as.character(boot_df[[exposure]]))

      cov_terms <- vapply(covariates, function(v) {
        if (v %in% names(outcome_poly)) {
          paste0("I(", v, "^", seq_len(outcome_poly[[v]]), ")",
                 collapse = " + ")
        } else {
          v
        }
      }, character(1))

      fml <- as.formula(paste(
        "dcs_overall ~", paste(c("t_num", cov_terms), collapse = " + ")
      ))

      fit <- tryCatch(
        lm(fml, data = boot_df, weights = w),
        error = function(e) NULL
      )
      if (is.null(fit)) next

      t_range <- as.numeric(levels(boot_df[[exposure]]))
      t_seq   <- seq(min(t_range), max(t_range), length.out = 101)
      exp_df  <- boot_df[rep(seq_len(nrow(boot_df)), each = length(t_seq)), ]
      exp_df$t_num <- rep(t_seq, nrow(boot_df))

      valid <- valid + 1
      store[[valid]] <- data.frame(
        treatment  = t_seq,
        prediction = tapply(predict(fit, newdata = exp_df), exp_df$t_num, mean)
      )
      if (valid %% 100 == 0) cat("  ", valid, "iterations\n")
    }

    results[[exposure]] <- store[seq_len(valid)]
  }
  results
}

################################################################################
# Run
################################################################################

set.seed(42)
boot_results <- bootstrap_dr(analysis_df, exposures, covariates, n_iter = 1000)

################################################################################
# Estimated Effects (25th-75th percentile contrast)
################################################################################

effect_table <- map_dfr(exposures, function(exposure) {
  t_num <- as.numeric(analysis_df[[exposure]])
  q25   <- quantile(t_num, 0.25, na.rm = TRUE)
  q75   <- quantile(t_num, 0.75, na.rm = TRUE)

  diffs <- sapply(boot_results[[exposure]], function(p) {
    p$prediction[which.min(abs(p$treatment - q75))] -
      p$prediction[which.min(abs(p$treatment - q25))]
  })
  diffs <- diffs[!is.na(diffs)]

  tibble(
    Exposure = exposure,
    ES       = round(mean(diffs), 3),
    CI_Lower = round(quantile(diffs, 0.025), 3),
    CI_Upper = round(quantile(diffs, 0.975), 3),
    P_Value  = round(2 * pnorm(-abs(mean(diffs) / sd(diffs))), 3)
  )
})

print(effect_table)