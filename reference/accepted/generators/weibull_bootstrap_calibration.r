###############################################################################
# weibull_bootstrap_calibration.r
# Round 24 Path B Test 1 — Weibull AFT Type-I correction via parametric bootstrap
#
# Per SR panel non-negotiable fix (b):
#   T1 power simulations show Weibull AFT Type-I rates of 7.8-13.1% at nominal
#   alpha=0.05 (n=91-98) — known small-sample bias of survreg's Wald inference.
#   The published H22/H3 results CSVs report nominal p-values that may overstate
#   significance.
#
# Procedure (parametric bootstrap):
#   For each of 8 Weibull main-panel specs (4 in H22, 4 in H3):
#     1. Fit a NULL model: same covariates BUT interaction term constrained to 0
#        (i.e., drop the eng_c × moderator_c term)
#     2. Generate B=999 bootstrap samples from the null model's predicted survival
#     3. For each bootstrap sample, refit the FULL model (with interaction)
#     4. Record the test statistic |z| = |beta_int / SE(beta_int)| under H0
#     5. Compute empirical critical value: 95th percentile of |z| under H0
#     6. Compute calibrated p-value for the OBSERVED beta_int as
#        2 * (1 - F_emp(|z_obs|)) where F_emp is the bootstrap CDF
#
# Output: results/weibull_bootstrap_calibration.csv
#   columns: family, spec, moderator, beta_int_obs, z_obs_nominal, p_nominal,
#            z_critical_emp_05, p_calibrated_05, p_calibrated_alpha_holm
#
# Per audit B5 residual (2026-05-07): now also calibrates:
#   H1 M1 Weibull AFT (focal: I(eng_prob_general^2) — curvature)
#   H2.1 M1, M2 Weibull AFT (focal: eng_c:moderator_c — interaction)
# Total: 11 specs (was 8: H22 M1-M4 + H3 M1-M4).
###############################################################################

rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
})

project_root <- Sys.getenv("ACCEPTED_ANALYSIS_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set ACCEPTED_ANALYSIS_ROOT to the original RR cowork directory.")
source_dir   <- file.path(project_root, "RR_datacode", "empirical clean")
results_dir  <- file.path(source_dir, "results")

# Source H22 and H3 into isolated child environments to avoid the rm(list=ls())
# that those scripts execute at top-of-file
h22_env <- new.env()
suppressMessages(invisible(capture.output(
  sys.source(file.path(source_dir, "H22_code.r"), envir = h22_env)
)))
H22_models <- mget(c("model1", "model2", "model3", "model4"), envir = h22_env)
H22_df <- h22_env$df

h3_env <- new.env()
suppressMessages(invisible(capture.output(
  sys.source(file.path(source_dir, "H3_code.r"), envir = h3_env)
)))
H3_models <- mget(c("model1", "model2", "model3", "model4"), envir = h3_env)
H3_df <- h3_env$df

# B5-residual: source H2.1 (model1, model2 = Weibull AFT) and H1 (model1 = Weibull AFT)
h21_env <- new.env()
suppressMessages(invisible(capture.output(
  sys.source(file.path(source_dir, "H21_code.r"), envir = h21_env)
)))
H21_models <- mget(c("model1", "model2"), envir = h21_env)
H21_df <- h21_env$df

h1_env <- new.env()
suppressMessages(invisible(capture.output(
  sys.source(file.path(source_dir, "H1_code.r"), envir = h1_env)
)))
H1_models <- list(model1 = h1_env$model1)  # H1 main-panel concession Weibull AFT (M1)
H1_df <- h1_env$df

cat("\n================================================================\n")
cat("WEIBULL AFT BOOTSTRAP CALIBRATION — Round 24 Path B SR fix (b)\n")
cat("================================================================\n")

read_env_int <- function(name, default) {
  v <- Sys.getenv(name, unset = "")
  if (!nzchar(v)) return(default)
  p <- suppressWarnings(as.integer(v))
  if (is.na(p) || p < 1L) return(default)
  p
}
detect_cores <- function() {
  if (.Platform$OS.type == "windows") return(1L)
  d <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (is.na(d) || d < 2L) return(1L)
  max(1L, as.integer(d) - 1L)
}

B              <- read_env_int("BOOT_REPS", 999L)
PARALLEL_CORES <- min(4L, detect_cores())
SEED           <- 24

cat(sprintf("Bootstrap reps: B=%d  Cores: %d\n", B, PARALLEL_CORES))

###############################################################################
# Spec configuration: 8 Weibull AFT main-panel specs
###############################################################################
# Each spec carries focal_coef + null_drop so the bootstrap routine handles
# both the H1 curvature case (drop squared term) and the H2.1/H2.2/H3
# interaction cases (drop the interaction term) uniformly.
WEIBULL_SPECS <- list(
  H1_M1   = list(family = "H1",   spec = "M1", moderator = NULL,
                 focal_coef = "I(eng_prob_general^2)",
                 null_drop  = "I(eng_prob_general^2)"),
  H21_M1  = list(family = "H2.1", spec = "M1", moderator = "domaut_c",
                 focal_coef = "eng_c:domaut_c",
                 null_drop  = "eng_c:domaut_c"),
  H21_M2  = list(family = "H2.1", spec = "M2", moderator = "stterr_c",
                 focal_coef = "eng_c:stterr_c",
                 null_drop  = "eng_c:stterr_c"),
  H22_M1 = list(family = "H2.2", spec = "M1", moderator = "fdiin_c",
                focal_coef = "eng_c:fdiin_c",   null_drop = "eng_c:fdiin_c"),
  H22_M2 = list(family = "H2.2", spec = "M2", moderator = "trade_c",
                focal_coef = "eng_c:trade_c",   null_drop = "eng_c:trade_c"),
  H22_M3 = list(family = "H2.2", spec = "M3", moderator = "lnaid_c",
                focal_coef = "eng_c:lnaid_c",   null_drop = "eng_c:lnaid_c"),
  H22_M4 = list(family = "H2.2", spec = "M4", moderator = "iodonor_c",
                focal_coef = "eng_c:iodonor_c", null_drop = "eng_c:iodonor_c"),
  H3_M1  = list(family = "H3",   spec = "M1", moderator = "repress_c",
                focal_coef = "eng_c:repress_c", null_drop = "eng_c:repress_c"),
  H3_M2  = list(family = "H3",   spec = "M2", moderator = "civlib_c",
                focal_coef = "eng_c:civlib_c", null_drop = "eng_c:civlib_c"),
  H3_M3  = list(family = "H3",   spec = "M3", moderator = "clpol_c",
                focal_coef = "eng_c:clpol_c",   null_drop = "eng_c:clpol_c"),
  H3_M4  = list(family = "H3",   spec = "M4", moderator = "polity_c",
                focal_coef = "eng_c:polity_c", null_drop = "eng_c:polity_c")
)

# Get fitted model + appropriate df for a spec
get_model_and_data <- function(key) {
  cfg <- WEIBULL_SPECS[[key]]
  spec_idx <- substr(cfg$spec, 2, 2)
  if (cfg$family == "H1") {
    list(model = H1_models[[paste0("model", spec_idx)]],
         df    = H1_df)
  } else if (cfg$family == "H2.1") {
    list(model = H21_models[[paste0("model", spec_idx)]],
         df    = H21_df)
  } else if (cfg$family == "H2.2") {
    list(model = H22_models[[paste0("model", spec_idx)]],
         df    = H22_df)
  } else {  # H3
    list(model = H3_models[[paste0("model", spec_idx)]],
         df    = H3_df)
  }
}

# Build the Weibull AFT formula (full model)
# - H1     : eng_prob_general + I(eng_prob_general^2) + NAVCO 1.3 controls (D-22 spec from H1_code.r M1)
# - H2.1   : eng_c * moderator + NAVCO 1.3 parsimonious controls (D-30; no lnpop/lngdp/colonized_english)
# - H2.2   : eng_c * moderator + NAVCO 1.3 controls + lngdp
# - H3     : eng_c * moderator + NAVCO 1.3 controls + lngdp; D-36 control swap when moderator=repress_c
build_full_formula <- function(family, moderator, response = "event_concession") {
  if (family == "H1") {
    rhs <- "eng_prob_general + I(eng_prob_general^2) + nonviolent_camp + REGCHANGE + v2csreprss + lnnum_image_sum + lnpop + colonized_english"
  } else if (family == "H2.1") {
    rhs <- sprintf("eng_c * %s + nonviolent_camp + REGCHANGE + v2csreprss + lnnum_image_sum",
                   moderator)
  } else if (family == "H2.2") {
    rhs <- sprintf("eng_c * %s + nonviolent_camp + REGCHANGE + v2csreprss + lnnum_image_sum + lngdp",
                   moderator)
  } else {  # H3
    dom_ctrl <- if (moderator == "repress_c") "p_polity2" else "v2csreprss"
    rhs <- sprintf("eng_c * %s + nonviolent_camp + REGCHANGE + %s + lnnum_image_sum + lngdp",
                   moderator, dom_ctrl)
  }
  as.formula(sprintf("Surv(duration, %s) ~ %s", response, rhs))
}

# Build NULL model formula (drop the focal term: squared for H1, interaction otherwise)
build_null_formula <- function(family, moderator, response = "event_concession") {
  if (family == "H1") {
    # Drop only the squared term — keep linear + controls
    rhs <- "eng_prob_general + nonviolent_camp + REGCHANGE + v2csreprss + lnnum_image_sum + lnpop + colonized_english"
  } else if (family == "H2.1") {
    rhs <- sprintf("eng_c + %s + nonviolent_camp + REGCHANGE + v2csreprss + lnnum_image_sum",
                   moderator)
  } else if (family == "H2.2") {
    rhs <- sprintf("eng_c + %s + nonviolent_camp + REGCHANGE + v2csreprss + lnnum_image_sum + lngdp",
                   moderator)
  } else {  # H3
    dom_ctrl <- if (moderator == "repress_c") "p_polity2" else "v2csreprss"
    rhs <- sprintf("eng_c + %s + nonviolent_camp + REGCHANGE + %s + lnnum_image_sum + lngdp",
                   moderator, dom_ctrl)
  }
  as.formula(sprintf("Surv(duration, %s) ~ %s", response, rhs))
}

# Simulate Weibull event times under H0 (no interaction).
# Manually compute X matrix from a formula object + base_df to avoid formula-env
# lookup issues inside parallel::mclapply child processes.
# Takes:
#   null_formula_obj  — formula like ~ eng_c + moderator + ... (no Surv() wrapper)
#   beta_null         — named coefficient vector from fitted null survreg
#   scale_sigma       — null_model$scale
sim_from_null <- function(null_formula_obj, beta_null, scale_sigma,
                          base_df, response = "event_concession",
                          duration_name = "duration") {
  # Identify variables needed and complete-case filter
  needed_vars <- unique(c(all.vars(null_formula_obj), duration_name, response))
  needed_vars <- needed_vars[needed_vars %in% names(base_df)]
  cc <- complete.cases(base_df[, needed_vars, drop = FALSE])
  base_cc <- base_df[cc, , drop = FALSE]
  # Build X matrix manually from formula
  mm <- model.matrix(null_formula_obj, data = base_cc)
  # Align columns with beta_null by name
  beta_aligned <- beta_null[colnames(mm)]
  if (any(is.na(beta_aligned))) {
    return(NULL)
  }
  eta <- as.vector(mm %*% beta_aligned)
  T_sim <- rweibull(length(eta), shape = 1 / scale_sigma, scale = exp(eta))
  obs_dur <- base_cc[[duration_name]]
  out <- base_cc
  out[[duration_name]] <- pmin(T_sim, obs_dur)
  out[[response]]      <- as.integer(T_sim <= obs_dur)
  out
}

# Refit full model on simulated data, return |z| for interaction term
refit_and_extract_z <- function(sim_data, full_formula, focal_coef) {
  res <- tryCatch({
    m <- do.call(survreg,
                 list(formula = full_formula, data = sim_data, dist = "weibull"))
    tbl <- summary(m)$table
    z   <- tbl[focal_coef, "z"]
    abs(z)
  }, error = function(e) NA_real_)
  res
}

###############################################################################
# Main bootstrap loop: per-spec parametric bootstrap of |z| under H0
###############################################################################

calibrate_one_spec <- function(key) {
  cfg <- WEIBULL_SPECS[[key]]
  mod_label <- if (is.null(cfg$moderator)) "curvilinear (eng^2)" else cfg$moderator
  cat(sprintf("\n--- %s [%s × %s] ---\n", key, cfg$family, mod_label))
  md  <- get_model_and_data(key)
  full_f <- build_full_formula(cfg$family, cfg$moderator)
  null_f <- build_null_formula(cfg$family, cfg$moderator)
  focal_coef <- cfg$focal_coef

  # Refit null model for the bootstrap DGP. Use do.call so the formula and
  # data are inlined in the model$call (rather than referenced by name), which
  # avoids "object 'null_f' not found" when summary() is called downstream.
  null_model <- do.call(survreg,
                        list(formula = null_f, data = md$df, dist = "weibull"))
  cat(sprintf("  null model converged. scale=%.3f, n=%d\n",
              null_model$scale, length(null_model$linear.predictors)))

  # Observed |z| from the published full model
  full_tbl <- summary(md$model)$table
  z_obs <- abs(full_tbl[focal_coef, "z"])
  beta_obs <- coef(md$model)[focal_coef]
  p_nominal <- full_tbl[focal_coef, "p"]
  cat(sprintf("  observed: beta_int=%.3f  |z|=%.3f  p_nominal=%.4f\n",
              beta_obs, z_obs, p_nominal))

  # Bootstrap |z| under H0
  cat(sprintf("  running B=%d bootstrap reps on %d cores...\n", B, PARALLEL_CORES))
  set.seed(SEED)
  rep_seeds <- sample.int(.Machine$integer.max, B)

  # Build the RHS-only version of null_f for sim_from_null
  # null_f is `Surv(duration, event_concession) ~ rhs`; extract RHS via deparse(.[[3]])
  null_rhs_f <- as.formula(paste("~", deparse(null_f[[3]], width.cutoff = 500L)))
  null_beta <- coef(null_model)
  null_scale <- null_model$scale

  z_boot <- parallel::mclapply(seq_len(B), function(i) {
    set.seed(rep_seeds[i])
    sim_d <- tryCatch(
      sim_from_null(null_rhs_f, null_beta, null_scale, md$df,
                    response = "event_concession", duration_name = "duration"),
      error = function(e) NULL)
    if (is.null(sim_d)) return(NA_real_)
    refit_and_extract_z(sim_d, full_f, focal_coef)
  }, mc.cores = PARALLEL_CORES)
  z_boot <- as.numeric(unlist(z_boot))
  z_boot <- z_boot[!is.na(z_boot)]
  fit_ok_rate <- length(z_boot) / B

  # Empirical critical values + calibrated p-value
  z_crit_05 <- quantile(z_boot, 0.95)   # 95th percentile of |z| under H0 = critical value at alpha=0.05
  z_crit_10 <- quantile(z_boot, 0.90)   # alpha=0.10
  # Calibrated p-value: proportion of bootstrap |z| at least as large as observed
  p_calib_raw <- mean(z_boot >= z_obs)
  # Calibrated p-value at Holm-corrected threshold (alpha/k where k=8)
  z_crit_holm <- quantile(z_boot, 1 - 0.05 / 8)
  # Empirical alpha at nominal z=1.96 (for the calibration-inflation report)
  emp_alpha_at_nominal <- mean(z_boot >= 1.96)

  cat(sprintf("  z_critical_emp(0.05)=%.3f  vs  z_nominal=1.96  (inflation factor=%.2f)\n",
              z_crit_05, z_crit_05 / 1.96))
  cat(sprintf("  empirical alpha at z=1.96 nominal = %.3f (vs nominal 0.05)\n",
              emp_alpha_at_nominal))
  cat(sprintf("  calibrated p-value for observed |z|=%.3f: p_calib=%.4f (vs nominal %.4f)\n",
              z_obs, p_calib_raw, p_nominal))

  data.frame(
    family = cfg$family, spec = cfg$spec,
    moderator = if (is.null(cfg$moderator)) NA_character_ else cfg$moderator,
    focal_coef = cfg$focal_coef,
    beta_int_obs = beta_obs, z_obs_nominal = z_obs, p_nominal = p_nominal,
    z_critical_emp_05 = z_crit_05, z_critical_emp_10 = z_crit_10,
    z_critical_holm   = z_crit_holm,
    inflation_factor  = z_crit_05 / 1.96,
    empirical_alpha_at_nominal_z196 = emp_alpha_at_nominal,
    p_calibrated_emp  = p_calib_raw,
    n_boot_reps_ok    = length(z_boot),
    fit_ok_rate       = fit_ok_rate,
    stringsAsFactors  = FALSE
  )
}

t0 <- Sys.time()
cat(sprintf("\nStart time: %s\n", t0))

results <- do.call(rbind, lapply(names(WEIBULL_SPECS), calibrate_one_spec))

t1 <- Sys.time()
cat(sprintf("\nTotal wall time: %.1f min\n", as.numeric(difftime(t1, t0, units = "mins"))))

cat("\n================================================================\n")
cat("WEIBULL BOOTSTRAP CALIBRATION — SUMMARY\n")
cat("================================================================\n")
print(results, row.names = FALSE, digits = 3)

out_path <- file.path(results_dir, "weibull_bootstrap_calibration.csv")
write.csv(results, out_path, row.names = FALSE)
cat(sprintf("\nResults written to: %s\n", out_path))
