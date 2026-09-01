###############################################################################
# power_h22.r
# H2.2 ECONOMIC VULNERABILITY — POWER DIAGNOSTICS (Round 24 / Path B Test 1)
#
# Per Round 24 panel-locked plan (D1=B, D4=A, D5=A; refined for H2.x interactions):
#   - Single-coefficient Wald test on the interaction term beta(eng_c:moderator_c)
#     (NOT joint Wald — H2.x interactions are not collinear like H1's eng+eng^2)
#   - Magnitude grid: {0, beta_int_hat/4, beta_int_hat/2, beta_int_hat, 2*beta_int_hat}
#   - Reps: 2000 for Weibull AFT (M1-M4) / Firth PML (M5-M8)
#   - Per-spec independent simulation; Holm-adjusted-power proxy = power at alpha/k
#     (Bonferroni worst-case proxy for Holm step-down across the k=8 H2.2 family)
#   - MDE framing per D5: "Minimum Detectable Effect under our design at alpha=0.05"
#
# H2.2 main panel (D-45):
#   M1 / M2 / M3 / M4: Weibull AFT × {FDI, Trade, Aid, IODonor}
#   M5 / M6 / M7 / M8: Firth PML × {FDI, Trade, Aid, IODonor}
# All on concession DV (event_concession), NAVCO 1.3 only. k_family=8.
#
# Outputs:
#   results/power_h22.csv         (or _smoke.csv when env-vars set)
#   results/power_h22_summary.csv
#
# Known limitations from H1 implementation:
#   - Weibull AFT survreg may segfault under simulation at extreme magnitudes
#   - safe_sim wrapper catches per-spec failures and continues to next spec
###############################################################################

rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(logistf)
})

project_root <- Sys.getenv("ACCEPTED_ANALYSIS_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set ACCEPTED_ANALYSIS_ROOT to the original RR cowork directory.")
source_dir   <- file.path(project_root, "RR_datacode", "empirical clean")
results_dir  <- file.path(source_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Source H22_code.r to get fitted models + datasets
suppressMessages(
  invisible(capture.output(
    sys.source(file.path(source_dir, "H22_code.r"), envir = environment())
  ))
)
source(file.path(clean_dir, "power_quarantine_helpers.r"))

cat("\n================================================================\n")
cat("H2.2 POWER DIAGNOSTICS — Round 24 / Path B Test 1\n")
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

N_REPS_FAST    <- read_env_int("POWER_H22_REPS_FAST", 2000L)
SMOKE_MODE     <- nzchar(Sys.getenv("POWER_H22_REPS_FAST", unset = ""))
MAG_GRID       <- c(0, 0.25, 0.5, 1.0, 2.0)
ALPHA_LEVELS   <- c(0.05, 0.10)
SEED           <- 24
K_FAMILY       <- 8   # H2.2 main panel size (per D-45 / H22_code.r:788)
PARALLEL_CORES <- min(4L, detect_cores())  # cap at 4 for stability per H1 lesson

cat(sprintf("Settings: nreps=%d  cores=%d  k_family=%d\n",
            N_REPS_FAST, PARALLEL_CORES, K_FAMILY))
if (SMOKE_MODE) cat("Output mode: smoke-test files\n")
cat(sprintf("Magnitude grid: %s (multipliers of fitted beta_int)\n",
            paste(MAG_GRID, collapse = ",")))
cat(sprintf("Holm-proxy power = power at alpha/k = alpha/%d\n", K_FAMILY))

ctrl_logistf <- logistf::logistf.control(maxit = 1000, maxstep = 0.5)
pl_ctrl      <- logistf::logistpl.control(maxit = 1000)

run_replications <- function(rep_ids, fn) {
  safe_fn <- function(rep_id) tryCatch(fn(rep_id), error = function(e) NA_real_)
  if (PARALLEL_CORES <= 1L) return(lapply(rep_ids, safe_fn))
  parallel::mclapply(rep_ids, safe_fn, mc.cores = PARALLEL_CORES)
}

###############################################################################
# Spec configuration: 8 main-panel specs in H2.2
# Each spec: estimator + moderator + interaction-term name + DGP + refit
###############################################################################
H22_SPECS <- list(
  M1 = list(estimator = "Weibull AFT", moderator = "fdiin_c",   focal = "eng_c:fdiin_c"),
  M2 = list(estimator = "Weibull AFT", moderator = "trade_c",   focal = "eng_c:trade_c"),
  M3 = list(estimator = "Weibull AFT", moderator = "lnaid_c",   focal = "eng_c:lnaid_c"),
  M4 = list(estimator = "Weibull AFT", moderator = "iodonor_c", focal = "eng_c:iodonor_c"),
  M5 = list(estimator = "Firth PML",   moderator = "fdiin_c",   focal = "eng_c:fdiin_c"),
  M6 = list(estimator = "Firth PML",   moderator = "trade_c",   focal = "eng_c:trade_c"),
  M7 = list(estimator = "Firth PML",   moderator = "lnaid_c",   focal = "eng_c:lnaid_c"),
  M8 = list(estimator = "Firth PML",   moderator = "iodonor_c", focal = "eng_c:iodonor_c")
)

# Fetch the fitted model for each spec (built by H22_code.r in the sourced env).
# The H22_code.r naming convention: model<spec_id> per the comment block at L48-74.
get_fitted <- function(spec_id) {
  if (spec_id == "M1") return(model1)
  if (spec_id == "M2") return(model2)
  if (spec_id == "M3") return(model3)
  if (spec_id == "M4") return(model4)
  if (spec_id == "M5") return(model5)
  if (spec_id == "M6") return(model6)
  if (spec_id == "M7") return(model7)
  if (spec_id == "M8") return(model8)
  stop(sprintf("Unknown spec %s", spec_id))
}

###############################################################################
# DGP simulators
###############################################################################
complete_model_data <- function(base_df, formula_obj, extra_vars = character()) {
  vars <- unique(c(all.vars(formula_obj), extra_vars))
  vars <- vars[vars %in% names(base_df)]
  base_df[complete.cases(base_df[, vars, drop = FALSE]), , drop = FALSE]
}

# Firth PML: simulate Bernoulli outcome from logit(eta)
sim_firth <- function(base_df, formula_obj, beta_sim, response_name) {
  sim_df <- complete_model_data(base_df, formula_obj)
  mm <- model.matrix(formula_obj, sim_df)
  beta_aligned <- beta_sim[colnames(mm)]
  if (any(is.na(beta_aligned))) {
    stop(sprintf("Coefficient name mismatch: %s",
                 paste(setdiff(colnames(mm), names(beta_sim)), collapse = ", ")))
  }
  eta <- as.vector(mm %*% beta_aligned)
  p   <- plogis(eta)
  out <- sim_df
  out[[response_name]] <- as.integer(rbinom(length(p), size = 1, prob = p))
  out
}

# Weibull AFT: simulate event times from Weibull, censor at observed duration
sim_weibull <- function(base_df, formula_obj, beta_sim, scale_sigma,
                        duration_name = "duration", event_name = "event_concession") {
  sim_df <- complete_model_data(base_df, formula_obj, c(duration_name, event_name))
  mm <- model.matrix(formula_obj, sim_df)
  beta_aligned <- beta_sim[colnames(mm)]
  if (any(is.na(beta_aligned))) {
    stop(sprintf("Coefficient name mismatch: %s",
                 paste(setdiff(colnames(mm), names(beta_sim)), collapse = ", ")))
  }
  eta <- as.vector(mm %*% beta_aligned)
  T_sim <- rweibull(length(eta), shape = 1 / scale_sigma, scale = exp(eta))
  obs_dur <- sim_df[[duration_name]]
  out <- sim_df
  out[[duration_name]] <- pmin(T_sim, obs_dur)
  out[[event_name]]    <- as.integer(T_sim <= obs_dur)
  out
}

###############################################################################
# Per-spec refit + p-value extraction
###############################################################################
refit_spec <- function(spec_id, sim_data) {
  cfg <- H22_SPECS[[spec_id]]
  res <- tryCatch({
    if (cfg$estimator == "Weibull AFT") {
      f <- as.formula(sprintf(
        "Surv(duration, event_concession) ~ eng_c * %s + nonviolent_camp + REGCHANGE + v2csreprss + lnnum_image_sum + lngdp",
        cfg$moderator))
      m <- survreg(f, data = sim_data, dist = "weibull")
      tbl <- summary(m)$table
      list(p_focal = unname(tbl[cfg$focal, "p"]),
           beta_focal = unname(coef(m)[cfg$focal]),
           fit_ok = TRUE)
    } else if (cfg$estimator == "Firth PML") {
      f <- as.formula(sprintf(
        "concession_bin ~ eng_c * %s + nonviolent_camp + REGCHANGE + v2csreprss + lnlengthofcam + lnnum_image_sum + lngdp",
        cfg$moderator))
      m <- logistf(f, data = sim_data, control = ctrl_logistf, plcontrol = pl_ctrl, pl = TRUE)
      list(p_focal = unname(m$prob[cfg$focal]),
           beta_focal = unname(coef(m)[cfg$focal]),
           fit_ok = TRUE)
    } else stop(sprintf("Unknown estimator for %s", spec_id))
  }, error = function(e) list(p_focal = NA_real_, beta_focal = NA_real_, fit_ok = FALSE))
  res
}

###############################################################################
# Per-spec power simulation
###############################################################################
simulate_power_for_spec <- function(spec_id, n_reps) {
  cfg <- H22_SPECS[[spec_id]]
  cat(sprintf("\n--- Simulating %s [%s × %s] (n_reps=%d) ---\n",
              spec_id, cfg$estimator, cfg$moderator, n_reps))

  fitted   <- get_fitted(spec_id)
  beta_hat <- coef(fitted)
  if (!(cfg$focal %in% names(beta_hat))) {
    cat(sprintf("  *** focal '%s' not found in coefficients; skipping ***\n", cfg$focal))
    return(NULL)
  }
  beta_focal_hat <- unname(beta_hat[cfg$focal])
  cat(sprintf("  beta_int (fitted) = %.4f  (focal = %s)\n", beta_focal_hat, cfg$focal))

  # Build simulation function for this spec
  if (cfg$estimator == "Weibull AFT") {
    f_sim <- as.formula(sprintf(
      "~ eng_c * %s + nonviolent_camp + REGCHANGE + v2csreprss + lnnum_image_sum + lngdp",
      cfg$moderator))
    sim_fn <- function(beta_sim) {
      sim_weibull(df, formula_obj = f_sim, beta_sim = beta_sim,
                  scale_sigma = fitted$scale,
                  duration_name = "duration", event_name = "event_concession")
    }
  } else {  # Firth PML
    f_sim <- as.formula(sprintf(
      "~ eng_c * %s + nonviolent_camp + REGCHANGE + v2csreprss + lnlengthofcam + lnnum_image_sum + lngdp",
      cfg$moderator))
    sim_fn <- function(beta_sim) {
      sim_firth(df, formula_obj = f_sim, beta_sim = beta_sim, response_name = "concession_bin")
    }
  }

  out <- data.frame()
  for (mag in MAG_GRID) {
    beta_sim <- beta_hat
    beta_sim[cfg$focal] <- beta_focal_hat * mag

    cat(sprintf("  magnitude = %.2f x fitted (beta_int=%.3f)  ... ",
                mag, beta_sim[cfg$focal]))

    p_values <- run_replications(seq_len(n_reps), function(rep_id) {
      set.seed(SEED + rep_id * 1000 + which(MAG_GRID == mag))
      sim_data <- sim_fn(beta_sim)
      refit_spec(spec_id, sim_data)$p_focal
    })
    p_values <- as.numeric(unlist(p_values))
    fit_ok_rate <- mean(!is.na(p_values))

    for (alpha in ALPHA_LEVELS) {
      valid_p <- p_values[!is.na(p_values)]
      # PE3 R25 convention: _overall = defensive (NA-fits as non-rejections);
      # _conditional = diagnostic. Backward-compatible aliases default to overall.
      raw_cond  <- if (length(valid_p) > 0) mean(valid_p < alpha) else NA_real_
      holm_cond <- if (length(valid_p) > 0) mean(valid_p < alpha / K_FAMILY) else NA_real_
      raw_over  <- sum(p_values < alpha,            na.rm = TRUE) / n_reps
      holm_over <- sum(p_values < alpha / K_FAMILY, na.rm = TRUE) / n_reps
      out <- rbind(out, data.frame(
        spec = spec_id, estimator = cfg$estimator, moderator = cfg$moderator,
        focal_coef = cfg$focal, beta_focal_hat = beta_focal_hat,
        magnitude = mag, beta_focal_sim = beta_sim[cfg$focal],
        n_reps = n_reps, fit_ok_rate = fit_ok_rate,
        alpha = alpha,
        raw_power_overall            = raw_over,
        holm_proxy_power_overall     = holm_over,
        raw_power_conditional        = raw_cond,
        holm_proxy_power_conditional = holm_cond,
        raw_power = raw_over, holm_proxy_power = holm_over,
        stringsAsFactors = FALSE
      ))
    }
    cat(sprintf("raw_power(0.05)=%.3f  holm_proxy(0.05)=%.3f  fit_ok=%.3f\n",
                out$raw_power[out$alpha == 0.05 & out$magnitude == mag][1],
                out$holm_proxy_power[out$alpha == 0.05 & out$magnitude == mag][1],
                fit_ok_rate))
  }
  out
}

###############################################################################
# Main loop with safe_sim wrapper (catches per-spec crashes per H1 lesson)
###############################################################################
safe_sim <- function(spec_id, n_reps) {
  tryCatch(simulate_power_for_spec(spec_id, n_reps),
           error = function(e) {
             cat(sprintf("\n*** %s SIMULATION FAILED: %s ***\n", spec_id, conditionMessage(e)))
             NULL
           })
}

t0 <- Sys.time()
results <- do.call(rbind, list(
  safe_sim("M5", N_REPS_FAST),  # Firth specs first (most likely to succeed)
  safe_sim("M6", N_REPS_FAST),
  safe_sim("M7", N_REPS_FAST),
  safe_sim("M8", N_REPS_FAST),
  safe_sim("M1", N_REPS_FAST),  # Weibull (segfault risk per H1)
  safe_sim("M2", N_REPS_FAST),
  safe_sim("M3", N_REPS_FAST),
  safe_sim("M4", N_REPS_FAST)
))
t1 <- Sys.time()
cat(sprintf("\nTotal sim wall time: %.1f min\n",
            as.numeric(difftime(t1, t0, units = "mins"))))

# MDE summary
compute_mde <- function(power_df, alpha_val, power_target = 0.80) {
  sub <- subset(power_df, alpha == alpha_val & magnitude > 0)
  if (nrow(sub) == 0) return(c(raw_multiplier = Inf, holm_multiplier = Inf))
  rc <- sub$magnitude[!is.na(sub$raw_power) & sub$raw_power >= power_target]
  hc <- sub$magnitude[!is.na(sub$holm_proxy_power) & sub$holm_proxy_power >= power_target]
  c(raw_multiplier = ifelse(length(rc) > 0, min(rc), Inf),
    holm_multiplier = ifelse(length(hc) > 0, min(hc), Inf))
}

if (!is.null(results) && nrow(results) > 0) {
  summary_df <- do.call(rbind, lapply(unique(results$spec), function(sp) {
    spec_results <- subset(results, spec == sp)
    sub05 <- subset(spec_results, alpha == 0.05)
    mde05 <- compute_mde(spec_results, 0.05)
    data.frame(
      spec = sp,
      estimator = sub05$estimator[1],
      moderator = sub05$moderator[1],
      beta_focal_hat = sub05$beta_focal_hat[1],
      type1_at_null_05       = sub05$raw_power[sub05$magnitude == 0.0],
      raw_power_at_half_05   = sub05$raw_power[sub05$magnitude == 0.5],
      raw_power_at_fitted_05 = sub05$raw_power[sub05$magnitude == 1.0],
      raw_power_at_double_05 = sub05$raw_power[sub05$magnitude == 2.0],
      holm_power_at_fitted_05 = sub05$holm_proxy_power[sub05$magnitude == 1.0],
      raw_mde_mult_05        = mde05["raw_multiplier"],
      holm_mde_mult_05       = mde05["holm_multiplier"],
      stringsAsFactors = FALSE
    )
  }))

  cat("\n================================================================\n")
  cat("H2.2 POWER — PER-SPEC MDE SUMMARY (alpha = 0.05)\n")
  cat("================================================================\n")
  print(summary_df, row.names = FALSE)

  prefix <- if (SMOKE_MODE) "power_h22_smoke" else "power_h22"
  out_csv     <- file.path(results_dir, sprintf("%s.csv", prefix))
  summary_csv <- file.path(results_dir, sprintf("%s_summary.csv", prefix))
  detail_diag <- annotate_weibull_power_diagnostics(
    results, "h22", basename(out_csv), results_dir
  )
  results <- detail_diag$data
  summary_diag <- annotate_weibull_power_diagnostics(
    summary_df, "h22", basename(summary_csv), results_dir
  )
  summary_df <- summary_diag$data
  write_weibull_diagnostic_audit(
    list(detail_diag$audit, summary_diag$audit), results_dir
  )
  write.csv(results,    out_csv,     row.names = FALSE)
  write.csv(summary_df, summary_csv, row.names = FALSE)
  cat(sprintf("\nResults  written to: %s\n", out_csv))
  cat(sprintf("Summary  written to: %s\n", summary_csv))
} else {
  cat("\n*** NO SUCCESSFUL SIMULATIONS — all specs crashed ***\n")
}
