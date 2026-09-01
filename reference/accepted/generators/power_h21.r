###############################################################################
# power_h21.r
# H2.1 POLITICAL VULNERABILITY — POWER DIAGNOSTICS (R24 / Path B Test 1)
#
# Per panel-locked plan (D1=B, D4=A, D5=A; H2.x = single-coef Wald on beta_int):
#   - Magnitude grid: {0, 0.25, 0.5, 1.0, 2.0} multipliers of fitted beta_int
#   - Reps: 2000 fast (Weibull/Firth); 500 Cox (env-overridable for smoke)
#   - Holm-proxy power = power at alpha/k where k = 8 (H2.1 main per D-63)
#
# H2.1 main panel (D-63) — 8 specs:
#   M1 / M2: Weibull AFT × {domaut_c, stterr_c} on NAVCO 1.3
#   M3 / M4: Firth PML   × {domaut_c, stterr_c} on NAVCO 1.3
#   M7 / M8: Cox PH counting-process × {domaut_c, stterr_c} on NAVCO 2.1
#            (first-event rows keyed by navco21_id after the duplicate id-year
#             merge; SEs clustered on navco21_id)
#   M9 / M10: Firth PML yearly concession × {domaut_c, stterr_c} on NAVCO 2.1
#
# Known limitations from H1:
#   - Weibull AFT survreg may segfault at extreme magnitudes
#   - Cox PH counting-process baseline simulator is unreliable at small N=59
#     (degenerate event/censoring patterns at low magnitudes); Cox specs may
#     fall back to "did not converge" with descriptive EPV from H21_results.csv
###############################################################################

rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(logistf)
  library(simsurv)  # PE2 R25 fix: M7/M8 Cox now use simsurv exponential baseline
})

project_root <- Sys.getenv("ACCEPTED_ANALYSIS_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set ACCEPTED_ANALYSIS_ROOT to the original RR cowork directory.")
source_dir   <- file.path(project_root, "RR_datacode", "empirical clean")
results_dir  <- file.path(source_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

suppressMessages(
  invisible(capture.output(
    sys.source(file.path(source_dir, "H21_code.r"), envir = environment())
  ))
)

cat("\n================================================================\n")
cat("H2.1 POWER DIAGNOSTICS — Round 24 / Path B Test 1\n")
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

N_REPS_FAST    <- read_env_int("POWER_H21_REPS_FAST", 2000L)
N_REPS_COX     <- read_env_int("POWER_H21_REPS_COX",  500L)
SMOKE_MODE     <- nzchar(Sys.getenv("POWER_H21_REPS_FAST", unset = "")) ||
                  nzchar(Sys.getenv("POWER_H21_REPS_COX",  unset = ""))
MAG_GRID       <- c(0, 0.25, 0.5, 1.0, 2.0)
ALPHA_LEVELS   <- c(0.05, 0.10)
SEED           <- 24
K_FAMILY       <- 8   # H2.1 main panel size per D-63
PARALLEL_CORES <- min(4L, detect_cores())

cat(sprintf("Settings: nreps_fast=%d  nreps_cox=%d  cores=%d  k_family=%d\n",
            N_REPS_FAST, N_REPS_COX, PARALLEL_CORES, K_FAMILY))
if (SMOKE_MODE) cat("Output mode: smoke-test files\n")

ctrl_logistf <- logistf::logistf.control(maxit = 1000, maxstep = 0.5)
pl_ctrl      <- logistf::logistpl.control(maxit = 1000)

run_replications <- function(rep_ids, fn) {
  safe_fn <- function(rep_id) tryCatch(fn(rep_id), error = function(e) NA_real_)
  if (PARALLEL_CORES <= 1L) return(lapply(rep_ids, safe_fn))
  parallel::mclapply(rep_ids, safe_fn, mc.cores = PARALLEL_CORES)
}

H21_SPECS <- list(
  M1  = list(estimator = "Weibull AFT", moderator = "domaut_c", focal = "eng_c:domaut_c"),
  M2  = list(estimator = "Weibull AFT", moderator = "stterr_c", focal = "eng_c:stterr_c"),
  M3  = list(estimator = "Firth PML",   moderator = "domaut_c", focal = "eng_c:domaut_c",
             data_name = "df", response = "concession_bin"),
  M4  = list(estimator = "Firth PML",   moderator = "stterr_c", focal = "eng_c:stterr_c",
             data_name = "df", response = "concession_bin"),
  M7  = list(estimator = "Cox PH",      moderator = "domaut_c", focal = "eng_c:domaut_c"),
  M8  = list(estimator = "Cox PH",      moderator = "stterr_c", focal = "eng_c:stterr_c"),
  M9  = list(estimator = "Firth PML",   moderator = "domaut_c", focal = "eng_c:domaut_c",
             data_name = "d21", response = "concession_bin"),
  M10 = list(estimator = "Firth PML",   moderator = "stterr_c", focal = "eng_c:stterr_c",
             data_name = "d21", response = "concession_bin")
)

get_fitted <- function(spec_id) {
  if (spec_id == "M1") return(model1)
  if (spec_id == "M2") return(model2)
  if (spec_id == "M3") return(model3)
  if (spec_id == "M4") return(model4)
  if (spec_id == "M7") return(model7)
  if (spec_id == "M8") return(model8)
  if (spec_id == "M9") return(model9)
  if (spec_id == "M10") return(model10)
  stop(sprintf("Unknown spec %s", spec_id))
}

complete_model_data <- function(base_df, formula_obj, extra_vars = character()) {
  vars <- unique(c(all.vars(formula_obj), extra_vars))
  vars <- vars[vars %in% names(base_df)]
  base_df[complete.cases(base_df[, vars, drop = FALSE]), , drop = FALSE]
}

sim_firth <- function(base_df, formula_obj, beta_sim, response_name) {
  sim_df <- complete_model_data(base_df, formula_obj)
  mm <- model.matrix(formula_obj, sim_df)
  beta_aligned <- beta_sim[colnames(mm)]
  if (any(is.na(beta_aligned)))
    stop(sprintf("Coef mismatch: %s",
                 paste(setdiff(colnames(mm), names(beta_sim)), collapse = ", ")))
  eta <- as.vector(mm %*% beta_aligned)
  out <- sim_df
  out[[response_name]] <- as.integer(rbinom(length(eta), size = 1, prob = plogis(eta)))
  out
}

sim_weibull <- function(base_df, formula_obj, beta_sim, scale_sigma,
                        duration_name = "duration", event_name = "event_concession") {
  sim_df <- complete_model_data(base_df, formula_obj, c(duration_name, event_name))
  mm <- model.matrix(formula_obj, sim_df)
  beta_aligned <- beta_sim[colnames(mm)]
  if (any(is.na(beta_aligned)))
    stop(sprintf("Coef mismatch: %s",
                 paste(setdiff(colnames(mm), names(beta_sim)), collapse = ", ")))
  eta <- as.vector(mm %*% beta_aligned)
  T_sim <- rweibull(length(eta), shape = 1 / scale_sigma, scale = exp(eta))
  obs_dur <- sim_df[[duration_name]]
  out <- sim_df
  out[[duration_name]] <- pmin(T_sim, obs_dur)
  out[[event_name]]    <- as.integer(T_sim <= obs_dur)
  out
}

# --- Cox PH counting-process simulation: simsurv exponential baseline ---
# PE2 R25 fix: replace AFT-proxy approach (which crashed with
# `'names' attribute [9] must be the same length as the vector [2]` due to
# misalignment between coxph coef vector and survreg AFT coef vector at the
# `(Intercept)` term) with simsurv exponential baseline calibrated to the
# observed event rate. Mirrors power_cox_specs.r for H1 M4 and the H1 M4
# replacement in power_h1.r.
build_cox_simsurv <- function(spec_id) {
  cfg <- H21_SPECS[[spec_id]]
  cox_vars <- c("tstart", "tstop", "event_concession", "eng_c", cfg$moderator,
                "nonviolent", "REGCHANGE", "v2csreprss", "lnnum_image_sum",
                "navco21_id")
  cc <- complete.cases(d21_first_conc[, cox_vars, drop = FALSE])
  base_cc <- d21_first_conc[cc, , drop = FALSE]
  rhs <- as.formula(sprintf("~ eng_c * %s + nonviolent + REGCHANGE + v2csreprss + lnnum_image_sum",
                            cfg$moderator))
  X_mm <- model.matrix(rhs, data = base_cc)
  X_use <- X_mm[, setdiff(colnames(X_mm), "(Intercept)"), drop = FALSE]
  obs_event_rate <- mean(base_cc$event_concession, na.rm = TRUE)
  max_followup <- max(base_cc$tstop, na.rm = TRUE)
  lambda_base <- -log(max(1 - obs_event_rate, 0.01)) / max_followup
  list(base_df = base_cc, X_use = X_use,
       lambda = lambda_base, gamma = 1.0, maxt = max_followup)
}

build_formula <- function(spec_id, side = c("refit", "sim")) {
  side <- match.arg(side)
  cfg <- H21_SPECS[[spec_id]]
  if (cfg$estimator == "Weibull AFT") {
    rhs <- sprintf("eng_c * %s + nonviolent_camp + REGCHANGE + v2csreprss + lnnum_image_sum",
                   cfg$moderator)
    if (side == "refit") return(as.formula(sprintf("Surv(duration, event_concession) ~ %s", rhs)))
    return(as.formula(sprintf("~ %s", rhs)))
  } else if (cfg$estimator == "Firth PML") {
    if (identical(cfg$data_name, "d21")) {
      rhs <- sprintf("eng_c * %s + nonviolent + REGCHANGE + v2csreprss + time_in_campaign + lnnum_image_sum",
                     cfg$moderator)
    } else {
      rhs <- sprintf("eng_c * %s + nonviolent_camp + REGCHANGE + v2csreprss + lnlengthofcam + lnnum_image_sum",
                     cfg$moderator)
    }
    if (side == "refit") return(as.formula(sprintf("%s ~ %s", cfg$response, rhs)))
    return(as.formula(sprintf("~ %s", rhs)))
  } else {  # Cox PH
    rhs <- sprintf("eng_c * %s + nonviolent + REGCHANGE + v2csreprss + lnnum_image_sum",
                   cfg$moderator)
    if (side == "refit") return(as.formula(sprintf(
      "Surv(tstart, tstop, event_concession) ~ %s + cluster(navco21_id)", rhs)))
    return(as.formula(sprintf("~ %s", rhs)))
  }
}

refit_spec <- function(spec_id, sim_data) {
  cfg <- H21_SPECS[[spec_id]]
  res <- tryCatch({
    f <- build_formula(spec_id, "refit")
    if (cfg$estimator == "Weibull AFT") {
      m <- survreg(f, data = sim_data, dist = "weibull")
      tbl <- summary(m)$table
      list(p_focal = unname(tbl[cfg$focal, "p"]),
           beta_focal = unname(coef(m)[cfg$focal]), fit_ok = TRUE)
    } else if (cfg$estimator == "Firth PML") {
      m <- logistf(f, data = sim_data, control = ctrl_logistf, plcontrol = pl_ctrl, pl = TRUE)
      list(p_focal = unname(m$prob[cfg$focal]),
           beta_focal = unname(coef(m)[cfg$focal]), fit_ok = TRUE)
    } else {  # Cox PH
      m <- coxph(f, data = sim_data)
      tbl <- summary(m)$coefficients
      list(p_focal = unname(tbl[cfg$focal, "Pr(>|z|)"]),
           beta_focal = unname(coef(m)[cfg$focal]), fit_ok = TRUE)
    }
  }, error = function(e) list(p_focal = NA_real_, beta_focal = NA_real_, fit_ok = FALSE))
  res
}

simulate_power_for_spec <- function(spec_id, n_reps) {
  cfg <- H21_SPECS[[spec_id]]
  cat(sprintf("\n--- Simulating %s [%s × %s] (n_reps=%d) ---\n",
              spec_id, cfg$estimator, cfg$moderator, n_reps))
  fitted <- get_fitted(spec_id)
  beta_hat <- coef(fitted)
  if (!(cfg$focal %in% names(beta_hat))) {
    cat(sprintf("  *** focal '%s' not found; skipping ***\n", cfg$focal))
    return(NULL)
  }
  beta_focal_hat <- unname(beta_hat[cfg$focal])
  cat(sprintf("  beta_int (fitted) = %.4f  (focal = %s)\n", beta_focal_hat, cfg$focal))

  if (cfg$estimator == "Weibull AFT") {
    f_sim <- build_formula(spec_id, "sim")
    sim_fn <- function(beta_sim) sim_weibull(df, f_sim, beta_sim, fitted$scale,
                                              "duration", "event_concession")
  } else if (cfg$estimator == "Firth PML") {
    f_sim <- build_formula(spec_id, "sim")
    sim_base <- if (identical(cfg$data_name, "d21")) d21 else df
    sim_fn <- function(beta_sim) sim_firth(sim_base, f_sim, beta_sim, cfg$response)
  } else {  # Cox PH — simsurv exponential baseline (PE2 R25 fix)
    cat("  Building simsurv exponential-baseline DGP for Cox PH spec...\n")
    proxy <- tryCatch(build_cox_simsurv(spec_id),
                      error = function(e) {
                        cat(sprintf("  *** simsurv DGP build failed: %s ***\n", conditionMessage(e)))
                        NULL
                      })
    if (is.null(proxy)) return(NULL)
    sim_fn <- function(beta_sim) {
      beta_named <- beta_sim[setdiff(names(beta_sim), "(Intercept)")]
      common <- intersect(colnames(proxy$X_use), names(beta_named))
      sim <- tryCatch(simsurv(dist = "weibull",
                              lambdas = proxy$lambda, gammas = proxy$gamma,
                              x = as.data.frame(proxy$X_use[, common, drop = FALSE]),
                              betas = beta_named[common],
                              maxt = proxy$maxt),
                      error = function(e) NULL)
      if (is.null(sim)) stop("simsurv returned NULL")
      out <- proxy$base_df
      out$tstart <- 0
      out$tstop  <- pmin(sim$eventtime, proxy$maxt)
      out$event_concession <- as.integer(sim$status == 1 & sim$eventtime <= proxy$maxt)
      out
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
      vp <- p_values[!is.na(p_values)]
      # PE3 R25 convention: _overall = defensive (NA-fits as non-rejections);
      # _conditional = diagnostic. Backward-compatible aliases default to overall.
      raw_cond  <- if (length(vp) > 0) mean(vp < alpha) else NA_real_
      holm_cond <- if (length(vp) > 0) mean(vp < alpha / K_FAMILY) else NA_real_
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
    cat(sprintf("raw(0.05)=%.3f  holm(0.05)=%.3f  fit_ok=%.3f\n",
                out$raw_power[out$alpha == 0.05 & out$magnitude == mag][1],
                out$holm_proxy_power[out$alpha == 0.05 & out$magnitude == mag][1],
                fit_ok_rate))
  }
  out
}

safe_sim <- function(spec_id, n_reps) {
  tryCatch(simulate_power_for_spec(spec_id, n_reps),
           error = function(e) {
             cat(sprintf("\n*** %s SIMULATION FAILED: %s ***\n", spec_id, conditionMessage(e)))
             NULL
           })
}

t0 <- Sys.time()
results <- do.call(rbind, list(
  safe_sim("M3", N_REPS_FAST), safe_sim("M4", N_REPS_FAST),  # Firth first
  safe_sim("M9", N_REPS_FAST), safe_sim("M10", N_REPS_FAST), # CY Firth
  safe_sim("M1", N_REPS_FAST), safe_sim("M2", N_REPS_FAST),  # Weibull
  safe_sim("M7", N_REPS_COX),  safe_sim("M8", N_REPS_COX)    # Cox last
))
t1 <- Sys.time()
cat(sprintf("\nTotal sim wall time: %.1f min\n",
            as.numeric(difftime(t1, t0, units = "mins"))))

compute_mde <- function(power_df, alpha_val, target = 0.80) {
  sub <- subset(power_df, alpha == alpha_val & magnitude > 0)
  if (nrow(sub) == 0) return(c(raw_multiplier = Inf, holm_multiplier = Inf))
  rc <- sub$magnitude[!is.na(sub$raw_power) & sub$raw_power >= target]
  hc <- sub$magnitude[!is.na(sub$holm_proxy_power) & sub$holm_proxy_power >= target]
  c(raw_multiplier = ifelse(length(rc) > 0, min(rc), Inf),
    holm_multiplier = ifelse(length(hc) > 0, min(hc), Inf))
}

if (!is.null(results) && nrow(results) > 0) {
  summary_df <- do.call(rbind, lapply(unique(results$spec), function(sp) {
    sr <- subset(results, spec == sp); s05 <- subset(sr, alpha == 0.05)
    mde05 <- compute_mde(sr, 0.05)
    data.frame(
      spec = sp, estimator = s05$estimator[1], moderator = s05$moderator[1],
      beta_focal_hat = s05$beta_focal_hat[1],
      type1_at_null_05       = s05$raw_power[s05$magnitude == 0.0],
      raw_power_at_half_05   = s05$raw_power[s05$magnitude == 0.5],
      raw_power_at_fitted_05 = s05$raw_power[s05$magnitude == 1.0],
      raw_power_at_double_05 = s05$raw_power[s05$magnitude == 2.0],
      holm_power_at_fitted_05 = s05$holm_proxy_power[s05$magnitude == 1.0],
      raw_mde_mult_05        = mde05["raw_multiplier"],
      holm_mde_mult_05       = mde05["holm_multiplier"],
      stringsAsFactors = FALSE
    )
  }))
  cat("\n================================================================\n")
  cat("H2.1 POWER — PER-SPEC MDE SUMMARY (alpha = 0.05)\n")
  cat("================================================================\n")
  print(summary_df, row.names = FALSE)

  prefix <- if (SMOKE_MODE) "power_h21_smoke" else "power_h21"
  write.csv(results,    file.path(results_dir, sprintf("%s.csv", prefix)),
            row.names = FALSE)
  write.csv(summary_df, file.path(results_dir, sprintf("%s_summary.csv", prefix)),
            row.names = FALSE)
  cat(sprintf("\nResults written to: %s.csv (+ _summary.csv)\n", prefix))
}
