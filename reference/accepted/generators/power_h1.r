###############################################################################
# power_h1.r
# H1 CURVILINEAR — POWER DIAGNOSTICS (Round 24 / Path B Test 1 / D-58)
#
# Per Round 24 panel-locked plan (D1=B, D4=A, A2, D5=A):
#   - Coefficient MDE only: power for beta(eng^2)
#   - Magnitude grid: {0, beta_hat/4, beta_hat/2, beta_hat, 2 * beta_hat}
#   - Reps: 2000 for Weibull AFT (M1) / Firth PML (M2, M5); 500 for Cox PH (M4)
#   - Scope: H1 main panel only — M1, M2, M4, M5 (D-39 main panel)
#   - Per-spec independent simulation; Holm-adjusted-power proxy = power at alpha/k
#     (Bonferroni worst-case proxy for Holm step-down across the k=10 H1 family)
#   - MDE framing per D5: "Minimum Detectable Effect under our design at alpha=0.05"
#     NEVER use the words "post-hoc power" in the §5.2 prose
#
# Outputs:
#   results/power_h1.csv         — per-spec × per-magnitude power
#   results/power_h1_summary.csv — per-spec MDE estimates
#
# Supersedes legacy RR_datacode/pml/power_diagnostics.r (D-21).
###############################################################################

rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(logistf)
  library(simsurv)  # PE1 fix: M4 now uses simsurv exponential baseline
})

project_root <- Sys.getenv("ACCEPTED_ANALYSIS_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set ACCEPTED_ANALYSIS_ROOT to the original RR cowork directory.")
source_dir   <- file.path(project_root, "RR_datacode", "empirical clean")
results_dir  <- file.path(source_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Source H1_code.r to get fitted models + datasets in scope.
# We capture stdout to suppress the H1 results-table print.
suppressMessages(
  invisible(capture.output(
    sys.source(file.path(source_dir, "H1_code.r"), envir = environment())
  ))
)

cat("\n================================================================\n")
cat("H1 POWER DIAGNOSTICS — Round 24 / Path B Test 1\n")
cat("================================================================\n")

# Settings (panel-locked)
read_env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }
  parsed <- suppressWarnings(as.integer(value))
  if (length(parsed) != 1L || is.na(parsed) || parsed < 1L) {
    warning(sprintf("Ignoring invalid %s=%s; using default %d", name, value, default))
    return(default)
  }
  parsed
}

detect_parallel_cores <- function() {
  if (.Platform$OS.type == "windows") {
    return(1L)
  }
  detected <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (length(detected) != 1L || is.na(detected) || detected < 2L) {
    return(1L)
  }
  max(1L, as.integer(detected) - 1L)
}

N_REPS_FAST      <- read_env_int("POWER_H1_REPS_FAST", 2000L)   # Weibull AFT, Firth PML
N_REPS_COX       <- read_env_int("POWER_H1_REPS_COX", 500L)     # Cox PH counting-process (D1=B)
SMOKE_MODE       <- nzchar(Sys.getenv("POWER_H1_REPS_FAST", unset = "")) ||
                    nzchar(Sys.getenv("POWER_H1_REPS_COX", unset = ""))
MAG_GRID         <- c(0, 0.25, 0.5, 1.0, 2.0)  # multipliers of beta_hat (D4=A)
ALPHA_LEVELS     <- c(0.05, 0.10)
SEED             <- 24
K_FAMILY         <- 10  # H1 Holm family size (per H1_code.r main panel)
PARALLEL_CORES   <- detect_parallel_cores()

cat(sprintf("Settings: nreps_fast=%d  nreps_cox=%d  cores=%d\n",
            N_REPS_FAST, N_REPS_COX, PARALLEL_CORES))
if (SMOKE_MODE) {
  cat("Output mode: smoke-test files (replication count override detected)\n")
}
cat(sprintf("Magnitude grid: %s (multipliers of fitted beta_eng^2)\n",
            paste(MAG_GRID, collapse = ",")))
cat(sprintf("Holm family k=%d; Holm-adjusted-power proxy = power at alpha/k = alpha/%d\n",
            K_FAMILY, K_FAMILY))

ctrl_logistf <- logistf::logistf.control(maxit = 1000, maxstep = 0.5)
pl_ctrl      <- logistf::logistpl.control(maxit = 1000)

run_replications <- function(rep_ids, fn, force_serial = FALSE) {
  safe_fn <- function(rep_id) {
    tryCatch(fn(rep_id), error = function(e) NA_real_)
  }
  if (force_serial || PARALLEL_CORES <= 1L) {
    return(lapply(rep_ids, safe_fn))
  }
  parallel::mclapply(rep_ids, safe_fn, mc.cores = PARALLEL_CORES)
}

###############################################################################
# DGP simulators per H1 main-panel spec
# Each returns a fresh data frame with the simulated outcome column populated.
# Other covariates are held at their observed values (X fixed by design).
###############################################################################

# --- M2 / M5 Firth PML: simulate Bernoulli outcome from logit(eta) ---
complete_model_data <- function(base_df, formula_obj, extra_vars = character()) {
  vars <- unique(c(all.vars(formula_obj), extra_vars))
  vars <- vars[vars %in% names(base_df)]
  base_df[complete.cases(base_df[, vars, drop = FALSE]), , drop = FALSE]
}

sim_firth <- function(base_df, formula_obj, beta_sim, response_name) {
  sim_df <- complete_model_data(base_df, formula_obj)
  mm <- model.matrix(formula_obj, sim_df)
  # Align column order between mm and beta_sim
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

# --- M1 Weibull AFT: simulate event times from Weibull, censor at observed duration ---
sim_weibull <- function(base_df, formula_obj, beta_sim, scale_sigma,
                        duration_name = "duration",
                        event_name    = "event_concession") {
  sim_df <- complete_model_data(base_df, formula_obj, duration_name)
  mm <- model.matrix(formula_obj, sim_df)
  beta_aligned <- beta_sim[colnames(mm)]
  if (any(is.na(beta_aligned))) {
    stop(sprintf("Coefficient name mismatch: %s",
                 paste(setdiff(colnames(mm), names(beta_sim)), collapse = ", ")))
  }
  eta <- as.vector(mm %*% beta_aligned)
  # survreg Weibull: log T = eta + scale * Z where Z ~ standard Gumbel (extreme value)
  # Equivalently T ~ Weibull(shape = 1/scale, scale = exp(eta))
  T_sim <- rweibull(length(eta), shape = 1 / scale_sigma, scale = exp(eta))
  # Administrative censoring at original campaign duration
  obs_dur <- sim_df[[duration_name]]
  T_sim[!is.finite(T_sim)] <- obs_dur[!is.finite(T_sim)]
  min_duration <- min(obs_dur[is.finite(obs_dur) & obs_dur > 0], na.rm = TRUE) / 1000
  if (!is.finite(min_duration) || min_duration <= 0) {
    min_duration <- 1e-6
  }
  sim_duration <- pmax(pmin(T_sim, obs_dur), min_duration)
  out <- sim_df
  out[[duration_name]] <- sim_duration
  out[[event_name]]    <- as.integer(T_sim <= obs_dur)
  out
}

# --- DEPRECATED (PE1 R25 fix): fitted-baseline Cox CP simulator was degenerate
# at N=34 campaigns / 65 obs (Type-I -> 1.0 with cluster() vcov). Retained as
# dead code for archival reference only. M4 now uses simsurv exponential
# baseline (mirrors power_cox_specs.r). Safe to delete in a follow-up refactor
# along with sim_cox_from_baseline() and cumhaz_at(). ---
cumhaz_at <- function(basehaz_df, t) {
  idx <- findInterval(t, basehaz_df$time)
  ifelse(idx == 0L, 0, basehaz_df$hazard[idx])
}

sim_cox_from_baseline <- function(base_df, formula_obj, beta_sim, basehaz_df,
                                  id_name = "navco21_id",
                                  event_name = "event_concession") {
  mm <- model.matrix(formula_obj, base_df)
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  beta_aligned <- beta_sim[colnames(mm)]
  if (any(is.na(beta_aligned))) {
    stop(sprintf("Coefficient name mismatch: %s",
                 paste(setdiff(colnames(mm), names(beta_sim)), collapse = ", ")))
  }

  eta <- as.vector(mm %*% beta_aligned)
  h0_delta <- pmax(cumhaz_at(basehaz_df, base_df$tstop) -
                     cumhaz_at(basehaz_df, base_df$tstart), 0)
  interval_hazard <- h0_delta * exp(pmin(eta, 700))
  p_event <- ifelse(interval_hazard > 700, 1, 1 - exp(-interval_hazard))
  p_event <- pmin(pmax(p_event, 0), 1)

  out <- base_df
  out$.p_event <- p_event
  out <- out[order(out[[id_name]], out$tstart, out$tstop), ]
  keep <- rep(TRUE, nrow(out))
  out[[event_name]] <- 0L

  for (id in unique(out[[id_name]])) {
    idx <- which(out[[id_name]] == id)
    for (pos in seq_along(idx)) {
      row_idx <- idx[pos]
      if (runif(1) < out$.p_event[row_idx]) {
        out[[event_name]][row_idx] <- 1L
        if (pos < length(idx)) {
          keep[idx[(pos + 1L):length(idx)]] <- FALSE
        }
        break
      }
    }
  }

  out <- out[keep, , drop = FALSE]
  out$.p_event <- NULL
  out
}

###############################################################################
# Spec configurations
# Each entry holds: focal coefficient name, fitted model, base data, formula(s),
# DGP function, refit function, p-value extractor, n_reps.
###############################################################################

# --- "Any English effect" Wald test on (beta_eng, beta_eng^2): 2-df chi-sq.
#     BUG B3 fix (Round 25): this is an OMNIBUS test of H0: beta_eng =
#     beta_eng^2 = 0 — it answers "is there ANY English-share signal at all
#     (linear, curved, or both)", NOT specifically a curvature test.
#     The H1 main-panel CSV column previously named joint_wald_p was renamed
#     to any_eng_effect_p in H1_code.r for the same reason; we keep the
#     internal R helper name joint_wald_pair for backward continuity but the
#     output `test` column below is renamed to "any_eng_effect_eng_eng2".
#
#     Originally added to close collinearity-induced Type-I inflation
#     discovered in smoke testing (M1 60% Type-I at alpha=0.05) when the
#     simulation tested only the single-coefficient Wald on beta_eng^2.  The
#     2-df joint test has honest size; a single-coefficient companion is now
#     also reported (single_coef_p_eng2 below) so that the manuscript power
#     table covers BOTH the omnibus and the published-table 1-df test.
joint_wald_pair <- function(model, linear_name, quad_name,
                            vcov_override = NULL) {
  cf <- coef(model)
  if (!(linear_name %in% names(cf)) || !(quad_name %in% names(cf))) {
    return(NA_real_)
  }
  bvec <- c(cf[linear_name], cf[quad_name])
  V <- if (!is.null(vcov_override)) vcov_override else
       tryCatch(vcov(model), error = function(e) NULL)
  if (is.null(V) && inherits(model, "logistf")) V <- model$var
  if (is.null(V)) return(NA_real_)
  i1 <- match(linear_name, names(cf))
  i2 <- match(quad_name,   names(cf))
  V2 <- V[c(i1, i2), c(i1, i2)]
  if (any(!is.finite(c(bvec, V2)))) return(NA_real_)
  W <- tryCatch(as.numeric(t(bvec) %*% solve(V2) %*% bvec),
                error = function(e) NA_real_)
  if (!is.finite(W)) return(NA_real_)
  pchisq(W, df = 2, lower.tail = FALSE)
}

# --- Single-coefficient Wald test on beta_eng^2 (BUG B3 Part C addition).
#     This 1-df test mirrors what the published H1 table reports as p_eng2 in
#     the per-coefficient stars column.  Adding it gives the response memo
#     power coverage for the actual published test.  Note: under the joint
#     simulation DGP (both eng and eng^2 scaled together), the single-
#     coefficient test can have inflated Type-I under collinearity — that is
#     a property of the test, not a bug.  The omnibus joint Wald above is
#     the honest-size test; this one is the published-table-aligned test.
#     We do NOT use the LR test in simulation because it would require a
#     refit per replicate (~ 2x the runtime) for marginal additional info
#     beyond the Wald p_eng2 the manuscript already reports.
single_coef_p_eng2 <- function(model, quad_name, vcov_override = NULL) {
  cf <- coef(model)
  if (!(quad_name %in% names(cf))) return(NA_real_)
  b2 <- unname(cf[quad_name])
  V  <- if (!is.null(vcov_override)) vcov_override else
        tryCatch(vcov(model), error = function(e) NULL)
  if (is.null(V) && inherits(model, "logistf")) V <- model$var
  if (is.null(V)) return(NA_real_)
  i2 <- match(quad_name, names(cf))
  v2 <- V[i2, i2]
  if (any(!is.finite(c(b2, v2))) || v2 <= 0) return(NA_real_)
  z  <- b2 / sqrt(v2)
  2 * pnorm(-abs(z))
}

# Helper: linear and quadratic coef names per spec
focal_pair_names <- function(spec_id) {
  if (spec_id %in% c("M1", "M2"))      c("eng_prob_general", "I(eng_prob_general^2)")
  else if (spec_id %in% c("M4", "M5")) c("eng_prop_year",    "I(eng_prop_year^2)")
  else stop(sprintf("Unknown spec %s", spec_id))
}

# --- Refit functions return list(p_joint, p_eng2_single, beta_eng, beta_eng2,
#     fit_ok).  BUG B3 Part C: p_eng2_single added so the power table covers
#     BOTH the omnibus joint Wald and the published-table single-coefficient
#     test on beta_eng^2.
refit_M1 <- function(sim_data) {
  n_events <- sum(sim_data$event_concession == 1, na.rm = TRUE)
  n_censored <- sum(sim_data$event_concession == 0, na.rm = TRUE)
  if (n_events < 3L || n_censored < 3L ||
      any(!is.finite(sim_data$duration)) || any(sim_data$duration <= 0)) {
    return(list(p_joint = NA_real_, p_eng2_single = NA_real_,
                beta_eng = NA_real_, beta_eng2 = NA_real_, fit_ok = FALSE))
  }
  res <- tryCatch({
    m <- survreg(Surv(duration, event_concession) ~
                   eng_prob_general + I(eng_prob_general^2) +
                   nonviolent_camp + REGCHANGE + v2csreprss +
                   lnnum_image_sum + lnpop + colonized_english,
                 data = sim_data, dist = "weibull")
    pn <- focal_pair_names("M1")
    list(p_joint = joint_wald_pair(m, pn[1], pn[2]),
         p_eng2_single = single_coef_p_eng2(m, pn[2]),
         beta_eng  = unname(coef(m)[pn[1]]),
         beta_eng2 = unname(coef(m)[pn[2]]),
         fit_ok = TRUE)
  }, error = function(e) list(p_joint = NA_real_, p_eng2_single = NA_real_,
                              beta_eng = NA_real_, beta_eng2 = NA_real_,
                              fit_ok = FALSE))
  res
}
refit_M2 <- function(sim_data) {
  res <- tryCatch({
    m <- logistf(concession ~
                   eng_prob_general + I(eng_prob_general^2) +
                   nonviolent_camp + REGCHANGE + v2csreprss +
                   lnlengthofcam + lnnum_image_sum + lnpop + colonized_english,
                 data = sim_data, control = ctrl_logistf, plcontrol = pl_ctrl, pl = TRUE)
    pn <- focal_pair_names("M2")
    list(p_joint = joint_wald_pair(m, pn[1], pn[2]),
         p_eng2_single = single_coef_p_eng2(m, pn[2]),
         beta_eng  = unname(coef(m)[pn[1]]),
         beta_eng2 = unname(coef(m)[pn[2]]),
         fit_ok = TRUE)
  }, error = function(e) list(p_joint = NA_real_, p_eng2_single = NA_real_,
                              beta_eng = NA_real_, beta_eng2 = NA_real_,
                              fit_ok = FALSE))
  res
}
refit_M4 <- function(sim_data) {
  res <- tryCatch({
    m <- coxph(Surv(tstart, tstop, event_concession) ~
                 eng_prop_year + I(eng_prop_year^2) +
                 nonviolent + REGCHANGE + v2csreprss +
                 lnnum_image_sum + lnpop_yr + colonized_english +
                 cluster(navco21_id),
               data = sim_data)
    pn <- focal_pair_names("M4")
    # coxph with cluster() returns a robust vcov that should be used for Wald
    list(p_joint = joint_wald_pair(m, pn[1], pn[2], vcov_override = m$var),
         p_eng2_single = single_coef_p_eng2(m, pn[2], vcov_override = m$var),
         beta_eng  = unname(coef(m)[pn[1]]),
         beta_eng2 = unname(coef(m)[pn[2]]),
         fit_ok = TRUE)
  }, error = function(e) list(p_joint = NA_real_, p_eng2_single = NA_real_,
                              beta_eng = NA_real_, beta_eng2 = NA_real_,
                              fit_ok = FALSE))
  res
}
refit_M5 <- function(sim_data) {
  res <- tryCatch({
    m <- logistf(concession ~
                   eng_prop_year + I(eng_prop_year^2) +
                   nonviolent + REGCHANGE + v2csreprss +
                   time_in_campaign + lnnum_image_sum + lnpop_yr + colonized_english,
                 data = sim_data, control = ctrl_logistf, plcontrol = pl_ctrl, pl = TRUE)
    pn <- focal_pair_names("M5")
    list(p_joint = joint_wald_pair(m, pn[1], pn[2]),
         p_eng2_single = single_coef_p_eng2(m, pn[2]),
         beta_eng  = unname(coef(m)[pn[1]]),
         beta_eng2 = unname(coef(m)[pn[2]]),
         fit_ok = TRUE)
  }, error = function(e) list(p_joint = NA_real_, p_eng2_single = NA_real_,
                              beta_eng = NA_real_, beta_eng2 = NA_real_,
                              fit_ok = FALSE))
  res
}

###############################################################################
# Per-spec power simulation
###############################################################################

simulate_power_for_spec <- function(spec_id, n_reps) {
  cat(sprintf("\n--- Simulating %s (n_reps=%d) ---\n", spec_id, n_reps))

  # Spec-specific configuration
  cfg <- switch(spec_id,
    "M1" = list(
      model = model1, base_df = df, focal = "I(eng_prob_general^2)",
      refit = refit_M1,
      sim_fn = function(beta_sim) {
        sim_weibull(df,
                    formula_obj = ~ eng_prob_general + I(eng_prob_general^2) +
                               nonviolent_camp + REGCHANGE + v2csreprss +
                               lnnum_image_sum + lnpop + colonized_english,
                    beta_sim = beta_sim, scale_sigma = model1$scale,
                    duration_name = "duration", event_name = "event_concession")
      }
    ),
    "M2" = list(
      model = model2, base_df = df, focal = "I(eng_prob_general^2)",
      refit = refit_M2,
      sim_fn = function(beta_sim) {
        sim_firth(df,
                  formula_obj = ~ eng_prob_general + I(eng_prob_general^2) +
                                  nonviolent_camp + REGCHANGE + v2csreprss +
                                  lnlengthofcam + lnnum_image_sum + lnpop + colonized_english,
                  beta_sim = beta_sim, response_name = "concession")
      }
    ),
    "M4" = list(
      # M4: Cox PH on NAVCO 2.1 first-event panel d21_first_conc.
      # Per PE1/PE4 R25 fix: use simsurv exponential-baseline DGP (mirrors
      # power_cox_specs.r) instead of fitted-baseline sim_cox_from_baseline,
      # which was degenerate at small N (Type-I -> 1.0 with cluster() vcov).
      model = NULL, base_df = NULL, focal = "I(eng_prop_year^2)",
      refit = refit_M4,
      sim_fn = NULL  # built below after we extract calibration parameters
    ),
    "M5" = list(
      model = model5, base_df = d21, focal = "I(eng_prop_year^2)",
      refit = refit_M5,
      sim_fn = function(beta_sim) {
        sim_firth(d21,
                  formula_obj = ~ eng_prop_year + I(eng_prop_year^2) +
                                  nonviolent + REGCHANGE + v2csreprss +
                                  time_in_campaign + lnnum_image_sum + lnpop_yr + colonized_english,
                  beta_sim = beta_sim, response_name = "concession")
      }
    )
  )

  # PE1/PE4 R25 fix: M4 uses simsurv exponential-baseline DGP, NOT fitted-baseline.
  # The fitted-baseline approach (sim_cox_from_baseline above, retained as dead
  # code) was degenerate at N=34 / 65 obs producing Type-I -> 1.0 with cluster()
  # vcov. simsurv with exponential baseline calibrated to the observed event rate
  # (lambda = -log(1 - obs_event_rate) / max_followup, gamma = 1) is the same
  # approach that works in power_cox_specs.r and yields honest Type-I.
  if (spec_id == "M4") {
    cat("  Building simsurv exponential-baseline DGP for Cox PH spec M4...\n")
    cox_vars <- c("tstart", "tstop", "event_concession", "eng_prop_year",
                  "nonviolent", "REGCHANGE", "v2csreprss", "lnnum_image_sum",
                  "lnpop_yr", "colonized_english", "navco21_id")
    cox_base <- d21_first_conc[complete.cases(d21_first_conc[, cox_vars]), , drop = FALSE]
    rhs_proxy <- ~ eng_prop_year + I(eng_prop_year^2) +
                   nonviolent + REGCHANGE + v2csreprss +
                   lnnum_image_sum + lnpop_yr + colonized_english
    X_mm <- model.matrix(rhs_proxy, data = cox_base)
    X_use_m4 <- X_mm[, setdiff(colnames(X_mm), "(Intercept)"), drop = FALSE]
    obs_event_rate <- mean(cox_base$event_concession, na.rm = TRUE)
    max_followup <- max(cox_base$tstop, na.rm = TRUE)
    lambda_base <- -log(max(1 - obs_event_rate, 0.01)) / max_followup
    cfg$model   <- model4
    cfg$base_df <- cox_base
    cfg$sim_fn  <- function(beta_sim) {
      beta_named <- beta_sim[setdiff(names(beta_sim), "(Intercept)")]
      common <- intersect(colnames(X_use_m4), names(beta_named))
      sim <- tryCatch(simsurv(dist = "weibull",
                              lambdas = lambda_base, gammas = 1.0,
                              x = as.data.frame(X_use_m4[, common, drop = FALSE]),
                              betas = beta_named[common],
                              maxt = max_followup),
                      error = function(e) NULL)
      if (is.null(sim)) stop("simsurv returned NULL")
      out <- cox_base
      out$tstart <- 0
      out$tstop  <- pmin(sim$eventtime, max_followup)
      out$event_concession <- as.integer(sim$status == 1 & sim$eventtime <= max_followup)
      out
    }
  }

  # Joint-magnitude scaling per refined D4: scale BOTH beta_eng AND beta_eng^2
  # by the same multiplier (effect-size scaling of the entire H1 inverted-U
  # pattern), then test the joint H0: beta_eng = beta_eng^2 = 0 via Wald chi-sq
  # df=2. Fixes M1 collinearity-induced 60% Type-I from the single-coefficient
  # approach; substantively answers "is there an English-share signal at all".
  pn <- focal_pair_names(spec_id)
  beta_eng_hat   <- unname(coef(cfg$model)[pn[1]])
  beta_eng2_hat  <- unname(coef(cfg$model)[pn[2]])
  beta_hat       <- coef(cfg$model)
  cat(sprintf("  fitted: beta(eng) = %.4f, beta(eng^2) = %.4f\n",
              beta_eng_hat, beta_eng2_hat))

  # BUG B3 Part C: each replicate now returns BOTH the joint Wald p-value
  # (omnibus "any English effect") AND the single-coefficient Wald p on
  # beta_eng^2 (mirrors the published H1 table p_eng2).  Two test rows are
  # written per (spec, magnitude, alpha):
  #   test = "any_eng_effect_eng_eng2" : 2-df joint Wald (renamed from
  #                                      "joint_wald_eng_eng2" in BUG B3 fix
  #                                      to honestly label the omnibus null)
  #   test = "p_eng2_single_coef"      : 1-df Wald on beta_eng^2 only
  #                                      (matches published-table star)
  out <- data.frame()
  test_specs <- list(
    list(name = "any_eng_effect_eng_eng2", extractor = function(r) r$p_joint),
    list(name = "p_eng2_single_coef",      extractor = function(r) r$p_eng2_single)
  )

  for (mag in MAG_GRID) {
    beta_sim <- beta_hat
    beta_sim[pn[1]] <- beta_eng_hat  * mag
    beta_sim[pn[2]] <- beta_eng2_hat * mag

    cat(sprintf("  magnitude = %.2f x fitted (eng=%.3f, eng^2=%.3f)  ... ",
                mag, beta_sim[pn[1]], beta_sim[pn[2]]))

    # Run replications: each returns list(p_joint, p_eng2_single).
    rep_results <- run_replications(seq_len(n_reps), function(rep_id) {
      set.seed(SEED + rep_id * 1000 + which(MAG_GRID == mag))
      sim_data <- cfg$sim_fn(beta_sim)
      r <- cfg$refit(sim_data)
      list(p_joint = r$p_joint, p_eng2_single = r$p_eng2_single)
    }, force_serial = spec_id %in% c("M1", "M4"))

    p_joint_vec        <- as.numeric(vapply(rep_results,
                                            function(r) if (is.list(r)) r$p_joint else NA_real_,
                                            numeric(1)))
    p_eng2_single_vec  <- as.numeric(vapply(rep_results,
                                            function(r) if (is.list(r)) r$p_eng2_single else NA_real_,
                                            numeric(1)))

    for (ts in test_specs) {
      p_values <- if (ts$name == "any_eng_effect_eng_eng2") p_joint_vec else p_eng2_single_vec
      fit_ok_rate <- mean(!is.na(p_values))
      for (alpha in ALPHA_LEVELS) {
        valid_p <- p_values[!is.na(p_values)]
        # PE3 R25 convention preserved.
        raw_power_conditional   <- if (length(valid_p) > 0) mean(valid_p < alpha) else NA_real_
        holm_power_conditional  <- if (length(valid_p) > 0) mean(valid_p < alpha / K_FAMILY) else NA_real_
        raw_power_overall       <- sum(p_values < alpha,            na.rm = TRUE) / n_reps
        holm_power_overall      <- sum(p_values < alpha / K_FAMILY, na.rm = TRUE) / n_reps
        out <- rbind(out, data.frame(
          spec = spec_id,
          test = ts$name,
          beta_eng_hat   = beta_eng_hat,
          beta_eng2_hat  = beta_eng2_hat,
          magnitude = mag,
          beta_eng_sim  = beta_sim[pn[1]],
          beta_eng2_sim = beta_sim[pn[2]],
          n_reps = n_reps, fit_ok_rate = fit_ok_rate,
          alpha = alpha,
          raw_power_overall            = raw_power_overall,
          holm_proxy_power_overall     = holm_power_overall,
          raw_power_conditional        = raw_power_conditional,
          holm_proxy_power_conditional = holm_power_conditional,
          raw_power = raw_power_overall, holm_proxy_power = holm_power_overall,
          stringsAsFactors = FALSE
        ))
      }
    }
    # Print a per-magnitude line covering BOTH tests at alpha = 0.05.
    rp_any <- out$raw_power[out$test == "any_eng_effect_eng_eng2" &
                            out$alpha == 0.05 & out$magnitude == mag][1]
    rp_p2  <- out$raw_power[out$test == "p_eng2_single_coef" &
                            out$alpha == 0.05 & out$magnitude == mag][1]
    cat(sprintf("any_eff_pwr(0.05)=%.3f  p_eng2_pwr(0.05)=%.3f  fit_ok=%.3f\n",
                rp_any, rp_p2, mean(!is.na(p_joint_vec))))
  }
  out
}

###############################################################################
# Main loop
###############################################################################

safe_sim <- function(spec_id, n_reps) {
  res <- tryCatch(simulate_power_for_spec(spec_id, n_reps),
                  error = function(e) {
                    msg <- conditionMessage(e)
                    cat(sprintf("\n*** %s SIMULATION FAILED: %s ***\n", spec_id, msg))
                    cat(sprintf("    Writing NA placeholder rows for %s with failure_reason.\n", spec_id))
                    # PE1 R25 fix: defensive placeholder so the CSV preserves the
                    # panel scope (M1/M2/M4/M5) even when a spec crashes — prior
                    # behavior returned NULL and the spec silently disappeared
                    # from the production CSV.
                    # BUG B3 Part C: emit BOTH test rows per (mag, alpha).
                    test_names <- c("any_eng_effect_eng_eng2", "p_eng2_single_coef")
                    do.call(rbind, lapply(MAG_GRID, function(mag) {
                      do.call(rbind, lapply(ALPHA_LEVELS, function(a) {
                        do.call(rbind, lapply(test_names, function(tn) {
                          data.frame(
                            spec = spec_id, test = tn,
                            beta_eng_hat = NA_real_, beta_eng2_hat = NA_real_,
                            magnitude = mag,
                            beta_eng_sim = NA_real_, beta_eng2_sim = NA_real_,
                            n_reps = n_reps, fit_ok_rate = 0,
                            alpha = a,
                            raw_power_overall = NA_real_,
                            holm_proxy_power_overall = NA_real_,
                            raw_power_conditional = NA_real_,
                            holm_proxy_power_conditional = NA_real_,
                            raw_power = NA_real_, holm_proxy_power = NA_real_,
                            failure_reason = msg,
                            stringsAsFactors = FALSE)
                        }))
                      }))
                    }))
                  })
  if (!is.null(res) && !"failure_reason" %in% names(res)) res$failure_reason <- NA_character_
  res
}

t0 <- Sys.time()
results <- do.call(rbind, list(
  safe_sim("M2", N_REPS_FAST),  # Firth PML — cleanest
  safe_sim("M1", N_REPS_FAST),  # Weibull AFT
  safe_sim("M5", N_REPS_FAST),  # Firth PML yearly
  safe_sim("M4", N_REPS_COX)    # Cox PH on NAVCO 2.1 first-event panel.
                                 # Per PE1 R25 fix: simulates from simsurv with
                                 # an exponential baseline calibrated to the
                                 # observed event rate (mirrors
                                 # power_cox_specs.r), superseding the degenerate
                                 # fitted-baseline DGP that produced Type-I = 1.0
                                 # in earlier rounds. Power numbers should now
                                 # align with power_cox_specs.csv H1 M4 row.
))
t1 <- Sys.time()
cat(sprintf("\nTotal sim wall time: %.1f min\n",
            as.numeric(difftime(t1, t0, units = "mins"))))

# Compute MDE per spec at alpha = 0.05 (raw and Holm-adjusted-proxy)
compute_mde <- function(power_df, alpha_val, power_target = 0.80) {
  # MDE = smallest absolute magnitude multiplier at which power >= power_target
  sub <- subset(power_df, alpha == alpha_val & magnitude > 0)
  if (nrow(sub) == 0) return(c(raw = NA_real_, holm = NA_real_))
  raw_candidates <- sub$magnitude[!is.na(sub$raw_power) & sub$raw_power >= power_target]
  holm_candidates <- sub$magnitude[!is.na(sub$holm_proxy_power) &
                                     sub$holm_proxy_power >= power_target]
  raw_mde  <- if (length(raw_candidates) > 0) min(raw_candidates, na.rm = TRUE) else Inf
  holm_mde <- if (length(holm_candidates) > 0) min(holm_candidates, na.rm = TRUE) else Inf
  c(raw_multiplier = ifelse(is.finite(raw_mde),  raw_mde,  Inf),
    holm_multiplier = ifelse(is.finite(holm_mde), holm_mde, Inf))
}

# BUG B3 Part C: emit one summary row per (spec, test) so the response memo
# can read off MDE for both the omnibus and the published-table 1-df test.
summary_df <- do.call(rbind, lapply(unique(results$spec), function(sp) {
  spec_results <- subset(results, spec == sp)
  do.call(rbind, lapply(unique(spec_results$test), function(tn) {
    test_results <- subset(spec_results, test == tn)
    sub05 <- subset(test_results, alpha == 0.05)
    mde05 <- compute_mde(test_results, 0.05)
    data.frame(
      spec = sp,
      test = tn,
      beta_eng_hat   = sub05$beta_eng_hat[1],
      beta_eng2_hat  = sub05$beta_eng2_hat[1],
      type1_at_null_05         = sub05$raw_power[sub05$magnitude == 0.0],
      raw_power_at_half_05     = sub05$raw_power[sub05$magnitude == 0.5],
      raw_power_at_fitted_05   = sub05$raw_power[sub05$magnitude == 1.0],
      raw_power_at_double_05   = sub05$raw_power[sub05$magnitude == 2.0],
      holm_power_at_fitted_05  = sub05$holm_proxy_power[sub05$magnitude == 1.0],
      raw_mde_mult_05          = mde05["raw_multiplier"],
      holm_mde_mult_05         = mde05["holm_multiplier"],
      stringsAsFactors = FALSE
    )
  }))
}))

cat("\n================================================================\n")
cat("H1 POWER — PER-SPEC MDE SUMMARY (alpha = 0.05)\n")
cat("================================================================\n")
print(summary_df, row.names = FALSE)

# Write outputs
output_prefix <- if (SMOKE_MODE) "power_h1_smoke" else "power_h1"
out_csv      <- file.path(results_dir, sprintf("%s.csv", output_prefix))
summary_csv  <- file.path(results_dir, sprintf("%s_summary.csv", output_prefix))
write.csv(results,    out_csv,     row.names = FALSE)
write.csv(summary_df, summary_csv, row.names = FALSE)
cat(sprintf("\nResults  written to: %s\n", out_csv))
cat(sprintf("Summary  written to: %s\n", summary_csv))
