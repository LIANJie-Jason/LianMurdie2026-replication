###############################################################################
# H1_code.r — Curvilinear Hypothesis (H1)
#
# Tests the inverted-U relationship between English-language protest signs and
# campaign outcomes. The main-panel DV is concession (success or limited
# success). Strict-success twins and ordered DVs appear in the appendix.
#
# Data:
#   df_final.csv         (NAVCO 1.3 campaign-level, N=101, from data_clean.r)
#   df_navco21_panel.csv (NAVCO 2.1 source-campaign-year panel,
#                         built by data_clean_navco21.r)
#
# Specifications
# --------------
# MAIN PANEL — 4 concession-outcome models:
#   M1   Weibull AFT (NAVCO 1.3)  time to concession
#   M2   Firth PML   (NAVCO 1.3)  concession binary
#   M4   Cox PH      (NAVCO 2.1)  hazard of first concession (CP format)
#   M5   Firth PML   (NAVCO 2.1)  concession yearly binary
#
# APPENDIX — 6 strict-success or ordered-outcome models:
#   M1s  Weibull AFT (NAVCO 1.3)  time to SUCCESS (strict)
#   M2s  Firth PML   (NAVCO 1.3)  SUCCESS binary (strict)
#   M4s  Cox PH      (NAVCO 2.1)  hazard of SUCCESS (CP format, strict)
#   M5s  Firth PML   (NAVCO 2.1)  SUCCESS yearly binary (strict)
#   M3   bracl MPL_Jeffreys (NAVCO 1.3)  success3 (0/1/2)
#   M6   bracl MPL_Jeffreys (NAVCO 2.1)  progress_ord (6-level; failure separate)
#
# SENSITIVITY — outside the Holm/BH family:
#   M1_ln  log-normal AFT (NAVCO 1.3)  time to concession
#   M1_cox Cox PH         (NAVCO 1.3)  hazard of concession (campaign-lvl)
#   M4_AG  Cox PH         (NAVCO 2.1)  recurrent-event concession hazard
#   M5_GEE GEE logit      (NAVCO 2.1)  clustered yearly concession binary
#
# INFERENCE AND ROBUSTNESS
#   #1: Holm (FWER) + BH (FDR) p-value adjustment across the 10-spec family
#       for p_eng, p_eng2, any_eng_effect_p, curvature_lr_p
#       (`any_eng_effect_p` is an omnibus test, not a curvature test).
#   #2: M1 sensitivity under log-normal AFT + Cox PH (above)
#
# Controls (NAVCO 1.3): nonviolent_camp (1=nonviolent), REGCHANGE, v2csreprss, lnnum_image_sum,
#   lnpop, colonized_english  (+ lnlengthofcam for binary/ordered logit)
#   - REGCHANGE only (drop SECESSION) to avoid VIF=1.77M collinearity
#     and FSELFDET separation
#
# Controls (NAVCO 2.1): nonviolent (from prim_meth), REGCHANGE,
#   v2csreprss, lnnum_image_sum, lnpop_yr, colonized_english
#   (+ time_in_campaign raw for binary/ordered logit)
#
# Cox PH: counting-process format Surv(tstart, tstop, event)
#   - tstart = time_in_campaign - 1; tstop = time_in_campaign
#   - Ensures correct risk sets (Andersen & Gill 1982)
#
# Estimator notes
# ---------------
# Firth PML (logistf): profile penalized-likelihood CIs (pl=TRUE)
# bracl MPL_Jeffreys (brglm2): Wald CIs (profile unsupported for multi-eq bracl)
#   - Kosmidis & Firth (2011) extension to ordered outcomes
# Cox PH: NAVCO 2.1 first-event risk sets use navco21_id as the
#   source-campaign subject after the duplicate id-year merge.
#
# Author: Jie (Jason) Lian
# Final accepted specification: concession main panel with strict-success,
# ordered-outcome, distributional, and clustered-inference checks.
###############################################################################

# --- Setup ---
rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(geepack)
  library(survival)
  library(logistf)
  library(brglm2)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
rep_root <- Sys.getenv("REPLICATION_ROOT", unset = "")
if (!nzchar(rep_root)) {
  script_file <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
  rep_root <- dirname(dirname(dirname(script_file)))
}
source(file.path(rep_root, "code", "R", "00_setup.R"), local = FALSE)
P <- init_replication(rep_root)

navco13_path <- file.path(P$data, "df_final.csv")
navco21_path <- file.path(P$data, "df_navco21_panel.csv")
out_path     <- file.path(P$estimates, "H1_curvilinear_results.csv")
cache_path   <- file.path(P$cache, "H1_models.rds")
frames_path  <- file.path(P$cache, "H1_analysis_frames.rds")

###############################################################################
# NAVCO 1.3 PREPARATION
###############################################################################

df <- read.csv(navco13_path)
rep_assert_columns(
  df,
  c(
    "EYEAR", "BYEAR", "success", "limited", "success3", "viol", "nonviolent_camp", "REGCHANGE",
    "v2csreprss", "lnnum_image_sum", "lnpop", "colonized_english",
    "eng_prob_general", "lnlengthofcam"
  ),
  "df_final.csv"
)

df$viol              <- as.factor(df$viol)
df$nonviolent_camp   <- as.factor(df$nonviolent_camp)  # harmonized polarity (1=nonviolent)
df$colonized_english <- as.factor(df$colonized_english)
df$success3          <- ordered(df$success3)
df$REGCHANGE         <- as.numeric(as.character(df$REGCHANGE))
#: SECESSION dropped to avoid VIF=1.77M collinearity and FSELFDET separation

df$duration          <- df$EYEAR - df$BYEAR + 1

# Concession = success OR limited (partial-or-full concession per NAVCO)
df$event_concession  <- as.integer(
  as.numeric(as.character(df$success)) == 1 |
  as.numeric(as.character(df$limited)) == 1
)
df$concession        <- df$event_concession

# Strict success (NAVCO 1.3 success=1)
df$event_success     <- as.numeric(as.character(df$success))
df$success_bin       <- df$event_success

###############################################################################
# NAVCO 2.1 PREPARATION
###############################################################################

d21 <- read.csv(navco21_path)
rep_assert_columns(
  d21,
  c(
    "CAMPAIGN", "LOCATION", "Year", "panel_campaign_id", "nonviolent",
    "colonized_english", "camp_goals", "Num_Image", "progress",
    "concession", "time_in_campaign", "eng_prop_year", "v2csreprss", "lnpop_yr",
    "navco21_id"
  ),
  "df_navco21_panel.csv"
)
rep_assert_unique(d21, c("navco21_id", "Year"), "df_navco21_panel NAVCO 2.1 source-id-Year")

d21$nonviolent        <- as.factor(d21$nonviolent)   # prim_meth: 1=nonviolent
d21$colonized_english <- as.factor(d21$colonized_english)

# NAVCO 2.1 campaign goals: 0=regime change, 3=secession (codebook)
#: SECESSION dropped from controls (kept only REGCHANGE)
d21$REGCHANGE         <- as.numeric(d21$camp_goals == 0)

d21$lnnum_image_sum   <- log(d21$Num_Image)

# progress: 0=status quo, 1=visible gains, 2=limited concession,
# 3=sig concession, 4=complete success, 5=failure ( fix)
make_progress_ord6 <- function(x) {
  ordered(x, levels = c(5, 0, 1, 2, 3, 4))
}
d21$progress_ord      <- make_progress_ord6(d21$progress)

# Concession event: progress 2-4 that year
d21$event_concession  <- as.numeric(d21$concession)

# Strict success event: progress == 4 that year
d21$event_success_yr  <- as.numeric(d21$progress == 4)
d21$success_bin_yr    <- d21$event_success_yr

# Counting-process intervals for Cox PH
d21$tstart <- d21$time_in_campaign - 1
d21$tstop  <- d21$time_in_campaign

build_first_event_panel <- function(dat, id_col, event_col, order_col) {
  out <- dat %>%
    arrange(across(all_of(c(id_col, order_col)))) %>%
    group_by(across(all_of(id_col))) %>%
    mutate(
      .cum_event = cumsum(.data[[event_col]]),
      .prior_event = lag(.cum_event, default = 0L),
      .keep_risk = .prior_event == 0L,
      .event_first = as.integer(.data[[event_col]] == 1 & .cum_event == 1L)
    ) %>%
    filter(.keep_risk) %>%
    ungroup()

  out[[event_col]] <- out$.event_first

  out %>%
    dplyr::select(-.cum_event, -.prior_event, -.keep_risk, -.event_first)
}

assert_first_event_panel <- function(dat, id_col, event_col) {
  violations <- dat %>%
    group_by(across(all_of(id_col))) %>%
    summarise(
      n_events = sum(.data[[event_col]] == 1, na.rm = TRUE),
      rows_after_event = if (any(.data[[event_col]] == 1, na.rm = TRUE)) {
        dplyr::n() - max(which(.data[[event_col]] == 1))
      } else {
        0L
      },
      .groups = "drop"
    ) %>%
    filter(n_events > 1 | rows_after_event > 0)

  if (nrow(violations) > 0) {
 stop("First-event truncation failed for NAVCO 2.1 concession panel.", call. = FALSE)
  }
}

d21_first_conc <- build_first_event_panel(d21, "navco21_id", "event_concession", "time_in_campaign")
assert_first_event_panel(d21_first_conc, "navco21_id", "event_concession")

###############################################################################
# MAIN PANEL — NAVCO 1.3 (M1, M1s, M2, M2s)
###############################################################################

# M1: Weibull AFT, time to concession
model1 <- survreg(
  Surv(duration, event_concession) ~
    eng_prob_general + I(eng_prob_general^2) +
    nonviolent_camp + REGCHANGE + v2csreprss +
    lnnum_image_sum + lnpop + colonized_english,
  data = df, dist = "weibull"
)

# M1s: Weibull AFT, time to STRICT success
model1s <- survreg(
  Surv(duration, event_success) ~
    eng_prob_general + I(eng_prob_general^2) +
    nonviolent_camp + REGCHANGE + v2csreprss +
    lnnum_image_sum + lnpop + colonized_english,
  data = df, dist = "weibull"
)

# M2: Firth PML, concession binary
model2 <- logistf(
  concession ~
    eng_prob_general + I(eng_prob_general^2) +
    nonviolent_camp + REGCHANGE + v2csreprss +
    lnlengthofcam + lnnum_image_sum + lnpop + colonized_english,
  data = df,
  control   = logistf.control(maxit = 1000, maxstep = 0.5),
  plcontrol = logistpl.control(maxit = 1000),
  pl = TRUE
)

# M2s: Firth PML, STRICT success binary
model2s <- logistf(
  success_bin ~
    eng_prob_general + I(eng_prob_general^2) +
    nonviolent_camp + REGCHANGE + v2csreprss +
    lnlengthofcam + lnnum_image_sum + lnpop + colonized_english,
  data = df,
  control   = logistf.control(maxit = 1000, maxstep = 0.5),
  plcontrol = logistpl.control(maxit = 1000),
  pl = TRUE
)

###############################################################################
# APPENDIX — NAVCO 1.3 (M3)
###############################################################################

# M3: bracl ordered logit, success3 (0/1/2)
model3 <- bracl(
  success3 ~
    eng_prob_general + I(eng_prob_general^2) +
    nonviolent_camp + REGCHANGE + v2csreprss +
    lnlengthofcam + lnnum_image_sum + lnpop + colonized_english,
  data     = df,
  type     = "MPL_Jeffreys",
  link     = "logit",
  parallel = TRUE,
  maxit    = 500
)
model3_ci <- confint(model3, level = 0.95, type = "Wald")

###############################################################################
# MAIN PANEL — NAVCO 2.1 (M4, M4s, M5, M5s)
###############################################################################

# M4: Cox PH, hazard of FIRST concession (counting-process)
model4 <- coxph(
  Surv(tstart, tstop, event_concession) ~
    eng_prop_year + I(eng_prop_year^2) +
    nonviolent + REGCHANGE + v2csreprss +
    lnnum_image_sum + lnpop_yr + colonized_english +
    cluster(navco21_id),
  data = d21_first_conc
)

# M4s: Cox PH, hazard of STRICT success (counting-process)
model4s <- coxph(
  Surv(tstart, tstop, event_success_yr) ~
    eng_prop_year + I(eng_prop_year^2) +
    nonviolent + REGCHANGE + v2csreprss +
    lnnum_image_sum + lnpop_yr + colonized_english +
    cluster(navco21_id),
  data = d21
)

# M5: Firth PML, yearly concession binary
model5 <- logistf(
  concession ~
    eng_prop_year + I(eng_prop_year^2) +
    nonviolent + REGCHANGE + v2csreprss +
    time_in_campaign + lnnum_image_sum + lnpop_yr + colonized_english,
  data      = d21,
  control   = logistf.control(maxit = 1000, maxstep = 0.5),
  plcontrol = logistpl.control(maxit = 1000),
  pl        = TRUE
)

fit_gee21_binary <- function(dv_col) {
  vars <- c(
    dv_col, "eng_prop_year", "nonviolent", "REGCHANGE", "v2csreprss",
    "time_in_campaign", "lnnum_image_sum", "lnpop_yr", "colonized_english",
    "navco21_id"
  )
  dd <- d21[complete.cases(d21[, vars]), vars, drop = FALSE]
  dd <- dd %>% arrange(navco21_id, time_in_campaign)
  geeglm(
    as.formula(sprintf(
      "%s ~ eng_prop_year + I(eng_prop_year^2) + nonviolent + REGCHANGE + v2csreprss + time_in_campaign + lnnum_image_sum + lnpop_yr + colonized_english",
      dv_col
    )),
    id     = navco21_id,
    data   = dd,
    family = binomial,
    corstr = "exchangeable"
  )
}

# M5s: Firth PML, yearly STRICT success binary
model5s <- logistf(
  success_bin_yr ~
    eng_prop_year + I(eng_prop_year^2) +
    nonviolent + REGCHANGE + v2csreprss +
    time_in_campaign + lnnum_image_sum + lnpop_yr + colonized_english,
  data      = d21,
  control   = logistf.control(maxit = 1000, maxstep = 0.5),
  plcontrol = logistpl.control(maxit = 1000),
  pl        = TRUE
)

model5_gee <- fit_gee21_binary("concession")

###############################################################################
# APPENDIX — NAVCO 2.1 (M6)
###############################################################################

# M6: bracl ordered logit, yearly progress_ord (6-level; failure separate)
model6 <- bracl(
  progress_ord ~
    eng_prop_year + I(eng_prop_year^2) +
    nonviolent + REGCHANGE + v2csreprss +
    time_in_campaign + lnnum_image_sum + lnpop_yr + colonized_english,
  data     = d21,
  type     = "MPL_Jeffreys",
  link     = "logit",
  parallel = TRUE,
  maxit    = 500
)
model6_ci <- confint(model6, level = 0.95, type = "Wald")

# M4_AG: Cox PH, recurrent-event concession hazard (Andersen-Gill robustness)
model4_ag <- coxph(
  Surv(tstart, tstop, event_concession) ~
    eng_prop_year + I(eng_prop_year^2) +
    nonviolent + REGCHANGE + v2csreprss +
    lnnum_image_sum + lnpop_yr + colonized_english +
    cluster(navco21_id),
  data = d21
)

###############################################################################
# M1 SENSITIVITY — NAVCO 1.3
###############################################################################
# Check potential Weibull distributional leverage by replicating M1 under
# (a) log-normal AFT (non-monotonic hazard family) and (b) semiparametric Cox
# PH (no parametric baseline assumption).
# NOTE: NAVCO 1.3 is campaign-level (one row per campaign), NOT panel. Cox PH
#      here uses standard Surv(duration, event), not counting-process — Cox
#      counting-process is only for M4/M4s on the NAVCO 2.1 panel.
# These sensitivity rows are OUTSIDE the 10-spec family for Holm/BH correction
# (they are robustness checks of a single model, not additional hypothesis tests).

# M1_ln: log-normal AFT, time to concession
model1_ln <- survreg(
  Surv(duration, event_concession) ~
    eng_prob_general + I(eng_prob_general^2) +
    nonviolent_camp + REGCHANGE + v2csreprss +
    lnnum_image_sum + lnpop + colonized_english,
  data = df, dist = "lognormal"
)

# M1_cox: Cox PH (campaign-level, not counting-process)
model1_cox <- coxph(
  Surv(duration, event_concession) ~
    eng_prob_general + I(eng_prob_general^2) +
    nonviolent_camp + REGCHANGE + v2csreprss +
    lnnum_image_sum + lnpop + colonized_english,
  data = df
)

###############################################################################
# COX PH DIAGNOSTICS
###############################################################################

cox_ph_diag <- function(tag, model, focal_name, focal_label = "focal term") {
  z <- tryCatch(cox.zph(model), error = function(e) NULL)
  if (is.null(z)) {
    cat(sprintf("[cox.zph] %-7s : test failed (skipped)\n", tag))
    return(data.frame(
      model = tag,
      ph_global_p = NA_real_,
      ph_int_p = NA_real_,
      ph_focal_p = NA_real_,
      ph_flag = "PH_TEST_FAILED",
      stringsAsFactors = FALSE
    ))
  }

  tab <- z$table
  global_p <- if ("GLOBAL" %in% rownames(tab)) tab["GLOBAL", "p"] else NA_real_
  focal_p <- if (focal_name %in% rownames(tab)) tab[focal_name, "p"] else NA_real_
  ph_flag <- if (is.finite(focal_p) && focal_p < 0.05) {
    "PH_FAIL_FOCAL"
  } else if (is.finite(global_p) && global_p < 0.05) {
    "PH_FAIL_GLOBAL"
  } else if (!is.finite(focal_p)) {
    "PH_FOCAL_MISSING"
  } else {
    ""
  }

  cat(sprintf("[cox.zph] %-7s | global p = %s | %s p = %s\n",
              tag,
              if (is.finite(global_p)) formatC(global_p, format = "f", digits = 4) else "NA",
              focal_label,
              if (is.finite(focal_p)) formatC(focal_p, format = "f", digits = 4) else "NA"))
  data.frame(
    model = tag,
    ph_global_p = global_p,
    ph_int_p = NA_real_,
    ph_focal_p = focal_p,
    ph_flag = ph_flag,
    stringsAsFactors = FALSE
  )
}

cox_recurrent_ph_placeholder <- function(tag) {
  data.frame(
    model = tag,
    ph_global_p = NA_real_,
    ph_int_p = NA_real_,
    ph_focal_p = NA_real_,
    ph_flag = "PH_NOT_REPORTED_RECURRENT",
    stringsAsFactors = FALSE
  )
}

cat("\n--- Cox PH diagnostic (cox.zph scaled Schoenfeld residuals) ---\n")
cox_diag_tbl <- bind_rows(
  cox_ph_diag("M4",     model4,     "I(eng_prop_year^2)",    "quadratic English term"),
  cox_ph_diag("M4s",    model4s,    "I(eng_prop_year^2)",    "quadratic English term"),
  cox_ph_diag("M1_cox", model1_cox, "I(eng_prob_general^2)", "quadratic English term"),
  cox_recurrent_ph_placeholder("M4_AG")
)
cat("--- end Cox PH diagnostic ---\n\n")

###############################################################################
# EXTRACTION HELPERS
###############################################################################

# Curvature likelihood-ratio test helper.
# The joint Wald test (renamed any_eng_effect_p) is an OMNIBUS "any English
# effect" test that rejects under a strong linear effect even with zero
# curvature. The conceptually correct curvature test is the 1-df likelihood-
# ratio (or deviance) test of H0: beta_eng^2 = 0, comparing the full quadratic
# model against a refitted linear-only nested model. Wald p_eng2 is also a
# valid 1-df curvature test but is sensitive to the parameterization (Hauck-
# Donner effect under near-separation in Firth PML). We report BOTH p_eng2
# (already computed above) AND the LR test for triangulation.
#
# Implementation: refit the same model with the squared term dropped, then
# extract the LR statistic. Per estimator class:
#   - survreg: 2 * (model$loglik[2] - linear$loglik[2]) ~ chi-sq df=1
#               (anova.survreg fails with cryptic "incorrect number of
#                dimensions" on this nested pair; manual loglik diff is
#                numerically identical to anova when anova works.)
#   - coxph: 2 * (model$loglik[2] - linear$loglik[2]) ~ chi-sq df=1
#               (anova.coxph refuses with cluster/robust variances:
#                "Can't do anova tables with robust variances".  The Cox
#                partial loglik does NOT depend on the robust-variance
#                computation, so the manual statistic is well-defined.)
#   - logistf: anova.logistf(linear, full) returns chi-sq + p directly.
#               Uses penalized-likelihood-ratio test per Heinze & Schemper
#               (2002).  Manual 2*(loglik[2] - lin$loglik[2]) is NOT
#               equivalent here — the penalty terms differ across nested
#               models.
#   - bracl: deviance difference, df = 1.
#   - geeglm: no true likelihood; return NA.
# Returns NA on any refit / numerical failure.
compute_curvature_lr_p <- function(model, data, linear_name, quad_name, ...) {
  if (is.null(data)) return(NA_real_)
  tryCatch({
    f_full    <- formula(model)
    drop_term <- paste0(". ~ . - ", quad_name)
    f_linear  <- update(f_full, drop_term)
    if (inherits(model, "survreg")) {
      m_lin <- survreg(f_linear, data = data, dist = model$dist)
      LR    <- 2 * (model$loglik[2] - m_lin$loglik[2])
      if (!is.finite(LR) || LR < 0) return(NA_real_)
      return(pchisq(LR, df = 1, lower.tail = FALSE))
    } else if (inherits(model, "coxph")) {
 # Re-run coxph on same data; cluster/strata preserved via formula.
      m_lin <- coxph(f_linear, data = data)
      LR    <- 2 * (model$loglik[2] - m_lin$loglik[2])
      if (!is.finite(LR) || LR < 0) return(NA_real_)
      return(pchisq(LR, df = 1, lower.tail = FALSE))
    } else if (inherits(model, "logistf")) {
      m_lin <- logistf(f_linear, data = data, pl = TRUE,
                       control   = logistf.control(maxit = 1000, maxstep = 0.5),
                       plcontrol = logistpl.control(maxit = 1000))
      anv <- anova(m_lin, model)
 # anova.logistf returns a list with $pval (named numeric).
      p <- tryCatch(unname(anv$pval), error = function(e) NA_real_)
      if (length(p) != 1L || !is.finite(p)) return(NA_real_)
      return(p)
    } else if (inherits(model, "bracl")) {
      m_lin <- bracl(f_linear, data = data, type = "MPL_Jeffreys",
                     link = "logit", parallel = TRUE, maxit = 500)
      LR <- m_lin$deviance - model$deviance
      if (!is.finite(LR) || LR < 0) return(NA_real_)
      return(pchisq(LR, df = 1, lower.tail = FALSE))
    } else if (inherits(model, "geeglm")) {
 # GEE has no true likelihood — return NA.
      return(NA_real_)
    } else {
      return(NA_real_)
    }
  }, error = function(e) NA_real_)
}

extract_result <- function(tag, estimator, data_label, response, model,
                            linear_name, quad_name, N, n_events, eng_range,
                            data = NULL) {
  b1 <- unname(coef(model)[linear_name])
  b2 <- unname(coef(model)[quad_name])
  tp <- -b1 / (2 * b2)

 # Helper: robustly fetch 2x2 vcov block for [linear, quadratic] by index.
 # Some estimators (logistf$var) return matrices without dimnames — use positional lookup.
  get_vcov_block <- function(mdl) {
    nm <- names(coef(mdl))
    i1 <- match(linear_name, nm)
    i2 <- match(quad_name,   nm)
    if (is.na(i1) || is.na(i2)) return(NULL)
    V_full <- tryCatch(vcov(mdl), error = function(e) NULL)
    if (is.null(V_full) && inherits(mdl, "logistf")) V_full <- mdl$var
    if (is.null(V_full)) return(NULL)
    V_full[c(i1, i2), c(i1, i2)]
  }

  if (inherits(model, "survreg")) {
    tbl <- summary(model)$table
    p1  <- tbl[linear_name, "p"]
    p2  <- tbl[quad_name, "p"]
  } else if (inherits(model, "coxph")) {
    tbl <- summary(model)$coefficients
    p1  <- tbl[linear_name, "Pr(>|z|)"]
    p2  <- tbl[quad_name, "Pr(>|z|)"]
  } else if (inherits(model, "logistf")) {
    p1  <- unname(model$prob[linear_name])
    p2  <- unname(model$prob[quad_name])
  } else if (inherits(model, "geeglm")) {
    tbl <- summary(model)$coefficients
    p1  <- tbl[linear_name, "Pr(>|W|)"]
    p2  <- tbl[quad_name, "Pr(>|W|)"]
  } else if (inherits(model, "bracl")) {
    tbl <- summary(model)$coefficients
    p1  <- tbl[linear_name, "Pr(>|z|)"]
    p2  <- tbl[quad_name, "Pr(>|z|)"]
  } else {
    p1 <- NA_real_; p2 <- NA_real_
  }
  V <- get_vcov_block(model)

 # "Any English effect" Wald test: H0: beta_eng = beta_eng2 = 0 (chi-sq, df=2)
# This is an omnibus test that rejects under either a linear-only effect,
# curvature, or both. It is not a curvature test on its own; see
# compute_curvature_lr_p for the 1-df curvature test.
  any_eng_effect_p <- NA_real_
  if (!is.null(V) && all(is.finite(c(b1, b2))) && all(is.finite(V))) {
    bvec <- c(b1, b2)
    W    <- tryCatch(as.numeric(t(bvec) %*% solve(V) %*% bvec), error = function(e) NA_real_)
    if (is.finite(W)) any_eng_effect_p <- pchisq(W, df = 2, lower.tail = FALSE)
  }

# The curvature LR test refits the model without I(...^2) and
 # compares via likelihood-ratio / deviance test. df=1. NA for GEE.
  curvature_lr_p <- compute_curvature_lr_p(model, data, linear_name, quad_name)

 # Concavity check: inverted-U signature by estimator parameterization.
  #   - Weibull AFT (survreg): log T = Xb.  Inverted-U on P(event) <=> U on log T
  #     <=> beta2 > 0.
  #   - Firth PML / Cox PH (logit/hazard, log-hazard): inverted-U <=> beta2 < 0.
  #   - bracl (brglm2) parallel adjacent-category logit: parameterizes
  #     log(pi_j / pi_{j+1}) = alpha_j + x'beta with outcomes coded low-to-high,
  #     so beta > 0 lowers P(higher category).  Inverted-U on outcome <=> beta2 > 0
#     so the expected quadratic sign is positive.
  expected_b2_sign <- if (estimator %in% c("Weibull AFT", "log-normal AFT", "bracl")) "+" else "-"
  concavity_consistent <- (expected_b2_sign == "+" && b2 > 0) ||
                          (expected_b2_sign == "-" && b2 < 0)

 # Turning point in sample range?
  tp_in_range <- is.finite(tp) && tp >= eng_range[1] && tp <= eng_range[2]

# Turning-point inferential support.
 # Compute delta-method SE + 95% CI and Fieller (1954) 95% CI for TP = -b1/(2*b2).
 # Both methods use the 2x2 [b1, b2] vcov block already extracted as V above.
# Delta-method gradient:
  #   d TP / d b1 = -1 / (2 b2)
  #   d TP / d b2 =  b1 / (2 b2^2)
 # Var(TP) = (1/(4 b2^2)) Var(b1) + (b1^2/(4 b2^4)) Var(b2)
  #            - (b1/(2 b2^3)) Cov(b1, b2)
 # Fieller's method solves the quadratic A g^2 + B g + C = 0 derived from
  #   (b1 + 2 g b2)^2 <= z^2 ( Var(b1) + 4 g^2 Var(b2) + 4 g Cov(b1, b2)),
 # with A = 4 b2^2 - 4 z^2 Var(b2),
  #      B = 4 b1 b2 - 4 z^2 Cov(b1, b2),
  #      C = b1^2  - z^2 Var(b1).
 # If A <= 0 (b2 not separated from zero at the chosen level), the
 # confidence set is unbounded -- itself informative for the inverted-U claim.
  z95         <- qnorm(0.975)
  tp_se       <- NA_real_
  tp_lo_delta <- NA_real_
  tp_up_delta <- NA_real_
  tp_lo_fieller    <- NA_real_
  tp_up_fieller    <- NA_real_
  tp_fieller_unbd  <- NA
  if (!is.null(V) && all(is.finite(c(b1, b2))) && all(is.finite(V)) && abs(b2) > 0) {
    var_b1 <- V[1, 1]; var_b2 <- V[2, 2]; cov_b12 <- V[1, 2]
    var_tp <- (1 / (4 * b2^2)) * var_b1 +
              (b1^2 / (4 * b2^4)) * var_b2 -
              (b1 / (2 * b2^3)) * cov_b12
    if (is.finite(var_tp) && var_tp >= 0) {
      tp_se       <- sqrt(var_tp)
      tp_lo_delta <- tp - z95 * tp_se
      tp_up_delta <- tp + z95 * tp_se
    }
    A_q <- 4 * b2^2  - 4 * z95^2 * var_b2
    B_q <- 4 * b1 * b2 - 4 * z95^2 * cov_b12
    C_q <- b1^2 - z95^2 * var_b1
    disc <- B_q^2 - 4 * A_q * C_q
    if (is.finite(A_q) && A_q > 0 && is.finite(disc) && disc >= 0) {
      r1 <- (-B_q - sqrt(disc)) / (2 * A_q)
      r2 <- (-B_q + sqrt(disc)) / (2 * A_q)
 # Fieller boundary roots are for g where (b1 + 2 g b2) = 0 i.e. g = -b1/(2b2) = TP.
      tp_lo_fieller   <- min(r1, r2)
      tp_up_fieller   <- max(r1, r2)
      tp_fieller_unbd <- FALSE
    } else {
      tp_fieller_unbd <- TRUE
    }
  }
  tp_ci_in_support_delta   <- is.finite(tp_lo_delta) && is.finite(tp_up_delta) &&
                              tp_lo_delta >= eng_range[1] && tp_up_delta <= eng_range[2]
  tp_ci_in_support_fieller <- isFALSE(tp_fieller_unbd) &&
                              is.finite(tp_lo_fieller) && is.finite(tp_up_fieller) &&
                              tp_lo_fieller >= eng_range[1] &&
                              tp_up_fieller <= eng_range[2]

 # Events per variable (Peduzzi >=10 rule of thumb); NA for ordered/AFT
  k_params <- length(coef(model))
  epv      <- ifelse(is.na(n_events), NA_real_, n_events / k_params)

  data.frame(
    model = tag, estimator = estimator, data = data_label, response = response,
    beta_eng = b1, p_eng = p1,
    beta_eng2 = b2, p_eng2 = p2,
    any_eng_effect_p = any_eng_effect_p,
    curvature_lr_p = curvature_lr_p,
    turning_point = tp, tp_in_range = tp_in_range,
    tp_se = tp_se,
    tp_lower_delta = tp_lo_delta, tp_upper_delta = tp_up_delta,
    tp_lower_fieller = tp_lo_fieller, tp_upper_fieller = tp_up_fieller,
    tp_fieller_unbounded = tp_fieller_unbd,
    tp_ci_in_support_delta = tp_ci_in_support_delta,
    tp_ci_in_support_fieller = tp_ci_in_support_fieller,
    expected_b2_sign = expected_b2_sign,
    concavity_consistent = concavity_consistent,
    N = N, n_events = n_events, epv = epv,
    sig_05_eng2 = as.logical(!is.na(p2) && p2 < 0.05),
    sig_10_eng2 = as.logical(!is.na(p2) && p2 < 0.10),
    stringsAsFactors = FALSE
  )
}

# Sample sizes and event counts from model-specific complete-case frames
cc_stats <- function(dat, vars, event_var = NA_character_) {
  cc_mask <- complete.cases(dat[, vars, drop = FALSE])
  N <- sum(cc_mask)

  if (is.na(event_var) || !(event_var %in% names(dat))) {
    n_events <- NA_integer_
  } else {
    ev_vec <- suppressWarnings(as.numeric(dat[[event_var]]))
    n_events <- as.integer(sum(ev_vec[cc_mask] == 1, na.rm = TRUE))
  }

  list(N = N, n_events = n_events)
}

vars_13_surv <- function(dv) {
  c(
    "duration", dv, "eng_prob_general", "nonviolent_camp", "REGCHANGE", "v2csreprss",
    "lnnum_image_sum", "lnpop", "colonized_english"
  )
}

vars_13_bin <- function(dv) {
  c(
    dv, "eng_prob_general", "nonviolent_camp", "REGCHANGE", "v2csreprss",
    "lnlengthofcam", "lnnum_image_sum", "lnpop", "colonized_english"
  )
}

vars_21_surv <- function(dv) {
  c(
    "tstart", "tstop", dv, "eng_prop_year", "nonviolent", "REGCHANGE",
    "v2csreprss", "lnnum_image_sum", "lnpop_yr", "colonized_english", "navco21_id"
  )
}

vars_21_bin <- function(dv) {
  c(
    dv, "eng_prop_year", "nonviolent", "REGCHANGE", "v2csreprss",
    "time_in_campaign", "lnnum_image_sum", "lnpop_yr", "colonized_english"
  )
}

s_m1     <- cc_stats(df,  vars_13_surv("event_concession"), "event_concession")
s_m1s    <- cc_stats(df,  vars_13_surv("event_success"),    "event_success")
s_m2     <- cc_stats(df,  vars_13_bin("concession"),        "concession")
s_m2s    <- cc_stats(df,  vars_13_bin("success_bin"),       "success_bin")
s_m3     <- cc_stats(df,  vars_13_bin("success3"))
s_m4     <- cc_stats(d21_first_conc, vars_21_surv("event_concession"), "event_concession")
s_m4s    <- cc_stats(d21, vars_21_surv("event_success_yr"), "event_success_yr")
s_m4_ag  <- cc_stats(d21, vars_21_surv("event_concession"), "event_concession")
s_m5     <- cc_stats(d21, vars_21_bin("concession"),        "concession")
s_m5_gee <- cc_stats(d21, vars_21_bin("concession"),        "concession")
s_m5s    <- cc_stats(d21, vars_21_bin("success_bin_yr"),    "success_bin_yr")
s_m6     <- cc_stats(d21, vars_21_bin("progress_ord"))
s_m1_ln  <- cc_stats(df,  vars_13_surv("event_concession"), "event_concession")
s_m1_cox <- cc_stats(df,  vars_13_surv("event_concession"), "event_concession")

# In-sample ranges for turning-point check
rng_eng13 <- range(df$eng_prob_general, na.rm = TRUE)
rng_eng21 <- range(d21$eng_prop_year,   na.rm = TRUE)

results <- rbind(
  extract_result("M1",  "Weibull AFT", "NAVCO 1.3", "Time to concession",      model1,  "eng_prob_general", "I(eng_prob_general^2)", s_m1$N,   s_m1$n_events,   rng_eng13, data = df),
  extract_result("M1s", "Weibull AFT", "NAVCO 1.3", "Time to success (strict)", model1s, "eng_prob_general", "I(eng_prob_general^2)", s_m1s$N,  s_m1s$n_events,  rng_eng13, data = df),
  extract_result("M2",  "Firth PML",   "NAVCO 1.3", "Concession (binary)",      model2,  "eng_prob_general", "I(eng_prob_general^2)", s_m2$N,   s_m2$n_events,   rng_eng13, data = df),
  extract_result("M2s", "Firth PML",   "NAVCO 1.3", "Success (strict, binary)", model2s, "eng_prob_general", "I(eng_prob_general^2)", s_m2s$N,  s_m2s$n_events,  rng_eng13, data = df),
  extract_result("M3",  "bracl",       "NAVCO 1.3", "success3 (0/1/2) APP",     model3,  "eng_prob_general", "I(eng_prob_general^2)", s_m3$N,   s_m3$n_events,   rng_eng13, data = df),
  extract_result("M4",  "Cox PH",      "NAVCO 2.1", "Hazard of first concession (CY)",  model4,  "eng_prop_year",    "I(eng_prop_year^2)",    s_m4$N,   s_m4$n_events,   rng_eng21, data = d21_first_conc),
  extract_result("M4s", "Cox PH",      "NAVCO 2.1", "Hazard of success (CY, strict)", model4s, "eng_prop_year", "I(eng_prop_year^2)", s_m4s$N,  s_m4s$n_events,  rng_eng21, data = d21),
  extract_result("M5",  "Firth PML",   "NAVCO 2.1", "Concession (CY, binary)",  model5,  "eng_prop_year",    "I(eng_prop_year^2)",    s_m5$N,   s_m5$n_events,   rng_eng21, data = d21),
  extract_result("M5s", "Firth PML",   "NAVCO 2.1", "Success (CY, strict, bin.)", model5s, "eng_prop_year",  "I(eng_prop_year^2)",    s_m5s$N,  s_m5s$n_events,  rng_eng21, data = d21),
  extract_result("M6",  "bracl",       "NAVCO 2.1", "progress_ord (6-level; failure separate) APP",   model6,  "eng_prop_year",    "I(eng_prop_year^2)",    s_m6$N,   s_m6$n_events,   rng_eng21, data = d21)
)

###############################################################################
# FAMILY-WISE MULTIPLE-TESTING CORRECTION
###############################################################################
# The H1 family has 10 curvilinear specifications tested on the same
# theoretical claim. Apply Holm (FWER, strong control) and BH
#      (FDR, less conservative) to four p-value families:
#        - p_eng: single-coefficient Wald on the LINEAR English term
#        - p_eng2: single-coefficient Wald on the SQUARED English term
#        - any_eng_effect_p: 2-df joint Wald (omnibus "any English effect";
#                             rejects under linear-only OR curvature OR both —
#                             not a curvature test on its own)
#        - curvature_lr_p: 1-df likelihood-ratio test of curvature (full
#                             vs linear-only nested model; NA for GEE)
# Only the 10-spec family is adjusted. Sensitivity rows (M1_ln, M1_cox, M4_AG,
# M5_GEE) appended below are robustness checks, not additional family members.

# Pass explicit n = nrow(results) so that any future NA
# p-value (e.g., a Firth/bracl convergence failure) does NOT silently shrink
# the family from 10 to <10. Mirrors the H21/H22/H3 explicit-n pattern set by
# the multiplicity family.
# At this point in the script `results` contains EXACTLY the 10-spec family
# (M1, M1s, M2, M2s, M3, M4, M4s, M5, M5s, M6). The four sensitivity rows
# (M1_ln, M1_cox, M4_AG, M5_GEE) are appended AFTER the p.adjust block below
# and explicitly receive NA in *_holm/*_bh columns (they are not part of the
# multiplicity family).
.k_h1_family <- nrow(results)
stopifnot(.k_h1_family == 10L)

results$p_eng_holm           <- p.adjust(results$p_eng,            method = "holm", n = .k_h1_family)
results$p_eng_bh             <- p.adjust(results$p_eng,            method = "BH",   n = .k_h1_family)
results$p_eng2_holm          <- p.adjust(results$p_eng2,           method = "holm", n = .k_h1_family)
results$p_eng2_bh            <- p.adjust(results$p_eng2,           method = "BH",   n = .k_h1_family)
# any_eng_effect_p is a 2-df omnibus test.
results$any_eng_effect_holm  <- p.adjust(results$any_eng_effect_p, method = "holm", n = .k_h1_family)
results$any_eng_effect_bh    <- p.adjust(results$any_eng_effect_p, method = "BH",   n = .k_h1_family)
# curvature_lr_p is a 1-df LR test of curvature.
results$curvature_lr_holm    <- p.adjust(results$curvature_lr_p,   method = "holm", n = .k_h1_family)
results$curvature_lr_bh      <- p.adjust(results$curvature_lr_p,   method = "BH",   n = .k_h1_family)

# Flags for adjusted significance (post-correction)
results$sig_05_eng2_holm <- !is.na(results$p_eng2_holm) & results$p_eng2_holm < 0.05
results$sig_10_eng2_holm <- !is.na(results$p_eng2_holm) & results$p_eng2_holm < 0.10
results$sig_05_eng2_bh   <- !is.na(results$p_eng2_bh)   & results$p_eng2_bh   < 0.05
results$sig_10_eng2_bh   <- !is.na(results$p_eng2_bh)   & results$p_eng2_bh   < 0.10

###############################################################################
# SENSITIVITY ROWS (outside family correction)
###############################################################################

sens_rows <- rbind(
  extract_result("M1_ln",  "log-normal AFT", "NAVCO 1.3",
                 "Time to concession (log-normal)", model1_ln,
                 "eng_prob_general", "I(eng_prob_general^2)",
                 s_m1_ln$N,  s_m1_ln$n_events, rng_eng13, data = df),
  extract_result("M1_cox", "Cox PH",         "NAVCO 1.3",
                 "Hazard of concession (Cox PH, campaign-lvl)", model1_cox,
                 "eng_prob_general", "I(eng_prob_general^2)",
                 s_m1_cox$N, s_m1_cox$n_events, rng_eng13, data = df),
  extract_result("M4_AG",  "Cox PH",         "NAVCO 2.1",
                 "Hazard of concession (CY, recurrent-event)", model4_ag,
                 "eng_prop_year", "I(eng_prop_year^2)",
                 s_m4_ag$N, s_m4_ag$n_events, rng_eng21, data = d21),
  extract_result("M5_GEE", "GEE logit",      "NAVCO 2.1",
                 "Concession (CY, binary, source-clustered GEE)", model5_gee,
                 "eng_prop_year", "I(eng_prop_year^2)",
                 s_m5_gee$N, s_m5_gee$n_events, rng_eng21, data = d21)
)
# Sensitivity rows don't receive family-wise adjustment (not part of test family)
sens_rows$p_eng_holm           <- NA_real_
sens_rows$p_eng_bh             <- NA_real_
sens_rows$p_eng2_holm          <- NA_real_
sens_rows$p_eng2_bh            <- NA_real_
# Omnibus any-English-effect adjustments.
sens_rows$any_eng_effect_holm  <- NA_real_
sens_rows$any_eng_effect_bh    <- NA_real_
sens_rows$curvature_lr_holm    <- NA_real_
sens_rows$curvature_lr_bh      <- NA_real_
sens_rows$sig_05_eng2_holm     <- NA
sens_rows$sig_10_eng2_holm     <- NA
sens_rows$sig_05_eng2_bh       <- NA
sens_rows$sig_10_eng2_bh       <- NA

results <- rbind(results, sens_rows)
results <- results %>%
  left_join(cox_diag_tbl, by = "model") %>%
  mutate(
    ph_global_p = ifelse(is.na(ph_global_p), NA_real_, ph_global_p),
    ph_int_p = ifelse(is.na(ph_int_p), NA_real_, ph_int_p),
    ph_focal_p = ifelse(is.na(ph_focal_p), NA_real_, ph_focal_p),
    ph_flag = ifelse(is.na(ph_flag), "", ph_flag),
    cox_estimand = case_when(
      estimator != "Cox PH" ~ "",
      grepl("_AG$", model) ~ "recurrent_event_concession",
      grepl("campaign-lvl|campaign-level", response, ignore.case = TRUE) ~ "campaign_level_concession",
      grepl("first concession", response, ignore.case = TRUE) ~ "first_concession",
      grepl("success", response, ignore.case = TRUE) ~ "strict_success",
      grepl("concession", response, ignore.case = TRUE) ~ "concession",
      TRUE ~ "cox_unspecified"
    ),
    ph_flag = case_when(
      estimator == "Cox PH" & grepl("_AG$", model) ~ "PH_NOT_REPORTED_RECURRENT",
      TRUE ~ ph_flag
    ),
    interpretation_note = case_when(
      ph_flag == "PH_NOT_REPORTED_RECURRENT" ~
        "Recurrent-event Cox robustness estimand; first-event PH diagnostic not reported",
      ph_flag == "PH_FAIL_FOCAL" ~
        "PH focal-term failure; interpret cautiously",
      ph_flag == "PH_FAIL_GLOBAL" ~
        "PH global failure; interpret cautiously",
      ph_flag == "PH_TEST_FAILED" ~
        "PH test failed; interpret cautiously",
      ph_flag == "PH_FOCAL_MISSING" ~
        "PH focal-term row missing from cox.zph output",
      TRUE ~ ""
    )
  )

###############################################################################
# PRINT FULL COEFFICIENT TABLES (: all terms, quadratic at top)
###############################################################################

cat("\n==== M1 Weibull AFT — Time to concession (NAVCO 1.3) ====\n")
print(summary(model1))
cat("\n==== M1s Weibull AFT — Time to STRICT success (NAVCO 1.3) ====\n")
print(summary(model1s))
cat("\n==== M2 Firth PML — Concession binary (NAVCO 1.3) ====\n")
print(summary(model2))
cat("\n==== M2s Firth PML — STRICT success binary (NAVCO 1.3) ====\n")
print(summary(model2s))
cat("\n==== M3 bracl — success3 (APPENDIX, NAVCO 1.3) ====\n")
print(summary(model3)); cat("95% Wald CIs:\n"); print(model3_ci)
cat("\n==== M4 Cox PH — Hazard of FIRST concession (NAVCO 2.1 CY, CP) ====\n")
print(summary(model4))
cat("\n==== M4s Cox PH — Hazard of STRICT success (NAVCO 2.1 CY, CP) ====\n")
print(summary(model4s))
cat("\n==== M4_AG Cox PH — Recurrent-event concession hazard (NAVCO 2.1 CY, CP) ====\n")
print(summary(model4_ag))
cat("\n==== M5 Firth PML — Concession yearly binary (NAVCO 2.1) ====\n")
print(summary(model5))
cat("\n==== M5_GEE GEE logit — Concession yearly binary, source-clustered (NAVCO 2.1) ====\n")
print(summary(model5_gee))
cat("\n==== M5s Firth PML — STRICT success yearly binary (NAVCO 2.1) ====\n")
print(summary(model5s))
cat("\n==== M6 bracl — progress_ord (6-level; APPENDIX, NAVCO 2.1) ====\n")
print(summary(model6)); cat("95% Wald CIs:\n"); print(model6_ci)

cat("\n==== M1_ln log-normal AFT — Time to concession (SENSITIVITY, NAVCO 1.3) ====\n")
print(summary(model1_ln))
cat("\n==== M1_cox Cox PH — Hazard of concession (SENSITIVITY, NAVCO 1.3) ====\n")
print(summary(model1_cox))

###############################################################################
# SUMMARY AND OUTPUT
###############################################################################

cat("\n================================================================\n")
cat("H1 CURVILINEAR — COMPACT HEADLINE (10-spec family + sensitivity)\n")
cat("================================================================\n")

fmt_p <- function(p) ifelse(is.na(p), "   NA ", sprintf("%6.3f", p))
fmt_b <- function(b) ifelse(is.na(b), "    NA ", sprintf("%7.2f", b))

headline <- data.frame(
  model     = results$model,
  estimator = substr(results$estimator, 1, 15),
  data      = substr(results$data, 1, 9),
  beta2     = fmt_b(results$beta_eng2),
  p2_raw    = fmt_p(results$p_eng2),
  p2_Holm   = fmt_p(results$p_eng2_holm),
  p2_BH     = fmt_p(results$p_eng2_bh),
# Joint Wald is an "any English effect" omnibus test, not curvature.
  pAny_raw     = fmt_p(results$any_eng_effect_p),
  pAny_Holm    = fmt_p(results$any_eng_effect_holm),
  pAny_BH      = fmt_p(results$any_eng_effect_bh),
  pCurvLR_raw  = fmt_p(results$curvature_lr_p),
  pCurvLR_Holm = fmt_p(results$curvature_lr_holm),
  pCurvLR_BH   = fmt_p(results$curvature_lr_bh),
  TP        = sprintf("%5.3f", results$turning_point),
  TPinRng   = results$tp_in_range,
  TP_CIdlt  = ifelse(is.na(results$tp_lower_delta) | is.na(results$tp_upper_delta), "[NA]",
                     sprintf("[%5.3f,%5.3f]", results$tp_lower_delta, results$tp_upper_delta)),
  TP_CIfie  = ifelse(isTRUE(results$tp_fieller_unbounded), "[unbnd]",
                     ifelse(is.na(results$tp_lower_fieller) | is.na(results$tp_upper_fieller), "[NA]",
                            sprintf("[%5.3f,%5.3f]", results$tp_lower_fieller, results$tp_upper_fieller))),
  concavOK  = results$concavity_consistent,
  N         = results$N,
  ev        = results$n_events,
  stringsAsFactors = FALSE
)
print(headline, row.names = FALSE)

cat("\n--- NOTES ---\n")
cat("p2_raw = raw two-sided p for beta(eng^2) (single-coefficient Wald)\n")
cat("p2_Holm = Holm-Bonferroni adjusted across 10-spec family (NA for sensitivity rows)\n")
cat("p2_BH = Benjamini-Hochberg FDR adjusted across 10-spec family\n")
cat("pAny_raw = raw joint Wald p for H0: beta(eng) = beta(eng^2) = 0 (df=2)\n")
cat(" NOTE: this tests 'any English effect (linear OR squared)' — NOT curvature.\n")
cat("pAny_Holm = Holm-Bonferroni adjusted any-English-effect p across family\n")
cat("pAny_BH = BH-adjusted any-English-effect p across family\n")
cat("pCurvLR_raw = raw 1-df LR test of curvature (full vs linear-only nested model)\n")
cat(" Conceptually correct curvature test; complements Wald p_eng2.\n")
cat(" NA for GEE (no true likelihood).\n")
cat(" NOTE: for Firth PML, pCurvLR_raw == p_eng2 by construction\n")
cat(" (logistf$prob already reports the per-coefficient PLR p,\n")
cat(" which equals the 1-df nested anova p; this is not a bug).\n")
cat("pCurvLR_Holm = Holm-adjusted curvature LR p across family\n")
cat("pCurvLR_BH = BH-adjusted curvature LR p across family\n")
cat("TP = turning point = -b1/(2*b2); TPinRng = TP in sample range\n")
cat("TP_CIdlt = 95% CI on TP via delta-method\n")
cat("TP_CIfie = 95% CI on TP via Fieller (1954); '[unbnd]' = beta_eng^2 not separated\n")
cat(" from zero at the 95% level (informative for inverted-U claim)\n")
cat("concavOK = sign(beta_eng^2) matches inverted-U for the estimator's parameterization\n")
cat(" (Weibull AFT & bracl: +; Firth PML & Cox PH: -; log-normal AFT: +)\n")
cat("M5_GEE = exchangeable-correlation GEE companion for the yearly concession model;\n")
cat(" use as the source-campaign clustered inference check for repeated source-years.\n")
cat("M5s/M6 remain pooled yearly appendix models and should be treated as weaker evidence.\n")

cat("\n================================================================\n")
cat("H1 CURVILINEAR — FULL RESULTS (all columns)\n")
cat("================================================================\n")
print(results, row.names = FALSE)

expected_result_cols <- c(
  "model", "estimator", "data", "response", "beta_eng", "p_eng", "beta_eng2",
# `any_eng_effect_p` is omnibus; `curvature_lr_p` is the curvature test.
  "p_eng2", "any_eng_effect_p", "curvature_lr_p",
  "turning_point", "tp_in_range",
  "tp_se", "tp_lower_delta", "tp_upper_delta",
  "tp_lower_fieller", "tp_upper_fieller", "tp_fieller_unbounded",
  "tp_ci_in_support_delta", "tp_ci_in_support_fieller",
  "expected_b2_sign",
  "concavity_consistent", "N", "n_events", "epv", "p_eng_holm", "p_eng_bh",
  "p_eng2_holm", "p_eng2_bh",
  "any_eng_effect_holm", "any_eng_effect_bh",
  "curvature_lr_holm", "curvature_lr_bh",
  "ph_global_p", "ph_int_p", "ph_focal_p", "ph_flag", "cox_estimand",
  "interpretation_note"
)
missing_result_cols <- setdiff(expected_result_cols, names(results))
if (length(missing_result_cols) > 0) {
  stop(
    sprintf("results is missing expected columns: %s", paste(missing_result_cols, collapse = ", ")),
    call. = FALSE
  )
}
if (nrow(results) != 14) {
  stop(sprintf("Expected 14 H1 result rows, found %d.", nrow(results)), call. = FALSE)
}
if (anyDuplicated(results$model) > 0) {
 stop("Duplicate model labels detected in H1 results.", call. = FALSE)
}

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
write.csv(results, out_path, row.names = FALSE)
models <- list(
  M1 = model1, M1s = model1s, M2 = model2, M2s = model2s,
  M3 = model3, M4 = model4, M4s = model4s, M5 = model5,
  M5s = model5s, M6 = model6, M1_ln = model1_ln,
  M1_cox = model1_cox, M4_AG = model4_ag, M5_GEE = model5_gee
)
stopifnot(identical(names(models), results$model))
saveRDS(models, cache_path, version = 3)
saveRDS(list(df = df, d21 = d21, d21_first_conc = d21_first_conc), frames_path, version = 3)
cat(sprintf("\nResults written to: %s\n", out_path))
cat(sprintf("Model bundle written to: %s\n", cache_path))
cat(sprintf("Analysis frames written to: %s\n", frames_path))
