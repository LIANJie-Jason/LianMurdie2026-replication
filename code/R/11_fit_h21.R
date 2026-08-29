###############################################################################
# H21_code.r
# H2.1: POLITICAL VULNERABILITY TO EXTERNAL PRESSURE
#
# Tests whether English-language protest signs have stronger positive effects
# on campaign outcomes when the state is more politically vulnerable to
# external pressure — operationalized by two V-Dem indicators:
#   v2svdomaut: Domestic Autonomy  (higher = less externally constrained)
#   v2svstterr: State Authority Over Territory (higher = stronger reach)
#
# Theoretical logic:
# English signs mobilize international audiences and leverage external
# pressure. That leverage should bite harder on states that are more
# dependent on external actors (lower domestic autonomy) or that have
# weaker internal authority (lower state-territory control). In centered
# form: the English marginal effect should weaken as the moderator rises.
#
# SPECIFICATION (centered and parsimonious):
#   - eng, moderators mean-centered (VIF drops from 19–60 to ~1)
#   - Controls (NAVCO 1.3): nonviolent_camp (1=nonviolent), REGCHANGE, v2csreprss, lnnum_image_sum
#       (+ lnlengthofcam for binary/ordered DV models)
#   - Controls (NAVCO 2.1): nonviolent, REGCHANGE, v2csreprss, lnnum_image_sum
#       (+ time_in_campaign for binary/ordered DV models)
#   - lnpop/colonized_english dropped to improve EPV
#   - Cox PH uses counting-process format (Andersen & Gill 1982)
#
# DUAL-DV STRUCTURE:
#   Report BOTH concession (success or limited, per NAVCO 1.3) AND strict
#   success across estimator families. Concession is the mechanism-relevant
#   DV (partial government response); strict success is the harder outcome.
#   Ordered-logit models (success3 / progress_ord) include both levels by
#   construction.
#
# MODELS (20 core + 2 robustness rows):
#   NAVCO 1.3 (N=101, campaign-level):
#     M1    Weibull AFT × domaut  — time to concession
#     M1s   Weibull AFT × domaut  — time to success (strict)
#     M2    Weibull AFT × stterr  — time to concession
#     M2s   Weibull AFT × stterr  — time to success (strict)
#     M3    Firth PML   × domaut  — concession (binary)
#     M3s   Firth PML   × domaut  — success (strict, binary)
#     M4    Firth PML   × stterr  — concession (binary)
#     M4s   Firth PML   × stterr  — success (strict, binary)
#     M5    bracl MPL   × domaut  — success3 (ordered 0/1/2)
#     M6    bracl MPL   × stterr  — success3 (ordered 0/1/2)
#   NAVCO 2.1 (source-campaign-year panel):
#     M7    Cox PH × domaut       — hazard of first concession (CP, clustered)
#     M7s   Cox PH × domaut       — hazard of success   (CP, clustered)
#     M8    Cox PH × stterr       — hazard of first concession (CP, clustered)
#     M8s   Cox PH × stterr       — hazard of success   (CP, clustered)
#     M9    Firth PML × domaut    — yearly concession (binary)
#     M9s   Firth PML × domaut    — yearly success    (binary)
#     M10   Firth PML × stterr    — yearly concession (binary)
#     M10s  Firth PML × stterr    — yearly success    (binary)
#     M11   bracl MPL × domaut    — progress_ord (ordered 0-4)
#     M12   bracl MPL × stterr    — progress_ord (ordered 0-4)
#
# EXPECTED INTERACTION SIGN (H2.1):
#   Benefit direction by estimator (which sign of eng_c means "English helps"):
#     - Weibull AFT: benefit = NEGATIVE eng (shorter time to event)
#                            ⇒ interaction sign under H2.1: POSITIVE
#     - Firth PML (logit): benefit = POSITIVE eng (more likely event)
#                            ⇒ interaction sign: NEGATIVE
#     - Cox PH (log-hazard): benefit = POSITIVE eng (higher hazard)
#                            ⇒ interaction sign: NEGATIVE
#     - bracl MPL (parallel adjacent-category logit, low-to-high ordered):
#                            β > 0 lowers P(higher category), so
#                            benefit = NEGATIVE eng
#                            ⇒ interaction sign: POSITIVE
#   The bracl sign convention follows its low-to-high outcome parameterization.
#
# OUTPUT:
#   RR_datacode/empirical clean/results/H21_results.csv — 24-row compact headline table.
#   Includes 20 core rows plus 4 robustness rows
#   (2 recurrent-event Cox + 2 clustered-GEE yearly concession checks).
#   Columns: model, estimator, data, response, moderator,
#            beta_eng_c, p_eng_c, beta_mod_c, p_mod_c,
#            beta_int,   p_int,   joint_wald_p (df=2, H0: β_eng=β_int=0),
#            expected_int_sign, sign_consistent,
#            N, n_events, epv, sig_05_int, sig_10_int
#
# Holm is applied within the 8-row main panel
# {M1,M2,M3,M4,M7,M8,M9,M10} and BH within the 12-row appendix core panel.
# M9 and M10 keep the main panel symmetric across survival and binary
# estimators in the two datasets.
# The 4 robustness rows (M7_AG, M8_AG, M9_GEE, M10_GEE) are excluded from
# the multiplicity families and reported descriptively only.
#
# Author: Jie (Jason) Lian
# Final accepted specification.
###############################################################################

# --- Setup ---
rm(list = ls())
library(dplyr)
library(geepack)
library(survival)   # survreg (Weibull AFT), coxph (Cox PH)
library(logistf)    # logistf (Firth PML with profile penalized CI)
library(brglm2)     # bracl (parallel adjacent-category logit, Jeffreys prior)

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
out_path     <- file.path(P$estimates, "H21_results.csv")
cache_path   <- file.path(P$cache, "H21_models.rds")
frames_path  <- file.path(P$cache, "H21_analysis_frames.rds")

assert_required_cols <- function(dat, required_cols, label) {
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop(
      sprintf("%s is missing required columns: %s", label, paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }
}

assert_unique_keys <- function(dat, keys, label) {
  dupes <- dat[duplicated(dat[keys]) | duplicated(dat[keys], fromLast = TRUE), keys, drop = FALSE]
  if (nrow(dupes) > 0) {
    stop(sprintf("Duplicate %s keys detected.", label), call. = FALSE)
  }
}

df  <- read.csv(navco13_path)
d21 <- read.csv(navco21_path)
assert_required_cols(
  df,
  c(
    "EYEAR", "BYEAR", "success", "limited", "success3", "viol", "nonviolent_camp", "REGCHANGE",
    "v2csreprss", "eng_prob_general", "lnnum_image_sum", "lnlengthofcam",
    "v2svdomaut", "v2svstterr"
  ),
  "df_final.csv"
)
assert_required_cols(
  d21,
  c(
    "CAMPAIGN", "LOCATION", "Year", "panel_campaign_id", "nonviolent",
    "camp_goals", "Num_Image", "progress", "concession", "time_in_campaign",
    "eng_prop_year", "v2svdomaut", "v2svstterr", "v2csreprss", "navco21_id"
  ),
  "df_navco21_panel.csv"
)
assert_unique_keys(df, c("CAMPAIGN", "LOCATION"), "df_final CAMPAIGN-LOCATION")
assert_unique_keys(d21, c("navco21_id", "Year"), "df_navco21_panel NAVCO 2.1 source-id-Year")

###############################################################################
# NAVCO 1.3 PREPARATION
###############################################################################

df$viol      <- as.factor(df$viol)
df$nonviolent_camp <- as.factor(df$nonviolent_camp)  # harmonized polarity (1=nonviolent)
df$success3  <- ordered(df$success3)
df$REGCHANGE <- as.numeric(as.character(df$REGCHANGE))

# Survival outcomes
df$duration          <- df$EYEAR - df$BYEAR + 1
df$event_concession  <- as.integer(
  as.numeric(as.character(df$success)) == 1 |
  as.numeric(as.character(df$limited)) == 1
)
df$event_success     <- as.numeric(as.character(df$success))

# Binary outcomes
df$concession_bin <- df$event_concession
df$success_bin    <- df$event_success

# Mean-center interaction constituents
df$eng_c    <- df$eng_prob_general - mean(df$eng_prob_general, na.rm = TRUE)
df$domaut_c <- df$v2svdomaut       - mean(df$v2svdomaut,       na.rm = TRUE)
df$stterr_c <- df$v2svstterr       - mean(df$v2svstterr,       na.rm = TRUE)

cat(sprintf("NAVCO 1.3 centering: mean(eng)=%.4f, mean(domaut)=%.4f, mean(stterr)=%.4f\n",
            mean(df$eng_prob_general, na.rm = TRUE),
            mean(df$v2svdomaut,       na.rm = TRUE),
            mean(df$v2svstterr,       na.rm = TRUE)))

# Full-sample event counts (diagnostic only; per-model EPV uses CC-mask counts
# computed in the sample-size section below)
cat(sprintf("NAVCO 1.3 events (full sample): concession=%d, success=%d (N=%d)\n",
            sum(df$event_concession, na.rm = TRUE),
            sum(df$event_success,    na.rm = TRUE),
            nrow(df)))

###############################################################################
# NAVCO 2.1 PREPARATION
###############################################################################

d21$nonviolent       <- as.factor(d21$nonviolent)
d21$REGCHANGE        <- as.numeric(d21$camp_goals == 0)
d21$lnnum_image_sum  <- log(d21$Num_Image)
# progress codes (raw NAVCO 2.1):
#   0 = status quo (no campaign-year change)
#   1 = visible gains short of concessions
#   2 = limited concessions achieved
#   3 = significant concessions achieved
#   4 = complete success
#   5 = FAILURE (campaign ended in defeat)
# Per data_clean_navco21.r: use a 6-level ordered scale with
# failure < status quo < visible gains < limited concession <
# significant concession < complete success.
make_progress_ord6 <- function(x) {
  ordered(x, levels = c(5, 0, 1, 2, 3, 4))
}
d21$progress_ord     <- make_progress_ord6(d21$progress)
d21$event_concession <- as.numeric(d21$concession)       # binary: progress 2-4
d21$event_success_yr <- as.numeric(d21$progress == 4)    # strict success
d21$success_bin_yr   <- d21$event_success_yr
d21$concession_bin   <- d21$event_concession

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

# Mean-center interaction constituents
d21$eng_c    <- d21$eng_prop_year - mean(d21$eng_prop_year, na.rm = TRUE)
d21$domaut_c <- d21$v2svdomaut    - mean(d21$v2svdomaut,    na.rm = TRUE)
d21$stterr_c <- d21$v2svstterr    - mean(d21$v2svstterr,    na.rm = TRUE)

d21_first_conc <- build_first_event_panel(d21, "navco21_id", "event_concession", "time_in_campaign")
assert_first_event_panel(d21_first_conc, "navco21_id", "event_concession")

cat(sprintf("NAVCO 2.1 centering: mean(eng)=%.4f, mean(domaut)=%.4f, mean(stterr)=%.4f\n",
            mean(d21$eng_prop_year, na.rm = TRUE),
            mean(d21$v2svdomaut,    na.rm = TRUE),
            mean(d21$v2svstterr,    na.rm = TRUE)))

# Full-sample event counts (diagnostic only; per-model EPV uses CC-mask counts)
cat(sprintf("NAVCO 2.1 events (full sample): concession=%d, success=%d (N=%d CY)\n",
            sum(d21$event_concession,  na.rm = TRUE),
            sum(d21$event_success_yr,  na.rm = TRUE),
            nrow(d21)))

###############################################################################
# FIT HELPERS (one per [dataset, estimator] combo)
###############################################################################
# Each helper returns the fitted model for a given (dv_name, moderator_name)
# pair. Variable names are injected via sprintf into the formula string, which
# is then parsed by as.formula — this avoids fragile non-standard evaluation.

fit_weibull13 <- function(dv_event, moderator) {
  f <- as.formula(sprintf(
    "Surv(duration, %s) ~ eng_c * %s + nonviolent_camp + REGCHANGE + v2csreprss + lnnum_image_sum",
    dv_event, moderator))
  survreg(f, data = df, dist = "weibull")
}

fit_firth13 <- function(dv_col, moderator) {
  f <- as.formula(sprintf(
    "%s ~ eng_c * %s + nonviolent_camp + REGCHANGE + v2csreprss + lnlengthofcam + lnnum_image_sum",
    dv_col, moderator))
  logistf(f, data = df,
          control   = logistf.control(maxit = 1000, maxstep = 0.5),
          plcontrol = logistpl.control(maxit = 1000),
          pl        = TRUE)
}

fit_bracl13 <- function(dv_col, moderator) {
  f <- as.formula(sprintf(
    "%s ~ eng_c * %s + nonviolent_camp + REGCHANGE + v2csreprss + lnlengthofcam + lnnum_image_sum",
    dv_col, moderator))
  bracl(f, data = df,
        type = "MPL_Jeffreys", link = "logit", parallel = TRUE, maxit = 500)
}

fit_cox21 <- function(dv_event, moderator, data = d21) {
  f <- as.formula(sprintf(
    "Surv(tstart, tstop, %s) ~ eng_c * %s + nonviolent + REGCHANGE + v2csreprss + lnnum_image_sum + cluster(navco21_id)",
    dv_event, moderator))
  coxph(f, data = data)
}

fit_firth21 <- function(dv_col, moderator) {
  f <- as.formula(sprintf(
    "%s ~ eng_c * %s + nonviolent + REGCHANGE + v2csreprss + time_in_campaign + lnnum_image_sum",
    dv_col, moderator))
  logistf(f, data = d21,
          control   = logistf.control(maxit = 1000, maxstep = 0.5),
          plcontrol = logistpl.control(maxit = 1000),
          pl        = TRUE)
}

fit_gee21 <- function(dv_col, moderator) {
  vars <- c(
    dv_col, "eng_c", moderator, "nonviolent", "REGCHANGE",
    "v2csreprss", "time_in_campaign", "lnnum_image_sum", "navco21_id"
  )
  dd <- d21[complete.cases(d21[, vars]), vars, drop = FALSE]
  dd <- dd %>% arrange(navco21_id, time_in_campaign)
  f <- as.formula(sprintf(
    "%s ~ eng_c * %s + nonviolent + REGCHANGE + v2csreprss + time_in_campaign + lnnum_image_sum",
    dv_col, moderator))
  geeglm(f, data = dd, id = navco21_id, family = binomial, corstr = "exchangeable")
}

fit_bracl21 <- function(dv_col, moderator) {
  f <- as.formula(sprintf(
    "%s ~ eng_c * %s + nonviolent + REGCHANGE + v2csreprss + time_in_campaign + lnnum_image_sum",
    dv_col, moderator))
  bracl(f, data = d21,
        type = "MPL_Jeffreys", link = "logit", parallel = TRUE, maxit = 500)
}

###############################################################################
# FIT ALL 20 MODELS
###############################################################################

# NAVCO 1.3 — domaut block
model1   <- fit_weibull13("event_concession", "domaut_c")
model1s  <- fit_weibull13("event_success",    "domaut_c")
model3   <- fit_firth13  ("concession_bin",   "domaut_c")
model3s  <- fit_firth13  ("success_bin",      "domaut_c")
model5   <- fit_bracl13  ("success3",         "domaut_c")

# NAVCO 1.3 — stterr block
model2   <- fit_weibull13("event_concession", "stterr_c")
model2s  <- fit_weibull13("event_success",    "stterr_c")
model4   <- fit_firth13  ("concession_bin",   "stterr_c")
model4s  <- fit_firth13  ("success_bin",      "stterr_c")
model6   <- fit_bracl13  ("success3",         "stterr_c")

# NAVCO 2.1 — domaut block
model7   <- fit_cox21    ("event_concession", "domaut_c", data = d21_first_conc)
model7s  <- fit_cox21    ("event_success_yr", "domaut_c")
model9   <- fit_firth21  ("concession_bin",   "domaut_c")
model9s  <- fit_firth21  ("success_bin_yr",   "domaut_c")
model11  <- fit_bracl21  ("progress_ord",     "domaut_c")
model9_gee <- fit_gee21  ("concession_bin",   "domaut_c")

# NAVCO 2.1 — stterr block
model8   <- fit_cox21    ("event_concession", "stterr_c", data = d21_first_conc)
model8s  <- fit_cox21    ("event_success_yr", "stterr_c")
model10  <- fit_firth21  ("concession_bin",   "stterr_c")
model10s <- fit_firth21  ("success_bin_yr",   "stterr_c")
model12  <- fit_bracl21  ("progress_ord",     "stterr_c")
model10_gee <- fit_gee21 ("concession_bin",   "stterr_c")

model7_ag <- fit_cox21   ("event_concession", "domaut_c", data = d21)
model8_ag <- fit_cox21   ("event_concession", "stterr_c", data = d21)

###############################################################################
# COX PH DIAGNOSTICS
###############################################################################

cox_ph_diag <- function(tag, model, focal_name) {
  z <- tryCatch(cox.zph(model), error = function(e) NULL)
  if (is.null(z)) {
    cat(sprintf("[cox.zph] %-6s : test failed (skipped)\n", tag))
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
    "PH_FAIL_INTERACTION"
  } else if (is.finite(global_p) && global_p < 0.05) {
    "PH_FAIL_GLOBAL"
  } else if (!is.finite(focal_p)) {
    "PH_FOCAL_MISSING"
  } else {
    ""
  }

  cat(sprintf("[cox.zph] %-6s | global p = %s | %s p = %s\n",
              tag,
              if (is.finite(global_p)) formatC(global_p, format = "f", digits = 4) else "NA",
              focal_name,
              if (is.finite(focal_p)) formatC(focal_p, format = "f", digits = 4) else "NA"))
  data.frame(
    model = tag,
    ph_global_p = global_p,
    ph_int_p = focal_p,
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
  cox_ph_diag("M7",   model7,  "eng_c:domaut_c"),
  cox_ph_diag("M7s",  model7s, "eng_c:domaut_c"),
  cox_ph_diag("M8",   model8,  "eng_c:stterr_c"),
  cox_ph_diag("M8s",  model8s, "eng_c:stterr_c"),
  cox_recurrent_ph_placeholder("M7_AG"),
  cox_recurrent_ph_placeholder("M8_AG")
)
cat("--- end Cox PH diagnostic ---\n\n")

###############################################################################
# EXTRACTION HELPER — interaction headline per model
###############################################################################

extract_h21 <- function(tag, estimator, data_label, response, model,
                        moderator_name, N, n_events) {
  eng_name <- "eng_c"
  mod_name <- moderator_name
  int_name <- paste0("eng_c:", moderator_name)

  cf    <- coef(model)
  b_eng <- unname(cf[eng_name])
  b_mod <- unname(cf[mod_name])
  b_int <- unname(cf[int_name])

 # Robust vcov block extraction by positional indexing
  get_vcov_block <- function(mdl, nms) {
    all_nms <- names(coef(mdl))
    ii <- match(nms, all_nms)
    if (any(is.na(ii))) return(NULL)
    V_full <- tryCatch(vcov(mdl), error = function(e) NULL)
    if (is.null(V_full) && inherits(mdl, "logistf")) V_full <- mdl$var
    if (is.null(V_full)) return(NULL)
    V_full[ii, ii, drop = FALSE]
  }

 # Per-estimator p-value extraction
  if (inherits(model, "survreg")) {
    tbl   <- summary(model)$table
    p_eng <- tbl[eng_name, "p"]
    p_mod <- tbl[mod_name, "p"]
    p_int <- tbl[int_name, "p"]
  } else if (inherits(model, "coxph")) {
    tbl   <- summary(model)$coefficients
 # Cox PH with cluster reports robust SE and its Pr(>|z|)
    p_eng <- tbl[eng_name, "Pr(>|z|)"]
    p_mod <- tbl[mod_name, "Pr(>|z|)"]
    p_int <- tbl[int_name, "Pr(>|z|)"]
  } else if (inherits(model, "logistf")) {
    p_eng <- unname(model$prob[eng_name])
    p_mod <- unname(model$prob[mod_name])
    p_int <- unname(model$prob[int_name])
  } else if (inherits(model, "geeglm")) {
    tbl   <- summary(model)$coefficients
    p_eng <- tbl[eng_name, "Pr(>|W|)"]
    p_mod <- tbl[mod_name, "Pr(>|W|)"]
    p_int <- tbl[int_name, "Pr(>|W|)"]
  } else if (inherits(model, "bracl")) {
    tbl   <- summary(model)$coefficients
    p_eng <- tbl[eng_name, "Pr(>|z|)"]
    p_mod <- tbl[mod_name, "Pr(>|z|)"]
    p_int <- tbl[int_name, "Pr(>|z|)"]
  } else {
    p_eng <- NA_real_; p_mod <- NA_real_; p_int <- NA_real_
  }

 # Joint Wald test on (eng, int): H0: β_eng = β_int = 0 (df=2)
 # Interpretation: "English has no effect on the outcome, anywhere in the
 # moderator range". The main effect of the moderator is not part of the
 # H2.1 test (it's a confounder/control for the interaction).
  joint_wald_p <- NA_real_
  V2 <- get_vcov_block(model, c(eng_name, int_name))
  if (!is.null(V2) && all(is.finite(c(b_eng, b_int))) && all(is.finite(V2))) {
    bvec <- c(b_eng, b_int)
    W    <- tryCatch(as.numeric(t(bvec) %*% solve(V2) %*% bvec),
                     error = function(e) NA_real_)
    if (is.finite(W)) joint_wald_p <- pchisq(W, df = 2, lower.tail = FALSE)
  }

 # Expected interaction sign under H2.1
 # (see header for derivation and bracl sign convention)
  expected_int_sign <- switch(estimator,
    "Weibull AFT" = "+",
    "Firth PML"   = "-",
    "GEE logit"   = "-",
    "Cox PH"      = "-",
    "bracl"       = "+",
    NA_character_
  )
  sign_consistent <- is.finite(b_int) && !is.na(expected_int_sign) && (
    (expected_int_sign == "+" && b_int > 0) ||
    (expected_int_sign == "-" && b_int < 0)
  )

  k_params <- length(cf)
  epv      <- if (is.na(n_events)) NA_real_ else n_events / k_params

  data.frame(
    model = tag, estimator = estimator, data = data_label,
    response = response, moderator = moderator_name,
    beta_eng_c   = b_eng, p_eng_c = p_eng,
    beta_mod_c   = b_mod, p_mod_c = p_mod,
    beta_int     = b_int, p_int   = p_int,
    joint_wald_p = joint_wald_p,
    expected_int_sign = expected_int_sign,
    sign_consistent   = sign_consistent,
    N = N, n_events = n_events, epv = epv,
    sig_05_int = as.logical(!is.na(p_int) && p_int < 0.05),
    sig_10_int = as.logical(!is.na(p_int) && p_int < 0.10),
    stringsAsFactors = FALSE
  )
}

###############################################################################
# PER-MODEL SAMPLE SIZES (complete-case, explicit to avoid class-specific
# nobs quirks for bracl)
#
# cc_stats returns BOTH N and the n_events count WITHIN the complete-case
# mask — necessary because missingness differs across specs/outcomes, so
# full-sample event counts would bias Peduzzi's EPV.
# For ordered outcomes (bracl), event_var is NA so n_events returns NA.
###############################################################################

cc_stats <- function(dat, vars, event_var = NA_character_) {
  cc_mask <- complete.cases(dat[, vars, drop = FALSE])
  N       <- sum(cc_mask)
  if (is.na(event_var) || !(event_var %in% names(dat))) {
    n_events <- NA_integer_
  } else {
    ev_vec   <- suppressWarnings(as.numeric(dat[[event_var]]))
    n_events <- as.integer(sum(ev_vec[cc_mask] == 1, na.rm = TRUE))
  }
  list(N = N, n_events = n_events)
}

# NAVCO 1.3 — survival covariates (no lnlengthofcam)
vars_13_surv <- function(dv, moderator) {
  c("duration", dv, "eng_c", moderator, "nonviolent_camp",
    "REGCHANGE", "v2csreprss", "lnnum_image_sum")
}
# NAVCO 1.3 — binary/ordered covariates (+ lnlengthofcam)
vars_13_bin  <- function(dv, moderator) {
  c(dv, "eng_c", moderator, "nonviolent_camp",
    "REGCHANGE", "v2csreprss", "lnlengthofcam", "lnnum_image_sum")
}
# NAVCO 2.1 — survival covariates (no time_in_campaign)
vars_21_surv <- function(dv, moderator) {
  c("tstart", "tstop", dv, "eng_c", moderator, "nonviolent",
    "REGCHANGE", "v2csreprss", "lnnum_image_sum", "navco21_id")
}
# NAVCO 2.1 — binary/ordered covariates (+ time_in_campaign)
vars_21_bin  <- function(dv, moderator) {
  c(dv, "eng_c", moderator, "nonviolent",
    "REGCHANGE", "v2csreprss", "time_in_campaign", "lnnum_image_sum")
}

s_m1   <- cc_stats(df,  vars_13_surv("event_concession", "domaut_c"), "event_concession")
s_m1s  <- cc_stats(df,  vars_13_surv("event_success",    "domaut_c"), "event_success")
s_m2   <- cc_stats(df,  vars_13_surv("event_concession", "stterr_c"), "event_concession")
s_m2s  <- cc_stats(df,  vars_13_surv("event_success",    "stterr_c"), "event_success")
s_m3   <- cc_stats(df,  vars_13_bin ("concession_bin",   "domaut_c"), "concession_bin")
s_m3s  <- cc_stats(df,  vars_13_bin ("success_bin",      "domaut_c"), "success_bin")
s_m4   <- cc_stats(df,  vars_13_bin ("concession_bin",   "stterr_c"), "concession_bin")
s_m4s  <- cc_stats(df,  vars_13_bin ("success_bin",      "stterr_c"), "success_bin")
s_m5   <- cc_stats(df,  vars_13_bin ("success3",         "domaut_c"))                      # ordered
s_m6   <- cc_stats(df,  vars_13_bin ("success3",         "stterr_c"))                      # ordered

s_m7   <- cc_stats(d21_first_conc, vars_21_surv("event_concession", "domaut_c"), "event_concession")
s_m7s  <- cc_stats(d21, vars_21_surv("event_success_yr", "domaut_c"), "event_success_yr")
s_m8   <- cc_stats(d21_first_conc, vars_21_surv("event_concession", "stterr_c"), "event_concession")
s_m8s  <- cc_stats(d21, vars_21_surv("event_success_yr", "stterr_c"), "event_success_yr")
s_m7_ag <- cc_stats(d21, vars_21_surv("event_concession", "domaut_c"), "event_concession")
s_m8_ag <- cc_stats(d21, vars_21_surv("event_concession", "stterr_c"), "event_concession")
s_m9   <- cc_stats(d21, vars_21_bin ("concession_bin",   "domaut_c"), "concession_bin")
s_m9_gee <- cc_stats(d21, vars_21_bin ("concession_bin",   "domaut_c"), "concession_bin")
s_m9s  <- cc_stats(d21, vars_21_bin ("success_bin_yr",   "domaut_c"), "success_bin_yr")
s_m10  <- cc_stats(d21, vars_21_bin ("concession_bin",   "stterr_c"), "concession_bin")
s_m10_gee <- cc_stats(d21, vars_21_bin ("concession_bin",   "stterr_c"), "concession_bin")
s_m10s <- cc_stats(d21, vars_21_bin ("success_bin_yr",   "stterr_c"), "success_bin_yr")
s_m11  <- cc_stats(d21, vars_21_bin ("progress_ord",     "domaut_c"))                      # ordered
s_m12  <- cc_stats(d21, vars_21_bin ("progress_ord",     "stterr_c"))                      # ordered

###############################################################################
# BUILD RESULTS TABLE (20 rows)
###############################################################################

core_results <- rbind(
 # NAVCO 1.3 — campaign-level
  extract_h21("M1",   "Weibull AFT", "NAVCO 1.3",
              "Time to concession",               model1,   "domaut_c", s_m1$N,   s_m1$n_events),
  extract_h21("M1s",  "Weibull AFT", "NAVCO 1.3",
              "Time to success (strict)",         model1s,  "domaut_c", s_m1s$N,  s_m1s$n_events),
  extract_h21("M2",   "Weibull AFT", "NAVCO 1.3",
              "Time to concession",               model2,   "stterr_c", s_m2$N,   s_m2$n_events),
  extract_h21("M2s",  "Weibull AFT", "NAVCO 1.3",
              "Time to success (strict)",         model2s,  "stterr_c", s_m2s$N,  s_m2s$n_events),
  extract_h21("M3",   "Firth PML",   "NAVCO 1.3",
              "Concession (binary)",              model3,   "domaut_c", s_m3$N,   s_m3$n_events),
  extract_h21("M3s",  "Firth PML",   "NAVCO 1.3",
              "Success (strict, binary)",         model3s,  "domaut_c", s_m3s$N,  s_m3s$n_events),
  extract_h21("M4",   "Firth PML",   "NAVCO 1.3",
              "Concession (binary)",              model4,   "stterr_c", s_m4$N,   s_m4$n_events),
  extract_h21("M4s",  "Firth PML",   "NAVCO 1.3",
              "Success (strict, binary)",         model4s,  "stterr_c", s_m4s$N,  s_m4s$n_events),
  extract_h21("M5",   "bracl",       "NAVCO 1.3",
              "success3 (0/1/2) ordered",         model5,   "domaut_c", s_m5$N,   s_m5$n_events),
  extract_h21("M6",   "bracl",       "NAVCO 1.3",
              "success3 (0/1/2) ordered",         model6,   "stterr_c", s_m6$N,   s_m6$n_events),
 # NAVCO 2.1 — campaign-year panel
  extract_h21("M7",   "Cox PH",      "NAVCO 2.1",
              "Hazard of first concession (CY)",  model7,   "domaut_c", s_m7$N,   s_m7$n_events),
  extract_h21("M7s",  "Cox PH",      "NAVCO 2.1",
              "Hazard of success (CY, strict)",   model7s,  "domaut_c", s_m7s$N,  s_m7s$n_events),
  extract_h21("M8",   "Cox PH",      "NAVCO 2.1",
              "Hazard of first concession (CY)",  model8,   "stterr_c", s_m8$N,   s_m8$n_events),
  extract_h21("M8s",  "Cox PH",      "NAVCO 2.1",
              "Hazard of success (CY, strict)",   model8s,  "stterr_c", s_m8s$N,  s_m8s$n_events),
  extract_h21("M9",   "Firth PML",   "NAVCO 2.1",
              "Yearly concession (CY, binary)",   model9,   "domaut_c", s_m9$N,   s_m9$n_events),
  extract_h21("M9s",  "Firth PML",   "NAVCO 2.1",
              "Yearly success (CY, strict bin.)", model9s,  "domaut_c", s_m9s$N,  s_m9s$n_events),
  extract_h21("M10",  "Firth PML",   "NAVCO 2.1",
              "Yearly concession (CY, binary)",   model10,  "stterr_c", s_m10$N,  s_m10$n_events),
  extract_h21("M10s", "Firth PML",   "NAVCO 2.1",
              "Yearly success (CY, strict bin.)", model10s, "stterr_c", s_m10s$N, s_m10s$n_events),
  extract_h21("M11",  "bracl",       "NAVCO 2.1",
              "progress_ord (6-level; failure separate)",   model11,  "domaut_c", s_m11$N,  s_m11$n_events),
  extract_h21("M12",  "bracl",       "NAVCO 2.1",
              "progress_ord (6-level; failure separate)",   model12,  "stterr_c", s_m12$N,  s_m12$n_events)
)

###############################################################################
# MAIN-PANEL SCOPE + MULTIPLICITY CORRECTION
# ---------------------------------------------------------------------------
# M9 and M10 (Firth PML, NAVCO 2.1, yearly concession_bin × domaut_c /
#   stterr_c) are in the main panel so each dataset
#   contributes BOTH a survival estimator AND a binary estimator, removing
#   the asymmetric structure of the prior 6-model panel (NAVCO 1.3 had
#   Weibull AFT + Firth PML; NAVCO 2.1 had Cox PH only). M9/M10 are
#   substantively distinct from M3/M4 — they test concession on a
#   yearly basis (campaign-year unit) rather than ever-during-campaign,
#   so they complement rather than duplicate the NAVCO 1.3 binary models.
#
# Main panel (Table 2, 8 models): M1, M2, M3, M4, M7, M8, M9, M10
#   — concession-only, V-Dem moderators (domaut_c, stterr_c) as-is
# Appendix (12 models): strict-success twins (M1s, M2s, M3s, M4s, M7s, M8s,
#   M9s, M10s), ordered bracl (M5, M6, M11, M12)
# Robustness rows (separate, no Holm/BH labels): M7_AG, M8_AG, M9_GEE, M10_GEE
#
# Multiplicity:
#   Holm within the 8 main-panel p_int (confirmatory)
#   Benjamini-Hochberg within the 12 appendix p_int (exploratory)
#   NO pooling across panels.
#
# Substantive interpretation:
#   Table 2 reports disconfirmation; wrong-direction interactions
#   across the full 20-model stack are themselves substantively meaningful
#   (consistent with Risse-Ropp-Sikkink 1999 spiral model / Hafner-Burton
#   2008 selection effects; see manuscript Results §6 / Discussion §7).
#
# EPV warning:
#   Models with epv < 5 (Peduzzi et al. 1996 threshold) flagged.
#   M7/M8 NAVCO 2.1 Cox concession at EPV = 1.86 carry low-power warning.
#   M9/M10 NAVCO 2.1 Firth PML at EPV = 1.67 are also flagged; these low-power
#   models remain in the main panel for estimator symmetry.
###############################################################################

main_models <- c("M1", "M2", "M3", "M4", "M7", "M8", "M9", "M10")
core_results$main_panel <- core_results$model %in% main_models

# EPV warning (< 5 per Peduzzi et al. 1996)
core_results$epv_warning <- !is.na(core_results$epv) & core_results$epv < 5

# Multiplicity: Holm within the main panel (k=8), BH within the appendix (k=12).
# Separate panels — no pooling across.
core_results$p_int_holm <- NA_real_
core_results$p_int_bh   <- NA_real_
core_results$joint_wald_holm <- NA_real_
core_results$joint_wald_bh   <- NA_real_

main_idx <- which(core_results$main_panel)
app_idx  <- which(!core_results$main_panel)

if (length(main_idx) > 0) {
# Explicit n guarantees k=8 even if a p-value is NA; base p.adjust otherwise
# drops NAs and resets n to the number of nonmissing values.
  core_results$p_int_holm[main_idx] <- p.adjust(
    core_results$p_int[main_idx], method = "holm", n = length(main_idx)
  )
  core_results$joint_wald_holm[main_idx] <- p.adjust(
    core_results$joint_wald_p[main_idx], method = "holm", n = length(main_idx)
  )
}
if (length(app_idx) > 0) {
  core_results$p_int_bh[app_idx] <- p.adjust(
    core_results$p_int[app_idx], method = "BH", n = length(app_idx)
  )
  core_results$joint_wald_bh[app_idx] <- p.adjust(
    core_results$joint_wald_p[app_idx], method = "BH", n = length(app_idx)
  )
}

# Significance on adjusted p-values
core_results$sig_05_int_holm <- !is.na(core_results$p_int_holm) & core_results$p_int_holm < 0.05
core_results$sig_10_int_holm <- !is.na(core_results$p_int_holm) & core_results$p_int_holm < 0.10
core_results$sig_05_int_bh   <- !is.na(core_results$p_int_bh)   & core_results$p_int_bh   < 0.05
core_results$sig_10_int_bh   <- !is.na(core_results$p_int_bh)   & core_results$p_int_bh   < 0.10
core_results$robustness_row  <- FALSE

robustness_rows <- rbind(
  extract_h21("M7_AG", "Cox PH", "NAVCO 2.1",
              "Hazard of concession (CY, recurrent-event)", model7_ag, "domaut_c", s_m7_ag$N, s_m7_ag$n_events),
  extract_h21("M8_AG", "Cox PH", "NAVCO 2.1",
              "Hazard of concession (CY, recurrent-event)", model8_ag, "stterr_c", s_m8_ag$N, s_m8_ag$n_events),
  extract_h21("M9_GEE", "GEE logit", "NAVCO 2.1",
              "Yearly concession (CY, source-clustered GEE)", model9_gee, "domaut_c", s_m9_gee$N, s_m9_gee$n_events),
  extract_h21("M10_GEE", "GEE logit", "NAVCO 2.1",
              "Yearly concession (CY, source-clustered GEE)", model10_gee, "stterr_c", s_m10_gee$N, s_m10_gee$n_events)
)
robustness_rows$main_panel      <- FALSE
robustness_rows$p_int_holm      <- NA_real_
robustness_rows$p_int_bh        <- NA_real_
robustness_rows$joint_wald_holm <- NA_real_
robustness_rows$joint_wald_bh   <- NA_real_
robustness_rows$sig_05_int_holm <- NA
robustness_rows$sig_10_int_holm <- NA
robustness_rows$sig_05_int_bh   <- NA
robustness_rows$sig_10_int_bh   <- NA
robustness_rows$epv_warning     <- !is.na(robustness_rows$epv) & robustness_rows$epv < 5
robustness_rows$robustness_row  <- TRUE

results <- rbind(core_results, robustness_rows[, names(core_results)])
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
      ph_flag == "PH_FAIL_INTERACTION" ~
        "PH interaction failure; interpret cautiously",
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
# COMPACT PRINT — 20 CORE MODELS + 4 ROBUSTNESS ROWS
###############################################################################

fmt_p <- function(p) ifelse(is.na(p), "   NA ", sprintf("%6.3f", p))
fmt_b <- function(b) ifelse(is.na(b), "    NA ", sprintf("%7.2f", b))

main_rows <- core_results[core_results$main_panel, ]
app_rows  <- core_results[!core_results$main_panel, ]
k_main <- nrow(main_rows)   # k=8
k_app  <- nrow(app_rows)    # k=12

cat(sprintf("\n=========== H2.1 MAIN PANEL (Table 2, %d models, Holm-adjusted) ===========\n", k_main))
cat(sprintf("%-5s %-12s %-10s %-9s %8s %7s %8s %9s %4s %4s %5s %4s %3s\n",
            "Model", "Estimator", "Data", "Moderator",
            "b_int", "p_int", "p_Holm", "p_JW_H",
            "Exp", "Sign", "N", "EPV", "OK?"))
cat(paste(rep("-", 110), collapse = ""), "\n")
for (i in seq_len(nrow(main_rows))) {
  r <- main_rows[i, ]
  cat(sprintf("%-5s %-12s %-10s %-9s %8s %7s %8s %9s %4s %4s %5d %4s %3s\n",
              r$model, r$estimator, r$data, r$moderator,
              fmt_b(r$beta_int),   fmt_p(r$p_int),
              fmt_p(r$p_int_holm),
              fmt_p(r$joint_wald_holm),
              r$expected_int_sign,
              ifelse(isTRUE(r$sign_consistent), "cons", "WRONG"),
              r$N,
              ifelse(is.na(r$epv), "  NA", sprintf("%4.1f", r$epv)),
              ifelse(r$epv_warning, "!", " ")))
}

cat(sprintf("\n=========== H2.1 APPENDIX (%d models, BH-adjusted) ==========================\n", k_app))
cat(sprintf("%-5s %-12s %-10s %-9s %8s %7s %8s %9s %4s %4s %5s %4s %3s\n",
            "Model", "Estimator", "Data", "Moderator",
            "b_int", "p_int", "p_BH", "p_JW_BH",
            "Exp", "Sign", "N", "EPV", "OK?"))
cat(paste(rep("-", 110), collapse = ""), "\n")
for (i in seq_len(nrow(app_rows))) {
  r <- app_rows[i, ]
  cat(sprintf("%-5s %-12s %-10s %-9s %8s %7s %8s %9s %4s %4s %5d %4s %3s\n",
              r$model, r$estimator, r$data, r$moderator,
              fmt_b(r$beta_int),   fmt_p(r$p_int),
              fmt_p(r$p_int_bh),
              fmt_p(r$joint_wald_bh),
              r$expected_int_sign,
              ifelse(isTRUE(r$sign_consistent), "cons", "WRONG"),
              r$N,
              ifelse(is.na(r$epv), "  NA", sprintf("%4.1f", r$epv)),
              ifelse(r$epv_warning, "!", " ")))
}

cat("\nNotes:\n")
cat(" b_int/p_int = coefficient and p-value on eng_c × moderator_c\n")
cat(sprintf("  p_Holm      = Holm-Bonferroni FWER-adjusted (main panel, k=%d; confirmatory)\n", k_main))
cat(sprintf("  p_BH        = Benjamini-Hochberg FDR-adjusted (appendix, k=%d; exploratory)\n", k_app))
cat(" p_JW_H = Holm-adjusted joint Wald p on (eng_c, eng_c:mod_c) in the main panel\n")
cat(" p_JW_BH = BH-adjusted joint Wald p on (eng_c, eng_c:mod_c) in the appendix\n")
cat(" Exp = expected sign of interaction under H2.1 (see header)\n")
cat(" Sign = cons/WRONG vs. expected sign\n")
cat(" EPV = events per variable; '!' marks EPV<5 (Peduzzi et al. 1996)\n")
cat(" M9_GEE/M10_GEE = exchangeable GEE companions for yearly concession binaries;\n")
cat(" use as the source-campaign clustered inference check for repeated source-years.\n")
cat(" M9s/M10s/M11/M12 remain appendix-only pooled yearly models and should be read as weaker evidence.\n")
cat(sprintf("\nSummary (main panel, %d models, Holm-adjusted):\n", k_main))
cat(sprintf("  sign-consistent rows: %d/%d\n", sum(main_rows$sign_consistent, na.rm = TRUE), k_main))
cat(sprintf("  p_Holm < 0.05       : %d/%d\n", sum(main_rows$sig_05_int_holm, na.rm = TRUE), k_main))
cat(sprintf("  p_Holm < 0.10       : %d/%d\n", sum(main_rows$sig_10_int_holm, na.rm = TRUE), k_main))
cat(sprintf("  EPV warning (<5)    : %d/%d\n", sum(main_rows$epv_warning, na.rm = TRUE), k_main))
cat(sprintf("\nSummary (appendix, %d models, BH-adjusted):\n", k_app))
cat(sprintf("  sign-consistent rows: %d/%d\n", sum(app_rows$sign_consistent, na.rm = TRUE), k_app))
cat(sprintf("  p_BH < 0.05         : %d/%d\n", sum(app_rows$sig_05_int_bh, na.rm = TRUE), k_app))
cat(sprintf("  p_BH < 0.10         : %d/%d\n", sum(app_rows$sig_10_int_bh, na.rm = TRUE), k_app))
cat(sprintf("  EPV warning (<5)    : %d/%d\n", sum(app_rows$epv_warning, na.rm = TRUE), k_app))

robust_rows <- results[results$robustness_row, ]
if (nrow(robust_rows) > 0) {
 cat("\n=========== H2.1 RECURRENT-EVENT ROBUSTNESS (outside multiplicity family) ===\n")
  cat(sprintf("%-6s %-10s %-9s %8s %7s %9s %4s %5s\n",
              "Model", "Data", "Moderator", "b_int", "p_int", "p_JW_raw", "N", "ev"))
  cat(paste(rep("-", 78), collapse = ""), "\n")
  for (i in seq_len(nrow(robust_rows))) {
    r <- robust_rows[i, ]
    cat(sprintf("%-6s %-10s %-9s %8s %7s %9s %4d %5d\n",
                r$model, r$data, r$moderator,
                fmt_b(r$beta_int), fmt_p(r$p_int), fmt_p(r$joint_wald_p),
                r$N, r$n_events))
  }
}

###############################################################################
# FULL MODEL SUMMARIES (for appendix / diagnostics)
###############################################################################

cat("\n\n==== FULL MODEL SUMMARIES ====\n")

cat("\n---- M1 Weibull AFT × domaut_c — time to concession ----\n"); print(summary(model1))
cat("\n---- M1s Weibull AFT × domaut_c — time to STRICT success ----\n"); print(summary(model1s))
cat("\n---- M2 Weibull AFT × stterr_c — time to concession ----\n"); print(summary(model2))
cat("\n---- M2s Weibull AFT × stterr_c — time to STRICT success ----\n"); print(summary(model2s))
cat("\n---- M3 Firth PML × domaut_c — concession (binary) ----\n"); print(summary(model3))
cat("\n---- M3s Firth PML × domaut_c — STRICT success (binary) ----\n"); print(summary(model3s))
cat("\n---- M4 Firth PML × stterr_c — concession (binary) ----\n"); print(summary(model4))
cat("\n---- M4s Firth PML × stterr_c — STRICT success (binary) ----\n"); print(summary(model4s))
cat("\n---- M5 bracl × domaut_c — success3 (ordered) ----\n"); print(summary(model5))
cat("\n---- M6 bracl × stterr_c — success3 (ordered) ----\n"); print(summary(model6))
cat("\n---- M7 Cox PH × domaut_c — hazard of FIRST concession (CY) ----\n"); print(summary(model7))
cat("\n---- M7s Cox PH × domaut_c — hazard of STRICT success (CY) ----\n"); print(summary(model7s))
cat("\n---- M8 Cox PH × stterr_c — hazard of FIRST concession (CY) ----\n"); print(summary(model8))
cat("\n---- M8s Cox PH × stterr_c — hazard of STRICT success (CY) ----\n"); print(summary(model8s))
cat("\n---- M7_AG Cox PH × domaut_c — recurrent-event concession hazard (CY) ----\n"); print(summary(model7_ag))
cat("\n---- M8_AG Cox PH × stterr_c — recurrent-event concession hazard (CY) ----\n"); print(summary(model8_ag))
cat("\n---- M9_GEE GEE logit × domaut_c — yearly concession, source-clustered ----\n"); print(summary(model9_gee))
cat("\n---- M10_GEE GEE logit × stterr_c — yearly concession, source-clustered ----\n"); print(summary(model10_gee))
cat("\n---- M9 Firth PML × domaut_c — yearly concession ----\n"); print(summary(model9))
cat("\n---- M9s Firth PML × domaut_c — yearly STRICT success ----\n"); print(summary(model9s))
cat("\n---- M10 Firth PML × stterr_c — yearly concession ----\n"); print(summary(model10))
cat("\n---- M10s Firth PML × stterr_c — yearly STRICT success ----\n"); print(summary(model10s))
cat("\n---- M11 bracl × domaut_c — progress_ord (6-level; failure separate) ----\n"); print(summary(model11))
cat("\n---- M12 bracl × stterr_c — progress_ord (6-level; failure separate) ----\n"); print(summary(model12))

expected_result_cols <- c(
  "model", "estimator", "data", "response", "moderator", "beta_eng_c",
  "p_eng_c", "beta_mod_c", "p_mod_c", "beta_int", "p_int", "joint_wald_p",
  "expected_int_sign", "sign_consistent", "N", "n_events", "epv", "main_panel",
  "p_int_holm", "p_int_bh", "joint_wald_holm", "joint_wald_bh", "robustness_row",
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
if (nrow(results) != 24) {
  stop(sprintf("Expected 24 H2.1 result rows, found %d.", nrow(results)), call. = FALSE)
}
if (anyDuplicated(results$model) > 0) {
 stop("Duplicate model labels detected in H2.1 results.", call. = FALSE)
}

###############################################################################
# EXPORT CSV
###############################################################################

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
write.csv(results, out_path, row.names = FALSE)
models <- list(
  M1 = model1, M1s = model1s, M2 = model2, M2s = model2s,
  M3 = model3, M3s = model3s, M4 = model4, M4s = model4s,
  M5 = model5, M6 = model6, M7 = model7, M7s = model7s,
  M8 = model8, M8s = model8s, M9 = model9, M9s = model9s,
  M10 = model10, M10s = model10s, M11 = model11, M12 = model12,
  M7_AG = model7_ag, M8_AG = model8_ag,
  M9_GEE = model9_gee, M10_GEE = model10_gee
)
stopifnot(identical(names(models), results$model))
saveRDS(models, cache_path, version = 3)
saveRDS(list(df = df, d21 = d21, d21_first_conc = d21_first_conc), frames_path, version = 3)
cat(sprintf("\nResults exported to: %s\n", out_path))
cat(sprintf("Model bundle exported to: %s\n", cache_path))
cat(sprintf("Analysis frames exported to: %s\n", frames_path))

cat("\n\n[H21_code.r complete]\n")
