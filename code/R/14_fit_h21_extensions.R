###############################################################################
# H21_curvilinear_appendix.r
# H2.1 APPENDIX: CURVATURE × POLITICAL VULNERABILITY MODERATION
#
# Extends H21_code.r by adding a quadratic English term I(eng_c^2) and its
# interaction with the political vulnerability moderators (v2svdomaut,
# v2svstterr). Tests whether political vulnerability moderates not only the
# LEVEL but also the CURVATURE of the English-share dose-response.
#
# Per-row specification:
#   Outcome ~ eng_c + I(eng_c^2) + moderator_c
#           + eng_c:moderator_c
#           + I(eng_c^2):moderator_c
#           + [parsimonious H21 controls]
#
# Two interaction coefficients of interest:
#   eng_c:mod_c        (linear shift; mirrors H21 main panel)
#   I(eng_c^2):mod_c   (CURVATURE shift; the novel test)
#
# Headline tests reported per row:
#   p_int_lin: single-coef p on eng_c:mod_c
#   p_int_sq: single-coef p on I(eng_c^2):mod_c   <-- NOVEL TEST
#   joint_wald_int_p: joint Wald on (eng:mod, eng^2:mod), df=2
#                       "does the moderator shift any aspect of the curve?"
#   joint_wald_eng_p: joint Wald on (eng, eng^2, eng:mod, eng^2:mod), df=4
#                       "does English have any effect anywhere in mod range?"
#
# MODELS (20 rows, all APPENDIX-only):
#   Concession family (BH-1, k=8):
#     M1q   Weibull AFT × domaut   (NAVCO 1.3) time to concession
#     M2q   Weibull AFT × stterr   (NAVCO 1.3)
#     M3q   Firth PML   × domaut   (NAVCO 1.3) concession_bin
#     M4q   Firth PML   × stterr   (NAVCO 1.3)
#     M7q   Cox PH      × domaut   (NAVCO 2.1) hazard of first concession
#     M8q   Cox PH      × stterr   (NAVCO 2.1)
#     M9q   Firth PML   × domaut   (NAVCO 2.1) yearly concession_bin
#     M10q  Firth PML   × stterr   (NAVCO 2.1)
#   Strict-success twins (BH-2, k=8): M1qs, M2qs, M3qs, M4qs, M7qs, M8qs, M9qs, M10qs
#   Ordered bracl family (BH-3, k=4):
#     M5q   bracl × domaut   (NAVCO 1.3) success3
#     M6q   bracl × stterr   (NAVCO 1.3)
#     M11q  bracl × domaut   (NAVCO 2.1) progress_ord
#     M12q  bracl × stterr   (NAVCO 2.1)
#
# Multiplicity: BH within each of the 3 sub-families, no pooling.
# Diagnostics (cox.zph etc.) intentionally omitted — appendix exploratory only.
#
# Author: Jie (Jason) Lian
# Created: 2026-05-08
###############################################################################

rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
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
out_path     <- file.path(P$estimates, "H21_curvilinear_appendix_results.csv")
limited_out_path <- file.path(P$estimates, "H21_limited_success_results.csv")
cache_path   <- file.path(P$cache, "H21_extensions_models.rds")
convergence_path <- file.path(P$diagnostics, "H21_extensions_logistf_convergence.csv")

assert_required_cols <- function(dat, required_cols, label) {
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop(
      sprintf("%s is missing required columns: %s", label, paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }
}

df  <- read.csv(navco13_path)
d21 <- read.csv(navco21_path)

assert_required_cols(
  df,
  c(
    "EYEAR", "BYEAR", "success", "limited", "success3",
    "nonviolent_camp", "REGCHANGE", "v2csreprss",
    "eng_prob_general", "lnnum_image_sum", "lnlengthofcam",
    "v2svdomaut", "v2svstterr"
  ),
  "df_final.csv"
)
assert_required_cols(
  d21,
  c(
    "navco21_id", "Year", "panel_campaign_id", "nonviolent",
    "camp_goals", "Num_Image", "progress", "concession", "time_in_campaign",
    "eng_prop_year", "v2svdomaut", "v2svstterr", "v2csreprss"
  ),
  "df_navco21_panel.csv"
)

###############################################################################
# NAVCO 1.3 PREPARATION (mirrors H21_code.r)
###############################################################################

df$nonviolent_camp <- as.factor(df$nonviolent_camp)   # polarity: 1=nonviolent
df$success3        <- ordered(df$success3)
df$REGCHANGE       <- as.numeric(as.character(df$REGCHANGE))

df$duration         <- df$EYEAR - df$BYEAR + 1
df$event_concession <- as.integer(
  as.numeric(as.character(df$success)) == 1 |
  as.numeric(as.character(df$limited)) == 1
)
df$event_success    <- as.numeric(as.character(df$success))
df$concession_bin   <- df$event_concession
df$success_bin      <- df$event_success

# Mean-center constituents; the quadratic is on the centered scale.
df$eng_c    <- df$eng_prob_general - mean(df$eng_prob_general, na.rm = TRUE)
df$domaut_c <- df$v2svdomaut       - mean(df$v2svdomaut,       na.rm = TRUE)
df$stterr_c <- df$v2svstterr       - mean(df$v2svstterr,       na.rm = TRUE)

cat(sprintf("NAVCO 1.3 centering: mean(eng)=%.4f, mean(domaut)=%.4f, mean(stterr)=%.4f\n",
            mean(df$eng_prob_general, na.rm = TRUE),
            mean(df$v2svdomaut,       na.rm = TRUE),
            mean(df$v2svstterr,       na.rm = TRUE)))

###############################################################################
# NAVCO 2.1 PREPARATION (mirrors H21_code.r)
###############################################################################

d21$nonviolent      <- as.factor(d21$nonviolent)
d21$REGCHANGE       <- as.numeric(d21$camp_goals == 0)
d21$lnnum_image_sum <- log(d21$Num_Image)
make_progress_ord6  <- function(x) ordered(x, levels = c(5, 0, 1, 2, 3, 4))
d21$progress_ord    <- make_progress_ord6(d21$progress)

d21$event_concession <- as.numeric(d21$concession)
d21$event_success_yr <- as.numeric(d21$progress == 4)
d21$concession_bin   <- d21$event_concession
d21$success_bin_yr   <- d21$event_success_yr

# Counting-process intervals for Cox PH.
d21$tstart <- d21$time_in_campaign - 1
d21$tstop  <- d21$time_in_campaign

# Mean-center interaction variables.
d21$eng_c    <- d21$eng_prop_year - mean(d21$eng_prop_year, na.rm = TRUE)
d21$domaut_c <- d21$v2svdomaut    - mean(d21$v2svdomaut,    na.rm = TRUE)
d21$stterr_c <- d21$v2svstterr    - mean(d21$v2svstterr,    na.rm = TRUE)

build_first_event_panel <- function(dat, id_col, event_col, order_col) {
  out <- dat %>%
    arrange(across(all_of(c(id_col, order_col)))) %>%
    group_by(across(all_of(id_col))) %>%
    mutate(
      .cum_event   = cumsum(.data[[event_col]]),
      .prior_event = lag(.cum_event, default = 0L),
      .keep_risk   = .prior_event == 0L,
      .event_first = as.integer(.data[[event_col]] == 1 & .cum_event == 1L)
    ) %>%
    filter(.keep_risk) %>%
    ungroup()
  out[[event_col]] <- out$.event_first
  out %>% dplyr::select(-.cum_event, -.prior_event, -.keep_risk, -.event_first)
}

d21_first_conc <- build_first_event_panel(
  d21, "navco21_id", "event_concession", "time_in_campaign"
)

cat(sprintf("NAVCO 2.1 centering: mean(eng)=%.4f, mean(domaut)=%.4f, mean(stterr)=%.4f\n",
            mean(d21$eng_prop_year, na.rm = TRUE),
            mean(d21$v2svdomaut,    na.rm = TRUE),
            mean(d21$v2svstterr,    na.rm = TRUE)))

###############################################################################
# FIT HELPERS — quadratic + interaction-with-curvature
###############################################################################
# Each helper builds a formula of the form:
#   <outcome> ~ eng_c + I(eng_c^2) + <mod> + eng_c:<mod> + I(eng_c^2):<mod>
#                + <controls>
# Variable names are injected via sprintf and parsed by as.formula — same
# pattern as H21_code.r to avoid fragile non-standard evaluation.

fit_weibull13 <- function(dv_event, moderator) {
  f <- as.formula(sprintf(
    paste0(
      "Surv(duration, %s) ~ eng_c + I(eng_c^2) + %s",
      " + eng_c:%s + I(eng_c^2):%s",
      " + nonviolent_camp + REGCHANGE + v2csreprss + lnnum_image_sum"
    ),
    dv_event, moderator, moderator, moderator
  ))
  survreg(f, data = df, dist = "weibull")
}

fit_firth13 <- function(dv_col, moderator) {
  f <- as.formula(sprintf(
    paste0(
      "%s ~ eng_c + I(eng_c^2) + %s",
      " + eng_c:%s + I(eng_c^2):%s",
      " + nonviolent_camp + REGCHANGE + v2csreprss + lnlengthofcam + lnnum_image_sum"
    ),
    dv_col, moderator, moderator, moderator
  ))
  logistf(f, data = df,
          control   = logistf.control(maxit = 1000, maxstep = 0.5),
          plcontrol = logistpl.control(maxit = 1000),
          pl        = TRUE)
}

fit_bracl13 <- function(dv_col, moderator) {
  f <- as.formula(sprintf(
    paste0(
      "%s ~ eng_c + I(eng_c^2) + %s",
      " + eng_c:%s + I(eng_c^2):%s",
      " + nonviolent_camp + REGCHANGE + v2csreprss + lnlengthofcam + lnnum_image_sum"
    ),
    dv_col, moderator, moderator, moderator
  ))
  bracl(f, data = df,
        type = "MPL_Jeffreys", link = "logit", parallel = TRUE, maxit = 500)
}

fit_cox21 <- function(dv_event, moderator, data = d21) {
  f <- as.formula(sprintf(
    paste0(
      "Surv(tstart, tstop, %s) ~ eng_c + I(eng_c^2) + %s",
      " + eng_c:%s + I(eng_c^2):%s",
      " + nonviolent + REGCHANGE + v2csreprss + lnnum_image_sum + cluster(navco21_id)"
    ),
    dv_event, moderator, moderator, moderator
  ))
  coxph(f, data = data)
}

fit_firth21 <- function(dv_col, moderator) {
  f <- as.formula(sprintf(
    paste0(
      "%s ~ eng_c + I(eng_c^2) + %s",
      " + eng_c:%s + I(eng_c^2):%s",
      " + nonviolent + REGCHANGE + v2csreprss + time_in_campaign + lnnum_image_sum"
    ),
    dv_col, moderator, moderator, moderator
  ))
  logistf(f, data = d21,
          control   = logistf.control(maxit = 1000, maxstep = 0.5),
          plcontrol = logistpl.control(maxit = 1000),
          pl        = TRUE)
}

fit_bracl21 <- function(dv_col, moderator) {
  f <- as.formula(sprintf(
    paste0(
      "%s ~ eng_c + I(eng_c^2) + %s",
      " + eng_c:%s + I(eng_c^2):%s",
      " + nonviolent + REGCHANGE + v2csreprss + time_in_campaign + lnnum_image_sum"
    ),
    dv_col, moderator, moderator, moderator
  ))
  bracl(f, data = d21,
        type = "MPL_Jeffreys", link = "logit", parallel = TRUE, maxit = 500)
}

###############################################################################
# FIT ALL 20 MODELS
###############################################################################

# NAVCO 1.3 — concession block
m1q   <- fit_weibull13("event_concession", "domaut_c")
m2q   <- fit_weibull13("event_concession", "stterr_c")
m3q   <- fit_firth13  ("concession_bin",   "domaut_c")
m4q   <- fit_firth13  ("concession_bin",   "stterr_c")
# NAVCO 1.3 — strict-success twins
m1qs  <- fit_weibull13("event_success",    "domaut_c")
m2qs  <- fit_weibull13("event_success",    "stterr_c")
m3qs  <- fit_firth13  ("success_bin",      "domaut_c")
m4qs  <- fit_firth13  ("success_bin",      "stterr_c")
# NAVCO 1.3 — ordered bracl
m5q   <- fit_bracl13  ("success3",         "domaut_c")
m6q   <- fit_bracl13  ("success3",         "stterr_c")

# NAVCO 2.1 — concession block
m7q   <- fit_cox21    ("event_concession", "domaut_c", data = d21_first_conc)
m8q   <- fit_cox21    ("event_concession", "stterr_c", data = d21_first_conc)
m9q   <- fit_firth21  ("concession_bin",   "domaut_c")
m10q  <- fit_firth21  ("concession_bin",   "stterr_c")
# NAVCO 2.1 — strict-success twins
m7qs  <- fit_cox21    ("event_success_yr", "domaut_c")
m8qs  <- fit_cox21    ("event_success_yr", "stterr_c")
m9qs  <- fit_firth21  ("success_bin_yr",   "domaut_c")
m10qs <- fit_firth21  ("success_bin_yr",   "stterr_c")
# NAVCO 2.1 — ordered bracl
m11q  <- fit_bracl21  ("progress_ord",     "domaut_c")
m12q  <- fit_bracl21  ("progress_ord",     "stterr_c")

###############################################################################
# EXTRACTION HELPER — pull the two interaction terms + joint Wald tests
###############################################################################

extract_curvi <- function(tag, estimator, data_label, response, model,
                          moderator_name, N, n_events) {
  eng_lin  <- "eng_c"
  eng_sq   <- "I(eng_c^2)"
  mod_name <- moderator_name
  int_lin_a <- paste0("eng_c:",       moderator_name)
  int_sq_a  <- paste0("I(eng_c^2):",  moderator_name)
  int_lin_b <- paste0(moderator_name, ":eng_c")        # name-order fallback
  int_sq_b  <- paste0(moderator_name, ":I(eng_c^2)")

  cf      <- coef(model)
  all_nms <- names(cf)

  resolve_name <- function(candidates) {
    hit <- candidates[candidates %in% all_nms]
    if (length(hit) == 0) NA_character_ else hit[1]
  }
  int_lin_name <- resolve_name(c(int_lin_a, int_lin_b))
  int_sq_name  <- resolve_name(c(int_sq_a,  int_sq_b))

  pull_b <- function(nm) if (is.na(nm)) NA_real_ else unname(cf[nm])

  b_eng_lin <- pull_b(eng_lin)
  b_eng_sq  <- pull_b(eng_sq)
  b_mod     <- pull_b(mod_name)
  b_int_lin <- pull_b(int_lin_name)
  b_int_sq  <- pull_b(int_sq_name)

  get_vcov_block <- function(mdl, nms) {
    if (any(is.na(nms))) return(NULL)
    ii <- match(nms, all_nms)
    if (any(is.na(ii))) return(NULL)
    V_full <- tryCatch(vcov(mdl), error = function(e) NULL)
    if (is.null(V_full) && inherits(mdl, "logistf")) V_full <- mdl$var
    if (is.null(V_full)) return(NULL)
    V_full[ii, ii, drop = FALSE]
  }

  get_p <- function(nm) {
    if (is.na(nm)) return(NA_real_)
    if (inherits(model, "survreg")) {
      tbl <- summary(model)$table
      if (nm %in% rownames(tbl)) return(tbl[nm, "p"])
    } else if (inherits(model, "coxph")) {
      tbl <- summary(model)$coefficients
      if (nm %in% rownames(tbl)) return(tbl[nm, "Pr(>|z|)"])
    } else if (inherits(model, "logistf")) {
      pv <- model$prob
      if (nm %in% names(pv)) return(unname(pv[nm]))
    } else if (inherits(model, "bracl")) {
      tbl <- summary(model)$coefficients
      if (nm %in% rownames(tbl)) return(tbl[nm, "Pr(>|z|)"])
    }
    NA_real_
  }

  p_eng_lin <- get_p(eng_lin)
  p_eng_sq  <- get_p(eng_sq)
  p_mod     <- get_p(mod_name)
  p_int_lin <- get_p(int_lin_name)
  p_int_sq  <- get_p(int_sq_name)

 # Joint Wald: (eng:mod, eng^2:mod), df=2 — does the moderator shift the curve?
  joint_wald_int_p <- NA_real_
  if (!is.na(int_lin_name) && !is.na(int_sq_name)) {
    V2   <- get_vcov_block(model, c(int_lin_name, int_sq_name))
    bvec <- c(b_int_lin, b_int_sq)
    if (!is.null(V2) && all(is.finite(bvec)) && all(is.finite(V2))) {
      W <- tryCatch(as.numeric(t(bvec) %*% solve(V2) %*% bvec),
                    error = function(e) NA_real_)
      if (is.finite(W) && W >= 0) joint_wald_int_p <- pchisq(W, df = 2, lower.tail = FALSE)
    }
  }

 # Joint Wald: (eng, eng^2, eng:mod, eng^2:mod), df=4 — any English effect?
  joint_wald_eng_p <- NA_real_
  if (!is.na(int_lin_name) && !is.na(int_sq_name)) {
    V4    <- get_vcov_block(model, c(eng_lin, eng_sq, int_lin_name, int_sq_name))
    bvec4 <- c(b_eng_lin, b_eng_sq, b_int_lin, b_int_sq)
    if (!is.null(V4) && all(is.finite(bvec4)) && all(is.finite(V4))) {
      W <- tryCatch(as.numeric(t(bvec4) %*% solve(V4) %*% bvec4),
                    error = function(e) NA_real_)
      if (is.finite(W) && W >= 0) joint_wald_eng_p <- pchisq(W, df = 4, lower.tail = FALSE)
    }
  }

  k_params <- length(cf)
  epv      <- if (is.na(n_events)) NA_real_ else n_events / k_params

  data.frame(
    model = tag, estimator = estimator, data = data_label,
    response = response, moderator = moderator_name,
    beta_eng_lin = b_eng_lin, p_eng_lin = p_eng_lin,
    beta_eng_sq  = b_eng_sq,  p_eng_sq  = p_eng_sq,
    beta_mod     = b_mod,     p_mod     = p_mod,
    beta_int_lin = b_int_lin, p_int_lin = p_int_lin,
    beta_int_sq  = b_int_sq,  p_int_sq  = p_int_sq,
    joint_wald_int_p = joint_wald_int_p,
    joint_wald_eng_p = joint_wald_eng_p,
    N = N, n_events = n_events, epv = epv,
    stringsAsFactors = FALSE
  )
}

###############################################################################
# PER-MODEL SAMPLE SIZES (complete-case; H21_code.r convention)
###############################################################################

cc_stats <- function(dat, vars, event_var = NA_character_) {
  cc_mask <- complete.cases(dat[, vars, drop = FALSE])
  N <- sum(cc_mask)
  if (is.na(event_var) || !(event_var %in% names(dat))) {
    n_events <- NA_integer_
  } else {
    ev <- suppressWarnings(as.numeric(dat[[event_var]]))
    n_events <- as.integer(sum(ev[cc_mask] == 1, na.rm = TRUE))
  }
  list(N = N, n_events = n_events)
}

vars_13_surv <- function(dv, mod) c("duration", dv, "eng_c", mod, "nonviolent_camp",
                                    "REGCHANGE", "v2csreprss", "lnnum_image_sum")
vars_13_bin  <- function(dv, mod) c(dv, "eng_c", mod, "nonviolent_camp",
                                    "REGCHANGE", "v2csreprss", "lnlengthofcam",
                                    "lnnum_image_sum")
vars_21_surv <- function(dv, mod) c("tstart", "tstop", dv, "eng_c", mod, "nonviolent",
                                    "REGCHANGE", "v2csreprss", "lnnum_image_sum",
                                    "navco21_id")
vars_21_bin  <- function(dv, mod) c(dv, "eng_c", mod, "nonviolent",
                                    "REGCHANGE", "v2csreprss", "time_in_campaign",
                                    "lnnum_image_sum")

s_m1q   <- cc_stats(df,             vars_13_surv("event_concession", "domaut_c"), "event_concession")
s_m2q   <- cc_stats(df,             vars_13_surv("event_concession", "stterr_c"), "event_concession")
s_m3q   <- cc_stats(df,             vars_13_bin ("concession_bin",   "domaut_c"), "concession_bin")
s_m4q   <- cc_stats(df,             vars_13_bin ("concession_bin",   "stterr_c"), "concession_bin")
s_m1qs  <- cc_stats(df,             vars_13_surv("event_success",    "domaut_c"), "event_success")
s_m2qs  <- cc_stats(df,             vars_13_surv("event_success",    "stterr_c"), "event_success")
s_m3qs  <- cc_stats(df,             vars_13_bin ("success_bin",      "domaut_c"), "success_bin")
s_m4qs  <- cc_stats(df,             vars_13_bin ("success_bin",      "stterr_c"), "success_bin")
s_m5q   <- cc_stats(df,             vars_13_bin ("success3",         "domaut_c"))
s_m6q   <- cc_stats(df,             vars_13_bin ("success3",         "stterr_c"))

s_m7q   <- cc_stats(d21_first_conc, vars_21_surv("event_concession", "domaut_c"), "event_concession")
s_m8q   <- cc_stats(d21_first_conc, vars_21_surv("event_concession", "stterr_c"), "event_concession")
s_m9q   <- cc_stats(d21,            vars_21_bin ("concession_bin",   "domaut_c"), "concession_bin")
s_m10q  <- cc_stats(d21,            vars_21_bin ("concession_bin",   "stterr_c"), "concession_bin")
s_m7qs  <- cc_stats(d21,            vars_21_surv("event_success_yr", "domaut_c"), "event_success_yr")
s_m8qs  <- cc_stats(d21,            vars_21_surv("event_success_yr", "stterr_c"), "event_success_yr")
s_m9qs  <- cc_stats(d21,            vars_21_bin ("success_bin_yr",   "domaut_c"), "success_bin_yr")
s_m10qs <- cc_stats(d21,            vars_21_bin ("success_bin_yr",   "stterr_c"), "success_bin_yr")
s_m11q  <- cc_stats(d21,            vars_21_bin ("progress_ord",     "domaut_c"))
s_m12q  <- cc_stats(d21,            vars_21_bin ("progress_ord",     "stterr_c"))

###############################################################################
# BUILD RESULTS TABLE (20 rows)
###############################################################################

results <- rbind(
 # Concession family
  extract_curvi("M1q",  "Weibull AFT", "NAVCO 1.3", "Time to concession",          m1q,  "domaut_c", s_m1q$N,  s_m1q$n_events),
  extract_curvi("M2q",  "Weibull AFT", "NAVCO 1.3", "Time to concession",          m2q,  "stterr_c", s_m2q$N,  s_m2q$n_events),
  extract_curvi("M3q",  "Firth PML",   "NAVCO 1.3", "Concession (binary)",         m3q,  "domaut_c", s_m3q$N,  s_m3q$n_events),
  extract_curvi("M4q",  "Firth PML",   "NAVCO 1.3", "Concession (binary)",         m4q,  "stterr_c", s_m4q$N,  s_m4q$n_events),
  extract_curvi("M7q",  "Cox PH",      "NAVCO 2.1", "Hazard of first concession",  m7q,  "domaut_c", s_m7q$N,  s_m7q$n_events),
  extract_curvi("M8q",  "Cox PH",      "NAVCO 2.1", "Hazard of first concession",  m8q,  "stterr_c", s_m8q$N,  s_m8q$n_events),
  extract_curvi("M9q",  "Firth PML",   "NAVCO 2.1", "Yearly concession (binary)",  m9q,  "domaut_c", s_m9q$N,  s_m9q$n_events),
  extract_curvi("M10q", "Firth PML",   "NAVCO 2.1", "Yearly concession (binary)",  m10q, "stterr_c", s_m10q$N, s_m10q$n_events),
 # Strict-success twins
  extract_curvi("M1qs",  "Weibull AFT", "NAVCO 1.3", "Time to success (strict)",        m1qs,  "domaut_c", s_m1qs$N,  s_m1qs$n_events),
  extract_curvi("M2qs",  "Weibull AFT", "NAVCO 1.3", "Time to success (strict)",        m2qs,  "stterr_c", s_m2qs$N,  s_m2qs$n_events),
  extract_curvi("M3qs",  "Firth PML",   "NAVCO 1.3", "Success (strict, binary)",        m3qs,  "domaut_c", s_m3qs$N,  s_m3qs$n_events),
  extract_curvi("M4qs",  "Firth PML",   "NAVCO 1.3", "Success (strict, binary)",        m4qs,  "stterr_c", s_m4qs$N,  s_m4qs$n_events),
  extract_curvi("M7qs",  "Cox PH",      "NAVCO 2.1", "Hazard of success (strict)",      m7qs,  "domaut_c", s_m7qs$N,  s_m7qs$n_events),
  extract_curvi("M8qs",  "Cox PH",      "NAVCO 2.1", "Hazard of success (strict)",      m8qs,  "stterr_c", s_m8qs$N,  s_m8qs$n_events),
  extract_curvi("M9qs",  "Firth PML",   "NAVCO 2.1", "Yearly success (strict, binary)", m9qs,  "domaut_c", s_m9qs$N,  s_m9qs$n_events),
  extract_curvi("M10qs", "Firth PML",   "NAVCO 2.1", "Yearly success (strict, binary)", m10qs, "stterr_c", s_m10qs$N, s_m10qs$n_events),
 # Ordered bracl
  extract_curvi("M5q",  "bracl", "NAVCO 1.3", "success3 (0/1/2) ordered",                m5q,  "domaut_c", s_m5q$N,  s_m5q$n_events),
  extract_curvi("M6q",  "bracl", "NAVCO 1.3", "success3 (0/1/2) ordered",                m6q,  "stterr_c", s_m6q$N,  s_m6q$n_events),
  extract_curvi("M11q", "bracl", "NAVCO 2.1", "progress_ord (6-level; failure separate)", m11q, "domaut_c", s_m11q$N, s_m11q$n_events),
  extract_curvi("M12q", "bracl", "NAVCO 2.1", "progress_ord (6-level; failure separate)", m12q, "stterr_c", s_m12q$N, s_m12q$n_events)
)

###############################################################################
# MULTIPLICITY: BH within 3 sub-families (no pooling)
###############################################################################

family_concession <- c("M1q","M2q","M3q","M4q","M7q","M8q","M9q","M10q")
family_strict     <- c("M1qs","M2qs","M3qs","M4qs","M7qs","M8qs","M9qs","M10qs")
family_ordered    <- c("M5q","M6q","M11q","M12q")

results$family <- dplyr::case_when(
  results$model %in% family_concession ~ "concession",
  results$model %in% family_strict     ~ "strict_success",
  results$model %in% family_ordered    ~ "ordered",
  TRUE                                  ~ NA_character_
)

idx_conc <- which(results$family == "concession")
idx_str  <- which(results$family == "strict_success")
idx_ord  <- which(results$family == "ordered")

# Adjust each of the four p-summaries within each sub-family.
# Use explicit n=length(idx) to prevent NA-shrinkage.
for (col in c("p_int_lin", "p_int_sq", "joint_wald_int_p", "joint_wald_eng_p")) {
  bh_col <- paste0(col, "_bh")
  results[[bh_col]] <- NA_real_
  if (length(idx_conc) > 0) {
    results[[bh_col]][idx_conc] <- p.adjust(results[[col]][idx_conc],
                                            method = "BH", n = length(idx_conc))
  }
  if (length(idx_str) > 0) {
    results[[bh_col]][idx_str ] <- p.adjust(results[[col]][idx_str ],
                                            method = "BH", n = length(idx_str ))
  }
  if (length(idx_ord) > 0) {
    results[[bh_col]][idx_ord ] <- p.adjust(results[[col]][idx_ord ],
                                            method = "BH", n = length(idx_ord ))
  }
}

results$epv_warning <- !is.na(results$epv) & results$epv < 5

###############################################################################
# COMPACT PRINTOUT
###############################################################################

fmt_p <- function(p) ifelse(is.na(p), "   NA ", sprintf("%6.3f", p))
fmt_b <- function(b) ifelse(is.na(b), "    NA ", sprintf("%7.2f", b))

print_family <- function(label, fam_models, k) {
  cat(sprintf("\n=========== %s (k=%d, BH-adjusted within sub-family) ===========\n", label, k))
  cat(sprintf("%-6s %-12s %-10s %-9s %8s %7s %9s %7s %8s %10s %5s %4s %3s\n",
              "Model","Estimator","Data","Moderator",
              "b_int_sq","p_isq","p_isq_BH",
              "p_iln","p_JWint","p_JWint_BH","N","EPV","OK?"))
  cat(paste(rep("-", 130), collapse = ""), "\n")
  rows <- results[results$model %in% fam_models, ]
  for (i in seq_len(nrow(rows))) {
    r <- rows[i, ]
    cat(sprintf("%-6s %-12s %-10s %-9s %8s %7s %9s %7s %8s %10s %5d %4s %3s\n",
                r$model, r$estimator, r$data, r$moderator,
                fmt_b(r$beta_int_sq), fmt_p(r$p_int_sq), fmt_p(r$p_int_sq_bh),
                fmt_p(r$p_int_lin), fmt_p(r$joint_wald_int_p),
                fmt_p(r$joint_wald_int_p_bh),
                r$N,
                ifelse(is.na(r$epv), "  NA", sprintf("%4.1f", r$epv)),
                ifelse(r$epv_warning, "!", " ")))
  }
}

print_family("CONCESSION FAMILY",    family_concession, length(idx_conc))
print_family("STRICT-SUCCESS TWINS", family_strict,     length(idx_str))
print_family("ORDERED bracl FAMILY", family_ordered,    length(idx_ord))

cat("\nNotes:\n")
cat(" b_int_sq = coefficient on I(eng_c^2):mod_c (CURVATURE-shift; novel)\n")
cat(" p_isq = raw p on I(eng_c^2):mod_c\n")
cat(" p_isq_BH = BH-adjusted p on I(eng_c^2):mod_c within sub-family\n")
cat(" p_iln = raw p on eng_c:mod_c (linear-shift, mirrors H21)\n")
cat(" p_JWint = joint Wald p on (eng:mod, eng^2:mod), df=2\n")
cat(" p_JWint_BH = BH-adjusted joint Wald p within sub-family\n")
cat(" EPV = events per variable (Peduzzi 1996; '!' marks EPV<5)\n")
cat(" Appendix-only ( spirit, parsimony); diagnostics omitted by design.\n")

###############################################################################
# FULL MODEL SUMMARIES
###############################################################################

cat("\n\n==== FULL MODEL SUMMARIES ====\n")
cat("\n---- M1q Weibull AFT × domaut — concession ----\n"); print(summary(m1q))
cat("\n---- M2q Weibull AFT × stterr — concession ----\n"); print(summary(m2q))
cat("\n---- M3q Firth PML × domaut — concession (binary) ----\n"); print(summary(m3q))
cat("\n---- M4q Firth PML × stterr — concession (binary) ----\n"); print(summary(m4q))
cat("\n---- M1qs Weibull AFT × domaut — STRICT success ----\n"); print(summary(m1qs))
cat("\n---- M2qs Weibull AFT × stterr — STRICT success ----\n"); print(summary(m2qs))
cat("\n---- M3qs Firth PML × domaut — STRICT success ----\n"); print(summary(m3qs))
cat("\n---- M4qs Firth PML × stterr — STRICT success ----\n"); print(summary(m4qs))
cat("\n---- M5q bracl × domaut — success3 (NAVCO 1.3) ----\n"); print(summary(m5q))
cat("\n---- M6q bracl × stterr — success3 (NAVCO 1.3) ----\n"); print(summary(m6q))
cat("\n---- M7q Cox PH × domaut — first concession ----\n"); print(summary(m7q))
cat("\n---- M8q Cox PH × stterr — first concession ----\n"); print(summary(m8q))
cat("\n---- M7qs Cox PH × domaut — STRICT success ----\n"); print(summary(m7qs))
cat("\n---- M8qs Cox PH × stterr — STRICT success ----\n"); print(summary(m8qs))
cat("\n---- M9q Firth PML × domaut — yearly concession ----\n"); print(summary(m9q))
cat("\n---- M10q Firth PML × stterr — yearly concession ----\n"); print(summary(m10q))
cat("\n---- M9qs Firth PML × domaut — yearly STRICT success ----\n"); print(summary(m9qs))
cat("\n---- M10qs Firth PML × stterr — yearly STRICT success ----\n"); print(summary(m10qs))
cat("\n---- M11q bracl × domaut — progress_ord (NAVCO 2.1) ----\n"); print(summary(m11q))
cat("\n---- M12q bracl × stterr — progress_ord (NAVCO 2.1) ----\n"); print(summary(m12q))

###############################################################################
# EXPORT
###############################################################################

if (nrow(results) != 20) {
  stop(sprintf("Expected 20 result rows, found %d.", nrow(results)), call. = FALSE)
}
if (anyDuplicated(results$model) > 0) {
 stop("Duplicate model labels detected.", call. = FALSE)
}

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
write.csv(results, out_path, row.names = FALSE)

curvilinear_models <- list(
  M1q = m1q, M2q = m2q, M3q = m3q, M4q = m4q,
  M7q = m7q, M8q = m8q, M9q = m9q, M10q = m10q,
  M1qs = m1qs, M2qs = m2qs, M3qs = m3qs, M4qs = m4qs,
  M7qs = m7qs, M8qs = m8qs, M9qs = m9qs, M10qs = m10qs,
  M5q = m5q, M6q = m6q, M11q = m11q, M12q = m12q
)
stopifnot(identical(names(curvilinear_models), results$model))

source(file.path(P$code, "R", "model_h21_limited.R"), local = FALSE)
limited <- fit_h21_limited_extension(df, d21)
stopifnot(identical(names(limited$models), limited$results$model))
write.csv(limited$results, limited_out_path, row.names = FALSE)

models <- c(curvilinear_models, limited$models)

logistf_convergence <- bind_rows(lapply(names(models), function(model_id) {
  model <- models[[model_id]]
  if (!inherits(model, "logistf")) return(NULL)
  conv <- abs(model$conv[c("LL change", "max abs score", "beta change")])
  thresholds <- c(model$control$lconv, model$control$gconv, model$control$xconv)
  data.frame(
    model = model_id,
    ll_change = unname(conv[[1L]]),
    max_abs_score = unname(conv[[2L]]),
    beta_change = unname(conv[[3L]]),
    ll_threshold = thresholds[[1L]],
    score_threshold = thresholds[[2L]],
    beta_threshold = thresholds[[3L]],
    converged = all(is.finite(conv)) && all(conv <= thresholds),
    stringsAsFactors = FALSE
  )
}))
logistf_convergence$accepted_values_preserved <- TRUE
logistf_convergence$reporting_class <- ifelse(
  logistf_convergence$converged,
  "accepted_replicated",
  "accepted_archival_nonconverged"
)
write.csv(logistf_convergence, convergence_path, row.names = FALSE)
if (any(!logistf_convergence$converged)) {
  warning(
    sprintf(
      "Accepted H2.1 extension estimates preserved, but non-converged logistf model(s) are quarantined in diagnostics: %s",
      paste(logistf_convergence$model[!logistf_convergence$converged], collapse = ", ")
    ),
    call. = FALSE
  )
}
saveRDS(models, cache_path, version = 3)
cat(sprintf("\nResults exported to: %s\n", out_path))
cat(sprintf("Limited-success results exported to: %s\n", limited_out_path))
cat(sprintf("Model bundle exported to: %s\n", cache_path))
cat(sprintf("Convergence diagnostics exported to: %s\n", convergence_path))

cat("\n[H21_curvilinear_appendix.r complete]\n")
