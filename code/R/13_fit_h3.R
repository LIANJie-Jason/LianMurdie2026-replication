###############################################################################
# H3_code.r
# H3: DOMESTIC / INTERNAL VULNERABILITY
#
# Tests whether English-language protest signs are associated with more
# favorable campaign outcomes when domestic institutions are more open.
#
# Moderators:
#   v2csreprss  : higher = less CSO repression
#   v2x_civlib  : higher = more civil liberties
#   v2x_clpol   : higher = more political liberty
#   p_polity2   : higher = more democratic
#
# Core branch conventions:
#   - local empirical-clean inputs/outputs only
#   - first-event Cox concession rows as main estimand
#   - recurrent-event Cox concession rows retained only as labeled robustness
#   - source-clustered GEE companions for yearly concession-binary NAVCO 2.1 rows
#   - 6-level ordered NAVCO 2.1 outcome with failure separate from status quo
#   - PH diagnostics exported for Cox rows
#   - main panel = M1-M8 (Holm k=8), appendix = remaining 32 core rows (BH k=32)
#
# Output:
#   RR_datacode/empirical clean/results/H3_results.csv
#   48 rows = 40 core + 4 recurrent-event Cox robustness + 4 clustered-GEE rows
###############################################################################

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
out_path     <- file.path(P$estimates, "H3_results.csv")
cache_path   <- file.path(P$cache, "H3_models.rds")
frames_path  <- file.path(P$cache, "H3_analysis_frames.rds")

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

assert_no_nonfinite <- function(x, label) {
  bad <- sum(!is.finite(x[!is.na(x)]))
  if (bad > 0) {
    stop(sprintf("%s contains %d non-finite non-missing values.", label, bad), call. = FALSE)
  }
}

make_progress_ord6 <- function(x) {
  ordered(x, levels = c(5, 0, 1, 2, 3, 4))
}

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

###############################################################################
# LOAD DATA
###############################################################################

df  <- read.csv(navco13_path)
d21 <- read.csv(navco21_path)

assert_required_cols(
  df,
  c(
    "CAMPAIGN", "LOCATION", "EYEAR", "BYEAR", "success", "limited", "success3",
    "viol", "nonviolent_camp", "REGCHANGE", "eng_prob_general", "lnnum_image_sum", "lnlengthofcam",
    "v2csreprss", "v2x_civlib", "v2x_clpol", "p_polity2", "lngdp"
  ),
  "df_final.csv"
)
assert_required_cols(
  d21,
  c(
    "navco21_id", "Year", "CAMPAIGN", "LOCATION", "panel_campaign_id",
    "nonviolent", "camp_goals", "Num_Image", "progress", "concession",
    "time_in_campaign", "eng_prop_year", "v2csreprss", "v2x_civlib",
    "v2x_clpol", "p_polity2", "lngdp_yr"
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

df$duration         <- df$EYEAR - df$BYEAR + 1
df$event_concession <- as.integer(
  as.numeric(as.character(df$success)) == 1 |
  as.numeric(as.character(df$limited)) == 1
)
df$event_success    <- as.numeric(as.character(df$success))
df$concession_bin   <- df$event_concession
df$success_bin      <- df$event_success

df$eng_c     <- df$eng_prob_general - mean(df$eng_prob_general, na.rm = TRUE)
df$repress_c <- df$v2csreprss       - mean(df$v2csreprss,       na.rm = TRUE)
df$civlib_c  <- df$v2x_civlib       - mean(df$v2x_civlib,       na.rm = TRUE)
df$clpol_c   <- df$v2x_clpol        - mean(df$v2x_clpol,        na.rm = TRUE)
df$polity_c  <- df$p_polity2        - mean(df$p_polity2,        na.rm = TRUE)

###############################################################################
# NAVCO 2.1 PREPARATION
###############################################################################

d21$nonviolent      <- as.factor(d21$nonviolent)
d21$REGCHANGE       <- as.numeric(d21$camp_goals == 0)
d21$lnnum_image_sum <- log(d21$Num_Image)
d21$progress_ord    <- make_progress_ord6(d21$progress)
d21$event_concession <- as.numeric(d21$concession)
d21$event_success_yr <- as.numeric(d21$progress == 4)
d21$success_bin_yr   <- d21$event_success_yr
d21$concession_bin   <- d21$event_concession

d21$tstart <- d21$time_in_campaign - 1
d21$tstop  <- d21$time_in_campaign
assert_no_nonfinite(d21$tstart, "d21$tstart")
assert_no_nonfinite(d21$tstop, "d21$tstop")

d21$eng_c     <- d21$eng_prop_year - mean(d21$eng_prop_year, na.rm = TRUE)
d21$repress_c <- d21$v2csreprss    - mean(d21$v2csreprss,    na.rm = TRUE)
d21$civlib_c  <- d21$v2x_civlib    - mean(d21$v2x_civlib,    na.rm = TRUE)
d21$clpol_c   <- d21$v2x_clpol     - mean(d21$v2x_clpol,     na.rm = TRUE)
d21$polity_c  <- d21$p_polity2     - mean(d21$p_polity2,     na.rm = TRUE)

d21_first_conc <- build_first_event_panel(d21, "navco21_id", "event_concession", "time_in_campaign")
assert_first_event_panel(d21_first_conc, "navco21_id", "event_concession")

###############################################################################
# CONTROL-SWAP HELPERS
###############################################################################

get_dom_control_13 <- function(moderator) {
  if (moderator == "repress_c") "p_polity2" else "v2csreprss"
}

get_dom_control_21 <- function(moderator) {
  if (moderator == "repress_c") NA_character_ else "v2csreprss"
}

control_slot <- function(term) {
  if (is.na(term)) "" else paste0(" + ", term)
}

###############################################################################
# FIT HELPERS
###############################################################################

fit_weibull13 <- function(dv_event, moderator) {
  dom <- get_dom_control_13(moderator)
  f <- as.formula(sprintf(
    "Surv(duration, %s) ~ eng_c * %s + nonviolent_camp + REGCHANGE%s + lnnum_image_sum + lngdp",
    dv_event, moderator, control_slot(dom)
  ))
  survreg(f, data = df, dist = "weibull")
}

fit_firth13 <- function(dv_col, moderator) {
  dom <- get_dom_control_13(moderator)
  f <- as.formula(sprintf(
    "%s ~ eng_c * %s + nonviolent_camp + REGCHANGE%s + lnlengthofcam + lnnum_image_sum + lngdp",
    dv_col, moderator, control_slot(dom)
  ))
  logistf(
    f, data = df,
    control   = logistf.control(maxit = 1000, maxstep = 0.5),
    plcontrol = logistpl.control(maxit = 1000),
    pl        = TRUE
  )
}

fit_bracl13 <- function(dv_col, moderator) {
  dom <- get_dom_control_13(moderator)
  f <- as.formula(sprintf(
    "%s ~ eng_c * %s + nonviolent_camp + REGCHANGE%s + lnlengthofcam + lnnum_image_sum + lngdp",
    dv_col, moderator, control_slot(dom)
  ))
  bracl(f, data = df, type = "MPL_Jeffreys", link = "logit", parallel = TRUE, maxit = 500)
}

fit_cox21 <- function(dv_event, moderator, data = d21) {
  dom <- get_dom_control_21(moderator)
  f <- as.formula(sprintf(
    "Surv(tstart, tstop, %s) ~ eng_c * %s + nonviolent + REGCHANGE%s + lnnum_image_sum + lngdp_yr + cluster(navco21_id)",
    dv_event, moderator, control_slot(dom)
  ))
  coxph(f, data = data)
}

fit_firth21 <- function(dv_col, moderator) {
  dom <- get_dom_control_21(moderator)
  f <- as.formula(sprintf(
    "%s ~ eng_c * %s + nonviolent + REGCHANGE%s + time_in_campaign + lnnum_image_sum + lngdp_yr",
    dv_col, moderator, control_slot(dom)
  ))
  logistf(
    f, data = d21,
    control   = logistf.control(maxit = 1000, maxstep = 0.5),
    plcontrol = logistpl.control(maxit = 1000),
    pl        = TRUE
  )
}

fit_bracl21 <- function(dv_col, moderator) {
  dom <- get_dom_control_21(moderator)
  f <- as.formula(sprintf(
    "%s ~ eng_c * %s + nonviolent + REGCHANGE%s + time_in_campaign + lnnum_image_sum + lngdp_yr",
    dv_col, moderator, control_slot(dom)
  ))
  bracl(f, data = d21, type = "MPL_Jeffreys", link = "logit", parallel = TRUE, maxit = 500)
}

fit_gee21 <- function(dv_col, moderator) {
  dom <- get_dom_control_21(moderator)
  vars <- c(
    dv_col, "eng_c", moderator, "nonviolent", "REGCHANGE",
    "time_in_campaign", "lnnum_image_sum", "lngdp_yr", "navco21_id"
  )
  if (!is.na(dom)) vars <- c(vars, dom)
  dd <- d21[complete.cases(d21[, vars]), vars, drop = FALSE]
  dd <- dd %>% arrange(navco21_id, time_in_campaign)
  f <- as.formula(sprintf(
    "%s ~ eng_c * %s + nonviolent + REGCHANGE%s + time_in_campaign + lnnum_image_sum + lngdp_yr",
    dv_col, moderator, control_slot(dom)
  ))
  geeglm(f, data = dd, id = navco21_id, family = binomial, corstr = "exchangeable")
}

###############################################################################
# FIT ALL MODELS
###############################################################################

# NAVCO 1.3
model1   <- fit_weibull13("event_concession", "repress_c")
model1s  <- fit_weibull13("event_success",    "repress_c")
model2   <- fit_weibull13("event_concession", "civlib_c")
model2s  <- fit_weibull13("event_success",    "civlib_c")
model3   <- fit_weibull13("event_concession", "clpol_c")
model3s  <- fit_weibull13("event_success",    "clpol_c")
model4   <- fit_weibull13("event_concession", "polity_c")
model4s  <- fit_weibull13("event_success",    "polity_c")

model5   <- fit_firth13("concession_bin", "repress_c")
model5s  <- fit_firth13("success_bin",    "repress_c")
model6   <- fit_firth13("concession_bin", "civlib_c")
model6s  <- fit_firth13("success_bin",    "civlib_c")
model7   <- fit_firth13("concession_bin", "clpol_c")
model7s  <- fit_firth13("success_bin",    "clpol_c")
model8   <- fit_firth13("concession_bin", "polity_c")
model8s  <- fit_firth13("success_bin",    "polity_c")

model9   <- fit_bracl13("success3", "repress_c")
model10  <- fit_bracl13("success3", "civlib_c")
model11  <- fit_bracl13("success3", "clpol_c")
model12  <- fit_bracl13("success3", "polity_c")

# NAVCO 2.1 core
model13  <- fit_cox21("event_concession", "repress_c", data = d21_first_conc)
model13s <- fit_cox21("event_success_yr", "repress_c")
model14  <- fit_cox21("event_concession", "civlib_c",  data = d21_first_conc)
model14s <- fit_cox21("event_success_yr", "civlib_c")
model15  <- fit_cox21("event_concession", "clpol_c",   data = d21_first_conc)
model15s <- fit_cox21("event_success_yr", "clpol_c")
model16  <- fit_cox21("event_concession", "polity_c",  data = d21_first_conc)
model16s <- fit_cox21("event_success_yr", "polity_c")

model17  <- fit_firth21("concession_bin", "repress_c")
model17s <- fit_firth21("success_bin_yr", "repress_c")
model18  <- fit_firth21("concession_bin", "civlib_c")
model18s <- fit_firth21("success_bin_yr", "civlib_c")
model19  <- fit_firth21("concession_bin", "clpol_c")
model19s <- fit_firth21("success_bin_yr", "clpol_c")
model20  <- fit_firth21("concession_bin", "polity_c")
model20s <- fit_firth21("success_bin_yr", "polity_c")

model21  <- fit_bracl21("progress_ord", "repress_c")
model22  <- fit_bracl21("progress_ord", "civlib_c")
model23  <- fit_bracl21("progress_ord", "clpol_c")
model24  <- fit_bracl21("progress_ord", "polity_c")

# NAVCO 2.1 robustness
model13_ag <- fit_cox21("event_concession", "repress_c", data = d21)
model14_ag <- fit_cox21("event_concession", "civlib_c",  data = d21)
model15_ag <- fit_cox21("event_concession", "clpol_c",   data = d21)
model16_ag <- fit_cox21("event_concession", "polity_c",  data = d21)

model17_gee <- fit_gee21("concession_bin", "repress_c")
model18_gee <- fit_gee21("concession_bin", "civlib_c")
model19_gee <- fit_gee21("concession_bin", "clpol_c")
model20_gee <- fit_gee21("concession_bin", "polity_c")

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
  tab      <- z$table
  global_p <- tab["GLOBAL", "p"]
  int_p    <- if (focal_name %in% rownames(tab)) tab[focal_name, "p"] else NA_real_
  ph_flag  <- if (is.finite(int_p) && int_p < 0.05) {
    "PH_FAIL_INTERACTION"
  } else if (is.finite(global_p) && global_p < 0.05) {
    "PH_FAIL_GLOBAL"
  } else if (!is.finite(int_p)) {
    "PH_FOCAL_MISSING"
  } else {
    ""
  }
  cat(sprintf("[cox.zph] %-6s | global p = %.4f | %s p = %s\n",
              tag, global_p,
              focal_name,
              if (is.finite(int_p)) formatC(int_p, format = "f", digits = 4) else "NA"))
  data.frame(
    model = tag,
    ph_global_p = global_p,
    ph_int_p = int_p,
    ph_focal_p = int_p,
    ph_flag = ph_flag,
    stringsAsFactors = FALSE
  )
}

cat("\n--- Cox PH diagnostic (cox.zph() scaled Schoenfeld residuals) ---\n")
cox_diag_tbl <- bind_rows(
  cox_ph_diag("M13",  model13,  "eng_c:repress_c"),
  cox_ph_diag("M13s", model13s, "eng_c:repress_c"),
  cox_ph_diag("M14",  model14,  "eng_c:civlib_c"),
  cox_ph_diag("M14s", model14s, "eng_c:civlib_c"),
  cox_ph_diag("M15",  model15,  "eng_c:clpol_c"),
  cox_ph_diag("M15s", model15s, "eng_c:clpol_c"),
  cox_ph_diag("M16",  model16,  "eng_c:polity_c"),
  cox_ph_diag("M16s", model16s, "eng_c:polity_c")
)
cat("--- end Cox PH diagnostic ---\n\n")

###############################################################################
# EXTRACTION HELPER
###############################################################################

extract_h3 <- function(tag, estimator, data_label, response, model,
                       moderator_name, N, n_events, mod_direction = "pos") {
  eng_name <- "eng_c"
  mod_name <- moderator_name
  int_name <- paste0("eng_c:", moderator_name)

  cf    <- coef(model)
  b_eng <- unname(cf[eng_name])
  b_mod <- unname(cf[mod_name])
  b_int <- unname(cf[int_name])

  get_vcov_block <- function(mdl, nms) {
    all_nms <- names(coef(mdl))
    ii <- match(nms, all_nms)
    if (any(is.na(ii))) return(NULL)
    V_full <- tryCatch(vcov(mdl), error = function(e) NULL)
    if (is.null(V_full) && inherits(mdl, "logistf")) V_full <- mdl$var
    if (is.null(V_full)) return(NULL)
    V_full[ii, ii, drop = FALSE]
  }

  if (inherits(model, "survreg")) {
    tbl   <- summary(model)$table
    p_eng <- tbl[eng_name, "p"]
    p_mod <- tbl[mod_name, "p"]
    p_int <- tbl[int_name, "p"]
  } else if (inherits(model, "coxph")) {
    tbl   <- summary(model)$coefficients
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

  joint_wald_p <- NA_real_
  V2 <- get_vcov_block(model, c(eng_name, int_name))
  if (!is.null(V2) && all(is.finite(c(b_eng, b_int))) && all(is.finite(V2))) {
    bvec <- c(b_eng, b_int)
    W    <- tryCatch(as.numeric(t(bvec) %*% solve(V2) %*% bvec),
                     error = function(e) NA_real_)
    if (is.finite(W)) joint_wald_p <- pchisq(W, df = 2, lower.tail = FALSE)
  }

  expected_int_sign <- if (mod_direction == "pos") {
    switch(estimator,
           "Weibull AFT" = "-",
           "Firth PML"   = "+",
           "GEE logit"   = "+",
           "Cox PH"      = "+",
           "bracl"       = "-",
           NA_character_)
  } else {
    switch(estimator,
           "Weibull AFT" = "+",
           "Firth PML"   = "-",
           "GEE logit"   = "-",
           "Cox PH"      = "-",
           "bracl"       = "+",
           NA_character_)
  }

  sign_consistent <- is.finite(b_int) && !is.na(expected_int_sign) && (
    (expected_int_sign == "+" && b_int > 0) ||
      (expected_int_sign == "-" && b_int < 0)
  )

  k_params <- length(cf)
  epv      <- if (is.na(n_events)) NA_real_ else n_events / k_params

  data.frame(
    model             = tag,
    estimator         = estimator,
    data              = data_label,
    response          = response,
    moderator         = moderator_name,
    mod_direction     = mod_direction,
    beta_eng_c        = b_eng,
    p_eng_c           = p_eng,
    beta_mod_c        = b_mod,
    p_mod_c           = p_mod,
    beta_int          = b_int,
    p_int             = p_int,
    joint_wald_p      = joint_wald_p,
    expected_int_sign = expected_int_sign,
    sign_consistent   = sign_consistent,
    N                 = N,
    n_events          = n_events,
    epv               = epv,
    sig_05_int        = as.logical(!is.na(p_int) && p_int < 0.05),
    sig_10_int        = as.logical(!is.na(p_int) && p_int < 0.10),
    epv_warning       = !is.na(epv) & epv < 5,
    stringsAsFactors  = FALSE
  )
}

###############################################################################
# SAMPLE-SIZE HELPERS
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

vars_13_surv <- function(dv, moderator) {
  dom <- get_dom_control_13(moderator)
  base <- c("duration", dv, "eng_c", moderator, "nonviolent_camp", "REGCHANGE", "lnnum_image_sum", "lngdp")
  if (!is.na(dom)) c(base, dom) else base
}

vars_13_bin <- function(dv, moderator) {
  dom <- get_dom_control_13(moderator)
  base <- c(dv, "eng_c", moderator, "nonviolent_camp", "REGCHANGE", "lnlengthofcam", "lnnum_image_sum", "lngdp")
  if (!is.na(dom)) c(base, dom) else base
}

vars_21_surv <- function(dv, moderator) {
  dom <- get_dom_control_21(moderator)
  base <- c("tstart", "tstop", dv, "eng_c", moderator, "nonviolent",
            "REGCHANGE", "lnnum_image_sum", "lngdp_yr", "navco21_id")
  if (!is.na(dom)) c(base, dom) else base
}

vars_21_bin <- function(dv, moderator) {
  dom <- get_dom_control_21(moderator)
  base <- c(dv, "eng_c", moderator, "nonviolent",
            "REGCHANGE", "time_in_campaign", "lnnum_image_sum", "lngdp_yr")
  if (!is.na(dom)) c(base, dom) else base
}

s_m1  <- cc_stats(df, vars_13_surv("event_concession", "repress_c"), "event_concession")
s_m1s <- cc_stats(df, vars_13_surv("event_success",    "repress_c"), "event_success")
s_m2  <- cc_stats(df, vars_13_surv("event_concession", "civlib_c"),  "event_concession")
s_m2s <- cc_stats(df, vars_13_surv("event_success",    "civlib_c"),  "event_success")
s_m3  <- cc_stats(df, vars_13_surv("event_concession", "clpol_c"),    "event_concession")
s_m3s <- cc_stats(df, vars_13_surv("event_success",    "clpol_c"),    "event_success")
s_m4  <- cc_stats(df, vars_13_surv("event_concession", "polity_c"),  "event_concession")
s_m4s <- cc_stats(df, vars_13_surv("event_success",    "polity_c"),  "event_success")

s_m5  <- cc_stats(df, vars_13_bin("concession_bin", "repress_c"), "concession_bin")
s_m5s <- cc_stats(df, vars_13_bin("success_bin",    "repress_c"), "success_bin")
s_m6  <- cc_stats(df, vars_13_bin("concession_bin", "civlib_c"),  "concession_bin")
s_m6s <- cc_stats(df, vars_13_bin("success_bin",    "civlib_c"),  "success_bin")
s_m7  <- cc_stats(df, vars_13_bin("concession_bin", "clpol_c"),    "concession_bin")
s_m7s <- cc_stats(df, vars_13_bin("success_bin",    "clpol_c"),    "success_bin")
s_m8  <- cc_stats(df, vars_13_bin("concession_bin", "polity_c"),  "concession_bin")
s_m8s <- cc_stats(df, vars_13_bin("success_bin",    "polity_c"),  "success_bin")

s_m9  <- cc_stats(df, vars_13_bin("success3", "repress_c"))
s_m10 <- cc_stats(df, vars_13_bin("success3", "civlib_c"))
s_m11 <- cc_stats(df, vars_13_bin("success3", "clpol_c"))
s_m12 <- cc_stats(df, vars_13_bin("success3", "polity_c"))

s_m13  <- cc_stats(d21_first_conc, vars_21_surv("event_concession", "repress_c"), "event_concession")
s_m13s <- cc_stats(d21, vars_21_surv("event_success_yr", "repress_c"), "event_success_yr")
s_m14  <- cc_stats(d21_first_conc, vars_21_surv("event_concession", "civlib_c"),  "event_concession")
s_m14s <- cc_stats(d21, vars_21_surv("event_success_yr", "civlib_c"),  "event_success_yr")
s_m15  <- cc_stats(d21_first_conc, vars_21_surv("event_concession", "clpol_c"),    "event_concession")
s_m15s <- cc_stats(d21, vars_21_surv("event_success_yr", "clpol_c"),    "event_success_yr")
s_m16  <- cc_stats(d21_first_conc, vars_21_surv("event_concession", "polity_c"),  "event_concession")
s_m16s <- cc_stats(d21, vars_21_surv("event_success_yr", "polity_c"),  "event_success_yr")

s_m13_ag <- cc_stats(d21, vars_21_surv("event_concession", "repress_c"), "event_concession")
s_m14_ag <- cc_stats(d21, vars_21_surv("event_concession", "civlib_c"),  "event_concession")
s_m15_ag <- cc_stats(d21, vars_21_surv("event_concession", "clpol_c"),    "event_concession")
s_m16_ag <- cc_stats(d21, vars_21_surv("event_concession", "polity_c"),  "event_concession")

s_m17  <- cc_stats(d21, vars_21_bin("concession_bin", "repress_c"), "concession_bin")
s_m17s <- cc_stats(d21, vars_21_bin("success_bin_yr", "repress_c"), "success_bin_yr")
s_m18  <- cc_stats(d21, vars_21_bin("concession_bin", "civlib_c"),  "concession_bin")
s_m18s <- cc_stats(d21, vars_21_bin("success_bin_yr", "civlib_c"),  "success_bin_yr")
s_m19  <- cc_stats(d21, vars_21_bin("concession_bin", "clpol_c"),    "concession_bin")
s_m19s <- cc_stats(d21, vars_21_bin("success_bin_yr", "clpol_c"),    "success_bin_yr")
s_m20  <- cc_stats(d21, vars_21_bin("concession_bin", "polity_c"),  "concession_bin")
s_m20s <- cc_stats(d21, vars_21_bin("success_bin_yr", "polity_c"),  "success_bin_yr")

s_m17_gee <- cc_stats(d21, vars_21_bin("concession_bin", "repress_c"), "concession_bin")
s_m18_gee <- cc_stats(d21, vars_21_bin("concession_bin", "civlib_c"),  "concession_bin")
s_m19_gee <- cc_stats(d21, vars_21_bin("concession_bin", "clpol_c"),    "concession_bin")
s_m20_gee <- cc_stats(d21, vars_21_bin("concession_bin", "polity_c"),  "concession_bin")

s_m21 <- cc_stats(d21, vars_21_bin("progress_ord", "repress_c"))
s_m22 <- cc_stats(d21, vars_21_bin("progress_ord", "civlib_c"))
s_m23 <- cc_stats(d21, vars_21_bin("progress_ord", "clpol_c"))
s_m24 <- cc_stats(d21, vars_21_bin("progress_ord", "polity_c"))

###############################################################################
# BUILD RESULTS
###############################################################################

core_results <- rbind(
  extract_h3("M1",   "Weibull AFT", "NAVCO 1.3", "Time to concession",          model1,   "repress_c", s_m1$N,   s_m1$n_events),
  extract_h3("M1s",  "Weibull AFT", "NAVCO 1.3", "Time to success (strict)",    model1s,  "repress_c", s_m1s$N,  s_m1s$n_events),
  extract_h3("M2",   "Weibull AFT", "NAVCO 1.3", "Time to concession",          model2,   "civlib_c",  s_m2$N,   s_m2$n_events),
  extract_h3("M2s",  "Weibull AFT", "NAVCO 1.3", "Time to success (strict)",    model2s,  "civlib_c",  s_m2s$N,  s_m2s$n_events),
  extract_h3("M3",   "Weibull AFT", "NAVCO 1.3", "Time to concession",          model3,   "clpol_c",    s_m3$N,   s_m3$n_events),
  extract_h3("M3s",  "Weibull AFT", "NAVCO 1.3", "Time to success (strict)",    model3s,  "clpol_c",    s_m3s$N,  s_m3s$n_events),
  extract_h3("M4",   "Weibull AFT", "NAVCO 1.3", "Time to concession",          model4,   "polity_c",  s_m4$N,   s_m4$n_events),
  extract_h3("M4s",  "Weibull AFT", "NAVCO 1.3", "Time to success (strict)",    model4s,  "polity_c",  s_m4s$N,  s_m4s$n_events),

  extract_h3("M5",   "Firth PML",   "NAVCO 1.3", "Concession (binary)",         model5,   "repress_c", s_m5$N,   s_m5$n_events),
  extract_h3("M5s",  "Firth PML",   "NAVCO 1.3", "Success (strict, binary)",    model5s,  "repress_c", s_m5s$N,  s_m5s$n_events),
  extract_h3("M6",   "Firth PML",   "NAVCO 1.3", "Concession (binary)",         model6,   "civlib_c",  s_m6$N,   s_m6$n_events),
  extract_h3("M6s",  "Firth PML",   "NAVCO 1.3", "Success (strict, binary)",    model6s,  "civlib_c",  s_m6s$N,  s_m6s$n_events),
  extract_h3("M7",   "Firth PML",   "NAVCO 1.3", "Concession (binary)",         model7,   "clpol_c",    s_m7$N,   s_m7$n_events),
  extract_h3("M7s",  "Firth PML",   "NAVCO 1.3", "Success (strict, binary)",    model7s,  "clpol_c",    s_m7s$N,  s_m7s$n_events),
  extract_h3("M8",   "Firth PML",   "NAVCO 1.3", "Concession (binary)",         model8,   "polity_c",  s_m8$N,   s_m8$n_events),
  extract_h3("M8s",  "Firth PML",   "NAVCO 1.3", "Success (strict, binary)",    model8s,  "polity_c",  s_m8s$N,  s_m8s$n_events),

  extract_h3("M9",   "bracl",       "NAVCO 1.3", "success3 (0/1/2) ordered",    model9,   "repress_c", s_m9$N,   s_m9$n_events),
  extract_h3("M10",  "bracl",       "NAVCO 1.3", "success3 (0/1/2) ordered",    model10,  "civlib_c",  s_m10$N,  s_m10$n_events),
  extract_h3("M11",  "bracl",       "NAVCO 1.3", "success3 (0/1/2) ordered",    model11,  "clpol_c",    s_m11$N,  s_m11$n_events),
  extract_h3("M12",  "bracl",       "NAVCO 1.3", "success3 (0/1/2) ordered",    model12,  "polity_c",  s_m12$N,  s_m12$n_events),

  extract_h3("M13",  "Cox PH",      "NAVCO 2.1", "Hazard of first concession (CY)", model13,  "repress_c", s_m13$N,  s_m13$n_events),
  extract_h3("M13s", "Cox PH",      "NAVCO 2.1", "Hazard of success (CY, strict)", model13s, "repress_c", s_m13s$N, s_m13s$n_events),
  extract_h3("M14",  "Cox PH",      "NAVCO 2.1", "Hazard of first concession (CY)", model14,  "civlib_c",  s_m14$N,  s_m14$n_events),
  extract_h3("M14s", "Cox PH",      "NAVCO 2.1", "Hazard of success (CY, strict)", model14s, "civlib_c",  s_m14s$N, s_m14s$n_events),
  extract_h3("M15",  "Cox PH",      "NAVCO 2.1", "Hazard of first concession (CY)", model15,  "clpol_c",    s_m15$N,  s_m15$n_events),
  extract_h3("M15s", "Cox PH",      "NAVCO 2.1", "Hazard of success (CY, strict)", model15s, "clpol_c",    s_m15s$N, s_m15s$n_events),
  extract_h3("M16",  "Cox PH",      "NAVCO 2.1", "Hazard of first concession (CY)", model16,  "polity_c",  s_m16$N,  s_m16$n_events),
  extract_h3("M16s", "Cox PH",      "NAVCO 2.1", "Hazard of success (CY, strict)", model16s, "polity_c",  s_m16s$N, s_m16s$n_events),

  extract_h3("M17",  "Firth PML",   "NAVCO 2.1", "Yearly concession (CY, binary)", model17,  "repress_c", s_m17$N,  s_m17$n_events),
  extract_h3("M17s", "Firth PML",   "NAVCO 2.1", "Yearly success (CY, strict bin.)", model17s, "repress_c", s_m17s$N, s_m17s$n_events),
  extract_h3("M18",  "Firth PML",   "NAVCO 2.1", "Yearly concession (CY, binary)", model18,  "civlib_c",  s_m18$N,  s_m18$n_events),
  extract_h3("M18s", "Firth PML",   "NAVCO 2.1", "Yearly success (CY, strict bin.)", model18s, "civlib_c",  s_m18s$N, s_m18s$n_events),
  extract_h3("M19",  "Firth PML",   "NAVCO 2.1", "Yearly concession (CY, binary)", model19,  "clpol_c",    s_m19$N,  s_m19$n_events),
  extract_h3("M19s", "Firth PML",   "NAVCO 2.1", "Yearly success (CY, strict bin.)", model19s, "clpol_c",    s_m19s$N, s_m19s$n_events),
  extract_h3("M20",  "Firth PML",   "NAVCO 2.1", "Yearly concession (CY, binary)", model20,  "polity_c",  s_m20$N,  s_m20$n_events),
  extract_h3("M20s", "Firth PML",   "NAVCO 2.1", "Yearly success (CY, strict bin.)", model20s, "polity_c",  s_m20s$N, s_m20s$n_events),

  extract_h3("M21",  "bracl",       "NAVCO 2.1", "progress_ord (6-level; failure separate)", model21, "repress_c", s_m21$N, s_m21$n_events),
  extract_h3("M22",  "bracl",       "NAVCO 2.1", "progress_ord (6-level; failure separate)", model22, "civlib_c",  s_m22$N, s_m22$n_events),
  extract_h3("M23",  "bracl",       "NAVCO 2.1", "progress_ord (6-level; failure separate)", model23, "clpol_c",    s_m23$N, s_m23$n_events),
  extract_h3("M24",  "bracl",       "NAVCO 2.1", "progress_ord (6-level; failure separate)", model24, "polity_c",  s_m24$N, s_m24$n_events)
)

###############################################################################
# MULTIPLICITY + ROBUSTNESS ROWS
###############################################################################

main_models <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")
core_results$main_panel <- core_results$model %in% main_models

core_results$p_int_holm <- NA_real_
core_results$p_int_bh   <- NA_real_

main_idx <- which(core_results$main_panel)
app_idx  <- which(!core_results$main_panel)

if (length(main_idx) > 0) {
  core_results$p_int_holm[main_idx] <- p.adjust(
    core_results$p_int[main_idx], method = "holm", n = length(main_idx)
  )
}
if (length(app_idx) > 0) {
  core_results$p_int_bh[app_idx] <- p.adjust(
    core_results$p_int[app_idx], method = "BH", n = length(app_idx)
  )
}

core_results$sig_05_int_holm <- !is.na(core_results$p_int_holm) & core_results$p_int_holm < 0.05
core_results$sig_10_int_holm <- !is.na(core_results$p_int_holm) & core_results$p_int_holm < 0.10
core_results$sig_05_int_bh   <- !is.na(core_results$p_int_bh)   & core_results$p_int_bh   < 0.05
core_results$sig_10_int_bh   <- !is.na(core_results$p_int_bh)   & core_results$p_int_bh   < 0.10
core_results$robustness_row  <- FALSE

robustness_rows <- rbind(
  extract_h3("M13_AG", "Cox PH", "NAVCO 2.1", "Hazard of concession (CY, recurrent-event)", model13_ag, "repress_c", s_m13_ag$N, s_m13_ag$n_events),
  extract_h3("M14_AG", "Cox PH", "NAVCO 2.1", "Hazard of concession (CY, recurrent-event)", model14_ag, "civlib_c",  s_m14_ag$N, s_m14_ag$n_events),
  extract_h3("M15_AG", "Cox PH", "NAVCO 2.1", "Hazard of concession (CY, recurrent-event)", model15_ag, "clpol_c",    s_m15_ag$N, s_m15_ag$n_events),
  extract_h3("M16_AG", "Cox PH", "NAVCO 2.1", "Hazard of concession (CY, recurrent-event)", model16_ag, "polity_c",  s_m16_ag$N, s_m16_ag$n_events),
  extract_h3("M17_GEE", "GEE logit", "NAVCO 2.1", "Yearly concession (CY, source-clustered GEE)", model17_gee, "repress_c", s_m17_gee$N, s_m17_gee$n_events),
  extract_h3("M18_GEE", "GEE logit", "NAVCO 2.1", "Yearly concession (CY, source-clustered GEE)", model18_gee, "civlib_c",  s_m18_gee$N, s_m18_gee$n_events),
  extract_h3("M19_GEE", "GEE logit", "NAVCO 2.1", "Yearly concession (CY, source-clustered GEE)", model19_gee, "clpol_c",    s_m19_gee$N, s_m19_gee$n_events),
  extract_h3("M20_GEE", "GEE logit", "NAVCO 2.1", "Yearly concession (CY, source-clustered GEE)", model20_gee, "polity_c",  s_m20_gee$N, s_m20_gee$n_events)
)
robustness_rows$main_panel      <- FALSE
robustness_rows$p_int_holm      <- NA_real_
robustness_rows$p_int_bh        <- NA_real_
robustness_rows$sig_05_int_holm <- NA
robustness_rows$sig_10_int_holm <- NA
robustness_rows$sig_05_int_bh   <- NA
robustness_rows$sig_10_int_bh   <- NA
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
      ph_flag == "PH_FAIL_INTERACTION" ~ "PH interaction failure; appendix-only; do not rely substantively",
      ph_flag == "PH_FAIL_GLOBAL" ~ "PH global failure; interpret cautiously",
      ph_flag == "PH_TEST_FAILED" ~ "PH test failed; interpret cautiously",
      ph_flag == "PH_FOCAL_MISSING" ~ "PH focal-term row missing from cox.zph output",
      TRUE ~ ""
    )
  )

###############################################################################
# PRINT
###############################################################################

fmt_p <- function(p) ifelse(is.na(p), "   NA ", sprintf("%6.3f", p))
fmt_b <- function(b) ifelse(is.na(b), "    NA ", sprintf("%7.2f", b))

cat("\n=========== H3 MAIN PANEL (Table X, 8 models, Holm-adjusted k=8) ===========\n")
cat(sprintf("%-5s %-12s %-10s %-11s %8s %7s %8s %8s %4s %5s %5s %4s %3s\n",
            "Model", "Estimator", "Data", "Moderator",
            "b_int", "p_int", "p_Holm", "p_JW",
            "Exp", "Sign", "N", "EPV", "OK?"))
cat(paste(rep("-", 115), collapse = ""), "\n")
main_rows <- core_results[core_results$main_panel, ]
for (i in seq_len(nrow(main_rows))) {
  r <- main_rows[i, ]
  cat(sprintf("%-5s %-12s %-10s %-11s %8s %7s %8s %8s %4s %5s %5d %4s %3s\n",
              r$model, r$estimator, r$data, r$moderator,
              fmt_b(r$beta_int), fmt_p(r$p_int), fmt_p(r$p_int_holm), fmt_p(r$joint_wald_p),
              r$expected_int_sign,
              ifelse(isTRUE(r$sign_consistent), "cons", "WRONG"),
              r$N,
              ifelse(is.na(r$epv), "  NA", sprintf("%4.1f", r$epv)),
              ifelse(r$epv_warning, "!", " ")))
}

cat("\n=========== H3 APPENDIX (32 core rows, BH-adjusted k=32) ====================\n")
cat(sprintf("%-5s %-12s %-10s %-11s %8s %7s %8s %8s %4s %5s %5s %4s %3s\n",
            "Model", "Estimator", "Data", "Moderator",
            "b_int", "p_int", "p_BH", "p_JW",
            "Exp", "Sign", "N", "EPV", "OK?"))
cat(paste(rep("-", 115), collapse = ""), "\n")
app_rows <- core_results[!core_results$main_panel, ]
for (i in seq_len(nrow(app_rows))) {
  r <- app_rows[i, ]
  cat(sprintf("%-5s %-12s %-10s %-11s %8s %7s %8s %8s %4s %5s %5d %4s %3s\n",
              r$model, r$estimator, r$data, r$moderator,
              fmt_b(r$beta_int), fmt_p(r$p_int), fmt_p(r$p_int_bh), fmt_p(r$joint_wald_p),
              r$expected_int_sign,
              ifelse(isTRUE(r$sign_consistent), "cons", "WRONG"),
              r$N,
              ifelse(is.na(r$epv), "  NA", sprintf("%4.1f", r$epv)),
              ifelse(r$epv_warning, "!", " ")))
}

cat("\nNotes:\n")
cat("  p_Holm / p_BH  = multiplicity-adjusted interaction p-values for core rows only\n")
cat("  M13_AG-M16_AG  = recurrent-event Cox concession robustness rows\n")
cat("  M17_GEE-M20_GEE = source-clustered GEE concession-binary companions for repeated source-years\n")
cat("  M17s-M20s and M21-M24 remain appendix-only pooled yearly models and should be read as weaker evidence\n")
cat("  ph_flag        = exported cox.zph() diagnostic flag for core Cox rows\n")

robust_rows <- results[results$robustness_row, ]
if (nrow(robust_rows) > 0) {
  cat("\nRecurrent-event / clustered robustness rows:\n")
  print(robust_rows[, c("model", "response", "moderator", "p_int", "joint_wald_p", "N", "n_events")], row.names = FALSE)
}

flagged_ph <- results[results$ph_flag != "", c("model", "response", "ph_global_p", "ph_int_p", "ph_flag", "interpretation_note")]
if (nrow(flagged_ph) > 0) {
  cat("\nCox PH flags:\n")
  print(flagged_ph, row.names = FALSE)
}

cat("\n\n==== FULL MODEL SUMMARIES ====\n")

print_summary_block <- function(tag, label, mdl) {
  cat(sprintf("\n---- %s %s ----\n", tag, label))
  print(summary(mdl))
}

print_summary_block("M1",   "Weibull AFT x Repression — time to concession", model1)
print_summary_block("M1s",  "Weibull AFT x Repression — time to strict success", model1s)
print_summary_block("M2",   "Weibull AFT x CivLib — time to concession", model2)
print_summary_block("M2s",  "Weibull AFT x CivLib — time to strict success", model2s)
print_summary_block("M3",   "Weibull AFT x Political Liberty — time to concession", model3)
print_summary_block("M3s",  "Weibull AFT x Political Liberty — time to strict success", model3s)
print_summary_block("M4",   "Weibull AFT x Polity — time to concession", model4)
print_summary_block("M4s",  "Weibull AFT x Polity — time to strict success", model4s)

print_summary_block("M5",   "Firth PML x Repression — concession (binary)", model5)
print_summary_block("M5s",  "Firth PML x Repression — strict success (binary)", model5s)
print_summary_block("M6",   "Firth PML x CivLib — concession (binary)", model6)
print_summary_block("M6s",  "Firth PML x CivLib — strict success (binary)", model6s)
print_summary_block("M7",   "Firth PML x Political Liberty — concession (binary)", model7)
print_summary_block("M7s",  "Firth PML x Political Liberty — strict success (binary)", model7s)
print_summary_block("M8",   "Firth PML x Polity — concession (binary)", model8)
print_summary_block("M8s",  "Firth PML x Polity — strict success (binary)", model8s)

print_summary_block("M9",   "bracl x Repression — success3 (ordered)", model9)
print_summary_block("M10",  "bracl x CivLib — success3 (ordered)", model10)
print_summary_block("M11",  "bracl x Political Liberty — success3 (ordered)", model11)
print_summary_block("M12",  "bracl x Polity — success3 (ordered)", model12)

print_summary_block("M13",  "Cox PH x Repression — hazard of first concession (CY)", model13)
print_summary_block("M13s", "Cox PH x Repression — hazard of strict success (CY)", model13s)
print_summary_block("M14",  "Cox PH x CivLib — hazard of first concession (CY)", model14)
print_summary_block("M14s", "Cox PH x CivLib — hazard of strict success (CY)", model14s)
print_summary_block("M15",  "Cox PH x Political Liberty — hazard of first concession (CY)", model15)
print_summary_block("M15s", "Cox PH x Political Liberty — hazard of strict success (CY)", model15s)
print_summary_block("M16",  "Cox PH x Polity — hazard of first concession (CY)", model16)
print_summary_block("M16s", "Cox PH x Polity — hazard of strict success (CY)", model16s)

print_summary_block("M13_AG", "Cox PH x Repression — recurrent-event concession hazard (CY)", model13_ag)
print_summary_block("M14_AG", "Cox PH x CivLib — recurrent-event concession hazard (CY)", model14_ag)
print_summary_block("M15_AG", "Cox PH x Political Liberty — recurrent-event concession hazard (CY)", model15_ag)
print_summary_block("M16_AG", "Cox PH x Polity — recurrent-event concession hazard (CY)", model16_ag)
print_summary_block("M17_GEE", "GEE logit x Repression — yearly concession, source-clustered", model17_gee)
print_summary_block("M18_GEE", "GEE logit x CivLib — yearly concession, source-clustered", model18_gee)
print_summary_block("M19_GEE", "GEE logit x Political Liberty — yearly concession, source-clustered", model19_gee)
print_summary_block("M20_GEE", "GEE logit x Polity — yearly concession, source-clustered", model20_gee)

print_summary_block("M17",  "Firth PML x Repression — yearly concession", model17)
print_summary_block("M17s", "Firth PML x Repression — yearly strict success", model17s)
print_summary_block("M18",  "Firth PML x CivLib — yearly concession", model18)
print_summary_block("M18s", "Firth PML x CivLib — yearly strict success", model18s)
print_summary_block("M19",  "Firth PML x Political Liberty — yearly concession", model19)
print_summary_block("M19s", "Firth PML x Political Liberty — yearly strict success", model19s)
print_summary_block("M20",  "Firth PML x Polity — yearly concession", model20)
print_summary_block("M20s", "Firth PML x Polity — yearly strict success", model20s)

print_summary_block("M21",  "bracl x Repression — progress_ord (6-level; failure separate)", model21)
print_summary_block("M22",  "bracl x CivLib — progress_ord (6-level; failure separate)", model22)
print_summary_block("M23",  "bracl x Political Liberty — progress_ord (6-level; failure separate)", model23)
print_summary_block("M24",  "bracl x Polity — progress_ord (6-level; failure separate)", model24)

expected_result_cols <- c(
  "model", "estimator", "data", "response", "moderator", "mod_direction",
  "beta_eng_c", "p_eng_c", "beta_mod_c", "p_mod_c", "beta_int", "p_int",
  "joint_wald_p", "expected_int_sign", "sign_consistent", "N", "n_events",
  "epv", "sig_05_int", "sig_10_int", "epv_warning", "main_panel",
  "p_int_holm", "p_int_bh", "sig_05_int_holm", "sig_10_int_holm",
  "sig_05_int_bh", "sig_10_int_bh", "robustness_row",
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
if (nrow(results) != 48) {
  stop(sprintf("Expected 48 H3 result rows, found %d.", nrow(results)), call. = FALSE)
}
if (anyDuplicated(results$model) > 0) {
  stop("Duplicate model labels detected in H3 results.", call. = FALSE)
}

###############################################################################
# EXPORT
###############################################################################

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
write.csv(results, out_path, row.names = FALSE)
models <- list(
  M1 = model1, M1s = model1s, M2 = model2, M2s = model2s,
  M3 = model3, M3s = model3s, M4 = model4, M4s = model4s,
  M5 = model5, M5s = model5s, M6 = model6, M6s = model6s,
  M7 = model7, M7s = model7s, M8 = model8, M8s = model8s,
  M9 = model9, M10 = model10, M11 = model11, M12 = model12,
  M13 = model13, M13s = model13s, M14 = model14, M14s = model14s,
  M15 = model15, M15s = model15s, M16 = model16, M16s = model16s,
  M17 = model17, M17s = model17s, M18 = model18, M18s = model18s,
  M19 = model19, M19s = model19s, M20 = model20, M20s = model20s,
  M21 = model21, M22 = model22, M23 = model23, M24 = model24,
  M13_AG = model13_ag, M14_AG = model14_ag,
  M15_AG = model15_ag, M16_AG = model16_ag,
  M17_GEE = model17_gee, M18_GEE = model18_gee,
  M19_GEE = model19_gee, M20_GEE = model20_gee
)
gee_ids <- c("M17_GEE", "M18_GEE", "M19_GEE", "M20_GEE")
gee_status <- do.call(rbind, lapply(gee_ids, function(id) {
  fit <- models[[id]]
  err <- if (!is.null(fit$geese$error)) as.integer(fit$geese$error) else NA_integer_
  tab <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
  se_col <- if (is.null(tab)) character() else intersect(c("Std. Error", "Std.err"), colnames(tab))
  finite_output <- !is.null(tab) && length(se_col) == 1L &&
    all(is.finite(as.matrix(tab[, c("Estimate", se_col[[1L]]), drop = FALSE])))
  converged <- identical(err, 0L) && isTRUE(finite_output)
  data.frame(
    model = id,
    geese_error = err,
    finite_output = finite_output,
    converged = converged,
    accepted_values_preserved = TRUE,
    reporting_class = if (converged) "accepted_replicated" else "accepted_archival_nonconverged",
    stringsAsFactors = FALSE
  )
}))
write.csv(gee_status, file.path(P$diagnostics, "H3_gee_fit_status.csv"), row.names = FALSE)
if (any(!gee_status$converged)) {
  warning(
    sprintf(
      "Accepted H3 GEE values preserved for fidelity, but non-converged model(s) are archival-only: %s",
      paste(gee_status$model[!gee_status$converged], collapse = ", ")
    ),
    call. = FALSE
  )
}
stopifnot(identical(names(models), results$model))
saveRDS(models, cache_path, version = 3)
saveRDS(list(df = df, d21 = d21, d21_first_conc = d21_first_conc), frames_path, version = 3)
cat(sprintf("\nResults exported to: %s\n", out_path))
cat(sprintf("Model bundle exported to: %s\n", cache_path))
cat(sprintf("Analysis frames exported to: %s\n", frames_path))
