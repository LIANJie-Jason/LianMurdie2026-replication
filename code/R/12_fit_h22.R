###############################################################################
# H22_code.r
# H2.2: ECONOMIC VULNERABILITY TO EXTERNAL PRESSURE
#
# Tests whether English-language protest signs have stronger positive effects
# on campaign outcomes when the state is more economically vulnerable to
# external pressure — operationalized by four economic-dependence moderators:
#   wdi_fdiin: FDI net inflows (% of GDP)
#   fi_ftradeint_pd: Fraser Institute trade freedom / openness index
#   lnaid: log(bilateral + IO aid + 1), NAVCO 1.3
#                        (constructed in data_clean.r with NA preservation;
#                        read from df_final.csv here)
#   lnaid_yr: pre-computed annual log-aid, NAVCO 2.1 panel
#                        (constructed with NA preservation)
#   aid_crnio: number of IO donors to the country
#
# Theoretical logic:
# English signs mobilize international audiences and leverage external
# economic pressure. Leverage should bite harder on states that depend more
# on foreign capital (FDI), open trade (trade freedom), foreign aid (lnaid),
# or multilateral aid relationships (aid_crnio). Higher dependence = greater
# vulnerability = English signs mobilize that pressure more effectively.
#
# MODERATOR DIRECTION (differs from H2.1 sovereignty):
#   For H2.2, HIGHER moderator = MORE vulnerable = English helps MORE.
#   For H2.1, LOWER sovereignty = more vulnerable.
#   Consequently H2.2 expected interaction signs FLIP relative to H2.1:
#     Weibull AFT (benefit = NEGATIVE eng) → interaction sign NEGATIVE
#     Firth PML   (benefit = POSITIVE eng) → interaction sign POSITIVE
#     Cox PH      (benefit = POSITIVE eng) → interaction sign POSITIVE
#     bracl       (benefit = NEGATIVE eng) → interaction sign NEGATIVE
#   This is encoded via extract_h22(..., mod_direction = "pos").
#
# SPECIFICATION (parsimonious with an economic-confound adjustment):
#   - eng and moderators mean-centered
#   - NAVCO 1.3 controls: nonviolent_camp (1=nonviolent) + REGCHANGE + v2csreprss + lnnum_image_sum + lngdp
#       (+ lnlengthofcam for binary/ordered DV models)
#   - NAVCO 2.1 controls: nonviolent + REGCHANGE + v2csreprss + lnnum_image_sum + lngdp_yr
#       (+ time_in_campaign for binary/ordered DV models)
#   - lnpop / colonized_english dropped for H2 family parsimony
#   - lngdp / lngdp_yr INCLUDED: FDI%, trade openness, aid/GDP mechanically vary
#       with GDP size; omitting lngdp confounds the moderator with state size
#   - Cox PH uses counting-process format (Andersen & Gill 1982)
#
# DUAL-DV STRUCTURE:
#   Main manuscript "campaign success" is broad NAVCO success: full success OR
#   concession. In NAVCO 1.3 this is success OR limited; in NAVCO 2.1 this is
#   progress 2-4 (limited concession, significant concession, complete
#   success). Internal variable names retain the older concession shorthand
#   (`event_concession`, `concession_bin`) for continuity. Strict success means
#   full/complete success only. Ordered-logit models (success3 / progress_ord)
#   include both outcome levels by construction.
#
# MODELS (40 core + 4 robustness rows):
#   NAVCO 1.3 (N=101 max, per-model CC may be smaller):
#     M1/M1s   Weibull AFT × FDI     — time to broad success / strict success
#     M2/M2s   Weibull AFT × Trade   — same
#     M3/M3s   Weibull AFT × Aid     — same
#     M4/M4s   Weibull AFT × IODonor — same
#     M5/M5s   Firth PML   × FDI     — broad success / strict success (binary)
#     M6/M6s   Firth PML   × Trade   — same
#     M7/M7s   Firth PML   × Aid     — same
#     M8/M8s   Firth PML   × IODonor — same
#     M9       bracl       × FDI     — success3 (ordered 0/1/2)
#     M10      bracl       × Trade   — same
#     M11      bracl       × Aid     — same
#     M12      bracl       × IODonor — same
#   NAVCO 2.1 (source-campaign-year panel):
#     M13/M13s Cox PH      × FDI     — hazard of first broad success / strict success
#     M14/M14s Cox PH      × Trade   — same
#     M15/M15s Cox PH      × Aid     — same
#     M16/M16s Cox PH      × IODonor — same
#     M17/M17s Firth PML   × FDI     — yearly broad success / strict (binary)
#     M18/M18s Firth PML   × Trade   — same
#     M19/M19s Firth PML   × Aid     — same
#     M20/M20s Firth PML   × IODonor — same
#     M21      bracl       × FDI     — progress_ord (ordered 0-4)
#     M22      bracl       × Trade   — same
#     M23      bracl       × Aid     — same
#     M24      bracl       × IODonor — same
#
# MAIN PANEL SCOPE / MULTIPLICITY:
# Table 3's main panel is M1–M8 (NAVCO 1.3
# broad-success DV across 4 economic moderators × {Weibull AFT, Firth PML}).
# Appendix = the remaining 32 rows. Holm is applied within the 8-row main
# panel and BH within the 32-row appendix for p_int only. joint_wald_p is
# still reported raw in this script and is labeled as such in the compact
# console output.
#
# OUTPUT:
#   RR_datacode/empirical clean/results/H22_results.csv — 48-row compact headline table
#   (40 core rows + 8 robustness rows:
#    4 recurrent-event Cox + 4 clustered-GEE yearly broad-success checks)
#
# Author: Jie (Jason) Lian
# Last updated: 2026-04-14
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
out_path     <- file.path(P$estimates, "H22_results.csv")
cache_path   <- file.path(P$cache, "H22_models.rds")
frames_path  <- file.path(P$cache, "H22_analysis_frames.rds")

df  <- read.csv(navco13_path)
d21 <- read.csv(navco21_path)
rep_assert_columns(
  df,
  c(
    "CAMPAIGN", "LOCATION", "EYEAR", "BYEAR", "success", "limited", "success3",
    "viol", "nonviolent_camp", "REGCHANGE", "v2csreprss", "eng_prob_general", "lnnum_image_sum",
    "lnlengthofcam", "v2svdomaut", "v2svstterr", "wdi_fdiin", "fi_ftradeint_pd",
    "aid_crsc", "aid_crsio", "aid_crnio", "lnaid", "lngdp"
  ),
  "df_final.csv"
)
rep_assert_columns(
  d21,
  c(
    "navco21_id", "Year", "CAMPAIGN", "LOCATION", "panel_campaign_id",
    "nonviolent", "camp_goals", "Num_Image", "progress", "concession",
    "time_in_campaign", "eng_prop_year", "v2csreprss", "v2svdomaut",
    "v2svstterr", "wdi_fdiin", "fi_ftradeint_pd", "lnaid_yr", "aid_crnio",
    "lngdp_yr"
  ),
  "df_navco21_panel.csv"
)
rep_assert_unique(df, c("CAMPAIGN", "LOCATION"), "df_final CAMPAIGN-LOCATION")
rep_assert_unique(d21, c("navco21_id", "Year"), "df_navco21_panel NAVCO 2.1 source-id-Year")

###############################################################################
# NAVCO 1.3 PREPARATION
###############################################################################

df$viol      <- as.factor(df$viol)
df$nonviolent_camp <- as.factor(df$nonviolent_camp)  # harmonized polarity (1=nonviolent)
df$success3  <- ordered(df$success3)
df$REGCHANGE <- as.numeric(as.character(df$REGCHANGE))

# Survival outcomes. Broad NAVCO success equals full success OR concession;
# the internal "concession" variable names are legacy shorthand.
df$duration         <- df$EYEAR - df$BYEAR + 1
df$event_concession <- as.integer(
  as.numeric(as.character(df$success)) == 1 |
  as.numeric(as.character(df$limited)) == 1
)
df$event_success    <- as.numeric(as.character(df$success))

# Binary outcomes
df$concession_bin <- df$event_concession
df$success_bin    <- df$event_success

# lnaid is constructed with an NA-preserving convention and read from
# df_final.csv.
assert_no_nonfinite(df$eng_prob_general, "df$eng_prob_general")
assert_no_nonfinite(df$lnnum_image_sum, "df$lnnum_image_sum")
assert_no_nonfinite(df$lnaid, "df$lnaid")

# Mean-center interaction constituents.
df$eng_c     <- df$eng_prob_general - mean(df$eng_prob_general, na.rm = TRUE)
df$fdiin_c   <- df$wdi_fdiin        - mean(df$wdi_fdiin,        na.rm = TRUE)
df$trade_c   <- df$fi_ftradeint_pd  - mean(df$fi_ftradeint_pd,  na.rm = TRUE)
df$lnaid_c   <- df$lnaid            - mean(df$lnaid,            na.rm = TRUE)
df$iodonor_c <- df$aid_crnio        - mean(df$aid_crnio,        na.rm = TRUE)

cat(sprintf(
  "NAVCO 1.3 centering: mean(eng)=%.4f, mean(fdi)=%.4f, mean(trade)=%.4f, mean(lnaid)=%.4f, mean(iodonor)=%.4f\n",
  mean(df$eng_prob_general, na.rm = TRUE),
  mean(df$wdi_fdiin,        na.rm = TRUE),
  mean(df$fi_ftradeint_pd,  na.rm = TRUE),
  mean(df$lnaid,            na.rm = TRUE),
  mean(df$aid_crnio,        na.rm = TRUE)))

cat(sprintf("NAVCO 1.3 events (full sample): concession=%d, success=%d (N=%d)\n",
            sum(df$event_concession, na.rm = TRUE),
            sum(df$event_success,    na.rm = TRUE),
            nrow(df)))

###############################################################################
# NAVCO 2.1 PREPARATION
###############################################################################

d21$nonviolent      <- as.factor(d21$nonviolent)
d21$REGCHANGE       <- as.numeric(d21$camp_goals == 0)
d21$lnnum_image_sum <- log(d21$Num_Image)

# Per data_clean_navco21.r: use a 6-level ordered scale with
# failure < status quo < visible gains < limited concession <
# significant concession < complete success.
make_progress_ord6 <- function(x) {
  ordered(x, levels = c(5, 0, 1, 2, 3, 4))
}
d21$progress_ord     <- make_progress_ord6(d21$progress)
d21$event_concession <- as.numeric(d21$concession)        # broad success: progress 2-4
d21$event_success_yr <- as.numeric(d21$progress == 4)     # strict success
d21$success_bin_yr   <- d21$event_success_yr
d21$concession_bin   <- d21$event_concession
assert_no_nonfinite(d21$eng_prop_year, "d21$eng_prop_year")
assert_no_nonfinite(d21$lnnum_image_sum, "d21$lnnum_image_sum")
assert_no_nonfinite(d21$lnaid_yr, "d21$lnaid_yr")

# Counting-process intervals for Cox PH
d21$tstart <- d21$time_in_campaign - 1
d21$tstop  <- d21$time_in_campaign
assert_no_nonfinite(d21$tstart, "d21$tstart")
assert_no_nonfinite(d21$tstop, "d21$tstop")

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
d21$eng_c     <- d21$eng_prop_year   - mean(d21$eng_prop_year,   na.rm = TRUE)
d21$fdiin_c   <- d21$wdi_fdiin       - mean(d21$wdi_fdiin,       na.rm = TRUE)
d21$trade_c   <- d21$fi_ftradeint_pd - mean(d21$fi_ftradeint_pd, na.rm = TRUE)
d21$lnaid_c   <- d21$lnaid_yr        - mean(d21$lnaid_yr,        na.rm = TRUE)
d21$iodonor_c <- d21$aid_crnio       - mean(d21$aid_crnio,       na.rm = TRUE)

d21_first_conc <- build_first_event_panel(d21, "navco21_id", "event_concession", "time_in_campaign")
assert_first_event_panel(d21_first_conc, "navco21_id", "event_concession")

cat(sprintf(
  "NAVCO 2.1 centering: mean(eng)=%.4f, mean(fdi)=%.4f, mean(trade)=%.4f, mean(lnaid)=%.4f, mean(iodonor)=%.4f\n",
  mean(d21$eng_prop_year,   na.rm = TRUE),
  mean(d21$wdi_fdiin,       na.rm = TRUE),
  mean(d21$fi_ftradeint_pd, na.rm = TRUE),
  mean(d21$lnaid_yr,        na.rm = TRUE),
  mean(d21$aid_crnio,       na.rm = TRUE)))

cat(sprintf("NAVCO 2.1 events (full sample): concession=%d, success=%d (N=%d CY)\n",
            sum(d21$event_concession, na.rm = TRUE),
            sum(d21$event_success_yr, na.rm = TRUE),
            nrow(d21)))

###############################################################################
# FIT HELPERS — one per [dataset, estimator] combination
# Variable names are injected via sprintf into formula strings (avoids fragile
# non-standard evaluation). Every helper uses lngdp / lngdp_yr as an
# economic-confound control.
###############################################################################

fit_weibull13 <- function(dv_event, moderator) {
  f <- as.formula(sprintf(
    "Surv(duration, %s) ~ eng_c * %s + nonviolent_camp + REGCHANGE + v2csreprss + lnnum_image_sum + lngdp",
    dv_event, moderator))
  survreg(f, data = df, dist = "weibull")
}

fit_firth13 <- function(dv_col, moderator) {
  f <- as.formula(sprintf(
    "%s ~ eng_c * %s + nonviolent_camp + REGCHANGE + v2csreprss + lnlengthofcam + lnnum_image_sum + lngdp",
    dv_col, moderator))
  logistf(f, data = df,
          control   = logistf.control(maxit = 1000, maxstep = 0.5),
          plcontrol = logistpl.control(maxit = 1000),
          pl        = TRUE)
}

fit_bracl13 <- function(dv_col, moderator) {
  f <- as.formula(sprintf(
    "%s ~ eng_c * %s + nonviolent_camp + REGCHANGE + v2csreprss + lnlengthofcam + lnnum_image_sum + lngdp",
    dv_col, moderator))
  bracl(f, data = df,
        type = "MPL_Jeffreys", link = "logit", parallel = TRUE, maxit = 500)
}

fit_cox21 <- function(dv_event, moderator, data = d21) {
  f <- as.formula(sprintf(
    "Surv(tstart, tstop, %s) ~ eng_c * %s + nonviolent + REGCHANGE + v2csreprss + lnnum_image_sum + lngdp_yr + cluster(navco21_id)",
    dv_event, moderator))
  coxph(f, data = data)
}

fit_firth21 <- function(dv_col, moderator) {
  f <- as.formula(sprintf(
    "%s ~ eng_c * %s + nonviolent + REGCHANGE + v2csreprss + time_in_campaign + lnnum_image_sum + lngdp_yr",
    dv_col, moderator))
  logistf(f, data = d21,
          control   = logistf.control(maxit = 1000, maxstep = 0.5),
          plcontrol = logistpl.control(maxit = 1000),
          pl        = TRUE)
}

fit_gee21 <- function(dv_col, moderator) {
  vars <- c(
    dv_col, "eng_c", moderator, "nonviolent", "REGCHANGE",
    "v2csreprss", "time_in_campaign", "lnnum_image_sum", "lngdp_yr",
    "navco21_id"
  )
  dd <- d21[complete.cases(d21[, vars]), vars, drop = FALSE]
  dd <- dd %>% arrange(navco21_id, time_in_campaign)
  f <- as.formula(sprintf(
    "%s ~ eng_c * %s + nonviolent + REGCHANGE + v2csreprss + time_in_campaign + lnnum_image_sum + lngdp_yr",
    dv_col, moderator))
  geeglm(f, data = dd, id = navco21_id, family = binomial, corstr = "exchangeable")
}

fit_bracl21 <- function(dv_col, moderator) {
  f <- as.formula(sprintf(
    "%s ~ eng_c * %s + nonviolent + REGCHANGE + v2csreprss + time_in_campaign + lnnum_image_sum + lngdp_yr",
    dv_col, moderator))
  bracl(f, data = d21,
        type = "MPL_Jeffreys", link = "logit", parallel = TRUE, maxit = 500)
}

###############################################################################
# FIT ALL 40 MODELS
###############################################################################

# NAVCO 1.3 — 4 moderators × 5 model slots = 20 models
# Weibull AFT (concession + strict success)
model1   <- fit_weibull13("event_concession", "fdiin_c")
model1s  <- fit_weibull13("event_success",    "fdiin_c")
model2   <- fit_weibull13("event_concession", "trade_c")
model2s  <- fit_weibull13("event_success",    "trade_c")
model3   <- fit_weibull13("event_concession", "lnaid_c")
model3s  <- fit_weibull13("event_success",    "lnaid_c")
model4   <- fit_weibull13("event_concession", "iodonor_c")
model4s  <- fit_weibull13("event_success",    "iodonor_c")

# Firth PML (concession + strict success)
model5   <- fit_firth13("concession_bin", "fdiin_c")
model5s  <- fit_firth13("success_bin",    "fdiin_c")
model6   <- fit_firth13("concession_bin", "trade_c")
model6s  <- fit_firth13("success_bin",    "trade_c")
model7   <- fit_firth13("concession_bin", "lnaid_c")
model7s  <- fit_firth13("success_bin",    "lnaid_c")
model8   <- fit_firth13("concession_bin", "iodonor_c")
model8s  <- fit_firth13("success_bin",    "iodonor_c")

# bracl ordered (success3, 3-level)
model9   <- fit_bracl13("success3", "fdiin_c")
model10  <- fit_bracl13("success3", "trade_c")
model11  <- fit_bracl13("success3", "lnaid_c")
model12  <- fit_bracl13("success3", "iodonor_c")

# NAVCO 2.1 — 4 moderators × 5 model slots = 20 models
# Cox PH counting-process (first concession + strict success)
model13  <- fit_cox21("event_concession", "fdiin_c", data = d21_first_conc)
model13s <- fit_cox21("event_success_yr", "fdiin_c")
model14  <- fit_cox21("event_concession", "trade_c", data = d21_first_conc)
model14s <- fit_cox21("event_success_yr", "trade_c")
model15  <- fit_cox21("event_concession", "lnaid_c", data = d21_first_conc)
model15s <- fit_cox21("event_success_yr", "lnaid_c")
model16  <- fit_cox21("event_concession", "iodonor_c", data = d21_first_conc)
model16s <- fit_cox21("event_success_yr", "iodonor_c")

# Firth PML (yearly concession + strict success)
model17  <- fit_firth21("concession_bin", "fdiin_c")
model17s <- fit_firth21("success_bin_yr", "fdiin_c")
model18  <- fit_firth21("concession_bin", "trade_c")
model18s <- fit_firth21("success_bin_yr", "trade_c")
model19  <- fit_firth21("concession_bin", "lnaid_c")
model19s <- fit_firth21("success_bin_yr", "lnaid_c")
model20  <- fit_firth21("concession_bin", "iodonor_c")
model20s <- fit_firth21("success_bin_yr", "iodonor_c")
model17_gee <- fit_gee21("concession_bin", "fdiin_c")
model18_gee <- fit_gee21("concession_bin", "trade_c")
model19_gee <- fit_gee21("concession_bin", "lnaid_c")
model20_gee <- fit_gee21("concession_bin", "iodonor_c")

# bracl ordered (progress_ord, 5-level)
model21  <- fit_bracl21("progress_ord", "fdiin_c")
model22  <- fit_bracl21("progress_ord", "trade_c")
model23  <- fit_bracl21("progress_ord", "lnaid_c")
model24  <- fit_bracl21("progress_ord", "iodonor_c")

model13_ag <- fit_cox21("event_concession", "fdiin_c", data = d21)
model14_ag <- fit_cox21("event_concession", "trade_c", data = d21)
model15_ag <- fit_cox21("event_concession", "lnaid_c", data = d21)
model16_ag <- fit_cox21("event_concession", "iodonor_c", data = d21)

###############################################################################
# COX PH DIAGNOSTIC — cox.zph proportional hazards test
#
# For each NAVCO 2.1 Cox model (M13–M16s), test the proportional hazards
# assumption using scaled Schoenfeld residuals. Reports global p-value and
# per-term p-values, with emphasis on the eng_c × moderator_c interaction
# term (the key causal parameter). The diagnostic is exported into the results
# table so appendix-level Cox rows can be flagged when PH fails. Under the
# a focal-term PH failure is treated as appendix-only weak evidence.
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

cat("\n--- Cox PH diagnostic (cox.zph scaled Schoenfeld residuals) ---\n")
cox_diag_tbl <- bind_rows(
  cox_ph_diag("M13",  model13,  "eng_c:fdiin_c"),
  cox_ph_diag("M13s", model13s, "eng_c:fdiin_c"),
  cox_ph_diag("M14",  model14,  "eng_c:trade_c"),
  cox_ph_diag("M14s", model14s, "eng_c:trade_c"),
  cox_ph_diag("M15",  model15,  "eng_c:lnaid_c"),
  cox_ph_diag("M15s", model15s, "eng_c:lnaid_c"),
  cox_ph_diag("M16",  model16,  "eng_c:iodonor_c"),
  cox_ph_diag("M16s", model16s, "eng_c:iodonor_c")
)
cat("--- end Cox PH diagnostic ---\n\n")

###############################################################################
# EXTRACTION HELPER — interaction headline per model
#
# mod_direction: "pos" = higher moderator → MORE vulnerable (H2.2 default)
#                "neg" = lower moderator  → more vulnerable (H2.1 convention)
# Expected-sign lookup flips accordingly.
###############################################################################

extract_h22 <- function(tag, estimator, data_label, response, model,
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
 # Cox PH with cluster reports robust SE and "Pr(>|z|)"
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

 # Joint Wald on (eng, int): H0: β_eng = β_int = 0 (df=2) — "no English effect
 # anywhere in the moderator range". The moderator main effect is a control
 # under the interaction hypothesis, not part of the joint test.
  joint_wald_p <- NA_real_
  V2 <- get_vcov_block(model, c(eng_name, int_name))
  if (!is.null(V2) && all(is.finite(c(b_eng, b_int))) && all(is.finite(V2))) {
    bvec <- c(b_eng, b_int)
    W    <- tryCatch(as.numeric(t(bvec) %*% solve(V2) %*% bvec),
                     error = function(e) NA_real_)
    if (is.finite(W)) joint_wald_p <- pchisq(W, df = 2, lower.tail = FALSE)
  }

 # Expected interaction sign under H2.2 (mod_direction = "pos"):
  #   Weibull AFT (benefit = negative eng) → interaction NEGATIVE
  #   Firth PML   (benefit = positive eng) → interaction POSITIVE
  #   Cox PH      (benefit = positive eng) → interaction POSITIVE
  #   bracl       (benefit = negative eng) → interaction NEGATIVE
 # Under H2.1 (mod_direction = "neg") signs flip — retained for parity with H21_code.r.
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
# PER-MODEL SAMPLE SIZES (complete-case, explicit to avoid class-specific
# nobs quirks for bracl). cc_stats returns BOTH N and n_events counted
# inside the complete-case mask — full-sample event counts would bias EPV
# when missingness differs across moderator/DV.
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

# NAVCO 1.3 covariate lists (lngdp included per B1)
vars_13_surv <- function(dv, moderator) {
  c("duration", dv, "eng_c", moderator, "nonviolent_camp",
    "REGCHANGE", "v2csreprss", "lnnum_image_sum", "lngdp")
}
vars_13_bin  <- function(dv, moderator) {
  c(dv, "eng_c", moderator, "nonviolent_camp",
    "REGCHANGE", "v2csreprss", "lnlengthofcam", "lnnum_image_sum", "lngdp")
}
# NAVCO 2.1 covariate lists (lngdp_yr included per B1)
vars_21_surv <- function(dv, moderator) {
  c("tstart", "tstop", dv, "eng_c", moderator, "nonviolent",
    "REGCHANGE", "v2csreprss", "lnnum_image_sum", "lngdp_yr", "navco21_id")
}
vars_21_bin  <- function(dv, moderator) {
  c(dv, "eng_c", moderator, "nonviolent",
    "REGCHANGE", "v2csreprss", "time_in_campaign", "lnnum_image_sum", "lngdp_yr")
}

# NAVCO 1.3 — Weibull AFT
s_m1  <- cc_stats(df, vars_13_surv("event_concession", "fdiin_c"),   "event_concession")
s_m1s <- cc_stats(df, vars_13_surv("event_success",    "fdiin_c"),   "event_success")
s_m2  <- cc_stats(df, vars_13_surv("event_concession", "trade_c"),   "event_concession")
s_m2s <- cc_stats(df, vars_13_surv("event_success",    "trade_c"),   "event_success")
s_m3  <- cc_stats(df, vars_13_surv("event_concession", "lnaid_c"),   "event_concession")
s_m3s <- cc_stats(df, vars_13_surv("event_success",    "lnaid_c"),   "event_success")
s_m4  <- cc_stats(df, vars_13_surv("event_concession", "iodonor_c"), "event_concession")
s_m4s <- cc_stats(df, vars_13_surv("event_success",    "iodonor_c"), "event_success")

# NAVCO 1.3 — Firth PML
s_m5  <- cc_stats(df, vars_13_bin("concession_bin", "fdiin_c"),   "concession_bin")
s_m5s <- cc_stats(df, vars_13_bin("success_bin",    "fdiin_c"),   "success_bin")
s_m6  <- cc_stats(df, vars_13_bin("concession_bin", "trade_c"),   "concession_bin")
s_m6s <- cc_stats(df, vars_13_bin("success_bin",    "trade_c"),   "success_bin")
s_m7  <- cc_stats(df, vars_13_bin("concession_bin", "lnaid_c"),   "concession_bin")
s_m7s <- cc_stats(df, vars_13_bin("success_bin",    "lnaid_c"),   "success_bin")
s_m8  <- cc_stats(df, vars_13_bin("concession_bin", "iodonor_c"), "concession_bin")
s_m8s <- cc_stats(df, vars_13_bin("success_bin",    "iodonor_c"), "success_bin")

# NAVCO 1.3 — bracl ordered (no event count)
s_m9  <- cc_stats(df, vars_13_bin("success3", "fdiin_c"))
s_m10 <- cc_stats(df, vars_13_bin("success3", "trade_c"))
s_m11 <- cc_stats(df, vars_13_bin("success3", "lnaid_c"))
s_m12 <- cc_stats(df, vars_13_bin("success3", "iodonor_c"))

# NAVCO 2.1 — Cox PH
s_m13  <- cc_stats(d21_first_conc, vars_21_surv("event_concession", "fdiin_c"),   "event_concession")
s_m13s <- cc_stats(d21, vars_21_surv("event_success_yr", "fdiin_c"),   "event_success_yr")
s_m14  <- cc_stats(d21_first_conc, vars_21_surv("event_concession", "trade_c"),   "event_concession")
s_m14s <- cc_stats(d21, vars_21_surv("event_success_yr", "trade_c"),   "event_success_yr")
s_m15  <- cc_stats(d21_first_conc, vars_21_surv("event_concession", "lnaid_c"),   "event_concession")
s_m15s <- cc_stats(d21, vars_21_surv("event_success_yr", "lnaid_c"),   "event_success_yr")
s_m16  <- cc_stats(d21_first_conc, vars_21_surv("event_concession", "iodonor_c"), "event_concession")
s_m16s <- cc_stats(d21, vars_21_surv("event_success_yr", "iodonor_c"), "event_success_yr")
s_m13_ag <- cc_stats(d21, vars_21_surv("event_concession", "fdiin_c"),   "event_concession")
s_m14_ag <- cc_stats(d21, vars_21_surv("event_concession", "trade_c"),   "event_concession")
s_m15_ag <- cc_stats(d21, vars_21_surv("event_concession", "lnaid_c"),   "event_concession")
s_m16_ag <- cc_stats(d21, vars_21_surv("event_concession", "iodonor_c"), "event_concession")

# NAVCO 2.1 — Firth PML
s_m17  <- cc_stats(d21, vars_21_bin("concession_bin", "fdiin_c"),   "concession_bin")
s_m17_gee <- cc_stats(d21, vars_21_bin("concession_bin", "fdiin_c"),   "concession_bin")
s_m17s <- cc_stats(d21, vars_21_bin("success_bin_yr", "fdiin_c"),   "success_bin_yr")
s_m18  <- cc_stats(d21, vars_21_bin("concession_bin", "trade_c"),   "concession_bin")
s_m18_gee <- cc_stats(d21, vars_21_bin("concession_bin", "trade_c"),   "concession_bin")
s_m18s <- cc_stats(d21, vars_21_bin("success_bin_yr", "trade_c"),   "success_bin_yr")
s_m19  <- cc_stats(d21, vars_21_bin("concession_bin", "lnaid_c"),   "concession_bin")
s_m19_gee <- cc_stats(d21, vars_21_bin("concession_bin", "lnaid_c"),   "concession_bin")
s_m19s <- cc_stats(d21, vars_21_bin("success_bin_yr", "lnaid_c"),   "success_bin_yr")
s_m20  <- cc_stats(d21, vars_21_bin("concession_bin", "iodonor_c"), "concession_bin")
s_m20_gee <- cc_stats(d21, vars_21_bin("concession_bin", "iodonor_c"), "concession_bin")
s_m20s <- cc_stats(d21, vars_21_bin("success_bin_yr", "iodonor_c"), "success_bin_yr")

# NAVCO 2.1 — bracl ordered
s_m21 <- cc_stats(d21, vars_21_bin("progress_ord", "fdiin_c"))
s_m22 <- cc_stats(d21, vars_21_bin("progress_ord", "trade_c"))
s_m23 <- cc_stats(d21, vars_21_bin("progress_ord", "lnaid_c"))
s_m24 <- cc_stats(d21, vars_21_bin("progress_ord", "iodonor_c"))

###############################################################################
# BUILD RESULTS TABLE (40 core rows, mod_direction = "pos" default)
###############################################################################

core_results <- rbind(
 # NAVCO 1.3 — Weibull AFT
  extract_h22("M1",   "Weibull AFT", "NAVCO 1.3", "Time to success",             model1,   "fdiin_c",   s_m1$N,   s_m1$n_events),
  extract_h22("M1s",  "Weibull AFT", "NAVCO 1.3", "Time to success (strict)",    model1s,  "fdiin_c",   s_m1s$N,  s_m1s$n_events),
  extract_h22("M2",   "Weibull AFT", "NAVCO 1.3", "Time to success",             model2,   "trade_c",   s_m2$N,   s_m2$n_events),
  extract_h22("M2s",  "Weibull AFT", "NAVCO 1.3", "Time to success (strict)",    model2s,  "trade_c",   s_m2s$N,  s_m2s$n_events),
  extract_h22("M3",   "Weibull AFT", "NAVCO 1.3", "Time to success",             model3,   "lnaid_c",   s_m3$N,   s_m3$n_events),
  extract_h22("M3s",  "Weibull AFT", "NAVCO 1.3", "Time to success (strict)",    model3s,  "lnaid_c",   s_m3s$N,  s_m3s$n_events),
  extract_h22("M4",   "Weibull AFT", "NAVCO 1.3", "Time to success",             model4,   "iodonor_c", s_m4$N,   s_m4$n_events),
  extract_h22("M4s",  "Weibull AFT", "NAVCO 1.3", "Time to success (strict)",    model4s,  "iodonor_c", s_m4s$N,  s_m4s$n_events),
 # NAVCO 1.3 — Firth PML
  extract_h22("M5",   "Firth PML",   "NAVCO 1.3", "Success (binary)",            model5,   "fdiin_c",   s_m5$N,   s_m5$n_events),
  extract_h22("M5s",  "Firth PML",   "NAVCO 1.3", "Success (strict, binary)",    model5s,  "fdiin_c",   s_m5s$N,  s_m5s$n_events),
  extract_h22("M6",   "Firth PML",   "NAVCO 1.3", "Success (binary)",            model6,   "trade_c",   s_m6$N,   s_m6$n_events),
  extract_h22("M6s",  "Firth PML",   "NAVCO 1.3", "Success (strict, binary)",    model6s,  "trade_c",   s_m6s$N,  s_m6s$n_events),
  extract_h22("M7",   "Firth PML",   "NAVCO 1.3", "Success (binary)",            model7,   "lnaid_c",   s_m7$N,   s_m7$n_events),
  extract_h22("M7s",  "Firth PML",   "NAVCO 1.3", "Success (strict, binary)",    model7s,  "lnaid_c",   s_m7s$N,  s_m7s$n_events),
  extract_h22("M8",   "Firth PML",   "NAVCO 1.3", "Success (binary)",            model8,   "iodonor_c", s_m8$N,   s_m8$n_events),
  extract_h22("M8s",  "Firth PML",   "NAVCO 1.3", "Success (strict, binary)",    model8s,  "iodonor_c", s_m8s$N,  s_m8s$n_events),
 # NAVCO 1.3 — bracl ordered
  extract_h22("M9",   "bracl",       "NAVCO 1.3", "success3 (0/1/2) ordered",    model9,   "fdiin_c",   s_m9$N,   s_m9$n_events),
  extract_h22("M10",  "bracl",       "NAVCO 1.3", "success3 (0/1/2) ordered",    model10,  "trade_c",   s_m10$N,  s_m10$n_events),
  extract_h22("M11",  "bracl",       "NAVCO 1.3", "success3 (0/1/2) ordered",    model11,  "lnaid_c",   s_m11$N,  s_m11$n_events),
  extract_h22("M12",  "bracl",       "NAVCO 1.3", "success3 (0/1/2) ordered",    model12,  "iodonor_c", s_m12$N,  s_m12$n_events),
 # NAVCO 2.1 — Cox PH
  extract_h22("M13",  "Cox PH",      "NAVCO 2.1", "Hazard of first success (CY)", model13,  "fdiin_c",   s_m13$N,  s_m13$n_events),
  extract_h22("M13s", "Cox PH",      "NAVCO 2.1", "Hazard of success (CY, strict)", model13s, "fdiin_c",   s_m13s$N, s_m13s$n_events),
  extract_h22("M14",  "Cox PH",      "NAVCO 2.1", "Hazard of first success (CY)", model14,  "trade_c",   s_m14$N,  s_m14$n_events),
  extract_h22("M14s", "Cox PH",      "NAVCO 2.1", "Hazard of success (CY, strict)", model14s, "trade_c",   s_m14s$N, s_m14s$n_events),
  extract_h22("M15",  "Cox PH",      "NAVCO 2.1", "Hazard of first success (CY)", model15,  "lnaid_c",   s_m15$N,  s_m15$n_events),
  extract_h22("M15s", "Cox PH",      "NAVCO 2.1", "Hazard of success (CY, strict)", model15s, "lnaid_c",   s_m15s$N, s_m15s$n_events),
  extract_h22("M16",  "Cox PH",      "NAVCO 2.1", "Hazard of first success (CY)", model16,  "iodonor_c", s_m16$N,  s_m16$n_events),
  extract_h22("M16s", "Cox PH",      "NAVCO 2.1", "Hazard of success (CY, strict)", model16s, "iodonor_c", s_m16s$N, s_m16s$n_events),
 # NAVCO 2.1 — Firth PML
  extract_h22("M17",  "Firth PML",   "NAVCO 2.1", "Yearly success (CY, binary)",      model17,  "fdiin_c",   s_m17$N,  s_m17$n_events),
  extract_h22("M17s", "Firth PML",   "NAVCO 2.1", "Yearly success (CY, strict bin.)", model17s, "fdiin_c",   s_m17s$N, s_m17s$n_events),
  extract_h22("M18",  "Firth PML",   "NAVCO 2.1", "Yearly success (CY, binary)",      model18,  "trade_c",   s_m18$N,  s_m18$n_events),
  extract_h22("M18s", "Firth PML",   "NAVCO 2.1", "Yearly success (CY, strict bin.)", model18s, "trade_c",   s_m18s$N, s_m18s$n_events),
  extract_h22("M19",  "Firth PML",   "NAVCO 2.1", "Yearly success (CY, binary)",      model19,  "lnaid_c",   s_m19$N,  s_m19$n_events),
  extract_h22("M19s", "Firth PML",   "NAVCO 2.1", "Yearly success (CY, strict bin.)", model19s, "lnaid_c",   s_m19s$N, s_m19s$n_events),
  extract_h22("M20",  "Firth PML",   "NAVCO 2.1", "Yearly success (CY, binary)",      model20,  "iodonor_c", s_m20$N,  s_m20$n_events),
  extract_h22("M20s", "Firth PML",   "NAVCO 2.1", "Yearly success (CY, strict bin.)", model20s, "iodonor_c", s_m20s$N, s_m20s$n_events),
 # NAVCO 2.1 — bracl ordered
  extract_h22("M21",  "bracl",       "NAVCO 2.1", "progress_ord (6-level; failure separate)", model21,  "fdiin_c",   s_m21$N,  s_m21$n_events),
  extract_h22("M22",  "bracl",       "NAVCO 2.1", "progress_ord (6-level; failure separate)", model22,  "trade_c",   s_m22$N,  s_m22$n_events),
  extract_h22("M23",  "bracl",       "NAVCO 2.1", "progress_ord (6-level; failure separate)", model23,  "lnaid_c",   s_m23$N,  s_m23$n_events),
  extract_h22("M24",  "bracl",       "NAVCO 2.1", "progress_ord (6-level; failure separate)", model24,  "iodonor_c", s_m24$N,  s_m24$n_events)
)

###############################################################################
# MAIN-PANEL SCOPE + MULTIPLICITY
#
# Table 3 of the main article contains M1–M8 (NAVCO 1.3
# broad-success DV across 4 economic moderators × {Weibull AFT, Firth PML}).
# Everything else (strict-success variants, success3/progress_ord ordered
# specs, full NAVCO 2.1 panel) → Appendix robustness.
#
# Correction convention:
#   Holm–Bonferroni within the main panel (k = 8, confirmatory family).
#   Benjamini–Hochberg within the appendix (k = 32, exploratory family).
#   Explicit n = length(main_idx)/length(app_idx) prevents p.adjust from
#   silently shrinking the family when any p_int is NA.
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
  extract_h22("M13_AG", "Cox PH", "NAVCO 2.1", "Hazard of recurrent campaign success (CY)", model13_ag, "fdiin_c",   s_m13_ag$N, s_m13_ag$n_events),
  extract_h22("M14_AG", "Cox PH", "NAVCO 2.1", "Hazard of recurrent campaign success (CY)", model14_ag, "trade_c",   s_m14_ag$N, s_m14_ag$n_events),
  extract_h22("M15_AG", "Cox PH", "NAVCO 2.1", "Hazard of recurrent campaign success (CY)", model15_ag, "lnaid_c",   s_m15_ag$N, s_m15_ag$n_events),
  extract_h22("M16_AG", "Cox PH", "NAVCO 2.1", "Hazard of recurrent campaign success (CY)", model16_ag, "iodonor_c", s_m16_ag$N, s_m16_ag$n_events),
  extract_h22("M17_GEE", "GEE logit", "NAVCO 2.1", "Yearly success (CY, source-clustered GEE)", model17_gee, "fdiin_c", s_m17_gee$N, s_m17_gee$n_events),
  extract_h22("M18_GEE", "GEE logit", "NAVCO 2.1", "Yearly success (CY, source-clustered GEE)", model18_gee, "trade_c", s_m18_gee$N, s_m18_gee$n_events),
  extract_h22("M19_GEE", "GEE logit", "NAVCO 2.1", "Yearly success (CY, source-clustered GEE)", model19_gee, "lnaid_c", s_m19_gee$N, s_m19_gee$n_events),
  extract_h22("M20_GEE", "GEE logit", "NAVCO 2.1", "Yearly success (CY, source-clustered GEE)", model20_gee, "iodonor_c", s_m20_gee$N, s_m20_gee$n_events)
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
    source_outcome = case_when(
      model %in% c("M1", "M2", "M3", "M4") ~
        "event_concession_broad_navco13_success_or_limited",
      model %in% c("M5", "M6", "M7", "M8") ~
        "concession_bin_broad_navco13_success_or_limited",
      model %in% c("M1s", "M2s", "M3s", "M4s") ~
        "event_success_strict_navco13_success",
      model %in% c("M5s", "M6s", "M7s", "M8s") ~
        "success_bin_strict_navco13_success",
      model %in% c("M9", "M10", "M11", "M12") ~
        "success3_ordered_navco13",
      model %in% c("M13", "M14", "M15", "M16",
                   "M13_AG", "M14_AG", "M15_AG", "M16_AG") ~
        "event_concession_broad_navco21_progress_2_4",
      model %in% c("M17", "M18", "M19", "M20",
                   "M17_GEE", "M18_GEE", "M19_GEE", "M20_GEE") ~
        "concession_bin_broad_navco21_progress_2_4",
      model %in% c("M13s", "M14s", "M15s", "M16s") ~
        "event_success_yr_strict_navco21_progress_4",
      model %in% c("M17s", "M18s", "M19s", "M20s") ~
        "success_bin_yr_strict_navco21_progress_4",
      model %in% c("M21", "M22", "M23", "M24") ~
        "progress_ord_ordered_navco21",
      TRUE ~ ""
    ),
    ph_global_p = ifelse(is.na(ph_global_p), NA_real_, ph_global_p),
    ph_int_p = ifelse(is.na(ph_int_p), NA_real_, ph_int_p),
    ph_focal_p = ifelse(is.na(ph_focal_p), NA_real_, ph_focal_p),
    ph_flag = ifelse(is.na(ph_flag), "", ph_flag),
    cox_estimand = case_when(
      estimator != "Cox PH" ~ "",
      grepl("_AG$", model) ~ "recurrent_event_success_broad",
      grepl("strict", response, ignore.case = TRUE) &
        grepl("success", response, ignore.case = TRUE) ~ "strict_success",
      grepl("first success", response, ignore.case = TRUE) ~ "first_success_broad",
      grepl("success", response, ignore.case = TRUE) ~ "success_broad",
      TRUE ~ "cox_unspecified"
    ),
    ph_flag = case_when(
      estimator == "Cox PH" & grepl("_AG$", model) ~ "PH_NOT_REPORTED_RECURRENT",
      TRUE ~ ph_flag
    ),
    interpretation_note = case_when(
      ph_flag == "PH_NOT_REPORTED_RECURRENT" ~
        "Recurrent-event broad-success Cox robustness estimand; first-event PH diagnostic not reported",
      ph_flag == "PH_FAIL_INTERACTION" ~ "PH interaction failure; appendix-only; do not rely substantively",
      ph_flag == "PH_FAIL_GLOBAL" ~ "PH global failure; interpret cautiously",
      ph_flag == "PH_TEST_FAILED" ~ "PH test failed; interpret cautiously",
      ph_flag == "PH_FOCAL_MISSING" ~ "PH focal-term row missing from cox.zph output",
      TRUE ~ ""
    )
  )

###############################################################################
# COMPACT PRINT — 40 CORE MODELS + 4 ROBUSTNESS ROWS
###############################################################################

fmt_p <- function(p) ifelse(is.na(p), "   NA ", sprintf("%6.3f", p))
fmt_b <- function(b) ifelse(is.na(b), "    NA ", sprintf("%7.2f", b))

cat("\n=========== H2.2 MAIN PANEL (Table 3, 8 models, Holm-adjusted k=8) ===========\n")
cat(sprintf("%-5s %-12s %-10s %-10s %8s %7s %8s %11s %4s %5s %5s %4s %3s\n",
            "Model", "Estimator", "Data", "Moderator",
            "b_int", "p_int", "p_Holm", "p_JW_raw",
            "Exp", "Sign", "N", "EPV", "OK?"))
cat(paste(rep("-", 115), collapse = ""), "\n")
main_rows <- core_results[core_results$main_panel, ]
for (i in seq_len(nrow(main_rows))) {
  r <- main_rows[i, ]
  cat(sprintf("%-5s %-12s %-10s %-10s %8s %7s %8s %11s %4s %5s %5d %4s %3s\n",
              r$model, r$estimator, r$data, r$moderator,
              fmt_b(r$beta_int), fmt_p(r$p_int),
              fmt_p(r$p_int_holm),
              fmt_p(r$joint_wald_p),
              r$expected_int_sign,
              ifelse(isTRUE(r$sign_consistent), "cons", "WRONG"),
              r$N,
              ifelse(is.na(r$epv), "  NA", sprintf("%4.1f", r$epv)),
              ifelse(r$epv_warning, "!", " ")))
}

cat("\n=========== H2.2 APPENDIX (32 models, BH-adjusted k=32) ======================\n")
cat(sprintf("%-5s %-12s %-10s %-10s %8s %7s %8s %11s %4s %5s %5s %4s %3s\n",
            "Model", "Estimator", "Data", "Moderator",
            "b_int", "p_int", "p_BH", "p_JW_raw",
            "Exp", "Sign", "N", "EPV", "OK?"))
cat(paste(rep("-", 115), collapse = ""), "\n")
app_rows <- core_results[!core_results$main_panel, ]
for (i in seq_len(nrow(app_rows))) {
  r <- app_rows[i, ]
  cat(sprintf("%-5s %-12s %-10s %-10s %8s %7s %8s %11s %4s %5s %5d %4s %3s\n",
              r$model, r$estimator, r$data, r$moderator,
              fmt_b(r$beta_int), fmt_p(r$p_int),
              fmt_p(r$p_int_bh),
              fmt_p(r$joint_wald_p),
              r$expected_int_sign,
              ifelse(isTRUE(r$sign_consistent), "cons", "WRONG"),
              r$N,
              ifelse(is.na(r$epv), "  NA", sprintf("%4.1f", r$epv)),
              ifelse(r$epv_warning, "!", " ")))
}

cat("\nNotes:\n")
cat(" b_int / p_int = coefficient and p-value on eng_c × moderator_c\n")
cat(" p_Holm = Holm-Bonferroni FWER-adjusted (main panel, k=8; confirmatory)\n")
cat(" p_BH = Benjamini-Hochberg FDR-adjusted (appendix, k=32; exploratory)\n")
cat(" p_JW_raw = raw joint Wald on (eng_c, eng_c:mod_c); df=2; H0: no English effect\n")
cat(" Exp = expected sign under H2.2 (higher moderator → English helps more)\n")
cat(" Sign = cons/WRONG vs. expected sign\n")
cat(" EPV = events per variable; '!' marks EPV<5 (low-power warning)\n")
cat(" M17_GEE-M20_GEE = exchangeable GEE companions for yearly broad-success binaries;\n")
cat(" use as the source-campaign clustered inference check for repeated source-years.\n")
cat(" M17s-M20s and M21-M24 remain appendix-only pooled yearly models and should be read as weaker evidence.\n")
cat(" ph_flag = Cox PH diagnostic flag from cox.zph; rows with PH interaction failure\n")
cat(" should remain appendix-only and not be used as substantive support.\n")

cat("\nCore-table summary (40 models):\n")
cat(sprintf("  sign-consistent     : %d/40\n", sum(core_results$sign_consistent, na.rm = TRUE)))
cat(sprintf("  raw p_int < 0.05    : %d/40\n", sum(core_results$sig_05_int, na.rm = TRUE)))
cat(sprintf("  raw p_int < 0.10    : %d/40\n", sum(core_results$sig_10_int, na.rm = TRUE)))
cat(sprintf("  EPV warning (<5)    : %d/40\n", sum(core_results$epv_warning, na.rm = TRUE)))
cat("\nMain panel (k=8, Holm-adjusted):\n")
cat(sprintf("  p_Holm < 0.05       : %d/8\n", sum(core_results$sig_05_int_holm, na.rm = TRUE)))
cat(sprintf("  p_Holm < 0.10       : %d/8\n", sum(core_results$sig_10_int_holm, na.rm = TRUE)))
cat("\nAppendix (k=32, BH-adjusted):\n")
cat(sprintf("  p_BH   < 0.05       : %d/32\n", sum(core_results$sig_05_int_bh,   na.rm = TRUE)))
cat(sprintf("  p_BH   < 0.10       : %d/32\n", sum(core_results$sig_10_int_bh,   na.rm = TRUE)))

robust_rows <- results[results$robustness_row, ]
if (nrow(robust_rows) > 0) {
 cat("\nRecurrent-event Cox robustness (outside multiplicity family):\n")
  cat(sprintf("  rows                 : %d\n", nrow(robust_rows)))
  cat(sprintf("  raw p_int < 0.05     : %d/%d\n", sum(robust_rows$sig_05_int, na.rm = TRUE), nrow(robust_rows)))
  cat(sprintf("  sign-consistent      : %d/%d\n", sum(robust_rows$sign_consistent, na.rm = TRUE), nrow(robust_rows)))
}

flagged_ph <- results[results$ph_flag != "", c("model", "response", "ph_global_p", "ph_int_p", "ph_flag", "interpretation_note")]
if (nrow(flagged_ph) > 0) {
 cat("\nCox PH flags:\n")
  print(flagged_ph, row.names = FALSE)
}

###############################################################################
# FULL MODEL SUMMARIES (for appendix / diagnostics)
###############################################################################

cat("\n\n==== FULL MODEL SUMMARIES ====\n")

print_summary_block <- function(tag, label, mdl) {
  cat(sprintf("\n---- %s %s ----\n", tag, label)); print(summary(mdl))
}

# NAVCO 1.3 — Weibull
print_summary_block("M1",   "Weibull AFT × FDI — time to success",          model1)
print_summary_block("M1s",  "Weibull AFT × FDI — time to strict success",   model1s)
print_summary_block("M2",   "Weibull AFT × Trade — time to success",        model2)
print_summary_block("M2s",  "Weibull AFT × Trade — time to strict success", model2s)
print_summary_block("M3",   "Weibull AFT × Aid — time to success",          model3)
print_summary_block("M3s",  "Weibull AFT × Aid — time to strict success",   model3s)
print_summary_block("M4",   "Weibull AFT × IODonor — time to success",          model4)
print_summary_block("M4s",  "Weibull AFT × IODonor — time to strict success",   model4s)

# NAVCO 1.3 — Firth
print_summary_block("M5",   "Firth PML × FDI — success (binary)",        model5)
print_summary_block("M5s",  "Firth PML × FDI — strict success (binary)", model5s)
print_summary_block("M6",   "Firth PML × Trade — success (binary)",        model6)
print_summary_block("M6s",  "Firth PML × Trade — strict success (binary)", model6s)
print_summary_block("M7",   "Firth PML × Aid — success (binary)",        model7)
print_summary_block("M7s",  "Firth PML × Aid — strict success (binary)", model7s)
print_summary_block("M8",   "Firth PML × IODonor — success (binary)",        model8)
print_summary_block("M8s",  "Firth PML × IODonor — strict success (binary)", model8s)

# NAVCO 1.3 — bracl
print_summary_block("M9",   "bracl × FDI — success3 (ordered)",     model9)
print_summary_block("M10",  "bracl × Trade — success3 (ordered)",   model10)
print_summary_block("M11",  "bracl × Aid — success3 (ordered)",     model11)
print_summary_block("M12",  "bracl × IODonor — success3 (ordered)", model12)

# NAVCO 2.1 — Cox
print_summary_block("M13",  "Cox PH × FDI — hazard of first success (CY)", model13)
print_summary_block("M13s", "Cox PH × FDI — hazard of strict success (CY)", model13s)
print_summary_block("M14",  "Cox PH × Trade — hazard of first success (CY)", model14)
print_summary_block("M14s", "Cox PH × Trade — hazard of strict success (CY)", model14s)
print_summary_block("M15",  "Cox PH × Aid — hazard of first success (CY)", model15)
print_summary_block("M15s", "Cox PH × Aid — hazard of strict success (CY)", model15s)
print_summary_block("M16",  "Cox PH × IODonor — hazard of first success (CY)", model16)
print_summary_block("M16s", "Cox PH × IODonor — hazard of strict success (CY)", model16s)
print_summary_block("M13_AG", "Cox PH × FDI — recurrent-event success hazard (CY)", model13_ag)
print_summary_block("M14_AG", "Cox PH × Trade — recurrent-event success hazard (CY)", model14_ag)
print_summary_block("M15_AG", "Cox PH × Aid — recurrent-event success hazard (CY)", model15_ag)
print_summary_block("M16_AG", "Cox PH × IODonor — recurrent-event success hazard (CY)", model16_ag)
print_summary_block("M17_GEE", "GEE logit × FDI — yearly success, source-clustered", model17_gee)
print_summary_block("M18_GEE", "GEE logit × Trade — yearly success, source-clustered", model18_gee)
print_summary_block("M19_GEE", "GEE logit × Aid — yearly success, source-clustered", model19_gee)
print_summary_block("M20_GEE", "GEE logit × IODonor — yearly success, source-clustered", model20_gee)

# NAVCO 2.1 — Firth
print_summary_block("M17",  "Firth PML × FDI — yearly success",     model17)
print_summary_block("M17s", "Firth PML × FDI — yearly strict success", model17s)
print_summary_block("M18",  "Firth PML × Trade — yearly success",     model18)
print_summary_block("M18s", "Firth PML × Trade — yearly strict success", model18s)
print_summary_block("M19",  "Firth PML × Aid — yearly success",     model19)
print_summary_block("M19s", "Firth PML × Aid — yearly strict success", model19s)
print_summary_block("M20",  "Firth PML × IODonor — yearly success",     model20)
print_summary_block("M20s", "Firth PML × IODonor — yearly strict success", model20s)

# NAVCO 2.1 — bracl
print_summary_block("M21",  "bracl × FDI — progress_ord (6-level; failure separate)",     model21)
print_summary_block("M22",  "bracl × Trade — progress_ord (6-level; failure separate)",   model22)
print_summary_block("M23",  "bracl × Aid — progress_ord (6-level; failure separate)",     model23)
print_summary_block("M24",  "bracl × IODonor — progress_ord (6-level; failure separate)", model24)

expected_result_cols <- c(
  "model", "estimator", "data", "response", "source_outcome",
  "moderator", "mod_direction",
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
  stop(sprintf("Expected 48 H2.2 result rows, found %d.", nrow(results)), call. = FALSE)
}
if (anyDuplicated(results$model) > 0) {
 stop("Duplicate model labels detected in H2.2 results.", call. = FALSE)
}

###############################################################################
# EXPORT CSV
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
stopifnot(identical(names(models), results$model))
saveRDS(models, cache_path, version = 3)
saveRDS(list(df = df, d21 = d21, d21_first_conc = d21_first_conc), frames_path, version = 3)
cat(sprintf("\nResults exported to: %s\n", out_path))
cat(sprintf("Model bundle exported to: %s\n", cache_path))
cat(sprintf("Analysis frames exported to: %s\n", frames_path))
