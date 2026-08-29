# Fit the eight appendix-only H2.1 limited-concession specifications.
# Inputs must already carry the centered English and moderator variables used
# by the accepted H2.1 analysis.
fit_h21_limited_extension <- function(df, d21) {
  df$event_limited <- as.integer(as.numeric(as.character(df$limited)) == 1)
  df$limited_bin <- df$event_limited

  d21$event_limited_yr <- as.integer(d21$progress == 2)
  d21$limited_bin_yr <- d21$event_limited_yr

  first_event_panel <- function(dat, event_col) {
    out <- dat %>%
      arrange(navco21_id, time_in_campaign) %>%
      group_by(navco21_id) %>%
      mutate(
        .cum_event = cumsum(.data[[event_col]]),
        .prior_event = lag(.cum_event, default = 0L),
        .keep_risk = .prior_event == 0L,
        .event_first = as.integer(.data[[event_col]] == 1 & .cum_event == 1L)
      ) %>%
      filter(.keep_risk) %>%
      ungroup()
    out[[event_col]] <- out$.event_first
    out %>% select(-.cum_event, -.prior_event, -.keep_risk, -.event_first)
  }
  d21_first_lim <- first_event_panel(d21, "event_limited_yr")

  fit_weibull13_limited <- function(moderator) {
    f <- as.formula(sprintf(
      paste0(
        "Surv(duration, event_limited) ~ eng_c * %s + nonviolent_camp + ",
        "REGCHANGE + v2csreprss + lnnum_image_sum"
      ), moderator
    ))
    survreg(f, data = df, dist = "weibull")
  }
  fit_firth13_limited <- function(moderator) {
    f <- as.formula(sprintf(
      paste0(
        "limited_bin ~ eng_c * %s + nonviolent_camp + REGCHANGE + ",
        "v2csreprss + lnlengthofcam + lnnum_image_sum"
      ), moderator
    ))
    logistf(
      f, data = df,
      control = logistf.control(maxit = 1000, maxstep = 0.5),
      plcontrol = logistpl.control(maxit = 1000), pl = TRUE
    )
  }
  fit_cox21_limited <- function(moderator) {
    f <- as.formula(sprintf(
      paste0(
        "Surv(tstart, tstop, event_limited_yr) ~ eng_c * %s + nonviolent + ",
        "REGCHANGE + v2csreprss + lnnum_image_sum + cluster(navco21_id)"
      ), moderator
    ))
    coxph(f, data = d21_first_lim)
  }
  fit_firth21_limited <- function(moderator) {
    f <- as.formula(sprintf(
      paste0(
        "limited_bin_yr ~ eng_c * %s + nonviolent + REGCHANGE + ",
        "v2csreprss + time_in_campaign + lnnum_image_sum"
      ), moderator
    ))
    logistf(
      f, data = d21,
      control = logistf.control(maxit = 1000, maxstep = 0.5),
      plcontrol = logistpl.control(maxit = 1000), pl = TRUE
    )
  }

  models <- list(
    M1l = fit_weibull13_limited("domaut_c"),
    M2l = fit_weibull13_limited("stterr_c"),
    M3l = fit_firth13_limited("domaut_c"),
    M4l = fit_firth13_limited("stterr_c"),
    M7l = fit_cox21_limited("domaut_c"),
    M8l = fit_cox21_limited("stterr_c"),
    M9l = fit_firth21_limited("domaut_c"),
    M10l = fit_firth21_limited("stterr_c")
  )

  cc_stats <- function(dat, vars, event_var) {
    keep <- complete.cases(dat[, vars, drop = FALSE])
    list(
      N = sum(keep),
      n_events = as.integer(sum(dat[[event_var]][keep] == 1, na.rm = TRUE))
    )
  }
  vars_13_surv <- function(moderator) c(
    "duration", "event_limited", "eng_c", moderator, "nonviolent_camp",
    "REGCHANGE", "v2csreprss", "lnnum_image_sum"
  )
  vars_13_bin <- function(moderator) c(
    "limited_bin", "eng_c", moderator, "nonviolent_camp", "REGCHANGE",
    "v2csreprss", "lnlengthofcam", "lnnum_image_sum"
  )
  vars_21_surv <- function(moderator) c(
    "tstart", "tstop", "event_limited_yr", "eng_c", moderator,
    "nonviolent", "REGCHANGE", "v2csreprss", "lnnum_image_sum", "navco21_id"
  )
  vars_21_bin <- function(moderator) c(
    "limited_bin_yr", "eng_c", moderator, "nonviolent", "REGCHANGE",
    "v2csreprss", "time_in_campaign", "lnnum_image_sum"
  )

  stats <- list(
    M1l = cc_stats(df, vars_13_surv("domaut_c"), "event_limited"),
    M2l = cc_stats(df, vars_13_surv("stterr_c"), "event_limited"),
    M3l = cc_stats(df, vars_13_bin("domaut_c"), "limited_bin"),
    M4l = cc_stats(df, vars_13_bin("stterr_c"), "limited_bin"),
    M7l = cc_stats(d21_first_lim, vars_21_surv("domaut_c"), "event_limited_yr"),
    M8l = cc_stats(d21_first_lim, vars_21_surv("stterr_c"), "event_limited_yr"),
    M9l = cc_stats(d21, vars_21_bin("domaut_c"), "limited_bin_yr"),
    M10l = cc_stats(d21, vars_21_bin("stterr_c"), "limited_bin_yr")
  )

  extract_row <- function(tag, estimator, data_label, response, moderator) {
    model <- models[[tag]]
    cf <- coef(model)
    int_term <- paste0("eng_c:", moderator)
    if (inherits(model, "survreg")) {
      tab <- summary(model)$table
      p <- tab[c("eng_c", moderator, int_term), "p"]
    } else if (inherits(model, "coxph")) {
      tab <- summary(model)$coefficients
      p <- tab[c("eng_c", moderator, int_term), "Pr(>|z|)"]
    } else {
      p <- unname(model$prob[c("eng_c", moderator, int_term)])
    }
    s <- stats[[tag]]
    data.frame(
      model = tag, estimator = estimator, data = data_label,
      response = response, moderator = moderator,
      beta_eng = unname(cf["eng_c"]), p_eng = unname(p[[1L]]),
      beta_mod = unname(cf[moderator]), p_mod = unname(p[[2L]]),
      beta_int = unname(cf[int_term]), p_int = unname(p[[3L]]),
      N = s$N, n_events = s$n_events, epv = s$n_events / length(cf),
      stringsAsFactors = FALSE
    )
  }

  results <- bind_rows(
    extract_row("M1l", "Weibull AFT", "NAVCO 1.3", "Time to limited success", "domaut_c"),
    extract_row("M2l", "Weibull AFT", "NAVCO 1.3", "Time to limited success", "stterr_c"),
    extract_row("M3l", "Firth PML", "NAVCO 1.3", "Limited success (binary)", "domaut_c"),
    extract_row("M4l", "Firth PML", "NAVCO 1.3", "Limited success (binary)", "stterr_c"),
    extract_row("M7l", "Cox PH", "NAVCO 2.1", "Hazard of first limited", "domaut_c"),
    extract_row("M8l", "Cox PH", "NAVCO 2.1", "Hazard of first limited", "stterr_c"),
    extract_row("M9l", "Firth PML", "NAVCO 2.1", "Yearly limited (binary)", "domaut_c"),
    extract_row("M10l", "Firth PML", "NAVCO 2.1", "Yearly limited (binary)", "stterr_c")
  )
  results$p_int_bh <- p.adjust(results$p_int, method = "BH", n = nrow(results))
  results$epv_warning <- results$epv < 5
  results$sig_05_int <- !is.na(results$p_int) & results$p_int < 0.05
  results$sig_05_int_bh <- !is.na(results$p_int_bh) & results$p_int_bh < 0.05

  list(results = results, models = models)
}
