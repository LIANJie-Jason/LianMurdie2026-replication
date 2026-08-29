root <- Sys.getenv("REPLICATION_ROOT", unset = "")
if (!nzchar(root)) {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  file_arg <- sub("^--file=", "", file_arg[[1L]])
  file_arg <- gsub("~+~", " ", file_arg, fixed = TRUE)
  root <- normalizePath(file.path(dirname(file_arg), "..", ".."), mustWork = TRUE)
}
source(file.path(root, "code", "R", "00_setup.R"))
source(file.path(root, "code", "R", "R", "power_figure_helpers.R"))
P <- init_replication(root)

families <- c(h1 = "H1", h21 = "H2.1", h22 = "H2.2", h3 = "H3")
suffix <- if (Sys.getenv("REPLICATION_SMOKE") == "1") "_smoke" else "_recomputed"
fresh_paths <- file.path(P$estimates, paste0("power_", names(families), suffix, ".csv"))
reference_paths <- file.path(P$reference, "model_results", paste0("power_", names(families), ".csv"))
require_files(c(fresh_paths, reference_paths))

h1_summary <- function(x) {
  keys <- unique(x[c("spec", "test")])
  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, ]
    z <- x[x$spec == key$spec & x$test == key$test & x$alpha == .05, ]
    at <- function(m, col = "raw_power") z[z$magnitude == m, col][[1]]
    data.frame(
      spec = key$spec,
      estimator = c(M1 = "Weibull AFT", M2 = "Firth PML", M5 = "Firth PML")[[key$spec]],
      moderator = "(none - H1 unconditional curvilinear)", test = key$test,
      beta_eng_hat = z$beta_eng_hat[[1]], beta_eng2_hat = z$beta_eng2_hat[[1]],
      type1_at_null_05 = at(0), raw_power_at_half_05 = at(.5),
      raw_power_at_fitted_05 = at(1), raw_power_at_double_05 = at(2),
      holm_power_at_fitted_05 = at(1, "holm_proxy_power"),
      raw_mde_mult_05 = min(z$magnitude[z$magnitude > 0 & z$raw_power >= .8], Inf),
      holm_mde_mult_05 = min(z$magnitude[z$magnitude > 0 & z$holm_proxy_power >= .8], Inf),
      stringsAsFactors = FALSE
    )
  }))
}

annotate_weibull <- function(data, calibration, family) {
  fields <- c("weibull_power_citable", "weibull_power_status", "weibull_calibration_source",
              "weibull_empirical_alpha_nominal_05", "weibull_z_critical_emp_05",
              "weibull_z_critical_holm", "weibull_p_calibrated_observed",
              "weibull_calibration_fit_ok_rate", "weibull_power_note")
  for (field in fields) data[[field]] <- NA
  for (spec in unique(data$spec[data$estimator == "Weibull AFT"])) {
    index <- data$spec == spec & data$estimator == "Weibull AFT"
    c_row <- calibration[calibration$family == family & calibration$spec == spec, ]
    if (nrow(c_row) != 1L) stop(sprintf("Missing calibration for %s %s", family, spec))
    data$weibull_power_citable[index] <- FALSE
    data$weibull_power_status[index] <- "nominal_wald_diagnostic_not_size_calibrated"
    data$weibull_calibration_source[index] <- "weibull_bootstrap_calibration.csv"
    data$weibull_empirical_alpha_nominal_05[index] <- c_row$empirical_alpha_at_nominal_z196
    data$weibull_z_critical_emp_05[index] <- c_row$z_critical_emp_05
    data$weibull_z_critical_holm[index] <- c_row$z_critical_holm
    data$weibull_p_calibrated_observed[index] <- c_row$p_calibrated_emp
    data$weibull_calibration_fit_ok_rate[index] <- c_row$fit_ok_rate
    data$weibull_power_note[index] <- paste0("Nominal survreg Wald power/MDE diagnostic; bootstrap calibration ",
      "shows anti-conservative size, so this row is not size-calibrated citation-ready power.")
  }
  data
}

comparisons <- list()
for (i in seq_along(families)) {
  short <- names(families)[i]
  fresh <- read.csv(fresh_paths[i], stringsAsFactors = FALSE, check.names = FALSE)
  reference <- read.csv(reference_paths[i], stringsAsFactors = FALSE, check.names = FALSE)
  summary <- if (short == "h1") h1_summary(fresh) else power_summary(fresh)
  if (short %in% c("h22", "h3")) {
    calibration <- read.csv(file.path(P$diagnostics, paste0("weibull_bootstrap_calibration", suffix, ".csv")),
                            stringsAsFactors = FALSE, check.names = FALSE)
    fresh <- annotate_weibull(fresh, calibration, families[i])
    summary <- annotate_weibull(summary, calibration, families[i])
  }
  if (short == "h21") summary <- summary[setdiff(names(summary),
    c("mde_at_80_raw_interp", "mde_at_80_holm_interp", "fit_ok_rate_min", "fit_ok_rate_mean"))]
  if (short == "h3") summary <- summary[setdiff(names(summary),
    c("mde_at_80_raw_interp", "mde_at_80_holm_interp", "fit_ok_rate_min", "fit_ok_rate_mean"))]
  summary_reference <- read.csv(
    file.path(P$reference, "model_results", paste0("power_", short, "_summary.csv")),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  if (!setequal(names(summary), names(summary_reference))) {
    stop(sprintf("Fresh %s power summary schema differs from accepted schema.", short))
  }
  summary <- summary[names(summary_reference)]
  write.csv(fresh, fresh_paths[i], row.names = FALSE, na = "")
  write.csv(summary, file.path(P$diagnostics, paste0("power_", short, suffix, "_summary.csv")),
            row.names = FALSE, na = "")

  keys <- c("spec", "magnitude", "alpha")
  if (short == "h1") keys <- c(keys, "test")
  ref_small <- reference[c(keys, "n_reps", "raw_power", "holm_proxy_power")]
  names(ref_small)[names(ref_small) == "n_reps"] <- "accepted_n_reps"
  names(ref_small)[names(ref_small) == "raw_power"] <- "accepted_raw_power"
  names(ref_small)[names(ref_small) == "holm_proxy_power"] <- "accepted_holm_power"
  rep_assert(!anyDuplicated(fresh[keys]), "Fresh %s power keys are not unique", short)
  rep_assert(!anyDuplicated(ref_small[keys]), "Accepted %s power keys are not unique", short)
  cmp <- merge(fresh[c(keys, "n_reps", "raw_power", "holm_proxy_power")], ref_small,
               by = keys, all.x = TRUE, sort = FALSE)
  unmatched <- is.na(cmp$accepted_n_reps) | is.na(cmp$accepted_raw_power) |
    is.na(cmp$accepted_holm_power)
  if (any(unmatched)) {
    bad_keys <- apply(cmp[unmatched, keys, drop = FALSE], 1L, paste, collapse = "|")
    stop(sprintf("Fresh %s power cells have no accepted match: %s", short,
                 paste(head(bad_keys, 10L), collapse = ", ")), call. = FALSE)
  }
  cmp$family <- families[i]
  cmp$raw_abs_diff <- abs(cmp$raw_power - cmp$accepted_raw_power)
  cmp$holm_abs_diff <- abs(cmp$holm_proxy_power - cmp$accepted_holm_power)
  comparisons[[i]] <- cmp
}
comparison_columns <- unique(unlist(lapply(comparisons, names)))
comparisons <- lapply(comparisons, function(x) {
  for (name in setdiff(comparison_columns, names(x))) x[[name]] <- NA
  x[comparison_columns]
})
comparison <- do.call(rbind, comparisons)
two_sample_mc_p <- function(fresh_rate, fresh_n, accepted_rate, accepted_n) {
  pooled <- (fresh_rate * fresh_n + accepted_rate * accepted_n) / (fresh_n + accepted_n)
  standard_error <- sqrt(pooled * (1 - pooled) * (1 / fresh_n + 1 / accepted_n))
  difference <- abs(fresh_rate - accepted_rate)
  out <- rep(NA_real_, length(difference))
  zero_se <- is.finite(standard_error) & standard_error == 0
  out[zero_se] <- ifelse(difference[zero_se] == 0, 1, 0)
  usable <- is.finite(standard_error) & standard_error > 0
  out[usable] <- 2 * stats::pnorm(-difference[usable] / standard_error[usable])
  out
}
comparison$raw_mc_p <- two_sample_mc_p(
  comparison$raw_power, comparison$n_reps,
  comparison$accepted_raw_power, comparison$accepted_n_reps
)
comparison$holm_mc_p <- two_sample_mc_p(
  comparison$holm_proxy_power, comparison$n_reps,
  comparison$accepted_holm_power, comparison$accepted_n_reps
)
adjusted <- stats::p.adjust(c(comparison$raw_mc_p, comparison$holm_mc_p), method = "holm")
comparison$raw_mc_p_holm <- adjusted[seq_len(nrow(comparison))]
comparison$holm_mc_p_holm <- adjusted[nrow(comparison) + seq_len(nrow(comparison))]
comparison$familywise_alpha <- .05
comparison$comparison_method <- "two-sample Monte Carlo proportion z tests; Holm across raw and adjusted rates"
comparison$comparison_status <- if (Sys.getenv("REPLICATION_SMOKE") == "1") {
  "SMOKE_ADVISORY"
} else ifelse(comparison$raw_mc_p_holm >= .05 & comparison$holm_mc_p_holm >= .05,
              "PASS", "FAIL")
write.csv(comparison, file.path(P$diagnostics, paste0("power_reference_comparison", suffix, ".csv")),
          row.names = FALSE, na = "")

# Cox power values are not recomputed: the accepted counting-process DGP was
# invalid. Preserve its three reader-facing quarantine rows as immutable audit
# evidence and never mingle them with fresh power estimates.
q_ref <- file.path(P$reference, "model_results", "cox_power_quarantine_audit.csv")
require_files(q_ref)
if (!file.copy(q_ref,
               file.path(P$diagnostics, paste0("cox_power_quarantine_audit", suffix, ".csv")),
               overwrite = TRUE)) {
  stop("Failed to preserve the accepted Cox quarantine audit.", call. = FALSE)
}

cal_path <- file.path(P$diagnostics, paste0("weibull_bootstrap_calibration", suffix, ".csv"))
require_files(cal_path)
cal <- read.csv(cal_path, stringsAsFactors = FALSE, check.names = FALSE)
audit_ref <- read.csv(file.path(P$reference, "model_results", "weibull_power_diagnostic_audit.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
audit <- audit_ref
for (i in seq_len(nrow(audit))) {
  if (grepl("_smoke", audit$source_file[i], fixed = TRUE)) next
  short <- if (audit$family[i] == "H2.2") "h22" else "h3"
  sum_path <- file.path(P$diagnostics, paste0("power_", short, suffix, "_summary.csv"))
  sum_data <- read.csv(sum_path, stringsAsFactors = FALSE, check.names = FALSE)
  s <- sum_data[sum_data$spec == audit$spec[i], ]
  c_row <- cal[cal$family == audit$family[i] & cal$spec == audit$spec[i], ]
  audit$type1_at_null_05[i] <- s$type1_at_null_05
  audit$raw_power_at_fitted_05[i] <- s$raw_power_at_fitted_05
  audit$empirical_alpha_at_nominal_z196[i] <- c_row$empirical_alpha_at_nominal_z196
  audit$z_critical_emp_05[i] <- c_row$z_critical_emp_05
  audit$z_critical_holm[i] <- c_row$z_critical_holm
  audit$p_calibrated_emp[i] <- c_row$p_calibrated_emp
}
write.csv(audit, file.path(P$diagnostics, paste0("weibull_power_diagnostic_audit", suffix, ".csv")),
          row.names = FALSE, na = "")

if (Sys.getenv("REPLICATION_SMOKE") != "1" && any(comparison$comparison_status == "FAIL")) {
  stop("At least one fresh power cell differs after familywise Holm-adjusted Monte Carlo comparison.",
       call. = FALSE)
}
cat(sprintf("Wrote four power summaries and %d reference-comparison rows.\n", nrow(comparison)))
