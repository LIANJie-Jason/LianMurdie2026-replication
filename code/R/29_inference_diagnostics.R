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
rep_require_packages("survival")
suppressPackageStartupMessages(library(survival))

bundle_paths <- file.path(P$cache, c("H1_models.rds", "H21_models.rds", "H22_models.rds", "H3_models.rds"))
frame_paths <- file.path(P$cache, c("H1_analysis_frames.rds", "H21_analysis_frames.rds",
                                    "H22_analysis_frames.rds", "H3_analysis_frames.rds"))
cox_spec <- file.path(P$code, "specs", "inference", "cox_ph_membership.csv")
suffix <- if (Sys.getenv("REPLICATION_SMOKE") == "1") "_smoke" else ""
cal_path <- file.path(P$diagnostics, paste0("weibull_bootstrap_calibration", suffix, ".csv"))
q_path <- file.path(P$reference, "appendix_tables", "inference_diagnostics_cox_power_quarantine_table_body.csv")
require_files(c(bundle_paths, frame_paths, cox_spec, cal_path, q_path))
bundles <- setNames(lapply(bundle_paths, readRDS), c("H1", "H2.1", "H2.2", "H3"))
frames <- setNames(lapply(frame_paths, readRDS), names(bundles))
membership <- read.csv(cox_spec, stringsAsFactors = FALSE, check.names = FALSE)
membership$N <- NA_integer_
membership$Events <- NA_integer_
membership$`Global PH p` <- NA_real_
membership$`Focal PH p` <- NA_real_
membership$`Diagnostic status` <- NA_character_

extract_zph <- function(model, family, moderator) {
  data_environment <- new.env(parent = environment(model$terms))
  list2env(frames[[family]], envir = data_environment)
  environment(model$terms) <- data_environment
  z <- tryCatch(survival::cox.zph(model), error = function(error) {
    stop(sprintf("cox.zph failed for %s: %s", family, conditionMessage(error)), call. = FALSE)
  })
  tab <- z$table
  p_col <- ncol(tab)
  global <- if ("GLOBAL" %in% rownames(tab)) tab["GLOBAL", p_col] else NA_real_
  candidates <- if (family == "H1") grep("eng.*\\^2", rownames(tab), value = TRUE) else
    grep(paste0("eng_c.*", moderator, "|", moderator, ".*eng_c"), rownames(tab), value = TRUE)
  focal <- if (length(candidates)) tab[candidates[[1]], p_col] else NA_real_
  c(global = unname(global), focal = unname(focal))
}

for (i in seq_len(nrow(membership))) {
  family <- membership$Family[i]
  id <- membership$Model[i]
  model <- bundles[[family]][[id]]
  if (is.null(model)) stop(sprintf("Missing Cox model %s %s", family, id), call. = FALSE)
  membership$N[i] <- model$n
  membership$Events[i] <- model$nevent
  if (grepl("Recurrent-event", membership$Role[i], fixed = TRUE)) {
    membership$`Diagnostic status`[i] <- "Not reported for recurrent-event robustness"
    next
  }
  p <- extract_zph(model, family, membership$Moderator[i])
  membership$`Global PH p`[i] <- p[["global"]]
  membership$`Focal PH p`[i] <- p[["focal"]]
  membership$`Diagnostic status`[i] <- if (all(!is.finite(p))) {
    "cox.zph unavailable"
  } else if (any(p < .05, na.rm = TRUE)) "cox.zph flag" else "No cox.zph flag"
}
membership <- membership[c(
  "Family", "Main table", "Model", "Role", "Moderator", "Response / estimand",
  "Cox estimand", "Outcome source", "N", "Events", "Global PH p", "Focal PH p",
  "Diagnostic status", "Source", "Note"
)]
write.csv(membership,
          file.path(P$diagnostics, paste0("inference_diagnostics_cox_ph", suffix, ".csv")),
          row.names = FALSE, na = "NA")
membership_display <- membership
membership_display$`Global PH p` <- round(membership_display$`Global PH p`, 3)
membership_display$`Focal PH p` <- round(membership_display$`Focal PH p`, 3)
write.csv(membership_display,
          file.path(P$diagnostics, paste0("inference_diagnostics_cox_ph_table_body", suffix, ".csv")),
          row.names = FALSE, na = "NA")

quarantine <- read.csv(q_path, stringsAsFactors = FALSE, check.names = FALSE)
write.csv(quarantine,
          file.path(P$diagnostics, paste0("inference_diagnostics_cox_power_quarantine_table_body", suffix, ".csv")),
          row.names = FALSE, na = "NA")

cal <- read.csv(cal_path, stringsAsFactors = FALSE, check.names = FALSE)
weibull <- data.frame(
  Family = cal$family, Model = cal$spec, Moderator = cal$moderator,
  `Focal term` = cal$focal_coef, `Nominal p` = round(cal$p_nominal, 3),
  `Bootstrap-calibrated p` = round(cal$p_calibrated_emp, 3),
  `Empirical alpha@1.96` = round(cal$empirical_alpha_at_nominal_z196, 3),
  `Empirical z(.05)` = round(cal$z_critical_emp_05, 3),
  `Holm z crit` = round(cal$z_critical_holm, 3),
  `Bootstrap reps` = cal$n_boot_reps_ok,
  `Diagnostic status` = ifelse(cal$family %in% c("H2.2", "H3"),
    "Power row diagnostic-only; not size-calibrated", "Calibration diagnostic"),
  `Citable status` = ifelse(cal$family %in% c("H2.2", "H3"),
    "Diagnostic-only power", "Diagnostic"),
  Source = ifelse(cal$family %in% c("H2.2", "H3"),
    paste0("weibull_bootstrap_calibration.csv; power_", tolower(gsub("\\.", "", cal$family)), "_summary.csv"),
    "weibull_bootstrap_calibration.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
write.csv(weibull, file.path(P$diagnostics, paste0("inference_diagnostics_weibull_table_body", suffix, ".csv")),
          row.names = FALSE, na = "NA")
cat(sprintf("Wrote Cox PH (%d), Cox quarantine (%d), and Weibull (%d) diagnostic rows.\n",
            nrow(membership), nrow(quarantine), nrow(weibull)))
