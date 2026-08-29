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
rep_require_packages(c("survival", "logistf"))
suppressPackageStartupMessages({ library(survival); library(logistf) })
load_analysis_frames(P$cache, "H1")
bundle_path <- file.path(P$cache, "H1_models.rds")
require_files(bundle_path)
models <- readRDS(bundle_path)
n_reps <- read_positive_integer("POWER_REPS", if (Sys.getenv("REPLICATION_SMOKE") == "1") 10L else 2000L)
specs <- data.frame(spec = c("M1", "M2", "M5"),
                    estimator = c("Weibull AFT", "Firth PML", "Firth PML"))
out <- do.call(rbind, lapply(seq_len(nrow(specs)), function(i) {
  run_h1_power(models[[specs$spec[i]]], specs$spec[i], specs$estimator[i], n_reps)
}))
out$failure_reason <- ifelse(out$fit_ok_rate == 0, "all_refits_failed", NA_character_)
out <- out[c("spec", "test", "beta_eng_hat", "beta_eng2_hat", "magnitude",
             "beta_eng_sim", "beta_eng2_sim", "n_reps", "fit_ok_rate", "alpha",
             "raw_power_overall", "holm_proxy_power_overall", "raw_power_conditional",
             "holm_proxy_power_conditional", "raw_power", "holm_proxy_power", "failure_reason")]
suffix <- if (Sys.getenv("REPLICATION_SMOKE") == "1") "_smoke" else "_recomputed"
write.csv(out, file.path(P$estimates, paste0("power_h1", suffix, ".csv")), row.names = FALSE, na = "")
cat(sprintf("Wrote %d H1 power cells (%d replications per design cell).\n", nrow(out), n_reps))
