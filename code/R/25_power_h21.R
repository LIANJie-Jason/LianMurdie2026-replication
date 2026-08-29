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
load_analysis_frames(P$cache, "H21")
bundle_path <- file.path(P$cache, "H21_models.rds")
require_files(bundle_path)
models <- readRDS(bundle_path)
n_reps <- read_positive_integer("POWER_REPS", if (Sys.getenv("REPLICATION_SMOKE") == "1") 10L else 2000L)
specs <- data.frame(
  spec = c("M1", "M2", "M3", "M4", "M9", "M10"),
  estimator = c("Weibull AFT", "Weibull AFT", rep("Firth PML", 4)),
  moderator = c("domaut_c", "stterr_c", "domaut_c", "stterr_c", "domaut_c", "stterr_c"),
  stringsAsFactors = FALSE
)
out <- do.call(rbind, lapply(seq_len(nrow(specs)), function(i) {
  focal <- paste0("eng_c:", specs$moderator[i])
  run_focal_power(models[[specs$spec[i]]], "H2.1", specs$spec[i], specs$estimator[i],
                  specs$moderator[i], focal, 8L, n_reps)
}))
suffix <- if (Sys.getenv("REPLICATION_SMOKE") == "1") "_smoke" else "_recomputed"
write.csv(out, file.path(P$estimates, paste0("power_h21", suffix, ".csv")), row.names = FALSE, na = "")
cat(sprintf("Wrote %d H2.1 power cells; Cox M7/M8 remain quarantined.\n", nrow(out)))
