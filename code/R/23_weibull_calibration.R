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
require_files(bundle_paths)
bundles <- setNames(lapply(bundle_paths, readRDS), c("H1", "H2.1", "H2.2", "H3"))
B <- read_positive_integer("BOOT_REPS", if (Sys.getenv("REPLICATION_SMOKE") == "1") 20L else 999L)

specs <- rbind(
  data.frame(family = "H1", spec = "M1", moderator = NA, focal = "I(eng_prob_general^2)", k = 10L),
  data.frame(family = "H2.1", spec = c("M1", "M2"), moderator = c("domaut_c", "stterr_c"),
             focal = c("eng_c:domaut_c", "eng_c:stterr_c"), k = 8L),
  data.frame(family = "H2.2", spec = paste0("M", 1:4),
             moderator = c("fdiin_c", "trade_c", "lnaid_c", "iodonor_c"),
             focal = paste0("eng_c:", c("fdiin_c", "trade_c", "lnaid_c", "iodonor_c")), k = 8L),
  data.frame(family = "H3", spec = paste0("M", 1:4),
             moderator = c("repress_c", "civlib_c", "clpol_c", "polity_c"),
             focal = paste0("eng_c:", c("repress_c", "civlib_c", "clpol_c", "polity_c")), k = 8L)
)

calibrate <- function(family, spec, moderator, focal, k) {
  frame_family <- c(H1 = "H1", `H2.1` = "H21", `H2.2` = "H22", H3 = "H3")[[family]]
  load_analysis_frames(P$cache, frame_family)
  model <- bundles[[family]][[spec]]
  model_formula <- rep_model_formula(model)
  labels <- attr(terms(model_formula), "term.labels")
  keep <- labels[labels != focal]
  response <- paste(deparse(model_formula[[2]], width.cutoff = 500L), collapse = "")
  null_formula <- reformulate(keep, response = response, env = environment(model_formula))
  base_data <- model_simulation_data(model)
  # Preserve the exact fitted null frame. The bootstrap simulator must not
  # depend on recovering the local `base_data` symbol from the serialized call.
  null_model <- survival::survreg(null_formula, data = base_data, dist = model$dist,
                                  model = TRUE, x = TRUE, y = TRUE)

  observed <- summary(model)$table
  z_obs <- abs(unname(observed[focal, "z"]))
  p_nominal <- unname(observed[focal, "p"])
  set.seed(named_seed(paste0("CAL-", family), spec, 0))
  rep_seeds <- sample.int(.Machine$integer.max, B)
  z_boot <- vapply(rep_seeds, function(seed) {
    set.seed(seed)
    tryCatch({
      sim <- simulate_weibull_model(null_model, names(coef(null_model))[[1]], 1)
      # simulate_weibull_model rescales the named term; restore the fitted-null
      # DGP by passing magnitude one, then fit the accepted full specification.
      refit <- survival::survreg(model_formula, data = sim, dist = model$dist)
      abs(unname(summary(refit)$table[focal, "z"]))
    }, error = function(e) NA_real_)
  }, numeric(1))
  valid <- z_boot[is.finite(z_boot)]
  if (!length(valid)) {
    stop(sprintf("All Weibull calibration refits failed for %s %s.", family, spec),
         call. = FALSE)
  }
  data.frame(
    family = family, spec = spec, moderator = moderator, focal_coef = focal,
    beta_int_obs = unname(coef(model)[focal]), z_obs_nominal = z_obs,
    p_nominal = p_nominal,
    z_critical_emp_05 = unname(quantile(valid, .95, names = FALSE)),
    z_critical_emp_10 = unname(quantile(valid, .90, names = FALSE)),
    z_critical_holm = unname(quantile(valid, 1 - .05 / k, names = FALSE)),
    inflation_factor = unname(quantile(valid, .95, names = FALSE)) / qnorm(.975),
    empirical_alpha_at_nominal_z196 = mean(valid >= qnorm(.975)),
    p_calibrated_emp = mean(valid >= z_obs), n_boot_reps_ok = length(valid),
    fit_ok_rate = length(valid) / B, stringsAsFactors = FALSE
  )
}

out <- do.call(rbind, lapply(seq_len(nrow(specs)), function(i) {
  calibrate(specs$family[i], specs$spec[i], specs$moderator[i], specs$focal[i], specs$k[i])
}))
suffix <- if (Sys.getenv("REPLICATION_SMOKE") == "1") "_smoke" else "_recomputed"
write.csv(out, file.path(P$diagnostics, paste0("weibull_bootstrap_calibration", suffix, ".csv")), row.names = FALSE, na = "")
cat(sprintf("Wrote %d Weibull calibration rows using %d bootstrap draws per row.\n", nrow(out), B))
