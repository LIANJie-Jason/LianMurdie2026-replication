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
rep_require_packages("ggplot2")
suppressPackageStartupMessages({
  library(survival)
  library(logistf)
  library(brglm2)
  library(geepack)
})

suffix <- if (Sys.getenv("REPLICATION_SMOKE") == "1") "_smoke" else ""

# C25 is an accepted codestack display artifact. Validate its schema contract and
# materialize it byte-exactly; optional recomputation never overwrites it.
c25_ref <- file.path(P$reference, "appendix_tables", "C25_power_final_body.csv")
require_files(c25_ref)
c25 <- read.csv(c25_ref, stringsAsFactors = FALSE, check.names = FALSE)
expected_c25 <- c("Family", "Main table", "Main-table model label", "Estimator", "Focal test",
                  "Power at fitted effect", "False-positive rate", "Status / note")
if (!identical(names(c25), expected_c25) || nrow(c25) != 28L || "Source" %in% names(c25)) {
 stop("Accepted C25 violates the schema/order/row-count contract.", call. = FALSE)
}
if (!nzchar(suffix) &&
    !file.copy(c25_ref, file.path(P$appendix_tables, "C25_power_table_body.csv"),
               overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)) {
 stop("Could not materialize accepted C25 display.", call. = FALSE)
}

expected_rows <- c(34L, 3L, 11L)
if (!nzchar(suffix)) {
  c26_refs <- file.path(P$reference, "appendix_tables", c(
    "inference_diagnostics_cox_ph_table_body.csv",
    "inference_diagnostics_cox_power_quarantine_table_body.csv",
    "inference_diagnostics_weibull_table_body.csv"
  ))
  require_files(c26_refs)
  for (i in seq_along(c26_refs)) {
    body <- read.csv(c26_refs[[i]], stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(body) != expected_rows[[i]]) {
      stop(sprintf("Accepted C26 panel %d has %d rows; expected %d.",
                   i, nrow(body), expected_rows[[i]]), call. = FALSE)
    }
    destination <- file.path(P$appendix_tables, paste0("C26_panel_", i, "_table_body.csv"))
    if (!file.copy(c26_refs[[i]], destination, overwrite = TRUE,
                   copy.mode = TRUE, copy.date = FALSE)) {
      stop(sprintf("Could not materialize accepted C26 panel %d.", i), call. = FALSE)
    }
  }
} else {
  c26_names <- paste0(c("inference_diagnostics_cox_ph_table_body",
                        "inference_diagnostics_cox_power_quarantine_table_body",
                        "inference_diagnostics_weibull_table_body"), suffix, ".csv")
  c26_paths <- file.path(P$diagnostics, c26_names)
  require_files(c26_paths)
  for (i in seq_along(c26_paths)) {
    body <- read.csv(c26_paths[[i]], stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(body) != expected_rows[[i]]) {
      stop(sprintf("C26 smoke panel %d has %d rows; expected %d.",
                   i, nrow(body), expected_rows[[i]]), call. = FALSE)
    }
    write.csv(body,
              file.path(P$appendix_tables, paste0("C26_panel_", i, suffix, "_table_body.csv")),
              row.names = FALSE, na = "NA")
  }
}

if (!nzchar(suffix)) {
bundle_paths <- file.path(P$cache, c("H1_models.rds", "H3_models.rds"))
estimate_paths <- file.path(P$estimates, c("H1_curvilinear_results.csv", "H3_results.csv"))
require_files(c(bundle_paths, estimate_paths))
h1_models <- readRDS(bundle_paths[1]); h3_models <- readRDS(bundle_paths[2])
h1_results <- read.csv(estimate_paths[1], stringsAsFactors = FALSE, check.names = FALSE)
h3_results <- read.csv(estimate_paths[2], stringsAsFactors = FALSE, check.names = FALSE)

h1_forest <- do.call(rbind, lapply(seq_len(nrow(h1_results)), function(i) {
  id <- h1_results$model[i]
  model <- h1_models[[id]]
  if (is.null(model)) stop(sprintf("Missing H1 model %s", id), call. = FALSE)
  term <- grep("^I\\(.*\\^2\\)$", names(coef(model)), value = TRUE)
  if (length(term) != 1L) stop(sprintf("H1 model %s lacks one quadratic term", id), call. = FALSE)
  ci <- wald_forest_row(model, term)
  data.frame(model = id, estimator = h1_results$estimator[i], term = term,
             beta = ci["beta"], se = ci["se"], ci_lo = ci["lo"], ci_hi = ci["hi"],
             panel = if (id %in% c("M1_ln", "M1_cox", "M4_AG", "M5_GEE")) "Sensitivity" else "Main",
             stringsAsFactors = FALSE)
}))
h1_forest$model <- factor(h1_forest$model, levels = rev(h1_results$model))
p_c1 <- ggplot2::ggplot(h1_forest, ggplot2::aes(beta, model, colour = estimator, shape = panel)) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_lo, xmax = ci_hi),
                         orientation = "y", width = .22) +
  ggplot2::geom_point(size = 2) + ggplot2::theme_bw(base_size = 9) +
  ggplot2::labs(x = expression(beta[English^2]~"(95% Wald CI)"), y = NULL,
                colour = "Estimator", shape = "Panel")
write.csv(h1_forest, file.path(P$figures, "Figure_C1_plot_data.csv"), row.names = FALSE)
ggplot2::ggsave(file.path(P$figures, "Figure_C1.pdf"), p_c1, width = 9, height = 6.5,
                units = "in", device = grDevices::cairo_pdf)
ggplot2::ggsave(file.path(P$figures, "Figure_C1.png"), p_c1, width = 9, height = 6.5,
                units = "in", dpi = 300)

main_ids <- paste0("M", 1:8)
sensitivity_ids <- c(paste0("M", 13:16, "_AG"), paste0("M", 17:20, "_GEE"))
h3_forest <- do.call(rbind, lapply(seq_len(nrow(h3_results)), function(i) {
  id <- h3_results$model[i]
  model <- h3_models[[id]]
  if (is.null(model)) stop(sprintf("Missing H3 model %s", id), call. = FALSE)
  moderator <- h3_results$moderator[i]
  candidates <- c(paste0("eng_c:", moderator), paste0(moderator, ":eng_c"))
  term <- intersect(candidates, names(coef(model)))
  if (length(term) != 1L) stop(sprintf("H3 model %s lacks interaction %s", id, moderator), call. = FALSE)
  ci <- wald_forest_row(model, term)
  data.frame(model = id, estimator = h3_results$estimator[i], moderator = moderator,
             term = term, beta = ci["beta"], se = ci["se"], ci_lo = ci["lo"], ci_hi = ci["hi"],
             panel = if (id %in% main_ids) "Main" else if (id %in% sensitivity_ids) "Sensitivity" else "Appendix",
             stringsAsFactors = FALSE)
}))
h3_forest$model <- factor(h3_forest$model, levels = rev(unique(h3_results$model)))
p_c2 <- ggplot2::ggplot(h3_forest, ggplot2::aes(beta, model, colour = estimator, shape = panel)) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_lo, xmax = ci_hi),
                         orientation = "y", width = .20) +
  ggplot2::geom_point(size = 1.8) +
  ggplot2::facet_wrap(~ moderator, scales = "free", ncol = 2) +
  ggplot2::theme_bw(base_size = 8.5) +
  ggplot2::labs(x = expression(beta[English%*%Moderator]~"(95% Wald CI)"), y = NULL,
                colour = "Estimator", shape = "Panel")
write.csv(h3_forest, file.path(P$figures, "Figure_C2_plot_data.csv"), row.names = FALSE)
ggplot2::ggsave(file.path(P$figures, "Figure_C2.pdf"), p_c2, width = 11, height = 9,
                units = "in", device = grDevices::cairo_pdf)
ggplot2::ggsave(file.path(P$figures, "Figure_C2.png"), p_c2, width = 11, height = 9,
                units = "in", dpi = 300)
}
cat(if (nzchar(suffix)) "Wrote all three C26 smoke panels.\n" else
      "Wrote C25, all three C26 panels, and Appendix Figures C1-C2.\n")
