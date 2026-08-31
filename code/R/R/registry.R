rep_load_reference_manifest <- function(paths) {
  manifest_path <- file.path(paths$reference, "reference_manifest.csv")
  rep_assert_file(manifest_path, "Reference manifest")
  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("reference_path", "output_path", "sha256", "comparison", "abs_tolerance", "rel_tolerance")
  rep_assert_columns(manifest, required, "reference manifest")
  rep_assert(!anyDuplicated(manifest$reference_path), "Reference manifest has duplicate reference_path values")
  outputs <- manifest$output_path[nzchar(manifest$output_path)]
  rep_assert(!anyDuplicated(outputs), "Reference manifest has duplicate non-empty output_path values")
  manifest
}

rep_script_registry <- function() {
  c(
    "01_preflight.R", "02_validate_inputs.R", "03_validate_coverage_negative.R",
    "10_fit_h1.R", "11_fit_h21.R", "12_fit_h22.R", "13_fit_h3.R",
    "14_fit_h21_extensions.R", "15_fit_h3_extensions.R", "16_audit_fit_status.R",
    "20_fit_nonviolent.R", "21_fit_ongoing.R", "22_fit_sensemakr.R",
    "23_materialize_accepted_stochastic.R", "29_inference_diagnostics.R",
    "30_build_main_tables.R", "31_build_main_figures.R", "32_build_appendix_h1_h21.R",
    "33_build_appendix_h22.R", "34_build_appendix_h3.R", "35_build_appendix_general.R",
    "36_build_appendix_power_inference.R", "37_materialize_accepted_displays.R",
    "40_validate_against_reference.R", "41_write_output_manifest.R"
  )
}

rep_expected_exhibit_ids <- function() {
  c(paste0("Table_", 1:4), paste0("Figure_", 7:9), paste0("C", 1:26),
    "Figure_C1", "Figure_C2")
}

rep_runtime_cache_registry <- function() {
  c(
    "output/cache/H1_models.rds", "output/cache/H1_analysis_frames.rds",
    "output/cache/H21_models.rds", "output/cache/H21_analysis_frames.rds",
    "output/cache/H22_models.rds", "output/cache/H22_analysis_frames.rds",
    "output/cache/H22_quadratic_models.rds",
    "output/cache/H3_models.rds", "output/cache/H3_analysis_frames.rds",
    "output/cache/H21_extensions_models.rds",
    "output/cache/H3_additional_moderators_models.rds",
    "output/cache/nonviolent_focal_models.rds", "output/cache/ongoing_models.rds",
    "output/cache/sensemakr_models.rds"
  )
}

rep_auxiliary_output_registry <- function() {
  c(
    "output/diagnostics/accepted_fit_status_all.csv",
    "output/diagnostics/H21_extensions_logistf_convergence.csv",
    "output/diagnostics/H3_gee_fit_status.csv",
    "output/diagnostics/display_materialization.csv",
    paste0("output/figures/Figure_", rep(c("7", "8", "9"), each = 2L),
           rep(c(".pdf", ".png"), 3L)),
    paste0("output/figures/Figure_", rep(c("C1", "C2"), each = 2L),
           rep(c(".pdf", ".png"), 2L)),
    "output/diagnostics/stochastic_materialization.csv",
    "output/diagnostics/inference_diagnostics_cox_ph_table_body.csv",
    "output/diagnostics/inference_diagnostics_cox_power_quarantine_table_body.csv",
    "output/diagnostics/inference_diagnostics_weibull_table_body.csv"
  )
}

rep_smoke_output_registry <- function() {
  c(
    paste0("output/estimates/power_", c("h1", "h21", "h22", "h3"), "_smoke.csv"),
    paste0("output/diagnostics/power_", c("h1", "h21", "h22", "h3"),
           "_smoke_summary.csv"),
    "output/diagnostics/power_reference_comparison_smoke.csv",
    "output/diagnostics/cox_power_quarantine_audit_smoke.csv",
    "output/diagnostics/weibull_bootstrap_calibration_smoke.csv",
    "output/diagnostics/weibull_power_diagnostic_audit_smoke.csv",
    "output/diagnostics/inference_diagnostics_cox_ph_smoke.csv",
    "output/diagnostics/inference_diagnostics_cox_ph_table_body_smoke.csv",
    "output/diagnostics/inference_diagnostics_cox_power_quarantine_table_body_smoke.csv",
    "output/diagnostics/inference_diagnostics_weibull_table_body_smoke.csv",
    paste0("output/tables/appendix/C26_panel_", 1:3, "_smoke_table_body.csv")
  )
}

rep_optional_recomputed_output_registry <- function() {
  c(
    paste0("output/estimates/power_", c("h1", "h21", "h22", "h3"),
           "_recomputed.csv"),
    paste0("output/diagnostics/power_", c("h1", "h21", "h22", "h3"),
           "_recomputed_summary.csv"),
    "output/diagnostics/power_reference_comparison_recomputed.csv",
    "output/diagnostics/weibull_bootstrap_calibration_recomputed.csv",
    "output/diagnostics/weibull_power_diagnostic_audit_recomputed.csv",
    "output/diagnostics/cox_power_quarantine_audit_recomputed.csv"
  )
}

rep_active_output_registry <- function(paths) {
  references <- rep_load_reference_manifest(paths)
  reference_active <- nzchar(references$output_path) & references$comparison != "inventory"
  paths_out <- sort(unique(c(references$output_path[reference_active],
                             rep_auxiliary_output_registry())))
  rep_assert(length(paths_out) > 0L, "Output registry is empty")
  paths_out
}

rep_all_known_outputs <- function(paths) {
  retired <- "output/diagnostics/power_reference_comparison.csv"
  sort(unique(c(rep_active_output_registry(paths), rep_runtime_cache_registry(),
                rep_smoke_output_registry(), rep_optional_recomputed_output_registry(), retired)))
}
