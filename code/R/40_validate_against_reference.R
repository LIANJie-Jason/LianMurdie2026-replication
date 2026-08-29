root <- Sys.getenv("REPLICATION_ROOT", unset = "")
if (!nzchar(root)) stop("REPLICATION_ROOT is not set; launch with run_all.R", call. = FALSE)
source(file.path(root, "code", "R", "00_setup.R"))
paths <- init_replication(root)
profile <- Sys.getenv("REPLICATION_PROFILE", unset = "full")
run_started <- as.numeric(Sys.getenv("REPLICATION_RUN_STARTED_EPOCH", unset = "0"))
run_marker <- Sys.getenv("REPLICATION_RUN_MARKER", unset = "")
rep_assert(length(run_started) == 1L && is.finite(run_started) && run_started > 0,
           "Stage 40 requires a current run marker; launch with run_all.R")
rep_assert(nzchar(run_marker) && file.exists(run_marker),
           "Stage 40 cannot verify the current run marker")
rep_assert(abs(as.numeric(file.info(run_marker)$mtime) - run_started) <= 1e-6,
           "Stage 40 run-marker timestamp does not match the execution context")
manifest <- rep_load_reference_manifest(paths)

reference_checks <- lapply(seq_len(nrow(manifest)), function(index) {
  expected <- file.path(paths$root, manifest$reference_path[[index]])
  rep_assert_file(expected, "Accepted reference")
  actual_hash <- rep_sha256_file(expected)
  rep_assert(identical(actual_hash, tolower(manifest$sha256[[index]])),
             "Accepted reference hash mismatch for %s", manifest$reference_path[[index]])
  TRUE
})

active <- manifest$required_profile == "fast" | (profile == "full" & manifest$required_profile == "full")
active <- active & manifest$comparison != "inventory"
targets <- manifest[active, , drop = FALSE]

checks <- lapply(seq_len(nrow(targets)), function(index) {
  entry <- targets[index, , drop = FALSE]
  actual <- file.path(paths$root, entry$output_path)
  expected <- file.path(paths$root, entry$reference_path)
  produced <- file.exists(actual) && as.numeric(file.info(actual)$mtime) >= run_started
  result <- rep_compare_output(
    actual, expected, entry$comparison,
    abs_tol = as.numeric(entry$abs_tolerance), rel_tol = as.numeric(entry$rel_tolerance)
  )
  if (!produced) {
    result$passed <- FALSE
    result$detail <- paste("not produced in this run;", result$detail)
  }
  data.frame(
    output_path = entry$output_path,
    reference_path = entry$reference_path,
    comparison = result$method,
    produced_this_run = produced,
    passed = result$passed,
    blocking = !entry$comparison %in% c("advisory_csv", "advisory_numeric", "display_copy"),
    max_abs_diff = result$max_abs_diff,
    max_rel_diff = result$max_rel_diff,
    detail = result$detail,
    stringsAsFactors = FALSE
  )
})
checks <- if (length(checks)) do.call(rbind, checks) else data.frame()

visual_path <- file.path(paths$audit, "figure_visual_comparison.csv")
rep_assert_file(visual_path, "Advisory figure raster audit")
visual <- read.csv(visual_path, stringsAsFactors = FALSE, check.names = FALSE)
rep_assert_columns(visual, c("figure", "same_dimensions", "mean_absolute_pixel_difference",
                             "predeclared_mean_abs_threshold", "verdict"), "figure raster audit")
visual_expected <- c("Figure_7", "Figure_8", "Figure_9")
rep_assert(identical(as.character(visual$figure), visual_expected),
           "Figure raster audit must contain exactly Figure_7, Figure_8, and Figure_9")
visual_fresh <- as.numeric(file.info(visual_path)$mtime) >= run_started
visual_pass <- visual$same_dimensions %in% TRUE & visual$verdict == "PASS" &
  is.finite(visual$mean_absolute_pixel_difference) &
  visual$mean_absolute_pixel_difference <= visual$predeclared_mean_abs_threshold &
  visual$predeclared_mean_abs_threshold == 0.05
visual_checks <- data.frame(
  output_path = file.path("output", "figures", paste0(visual$figure, ".pdf")),
  reference_path = file.path("reference", "accepted", "figures", c(
    "fig_marginal_h1.pdf", "fig_h22_table3_cols3_4_trade_interaction.pdf",
    "fig_h3_table4_cols1_2_5_6_interactions.pdf")),
  comparison = "pdf_raster_mae<=0.05",
  produced_this_run = rep(visual_fresh, nrow(visual)),
  passed = visual_pass & visual_fresh,
    blocking = FALSE,
  max_abs_diff = visual$max_absolute_pixel_difference,
  max_rel_diff = NA_real_,
  detail = sprintf("raster MAE %.7f <= %.2f; dimensions_equal=%s",
                   visual$mean_absolute_pixel_difference,
                   visual$predeclared_mean_abs_threshold, visual$same_dimensions),
  stringsAsFactors = FALSE
)
checks <- rbind(checks, visual_checks)

contract_row <- function(output_path, reference_path, comparison, fresh, passed,
                         detail, max_abs_diff = NA_real_) {
  data.frame(
    output_path = output_path, reference_path = reference_path,
    comparison = comparison, produced_this_run = fresh,
    passed = isTRUE(passed) && isTRUE(fresh), blocking = TRUE,
    max_abs_diff = max_abs_diff, max_rel_diff = NA_real_, detail = detail,
    stringsAsFactors = FALSE
  )
}
is_fresh <- function(path) file.exists(path) && as.numeric(file.info(path)$mtime) >= run_started

h21_status_path <- file.path(paths$diagnostics, "H21_extensions_logistf_convergence.csv")
h3_status_path <- file.path(paths$diagnostics, "H3_gee_fit_status.csv")
all_fit_status_path <- file.path(paths$diagnostics, "accepted_fit_status_all.csv")
rep_assert_file(h21_status_path, "H2.1 accepted-fit status diagnostic")
rep_assert_file(h3_status_path, "H3 accepted-fit status diagnostic")
rep_assert_file(all_fit_status_path, "Cross-family accepted-fit status diagnostic")
h21_status <- read.csv(h21_status_path, stringsAsFactors = FALSE, check.names = FALSE)
h3_status <- read.csv(h3_status_path, stringsAsFactors = FALSE, check.names = FALSE)
all_fit_status <- read.csv(all_fit_status_path, stringsAsFactors = FALSE, check.names = FALSE)
h21_bad <- h21_status[!(h21_status$converged %in% TRUE), , drop = FALSE]
h3_bad <- h3_status[!(h3_status$converged %in% TRUE), , drop = FALSE]
h21_contract <- nrow(h21_bad) == 1L && h21_bad$model[[1L]] == "M9q" &&
  h21_bad$accepted_values_preserved[[1L]] %in% TRUE &&
  h21_bad$reporting_class[[1L]] == "accepted_archival_nonconverged"
h3_contract <- nrow(h3_bad) == 1L && h3_bad$model[[1L]] == "M20_GEE" &&
  h3_bad$accepted_values_preserved[[1L]] %in% TRUE &&
  h3_bad$reporting_class[[1L]] == "accepted_archival_nonconverged"
flagged_fit_keys <- paste(all_fit_status$family[!all_fit_status$converged],
                          all_fit_status$model[!all_fit_status$converged], sep = "::")
expected_flagged_fit_keys <- c(
  "H1::M1s", "H1::M4s",
  "H21::M1s", "H21::M2s", "H21::M7s", "H21::M8s",
  paste0("H22::M", 1:4, "s"), paste0("H3::M", 1:4, "s"),
  paste0("H22::M", 13:16, "s"), paste0("H3::M", 13:16, "s"),
  "H3::M20_GEE", "H21_extensions::M9q", "H21_extensions::M1qs",
  "H21_extensions::M2qs", "H21_extensions::M7qs",
  "H21_extensions::M8qs", "H21_extensions::M7l", "H21_extensions::M8l"
)
all_fit_contract <- nrow(all_fit_status) == 140L &&
  setequal(flagged_fit_keys, expected_flagged_fit_keys) &&
  all(all_fit_status$accepted_values_preserved %in% TRUE) &&
  all(all_fit_status$reporting_class[all_fit_status$converged] == "accepted_replicated") &&
  all(all_fit_status$reporting_class[!all_fit_status$converged] ==
        "accepted_archival_nonconverged")
checks <- rbind(
  checks,
  contract_row("output/diagnostics/accepted_fit_status_all.csv",
               "README.md#accepted-fit-disclosures", "cross_family_fit_status_contract",
               is_fresh(all_fit_status_path), all_fit_contract,
               "140 supported fits audited; 30 accepted Weibull/Cox/logistf/GEE fits are explicitly archival/nonconverged"),
  contract_row("output/diagnostics/H21_extensions_logistf_convergence.csv",
               "README.md#accepted-fit-disclosures", "accepted_fit_status_contract",
               is_fresh(h21_status_path), h21_contract,
               "M9q is the sole nonconverged H2.1-extension logistf fit and remains archival-only"),
  contract_row("output/diagnostics/H3_gee_fit_status.csv",
               "README.md#accepted-fit-disclosures", "accepted_fit_status_contract",
               is_fresh(h3_status_path), h3_contract,
               "M20_GEE is the sole nonconverged H3 GEE fit and remains archival-only")
)

display_path <- file.path(paths$diagnostics, "display_materialization.csv")
rep_assert_file(display_path, "Accepted display materialization audit")
display_lineage <- read.csv(display_path, stringsAsFactors = FALSE, check.names = FALSE)
rep_assert_columns(display_lineage,
                   c("exhibit_id", "component_id", "reference_path", "output_path",
                     "precopy_origin", "rebuilt_sha256", "rebuilt_exact_string_match",
                     "rebuild_detail", "rebuild_mismatches", "rebuild_mismatch_fit_key",
                     "reference_sha256", "output_sha256", "lineage_class", "status"),
                   "accepted display materialization audit")
expected_display_rows <- if (profile == "full") 32L else 28L
expected_upstream_rows <- if (profile == "full") 4L else 0L
independent_rebuilds <- display_lineage$precopy_origin == "independent_table_builder"
rebuild_nonmatches <- independent_rebuilds & !display_lineage$rebuilt_exact_string_match
rebuild_mismatch_cells <- sum(vapply(
  display_lineage$rebuild_mismatches[rebuild_nonmatches],
  function(value) if (nzchar(value)) length(strsplit(value, " | ", fixed = TRUE)[[1L]]) else 0L,
  integer(1)
))
display_contract <- nrow(display_lineage) == expected_display_rows &&
  sum(independent_rebuilds) == 28L &&
  sum(display_lineage$precopy_origin == "accepted_source_materialized_upstream") ==
    expected_upstream_rows &&
  all(nzchar(display_lineage$rebuilt_sha256)) &&
  all(nzchar(display_lineage$rebuild_detail)) &&
  setequal(display_lineage$exhibit_id[rebuild_nonmatches], c("C4", "C8")) &&
  rebuild_mismatch_cells == 3L &&
  all(nzchar(display_lineage$rebuild_mismatches[rebuild_nonmatches])) &&
  all(display_lineage$rebuild_mismatch_fit_key[rebuild_nonmatches] %in% flagged_fit_keys) &&
  all(display_lineage$lineage_class == "accepted_display_copy") &&
  all(display_lineage$status == "PASS") &&
  all(display_lineage$reference_sha256 == display_lineage$output_sha256) &&
  all(vapply(seq_len(nrow(display_lineage)), function(index) {
    output <- file.path(paths$root, display_lineage$output_path[[index]])
    file.exists(output) &&
      identical(rep_sha256_file(output), display_lineage$output_sha256[[index]]) &&
      as.numeric(file.info(output)$mtime) >= run_started
  }, logical(1)))
checks <- rbind(
  checks,
  contract_row("output/diagnostics/display_materialization.csv",
               "reference/accepted/exhibit_manifest.csv",
               "accepted_display_materialization_contract",
               is_fresh(display_path), display_contract,
               sprintf(paste0("%d accepted displays copied byte-exactly; 28 independent ",
                              "builder bodies audited before copy (%d exact; 3 bounded nuisance-term ",
                              "cells in archival fits C4/C8 preserved from accepted display)"),
                       expected_display_rows,
                       sum(display_lineage$rebuilt_exact_string_match & independent_rebuilds)))
)

c22_output_path <- file.path(paths$appendix_tables, "C22_nonviolent_final_body.csv")
rep_assert_file(c22_output_path, "Fresh C22 output")
c22_output <- read.csv(c22_output_path, stringsAsFactors = FALSE, check.names = FALSE)
c22_contract <- nrow(c22_output) == 31L && "Model" %in% names(c22_output) &&
  sum(c22_output$Model == "NVH1", na.rm = TRUE) == 1L &&
  identical(as.character(c22_output$Model[[1L]]), "NVH1")
checks <- rbind(
  checks,
  contract_row("output/tables/appendix/C22_nonviolent_final_body.csv",
               "reference/accepted/exhibit_manifest.csv", "d69_accepted_display_contract",
               is_fresh(c22_output_path), c22_contract,
               "Accepted C22 display copy has 31 rows with exactly one leading NVH1; the pre-copy builder contract is recorded separately")
)

if (profile == "full") {
  c26_diagnostic_outputs <- file.path(paths$diagnostics, c(
    "inference_diagnostics_cox_ph_table_body.csv",
    "inference_diagnostics_cox_power_quarantine_table_body.csv",
    "inference_diagnostics_weibull_table_body.csv"
  ))
  c26_diagnostic_references <- file.path(paths$reference, "appendix_tables", c(
    "inference_diagnostics_cox_ph_table_body.csv",
    "inference_diagnostics_cox_power_quarantine_table_body.csv",
    "inference_diagnostics_weibull_table_body.csv"
  ))
  c26_contract_names <- c(
    "fresh_c26_cox_ph_numeric_contract",
    "accepted_c26_quarantine_copy_numeric_contract",
    "accepted_c26_bootstrap_display_numeric_contract"
  )
  for (index in seq_along(c26_diagnostic_outputs)) {
    comparison <- rep_compare_csv_numeric(c26_diagnostic_outputs[[index]],
                                          c26_diagnostic_references[[index]],
                                          abs_tol = 5e-4, rel_tol = 0)
    checks <- rbind(
      checks,
      contract_row(rep_relative_path(c26_diagnostic_outputs[[index]], paths$root),
                   rep_relative_path(c26_diagnostic_references[[index]], paths$root),
                   c26_contract_names[[index]],
                   is_fresh(c26_diagnostic_outputs[[index]]), comparison$passed,
                   comparison$detail, comparison$max_abs_diff)
    )
  }

  c25_output_path <- file.path(paths$appendix_tables, "C25_power_table_body.csv")
  rep_assert_file(c25_output_path, "Fresh C25 display output")
  c25_output <- read.csv(c25_output_path, stringsAsFactors = FALSE, check.names = FALSE)
  c25_columns <- c("Family", "Main table", "Main-table model label", "Estimator",
                   "Focal test", "Power at fitted effect", "False-positive rate",
                   "Status / note")
  c25_contract <- nrow(c25_output) == 28L && identical(names(c25_output), c25_columns) &&
    !"Source" %in% names(c25_output) &&
    match("Power at fitted effect", names(c25_output)) <
      match("False-positive rate", names(c25_output))
  checks <- rbind(
    checks,
    contract_row("output/tables/appendix/C25_power_table_body.csv",
                 "reference/accepted/exhibit_manifest.csv", "d70_accepted_display_contract",
                 is_fresh(c25_output_path), c25_contract,
                 "Accepted C25 display copy has 28 rows, accepted column order, and no Source column")
  )

  materialization_path <- file.path(paths$diagnostics, "stochastic_materialization.csv")
  rep_assert_file(materialization_path, "Accepted stochastic materialization audit")
  materialization <- read.csv(materialization_path, stringsAsFactors = FALSE,
                              check.names = FALSE)
  rep_assert_columns(materialization,
                     c("artifact", "reference_path", "output_path", "reference_sha256",
                       "output_sha256", "lineage_class", "status"),
                     "accepted stochastic materialization audit")
  expected_stochastic <- c(
    "power_h1.csv", "power_h1_summary.csv", "power_h21.csv", "power_h21_summary.csv",
    "power_h22.csv", "power_h22_summary.csv", "power_h3.csv", "power_h3_summary.csv",
    "weibull_bootstrap_calibration.csv", "weibull_power_diagnostic_audit.csv",
    "cox_power_quarantine_audit.csv"
  )
  materialization_contract <- nrow(materialization) == 11L &&
    identical(as.character(materialization$artifact), expected_stochastic) &&
    all(materialization$lineage_class == "accepted_codestack_copy") &&
    all(materialization$status == "PASS") &&
    all(materialization$reference_sha256 == materialization$output_sha256) &&
    all(vapply(seq_len(nrow(materialization)), function(index) {
      output <- file.path(paths$root, materialization$output_path[[index]])
      file.exists(output) &&
        identical(rep_sha256_file(output), materialization$output_sha256[[index]]) &&
        as.numeric(file.info(output)$mtime) >= run_started
    }, logical(1)))
  checks <- rbind(
    checks,
    contract_row("output/diagnostics/stochastic_materialization.csv",
                 "reference/accepted/provenance/accepted_stochastic_artifacts.csv",
                 "accepted_codestack_materialization_contract",
                 is_fresh(materialization_path), materialization_contract,
                 "11 authoritative accepted stochastic artifacts were copied byte-exactly with explicit lineage")
  )
}
report_path <- file.path(paths$audit, "reference_validation.csv")
write.csv(checks, report_path, row.names = FALSE, na = "")

all_passed <- nrow(checks) > 0 && all(checks$passed | !checks$blocking)
final_status <- if (!all_passed) "FAIL" else if (profile == "full") "PASS" else "INCOMPLETE"
summary <- data.frame(profile = profile, checked = nrow(checks), passed = sum(checks$passed),
                      failed = sum(!checks$passed), blocking_failed = sum(!checks$passed & checks$blocking),
                      advisory_failed = sum(!checks$passed & !checks$blocking), final_status = final_status,
                      stringsAsFactors = FALSE)
write.csv(summary, file.path(paths$audit, "validation_summary.csv"), row.names = FALSE, na = "")
rep_assert(all_passed, "Reference validation failed for %d blocking output(s); see %s",
           sum(!checks$passed & checks$blocking), report_path)
cat(sprintf("Reference validation: %s (%d/%d artifacts passed)\n", final_status, sum(checks$passed), nrow(checks)))
