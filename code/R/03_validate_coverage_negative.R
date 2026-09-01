root <- Sys.getenv("REPLICATION_ROOT", unset = "")
if (!nzchar(root)) stop("REPLICATION_ROOT is not set; launch with run_all.R", call. = FALSE)
source(file.path(root, "code", "R", "00_setup.R"))
source(file.path(root, "code", "R", "R", "power_figure_helpers.R"))
paths <- init_replication(root)

coverage_rows <- list()
record_coverage <- function(check, passed, detail) {
  coverage_rows[[length(coverage_rows) + 1L]] <<- data.frame(
    check_id = check, passed = isTRUE(passed), detail = detail, stringsAsFactors = FALSE)
}

exhibit_path <- file.path(paths$reference, "exhibit_manifest.csv")
rep_assert_file(exhibit_path, "Accepted exhibit manifest")
exhibits <- read.csv(exhibit_path, stringsAsFactors = FALSE, check.names = FALSE)
required_fields <- c("exhibit_id", "component_id", "artifact_kind", "reference_path",
                     "reference_sha256", "output_path", "validation", "abs_tolerance",
                     "rel_tolerance", "provenance_artifact", "provenance_reference")
rep_assert_columns(exhibits, required_fields, "accepted exhibit manifest")
record_coverage("exact_exhibit_ids",
                identical(sort(unique(exhibits$exhibit_id)), sort(rep_expected_exhibit_ids())),
                sprintf("found %d unique exhibits; expected %d",
                        length(unique(exhibits$exhibit_id)), length(rep_expected_exhibit_ids())))

expected_components <- setNames(rep(1L, length(rep_expected_exhibit_ids())), rep_expected_exhibit_ids())
expected_components[paste0("Figure_", c(7:9, "C1", "C2"))] <- 2L
expected_components["C26"] <- 3L
actual_components <- table(exhibits$exhibit_id)
component_ok <- identical(as.integer(actual_components[names(expected_components)]),
                          as.integer(expected_components))
record_coverage("exact_component_counts", component_ok,
                sprintf("manifest contains %d exhibit components", nrow(exhibits)))

component_keys <- paste(exhibits$exhibit_id, exhibits$component_id, sep = "::")
record_coverage("unique_components", !anyDuplicated(component_keys),
                "exhibit/component keys are unique")

hash_ok <- vapply(seq_len(nrow(exhibits)), function(index) {
  path <- file.path(paths$root, exhibits$reference_path[[index]])
  file.exists(path) && identical(rep_sha256_file(path), exhibits$reference_sha256[[index]])
}, logical(1))
record_coverage("component_hashes", all(hash_ok),
                sprintf("%d/%d exhibit component hashes match", sum(hash_ok), length(hash_ok)))

reference_manifest <- rep_load_reference_manifest(paths)
artifact_reference <- grepl("reference/accepted/(main_tables|appendix_tables|figures|figure_data)/",
                            reference_manifest$reference_path)
record_coverage("no_orphan_exhibit_artifacts",
                setequal(reference_manifest$reference_path[artifact_reference], exhibits$reference_path),
                "accepted exhibit artifact namespace equals exhibit manifest exactly")
record_coverage("output_mapping",
                all(exhibits$output_path == reference_manifest$output_path[
                  match(exhibits$reference_path, reference_manifest$reference_path)] |
                    exhibits$validation == "render"),
                "numerical exhibit output mappings agree with reference manifest")

c20 <- exhibits[exhibits$exhibit_id == "C20", , drop = FALSE]
c20_contract_path <- if (nrow(c20) == 1L) file.path(paths$root, c20$provenance_reference) else ""
c20_contract <- if (nzchar(c20_contract_path) && file.exists(c20_contract_path)) {
  read.csv(c20_contract_path, stringsAsFactors = FALSE, check.names = FALSE)
} else data.frame(field = character(), value = character())
c20_values <- if (all(c("field", "value") %in% names(c20_contract))) {
  setNames(c20_contract$value, c20_contract$field)
} else character()
c20_get <- function(field, default = "") {
  if (field %in% names(c20_values)) c20_values[[field]] else default
}
c20_required <- c(
  "exhibit", "artifact_path", "artifact_sha256", "artifact_shape", "variable",
  "source", "transformation", "country_mapping", "campaign_window",
  "aggregation", "join_key", "model_producer", "table_producer",
  "model_reference", "table_reference"
)
c20_contract_fields_ok <- identical(as.character(c20_contract$field), c20_required)
c20_artifact_path <- c20_get("artifact_path")
c20_artifact_file <- file.path(paths$root, c20_artifact_path)
c20_input_manifest_path <- file.path(paths$data, "input_manifest.csv")
c20_input_manifest <- if (file.exists(c20_input_manifest_path)) {
  read.csv(c20_input_manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
} else data.frame()
c20_input <- if ("file" %in% names(c20_input_manifest)) {
  c20_input_manifest[c20_input_manifest$file == basename(c20_artifact_path), , drop = FALSE]
} else data.frame()
c20_required_columns <- if (nrow(c20_input) == 1L &&
                                "required_columns" %in% names(c20_input)) {
  strsplit(c20_input$required_columns[[1L]], ";", fixed = TRUE)[[1L]]
} else character()
c20_model_producer <- c20_get("model_producer")
c20_table_producer <- c20_get("table_producer")
c20_model_reference <- c20_get("model_reference")
c20_table_reference <- c20_get("table_reference")
c20_model_manifest <- reference_manifest[
  reference_manifest$reference_path == c20_model_reference, , drop = FALSE
]
c20_ok <- nrow(c20) == 1L && identical(c20$component_id[[1L]], "body") &&
  nzchar(c20$provenance_artifact[[1L]]) && nzchar(c20$provenance_reference[[1L]]) &&
  c20_contract_fields_ok &&
  identical(c20_get("exhibit"), "Appendix C20") &&
  identical(c20_artifact_path, "data/h3_additional_moderators.csv") &&
  identical(c20_artifact_path, c20$provenance_artifact[[1L]]) &&
  file.exists(c20_artifact_file) &&
  identical(rep_sha256_file(c20_artifact_file), c20_get("artifact_sha256")) &&
  nrow(c20_input) == 1L &&
  identical(tolower(c20_input$sha256[[1L]]), c20_get("artifact_sha256")) &&
  identical(sprintf("%d rows x %d columns", c20_input$rows[[1L]],
                    c20_input$columns[[1L]]), c20_get("artifact_shape")) &&
  identical(c20_input$key_columns[[1L]], "CAMPAIGN;LOCATION") &&
  "v2x_frassoc_thick" %in% c20_required_columns &&
  identical(c20_get("variable"), "v2x_frassoc_thick") &&
  identical(c20_get("source"), "V-Dem v15") &&
  identical(c20_get("transformation"),
            "campaign-period mean frozen before publication replication") &&
  identical(c20_get("country_mapping"),
            "countrycode COW numeric to ISO3; COW 345=SRB; COW 347=XKX") &&
  identical(c20_get("campaign_window"),
            "max(2009,BYEAR) through min(2019,EYEAR), inclusive") &&
  identical(c20_get("aggregation"),
            "arithmetic mean over nonmissing campaign-years") &&
  identical(c20_get("join_key"), "CAMPAIGN;LOCATION;iso3;BYEAR;EYEAR") &&
  identical(c20_model_producer, "code/R/15_fit_h3_extensions.R") &&
  identical(c20_table_producer, "code/R/34_build_appendix_h3.R") &&
  identical(c20_model_reference,
            "reference/accepted/model_results/H3_table4_additional_moderators_results.csv") &&
  identical(c20_table_reference,
            "reference/accepted/appendix_tables/H3_additional_moderators_table_body.csv") &&
  file.exists(file.path(paths$root, c20_model_producer)) &&
  file.exists(file.path(paths$root, c20_table_producer)) &&
  file.exists(file.path(paths$root, c20_model_reference)) &&
  file.exists(file.path(paths$root, c20_table_reference)) &&
  identical(c20_table_reference, c20$reference_path[[1L]]) &&
  nrow(c20_model_manifest) == 1L &&
  identical(c20_model_manifest$output_path[[1L]],
            "output/estimates/H3_table4_additional_moderators_results.csv") &&
  identical(c20_model_manifest$comparison[[1L]], "numeric")
record_coverage(
  "c20_frozen_input_integrity", c20_ok,
  "C20 frozen auxiliary input is hash-locked in the input manifest and linked to its model and table producers/references"
)

lineage_path <- file.path(paths$reference, "provenance", "archival_passthroughs.csv")
lineage <- if (file.exists(lineage_path)) {
  read.csv(lineage_path, stringsAsFactors = FALSE, check.names = FALSE)
} else data.frame()
lineage_required <- c("artifact", "reference_path", "reference_sha256", "lineage_class",
                      "independent_evidence", "interpretation")
lineage_ids <- c("C25_display", "C26_cox_quarantine_panel", "Cox_quarantine_audit")
lineage_ok <- all(lineage_required %in% names(lineage)) &&
  identical(as.character(lineage$artifact), lineage_ids) &&
  all(lineage$lineage_class == "frozen_passthrough") &&
  all(vapply(seq_len(nrow(lineage)), function(index) {
    artifact_path <- file.path(paths$root, lineage$reference_path[[index]])
    file.exists(artifact_path) &&
      identical(rep_sha256_file(artifact_path), lineage$reference_sha256[[index]])
  }, logical(1)))
record_coverage("archival_passthrough_lineage", lineage_ok,
                "C25 and Cox quarantine archives are hash-linked and never labeled as fresh estimates")
lineage_manifest_rows <- match(lineage$reference_path, reference_manifest$reference_path)
passthrough_nonblocking <- lineage_ok && !anyNA(lineage_manifest_rows) &&
  all(reference_manifest$comparison[lineage_manifest_rows] %in%
        c("advisory_csv", "display_copy", "inventory"))
record_coverage("archival_passthrough_nonblocking", passthrough_nonblocking,
                "Every frozen_passthrough reference is registered as advisory or inventory-only")

derived_path <- file.path(paths$reference, "provenance", "derived_numeric_references.csv")
derived <- if (file.exists(derived_path)) {
  read.csv(derived_path, stringsAsFactors = FALSE, check.names = FALSE)
} else data.frame()
derived_required <- c("artifact_id", "reference_path", "reference_sha256", "lineage_class",
                      "producer_path", "producer_sha256", "accepted_anchor",
                      "accepted_anchor_sha256", "validation_role", "disclosure")
derived_ids <- c("Figure_7_plot_data", "Figure_8_plot_data", "Figure_9_plot_data",
                 "Figure_C1_plot_data", "Figure_C2_plot_data",
                 "C9_stable_id_results", "C12_stable_id_results",
                 "C26_cox_ph_diagnostics")
derived_ok <- all(derived_required %in% names(derived)) &&
  identical(as.character(derived$artifact_id), derived_ids) &&
  all(derived$lineage_class == "frozen_derived_reference") &&
  all(vapply(seq_len(nrow(derived)), function(index) {
    reference_path <- file.path(paths$root, derived$reference_path[[index]])
    producer_path <- file.path(paths$root, derived$producer_path[[index]])
    anchor_path <- file.path(paths$root, derived$accepted_anchor[[index]])
    all(file.exists(reference_path), file.exists(producer_path), file.exists(anchor_path)) &&
      identical(rep_sha256_file(reference_path), derived$reference_sha256[[index]]) &&
      identical(rep_sha256_file(producer_path), derived$producer_sha256[[index]]) &&
      identical(rep_sha256_file(anchor_path), derived$accepted_anchor_sha256[[index]])
  }, logical(1)))
figure_data <- exhibits[exhibits$artifact_kind == "figure_data", , drop = FALSE]
expected_figure_data_ids <- c(paste0("Figure_", 7:9), "Figure_C1", "Figure_C2")
figure_provenance_ok <- nrow(figure_data) == 5L &&
  setequal(figure_data$exhibit_id, expected_figure_data_ids) &&
  all(nzchar(figure_data$provenance_artifact)) &&
  all(figure_data$provenance_reference ==
        "reference/accepted/provenance/derived_numeric_references.csv")
record_coverage("derived_reference_lineage", derived_ok && figure_provenance_ok,
                "Figure 7-9/C1-C2 plot data and C9/C12/C26 derived oracles are hash-linked and disclosed")

seed_lineage_path <- file.path(paths$reference, "provenance", "power_seed_lineage.csv")
seed_lineage <- if (file.exists(seed_lineage_path)) {
  read.csv(seed_lineage_path, stringsAsFactors = FALSE, check.names = FALSE)
} else data.frame()
seed_lineage_required <- c("artifact", "accepted_generator_sha256", "accepted_seed_contract",
                           "fresh_generator", "fresh_generator_sha256",
                           "fresh_seed_contract", "relationship")
seed_lineage_ids <- c("power_h1.csv", "power_h21.csv", "power_h22.csv", "power_h3.csv",
                      "weibull_bootstrap_calibration.csv")
seed_lineage_ok <- all(seed_lineage_required %in% names(seed_lineage)) &&
  identical(as.character(seed_lineage$artifact), seed_lineage_ids) &&
  all(grepl("distinct pseudorandom streams", seed_lineage$relationship, fixed = TRUE)) &&
  all(vapply(seq_len(nrow(seed_lineage)), function(index) {
    producer <- file.path(paths$root, seed_lineage$fresh_generator[[index]])
    file.exists(producer) &&
      identical(rep_sha256_file(producer), seed_lineage$fresh_generator_sha256[[index]])
  }, logical(1)))
record_coverage("power_seed_lineage", seed_lineage_ok,
                "Accepted seed-24 streams and fresh cell-hash streams are explicitly distinct")

stochastic_path <- file.path(paths$reference, "provenance", "accepted_stochastic_artifacts.csv")
stochastic <- if (file.exists(stochastic_path)) {
  read.csv(stochastic_path, stringsAsFactors = FALSE, check.names = FALSE)
} else data.frame()
stochastic_required <- c("artifact", "reference_path", "sha256", "producer_path",
                         "producer_sha256", "accepted_producer_sha256", "lineage_class")
stochastic_ids <- c(
  "power_h1.csv", "power_h1_summary.csv", "power_h21.csv", "power_h21_summary.csv",
  "power_h22.csv", "power_h22_summary.csv", "power_h3.csv", "power_h3_summary.csv",
  "weibull_bootstrap_calibration.csv", "weibull_power_diagnostic_audit.csv",
  "cox_power_quarantine_audit.csv"
)
stochastic_producer_paths <- file.path(paths$root, stochastic$producer_path)
stochastic_ok <- all(stochastic_required %in% names(stochastic)) &&
  identical(as.character(stochastic$artifact), stochastic_ids) &&
  all(stochastic$lineage_class == "accepted_codestack_copy") &&
  all(grepl("^[0-9a-f]{64}$", stochastic$accepted_producer_sha256)) &&
  all(vapply(seq_len(nrow(stochastic)), function(index) {
    reference <- file.path(paths$root, stochastic$reference_path[[index]])
    producer <- stochastic_producer_paths[[index]]
    reference_ok <- file.exists(reference) &&
      identical(rep_sha256_file(reference), stochastic$sha256[[index]])
    producer_ok <- file.exists(producer) &&
      identical(rep_sha256_file(producer), stochastic$producer_sha256[[index]])
    reference_ok && producer_ok
  }, logical(1)))
record_coverage("accepted_stochastic_codestack", stochastic_ok,
                "11 accepted stochastic references and their archived producer scripts hash-match")

seed_registry <- assert_unique_power_seeds()
record_coverage("unique_rng_streams", !anyDuplicated(seed_registry$seed),
                sprintf("%d power/calibration design cells have distinct deterministic seeds",
                        nrow(seed_registry)))

coverage <- do.call(rbind, coverage_rows)
write.csv(coverage, file.path(paths$audit, "rep_coverage.csv"), row.names = FALSE, na = "")

negative_rows <- list()
record_negative <- function(check, passed, detail) {
  negative_rows[[length(negative_rows) + 1L]] <<- data.frame(
    check_id = check, passed = isTRUE(passed), detail = detail, stringsAsFactors = FALSE)
}

accepted_files <- list.files(paths$reference, recursive = TRUE, full.names = FALSE)
output_files <- list.files(file.path(paths$root, "output"), recursive = TRUE, full.names = FALSE)
deleted_pattern <- "identification_(nonviolent|sensemakr)"
deleted_hits <- c(grep(deleted_pattern, accepted_files, value = TRUE),
                  grep(deleted_pattern, output_files, value = TRUE))
record_negative("deleted_identification_absent", !length(deleted_hits),
                if (!length(deleted_hits)) "deleted identification artifacts absent" else
                  paste(deleted_hits, collapse = "; "))

h3_reference <- rep_read_csv(file.path(paths$reference, "model_results", "H3_results.csv"))
h3_main <- h3_reference[h3_reference$main_panel %in% TRUE, , drop = FALSE]
expected_h3 <- c(M1 = "repress_c", M2 = "civlib_c", M3 = "clpol_c", M4 = "polity_c",
                 M5 = "repress_c", M6 = "civlib_c", M7 = "clpol_c", M8 = "polity_c")
h3_ok <- nrow(h3_main) == 8L && identical(as.character(h3_main$model), names(expected_h3)) &&
  identical(as.character(h3_main$moderator), unname(expected_h3))
h3_table <- rep_read_csv(file.path(paths$reference, "main_tables", "H3_main_table_body.csv"))
h3_movement <- any(grepl("movement|dmove", unlist(h3_table, use.names = FALSE), ignore.case = TRUE))
record_negative("no_stale_h3_movement_slot", h3_ok && !h3_movement,
                "H3 main panel is exactly repression/civil-liberty/political-liberty/Polity")

h21_reference <- rep_read_csv(file.path(paths$reference, "model_results", "H21_results.csv"))
h21_main <- h21_reference[h21_reference$main_panel %in% TRUE, , drop = FALSE]
h21_appendix <- h21_reference[!(h21_reference$main_panel %in% TRUE) &
                                 !(h21_reference$robustness_row %in% TRUE), , drop = FALSE]
h21_ids <- c("M1", "M2", "M3", "M4", "M7", "M8", "M9", "M10")
h21_shape_ok <- nrow(h21_main) == 8L && identical(as.character(h21_main$model), h21_ids) &&
  nrow(h21_appendix) == 12L && sum(!is.na(h21_reference$p_int_holm)) == 8L &&
  sum(!is.na(h21_reference$p_int_bh)) == 12L

h21_code <- file.path(paths$code, "11_fit_h21.R")
parsed_text <- paste(vapply(parse(h21_code, keep.source = FALSE), function(expression) {
  paste(deparse(expression, width.cutoff = 500L), collapse = " ")
}, character(1)), collapse = "\n")
old_logic <- grepl("k_main *<- *6|k_app *<- *14|n *= *6|n *= *14", parsed_text)
record_negative("no_h21_k6_k14_logic", h21_shape_ok && !old_logic,
                "H2.1 uses exactly k=8 main and k=12 appendix multiplicity families")

c22_reference <- rep_read_csv(file.path(paths$reference, "appendix_tables",
                                        "C22_nonviolent_final_body.csv"))
c22_contract <- function(data) nrow(data) == 31L && "Model" %in% names(data) &&
  sum(data$Model == "NVH1", na.rm = TRUE) == 1L &&
  identical(as.character(data$Model[[1L]]), "NVH1")
c22_mutant <- c22_reference[c22_reference$Model != "NVH1", , drop = FALSE]
c22_ok <- c22_contract(c22_reference) && !c22_contract(c22_mutant)
record_negative("c22_nvh1_contract", c22_ok,
                "Accepted C22 requires 31 rows with exactly one leading NVH1 row; deletion mutant rejected")

c25_reference <- rep_read_csv(file.path(paths$reference, "appendix_tables",
                                        "C25_power_final_body.csv"))
c25_columns <- c("Family", "Main table", "Main-table model label", "Estimator",
                 "Focal test", "Power at fitted effect", "False-positive rate",
                 "Status / note")
c25_contract <- function(data) nrow(data) == 28L && identical(names(data), c25_columns) &&
  !"Source" %in% names(data) &&
  match("Power at fitted effect", names(data)) <
    match("False-positive rate", names(data))
c25_mutant <- c25_reference
c25_mutant$Source <- "obsolete"
c25_ok <- c25_contract(c25_reference) && !c25_contract(c25_mutant)
record_negative("c25_column_contract", c25_ok,
                "Accepted C25 requires 28 rows, power before false-positive rate, and no Source column; Source mutant rejected")

negative <- do.call(rbind, negative_rows)
write.csv(negative, file.path(paths$audit, "rep_negative.csv"), row.names = FALSE, na = "")
rep_assert(all(coverage$passed), "REP-COVERAGE failed: %s",
           paste(coverage$check_id[!coverage$passed], collapse = ", "))
rep_assert(all(negative$passed), "REP-NEGATIVE failed: %s",
           paste(negative$check_id[!negative$passed], collapse = ", "))
cat(sprintf("REP-COVERAGE and REP-NEGATIVE passed: %d coverage + %d negative checks\n",
            nrow(coverage), nrow(negative)))
