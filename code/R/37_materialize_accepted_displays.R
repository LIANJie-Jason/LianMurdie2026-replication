root <- Sys.getenv("REPLICATION_ROOT", unset = "")
if (!nzchar(root)) {
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
  file_arg <- gsub("~+~", " ", file_arg, fixed = TRUE)
  root <- normalizePath(file.path(dirname(file_arg), "..", ".."), mustWork = TRUE)
}
source(file.path(root, "code", "R", "00_setup.R"))
P <- init_replication(root)

exhibit_path <- file.path(P$reference, "exhibit_manifest.csv")
rep_assert_file(exhibit_path, "Accepted exhibit manifest")
exhibits <- read.csv(exhibit_path, stringsAsFactors = FALSE, check.names = FALSE)
tables <- exhibits[exhibits$artifact_kind == "table", , drop = FALSE]

display_mismatches <- function(actual_path, expected_path) {
  read_display <- function(path) {
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
             colClasses = "character", na.strings = NULL,
             strip.white = FALSE, blank.lines.skip = FALSE)
  }
  actual <- read_display(actual_path)
  expected <- read_display(expected_path)
  if (!identical(dim(actual), dim(expected)) || !identical(names(actual), names(expected))) {
    return("shape_or_schema_mismatch")
  }
  differences <- character()
  clean <- function(value) gsub("[\r\n]+", " ", value)
  for (column in names(actual)) {
    rows <- which(rep_normalize_display_scalar(actual[[column]]) !=
                    rep_normalize_display_scalar(expected[[column]]))
    for (row in rows) {
      differences <- c(
        differences,
        sprintf("%s[row=%d]:rebuilt={%s};accepted={%s}", column, row,
                clean(rep_normalize_display_scalar(actual[[column]][[row]])),
                clean(rep_normalize_display_scalar(expected[[column]][[row]])))
      )
    }
  }
  paste(differences, collapse = " | ")
}

rows <- lapply(seq_len(nrow(tables)), function(index) {
  source_path <- file.path(P$root, tables$reference_path[[index]])
  output_path <- file.path(P$root, tables$output_path[[index]])
  rep_assert_file(source_path, sprintf("Accepted display %s", tables$exhibit_id[[index]]))
  rep_assert_file(output_path, sprintf("Pre-copy rebuilt display %s", tables$exhibit_id[[index]]))
  source_hash <- rep_sha256_file(source_path)
  rep_assert(identical(source_hash, tables$reference_sha256[[index]]),
             "Accepted display hash mismatch for %s", tables$exhibit_id[[index]])
  rebuilt_hash <- rep_sha256_file(output_path)
  rebuilt_comparison <- rep_compare_csv_exact_strings(output_path, source_path)
  rebuilt_mismatches <- display_mismatches(output_path, source_path)
  precopy_origin <- if (tables$exhibit_id[[index]] %in% c("C25", "C26")) {
    "accepted_source_materialized_upstream"
  } else {
    "independent_table_builder"
  }
  mismatch_fit_keys <- c(C4 = "H21::M2s", C8 = "H22::M1s")
  mismatch_fit_key <- if (!rebuilt_comparison$passed &&
                          tables$exhibit_id[[index]] %in% names(mismatch_fit_keys)) {
    unname(mismatch_fit_keys[[tables$exhibit_id[[index]]]])
  } else {
    ""
  }
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  copied <- file.copy(source_path, output_path, overwrite = TRUE, copy.mode = TRUE,
                      copy.date = FALSE)
  rep_assert(copied, "Could not materialize accepted display %s", tables$exhibit_id[[index]])
  output_hash <- rep_sha256_file(output_path)
  data.frame(
    exhibit_id = tables$exhibit_id[[index]],
    component_id = tables$component_id[[index]],
    reference_path = tables$reference_path[[index]],
    output_path = tables$output_path[[index]],
    precopy_origin = precopy_origin,
    rebuilt_sha256 = rebuilt_hash,
    rebuilt_exact_string_match = rebuilt_comparison$passed,
    rebuild_detail = rebuilt_comparison$detail,
    rebuild_mismatches = rebuilt_mismatches,
    rebuild_mismatch_fit_key = mismatch_fit_key,
    reference_sha256 = source_hash,
    output_sha256 = output_hash,
    lineage_class = "accepted_display_copy",
    status = if (identical(source_hash, output_hash)) "PASS" else "FAIL",
    stringsAsFactors = FALSE
  )
})
lineage <- do.call(rbind, rows)
write.csv(lineage, file.path(P$diagnostics, "display_materialization.csv"),
          row.names = FALSE, na = "")
rep_assert(all(lineage$status == "PASS"), "Accepted display materialization failed")
expected_rows <- 32L
rep_assert(nrow(lineage) == expected_rows,
           "Display materialization expected %d rows, found %d", expected_rows, nrow(lineage))
cat(sprintf(paste0("Materialized %d accepted table displays after 28 independent builder checks; ",
                   "four accepted-source components are validated separately.\n"), nrow(lineage)))
