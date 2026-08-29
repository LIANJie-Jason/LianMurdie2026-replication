root <- Sys.getenv("REPLICATION_ROOT", unset = "")
if (!nzchar(root)) stop("REPLICATION_ROOT is not set; launch with run_all.R", call. = FALSE)
source(file.path(root, "code", "R", "00_setup.R"))
paths <- init_replication(root)

manifest_path <- file.path(paths$data, "input_manifest.csv")
manifest <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
required_manifest_columns <- c("file", "rows", "columns", "analytic_unit", "key_columns", "required_columns", "sha256", "provenance")
rep_assert_columns(manifest, required_manifest_columns, "input manifest")
rep_assert(!anyDuplicated(manifest$file), "Input manifest has duplicate file values")

checks <- lapply(seq_len(nrow(manifest)), function(index) {
  entry <- manifest[index, , drop = FALSE]
  path <- file.path(paths$data, entry$file)
  rep_assert_file(path, "Merged analysis input")
  actual_hash <- rep_sha256_file(path)
  rep_assert(identical(actual_hash, tolower(entry$sha256)),
             "Input hash mismatch for %s (actual %s; expected %s)", entry$file, actual_hash, entry$sha256)
  data <- rep_read_csv(path)
  rep_assert(nrow(data) == entry$rows, "%s row count mismatch: %d != %d", entry$file, nrow(data), entry$rows)
  rep_assert(ncol(data) == entry$columns, "%s column count mismatch: %d != %d", entry$file, ncol(data), entry$columns)
  required <- strsplit(entry$required_columns, ";", fixed = TRUE)[[1L]]
  keys <- strsplit(entry$key_columns, ";", fixed = TRUE)[[1L]]
  rep_assert_columns(data, required, entry$file)
  rep_assert_unique(data, keys, entry$file)
  data.frame(file = entry$file, rows = nrow(data), columns = ncol(data), sha256 = actual_hash,
             shape_ok = TRUE, schema_ok = TRUE, key_unique = TRUE, stringsAsFactors = FALSE)
})

checks <- do.call(rbind, checks)
write.csv(checks, file.path(paths$audit, "input_validation.csv"), row.names = FALSE, na = "")
cat(sprintf("Input validation passed: %d frozen merged files\n", nrow(checks)))
