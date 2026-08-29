root <- Sys.getenv("REPLICATION_ROOT", unset = "")
if (!nzchar(root)) stop("REPLICATION_ROOT is not set; launch with run_all.R", call. = FALSE)
source(file.path(root, "code", "R", "00_setup.R"))
paths <- init_replication(root)
profile <- Sys.getenv("REPLICATION_PROFILE", unset = "full")
run_started <- as.numeric(Sys.getenv("REPLICATION_RUN_STARTED_EPOCH", unset = "0"))

# Finder/Dropbox may create metadata while the pipeline is running. These files
# are never publication artifacts; remove them immediately before the exact
# package/output inventory so the delivered directory stays clean.
metadata <- list.files(paths$root, pattern = "^[.]DS_Store$", recursive = TRUE,
                       full.names = TRUE, all.files = TRUE, no.. = TRUE)
if (length(metadata)) unlink(metadata, recursive = FALSE, force = TRUE)
rep_assert(!any(file.exists(metadata)), "Could not remove package metadata: %s",
           paste(metadata[file.exists(metadata)], collapse = ", "))
run_marker <- Sys.getenv("REPLICATION_RUN_MARKER", unset = "")
rep_assert(length(run_started) == 1L && is.finite(run_started) && run_started > 0,
           "Stage 41 requires a current run marker; launch with run_all.R")
rep_assert(nzchar(run_marker) && file.exists(run_marker),
           "Stage 41 cannot verify the current run marker")
rep_assert(abs(as.numeric(file.info(run_marker)$mtime) - run_started) <= 1e-6,
           "Stage 41 run-marker timestamp does not match the execution context")

# Model RDS files are required only while downstream scripts run. Remove these
# runtime intermediates after validation so the publication folder contains
# code, frozen inputs, numerical/display outputs, and audit evidence only.
runtime_cache <- file.path(paths$root, rep_runtime_cache_registry())
unlink(runtime_cache[file.exists(runtime_cache)], recursive = FALSE, force = TRUE)
rep_assert(!any(file.exists(runtime_cache)), "Could not remove runtime model caches")
if (dir.exists(paths$cache) && !length(list.files(paths$cache, all.files = TRUE, no.. = TRUE))) {
  unlink(paths$cache, recursive = TRUE, force = TRUE)
}

code_files <- c(
  file.path(paths$root, "run_all.R"),
  list.files(paths$code, pattern = "[.]R$", recursive = TRUE, full.names = TRUE),
  list.files(file.path(paths$root, "provenance"), pattern = "[.]R$",
             recursive = TRUE, full.names = TRUE)
)
code_files <- sort(unique(code_files[file.exists(code_files)]))
code_manifest <- data.frame(
  path = vapply(code_files, rep_relative_path, character(1), root = paths$root),
  bytes = vapply(code_files, rep_file_size, numeric(1)),
  sha256 = vapply(code_files, rep_sha256_file, character(1)),
  run_id = Sys.getenv("REPLICATION_RUN_ID", unset = "manual"),
  stringsAsFactors = FALSE
)
write.csv(code_manifest, file.path(paths$audit, "code_manifest.csv"),
          row.names = FALSE, na = "")

expected <- rep_active_output_registry(paths, profile)
expected_paths <- file.path(paths$root, expected)
exists <- file.exists(expected_paths)
fresh <- exists & as.numeric(file.info(expected_paths)$mtime) >= run_started
completeness <- data.frame(path = expected, exists = exists, produced_this_run = fresh, stringsAsFactors = FALSE)
write.csv(completeness, file.path(paths$audit, "output_completeness.csv"), row.names = FALSE, na = "")
rep_assert(all(exists), "Expected output files are missing: %s", paste(expected[!exists], collapse = ", "))
rep_assert(all(fresh), "Expected outputs were not produced in this run: %s", paste(expected[!fresh], collapse = ", "))

output_root <- file.path(paths$root, "output")
files <- list.files(output_root, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
files <- files[!dir.exists(files)]
rep_assert(length(files) > 0, "No replication outputs exist under %s", output_root)
files <- sort(files)
relative <- vapply(files, rep_relative_path, character(1), root = paths$root)
unexpected <- setdiff(relative, expected)
missing_from_tree <- setdiff(expected, relative)
rep_assert(!length(unexpected), "Unexpected output files are outside the allowlist: %s",
           paste(unexpected, collapse = ", "))
rep_assert(!length(missing_from_tree), "Allowlisted outputs are missing: %s",
           paste(missing_from_tree, collapse = ", "))
artifact_class <- sub("^output/([^/]+)(?:/.*)?$", "\\1", relative)
manifest <- data.frame(
  path = relative,
  artifact_class = artifact_class,
  bytes = vapply(files, rep_file_size, numeric(1)),
  sha256 = vapply(files, rep_sha256_file, character(1)),
  modified_utc = vapply(file.info(files)$mtime, format, character(1), format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  produced_this_run = as.numeric(file.info(files)$mtime) >= run_started,
  run_id = Sys.getenv("REPLICATION_RUN_ID", unset = "manual"),
  profile = profile,
  stringsAsFactors = FALSE
)
rep_assert(all(manifest$produced_this_run), "Output manifest contains stale artifacts")
write.csv(manifest, file.path(paths$audit, "output_manifest.csv"), row.names = FALSE, na = "")
cat(sprintf("Output manifest written: %d files\n", nrow(manifest)))
