# Package-local fallback: unlike base R's operator, this intentionally treats
# NULL, length-zero, and a leading empty string as absent path/config values.
`%||%` <- function(x, y) if (is.null(x) || !length(x) || !nzchar(x[[1L]])) y else x

# Shared package initialization. This file must remain base-R only.
.rep_env_root <- Sys.getenv("REPLICATION_ROOT", unset = "")
.rep_frame_file <- sys.frame(1)$ofile %||% ""
.rep_candidates <- c(
  if (nzchar(.rep_env_root)) file.path(.rep_env_root, "code", "R", "00_setup.R") else character(),
  .rep_frame_file,
  file.path(getwd(), "code", "R", "00_setup.R")
)
.rep_setup_file <- .rep_candidates[file.exists(.rep_candidates)][1L]
if (is.na(.rep_setup_file) || !nzchar(.rep_setup_file)) {
  stop("Cannot locate code/R/00_setup.R; set REPLICATION_ROOT before sourcing", call. = FALSE)
}
.rep_setup_file <- normalizePath(.rep_setup_file, mustWork = TRUE)
.rep_setup_dir <- dirname(.rep_setup_file)

for (.rep_helper in c("paths.R", "assertions.R", "hash.R", "comparison.R", "registry.R", "process.R")) {
  .rep_helper_path <- file.path(.rep_setup_dir, "R", .rep_helper)
  if (!file.exists(.rep_helper_path)) stop(sprintf("Missing shared helper: %s", .rep_helper_path), call. = FALSE)
  source(.rep_helper_path, local = FALSE, chdir = FALSE)
}
rm(.rep_helper, .rep_helper_path)

init_replication <- function(root = NULL) {
  root <- normalizePath(root %||% rep_find_root(), mustWork = TRUE)
  paths <- list(
    root = root,
    data = file.path(root, "data"),
    code = file.path(root, "code", "R"),
    estimates = file.path(root, "output", "estimates"),
    cache = file.path(root, "output", "cache"),
    main_tables = file.path(root, "output", "tables", "main"),
    appendix_tables = file.path(root, "output", "tables", "appendix"),
    figures = file.path(root, "output", "figures"),
    diagnostics = file.path(root, "output", "diagnostics"),
    audit = file.path(root, "audit"),
    reference = file.path(root, "reference", "accepted")
  )
  rep_assert_file(file.path(paths$data, "input_manifest.csv"), "Input manifest")
  rep_assert(dir.exists(paths$reference), "Accepted-reference directory is missing: %s", paths$reference)
  rep_ensure_dirs(paths[c("estimates", "cache", "main_tables", "appendix_tables", "figures", "diagnostics", "audit")])
  options(replication.root = root, stringsAsFactors = FALSE)
  paths
}
