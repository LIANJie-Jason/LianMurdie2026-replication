root <- Sys.getenv("REPLICATION_ROOT", unset = "")
if (!nzchar(root)) stop("REPLICATION_ROOT is not set; launch with run_all.R", call. = FALSE)
source(file.path(root, "code", "R", "00_setup.R"))
paths <- init_replication(root)
profile <- Sys.getenv("REPLICATION_PROFILE", unset = "full")

rep_assert(getRversion() >= "4.3.0", "R >= 4.3.0 is required; found %s", getRversion())
rep_assert(rep_command_exists("shasum") || rep_command_exists("sha256sum"),
           "A SHA-256 utility ('shasum' or 'sha256sum') is required")
rep_assert(rep_command_exists("pdftoppm"),
           "Poppler 'pdftoppm' is required for the advisory figure raster audit")

requirements_path <- file.path(paths$root, "environment", "required_packages.csv")
rep_assert_file(requirements_path, "Package requirements manifest")
requirements <- read.csv(requirements_path, stringsAsFactors = FALSE, check.names = FALSE)
rep_assert_columns(requirements, c("package", "minimum_version", "used_for"), "package requirements")
rep_require_packages(requirements$package)
too_old <- vapply(seq_len(nrow(requirements)), function(index) {
  packageVersion(requirements$package[[index]]) < package_version(requirements$minimum_version[[index]])
}, logical(1))
rep_assert(!any(too_old), "Packages below minimum version: %s", paste(requirements$package[too_old], collapse = ", "))

registry <- rep_script_registry(profile)
missing_scripts <- registry$script[!file.exists(file.path(paths$code, registry$script))]
rep_assert(!length(missing_scripts), "Pipeline scripts are missing: %s", paste(missing_scripts, collapse = ", "))

portable_r <- c(
  list.files(paths$code, pattern = "[.]R$", recursive = TRUE, full.names = TRUE),
  file.path(paths$root, "run_all.R")
)
portable_r <- unique(portable_r[file.exists(portable_r)])
absolute_hits <- unlist(lapply(portable_r, function(path) {
  lines <- readLines(path, warn = FALSE)
  unix_prefixes <- paste0("/", c("Users", "home", "opt"), "/")
  unix_absolute <- Reduce(`|`, lapply(unix_prefixes, function(prefix) {
    grepl(prefix, lines, fixed = TRUE)
  }))
  drive_forward <- grepl("[A-Za-z]:/", lines)
  drive_back <- grepl(paste0(":", strrep("\\", 2L)), lines, fixed = TRUE)
  unc_prefixes <- paste0(c("\"", "'"), strrep("\\", 4L))
  unc <- Reduce(`|`, lapply(unc_prefixes, function(prefix) grepl(prefix, lines, fixed = TRUE)))
  hit <- which(unix_absolute | drive_forward | drive_back | unc)
  if (length(hit)) sprintf("%s:%s", rep_relative_path(path, paths$root), paste(hit, collapse = ",")) else character()
}), use.names = FALSE)
rep_assert(!length(absolute_hits), "Machine-specific absolute paths found: %s", paste(absolute_hits, collapse = "; "))

probe <- tempfile("write-probe-", tmpdir = paths$audit)
rep_assert(file.create(probe), "Audit directory is not writable: %s", paths$audit)
unlink(probe)

session_path <- file.path(paths$audit, "sessionInfo.txt")
capture.output(sessionInfo(), file = session_path)
versions <- data.frame(
  component = c("R", requirements$package, "pdftoppm"),
  version = c(as.character(getRversion()),
              vapply(requirements$package, function(package) as.character(packageVersion(package)), character(1)),
              paste(system2("pdftoppm", "-v", stdout = TRUE, stderr = TRUE), collapse = " ")),
  stringsAsFactors = FALSE
)
write.csv(versions, file.path(paths$audit, "runtime_versions.csv"), row.names = FALSE, na = "")
cat(sprintf("Preflight passed: profile=%s; R=%s; packages=%d\n", profile, getRversion(), nrow(requirements)))
