root <- Sys.getenv("REPLICATION_ROOT", unset = "")
if (!nzchar(root)) {
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  if (!length(file_arg)) stop("Cannot locate run_all script; set REPLICATION_ROOT", call. = FALSE)
  runner_path <- sub("^--file=", "", file_arg[[1L]])
  runner_path <- gsub("~+~", " ", runner_path, fixed = TRUE)
  root <- normalizePath(file.path(dirname(runner_path), "..", ".."), mustWork = TRUE)
}
source(file.path(root, "code", "R", "00_setup.R"))
paths <- init_replication(root)
setwd(paths$root)

args <- commandArgs(trailingOnly = TRUE)
rep_assert(!length(args) ||
             (length(args) == 1L && args[[1L]] %in% c("full", "--profile=full")),
           "This package has one authoritative run; use --profile=full")
Sys.unsetenv(c("REPLICATION_SMOKE", "POWER_REPS", "BOOT_REPS"))

run_started <- Sys.time()
run_id <- sprintf("%s-pid%d",
                  format(run_started, "%Y%m%dT%H%M%OS6Z", tz = "UTC"),
                  Sys.getpid())
lock_dir <- file.path(paths$audit, ".pipeline.lock")
recovered_lock <- ""
if (dir.exists(lock_dir)) {
  owner_path <- file.path(lock_dir, "owner.txt")
  owner_lines <- if (file.exists(owner_path)) readLines(owner_path, warn = FALSE) else character()
  pid_line <- grep("^pid=", owner_lines, value = TRUE)
  owner_pid <- if (length(pid_line) == 1L) suppressWarnings(as.integer(sub("^pid=", "", pid_line))) else NA_integer_
  owner_alive <- if (is.na(owner_pid)) NA else
    tryCatch(tools::pskill(owner_pid, 0L), error = function(error) NA)
  if (identical(owner_alive, FALSE)) {
    stale_name <- sprintf("%s.stale.%s", lock_dir,
                          gsub("[^A-Za-z0-9_.-]", "_", run_id))
    if (file.rename(lock_dir, stale_name)) {
      recovered_lock <- stale_name
      message(sprintf("Recovered stale replication lock owned by inactive PID %d", owner_pid))
    }
  }
}
if (!dir.create(lock_dir, showWarnings = FALSE)) {
  owner_path <- file.path(lock_dir, "owner.txt")
  owner <- if (file.exists(owner_path)) paste(readLines(owner_path, warn = FALSE), collapse = "; ") else "unknown"
  rep_abort("Another replication run holds %s (owner: %s)", lock_dir, owner)
}
if (nzchar(recovered_lock) && dir.exists(recovered_lock)) {
  unlink(recovered_lock, recursive = TRUE, force = TRUE)
}
writeLines(c(sprintf("run_id=%s", run_id), sprintf("pid=%d", Sys.getpid()),
             sprintf("started_utc=%s", format(run_started, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))),
           file.path(lock_dir, "owner.txt"), useBytes = TRUE)
lock_guard <- new.env(parent = emptyenv())
reg.finalizer(lock_guard, function(environment) unlink(lock_dir, recursive = TRUE, force = TRUE),
              onexit = TRUE)
logs_root <- file.path(paths$audit, "logs")
prior_logs <- if (dir.exists(logs_root)) {
  list.dirs(logs_root, recursive = FALSE, full.names = TRUE)
} else character()
if (length(prior_logs)) unlink(prior_logs, recursive = TRUE, force = TRUE)
dir.create(logs_root, recursive = TRUE, showWarnings = FALSE)
log_dir <- file.path(logs_root, run_id)
rep_assert(!dir.exists(log_dir), "Run directory already exists: %s", log_dir)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# A run owns a clean, declared output namespace.  Remove only known generated
# artifacts; unknown files are preserved so the final allowlist gate can flag
# them rather than silently deleting user material.
known_output_paths <- file.path(paths$root, rep_all_known_outputs(paths))
existing_known <- known_output_paths[file.exists(known_output_paths)]
if (length(existing_known)) {
  unlink(existing_known, recursive = FALSE, force = TRUE)
  rep_assert(!any(file.exists(existing_known)), "Could not clear prior generated outputs")
}
current_audits <- file.path(paths$audit, c(
  "run_context.csv", "pipeline_status.csv", "runtime_versions.csv", "sessionInfo.txt",
  "input_validation.csv", "rep_coverage.csv", "rep_negative.csv",
  "figure_visual_comparison.csv", "reference_validation.csv", "validation_summary.csv",
  "output_completeness.csv", "output_manifest.csv", "code_manifest.csv"
))
unlink(current_audits[file.exists(current_audits)], recursive = FALSE, force = TRUE)

marker_path <- file.path(log_dir, "RUN_MARKER")
writeLines(run_id, marker_path, useBytes = TRUE)
marker_epoch <- as.numeric(file.info(marker_path)$mtime)
rep_assert(is.finite(marker_epoch), "Could not create the run marker")
Sys.setenv(
  REPLICATION_ROOT = paths$root,
  REPLICATION_RUN_ID = run_id,
  REPLICATION_RUN_STARTED_EPOCH = sprintf("%.9f", marker_epoch),
  REPLICATION_RUN_MARKER = marker_path
)

context <- data.frame(
  run_id = run_id,
  profile = "full",
  started_utc = format(run_started, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  r_version = as.character(getRversion()),
  platform = R.version$platform,
  stringsAsFactors = FALSE
)
write.csv(context, file.path(paths$audit, "run_context.csv"), row.names = FALSE, na = "")

registry <- rep_script_registry()
execution_status <- tryCatch({
  for (script in registry) {
    message(sprintf("[%s] running %s", format(Sys.time(), "%H:%M:%S"), script))
    rep_run_script(script, paths, log_dir)
  }
  "PASS"
}, error = function(error) {
  writeLines(conditionMessage(error), file.path(log_dir, "FAILED.txt"))
  message(conditionMessage(error))
  "FAIL"
})

finished <- Sys.time()
pipeline_status <- transform(context,
  finished_utc = format(finished, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  elapsed_seconds = as.numeric(difftime(finished, run_started, units = "secs")),
  status = execution_status
)
write.csv(pipeline_status, file.path(paths$audit, "pipeline_status.csv"), row.names = FALSE, na = "")
if (!identical(execution_status, "PASS")) quit(save = "no", status = 1L)
message("Replication pipeline completed: profile=full, status=PASS")
