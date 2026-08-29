rep_run_script <- function(script, paths, log_dir, extra_args = character()) {
  script_path <- file.path(paths$code, script)
  rep_assert_file(script_path, "Pipeline script")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  log_path <- file.path(log_dir, sub("[.]R$", ".log", basename(script)))
  rscript <- file.path(R.home("bin"), "Rscript")
  env <- sprintf("REPLICATION_ROOT=%s", shQuote(paths$root))
  args <- c("--vanilla", shQuote(script_path), extra_args)
  started <- Sys.time()
  output <- system2(rscript, args, stdout = TRUE, stderr = TRUE, env = env)
  status <- as.integer(attr(output, "status") %||% 0L)
  writeLines(c(sprintf("script: %s", script), sprintf("started: %s", format(started, tz = "UTC")),
               sprintf("finished: %s", format(Sys.time(), tz = "UTC")), sprintf("exit_status: %d", status), "", output),
             log_path, useBytes = TRUE)
  if (status != 0L) rep_abort("Pipeline stage failed (%s). See %s", script, log_path)
  invisible(log_path)
}
