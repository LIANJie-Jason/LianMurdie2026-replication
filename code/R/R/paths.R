rep_find_root <- function(start = getwd()) {
  env_root <- Sys.getenv("REPLICATION_ROOT", unset = "")
  candidates <- c(env_root, start)

  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    script_path <- sub("^--file=", "", file_arg[[1L]])
    candidates <- c(candidates, dirname(normalizePath(script_path, mustWork = FALSE)))
  }

  candidates <- unique(candidates[nzchar(candidates)])
  for (candidate in candidates) {
    current <- normalizePath(candidate, mustWork = FALSE)
    repeat {
      marker <- file.path(current, "data", "input_manifest.csv")
      setup <- file.path(current, "code", "R", "00_setup.R")
      if (file.exists(marker) && file.exists(setup)) return(current)
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }
  stop("Cannot locate replication package root. Set REPLICATION_ROOT or run from inside the package.", call. = FALSE)
}

rep_relative_path <- function(path, root) {
  path <- normalizePath(path, mustWork = FALSE)
  root <- normalizePath(root, mustWork = TRUE)
  prefix <- paste0(root, .Platform$file.sep)
  if (identical(path, root)) return(".")
  if (!startsWith(path, prefix)) {
    stop(sprintf("Path is outside replication root: %s", path), call. = FALSE)
  }
  substring(path, nchar(prefix) + 1L)
}

rep_ensure_dirs <- function(paths) {
  failures <- character()
  for (path in unname(unlist(paths, use.names = FALSE))) {
    if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
      failures <- c(failures, path)
    }
  }
  if (length(failures)) {
    stop(sprintf("Could not create directories: %s", paste(failures, collapse = ", ")), call. = FALSE)
  }
  invisible(paths)
}
