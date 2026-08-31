rep_abort <- function(..., call. = FALSE) {
  stop(sprintf(...), call. = call.)
}

rep_assert <- function(condition, ..., call. = FALSE) {
  if (!isTRUE(condition)) rep_abort(..., call. = call.)
  invisible(TRUE)
}

rep_assert_file <- function(path, label = "Required file") {
  rep_assert(file.exists(path), "%s is missing: %s", label, path)
  rep_assert(!dir.exists(path), "%s is a directory, not a file: %s", label, path)
  invisible(path)
}

rep_assert_columns <- function(data, required, label = "data") {
  missing <- setdiff(required, names(data))
  rep_assert(!length(missing), "%s is missing columns: %s", label, paste(missing, collapse = ", "))
  invisible(data)
}

rep_assert_unique <- function(data, columns, label = "data") {
  rep_assert_columns(data, columns, label)
  key <- do.call(paste, c(data[columns], sep = "\r"))
  duplicated_rows <- which(duplicated(key) | duplicated(key, fromLast = TRUE))
  rep_assert(!length(duplicated_rows), "%s has duplicate key rows for [%s]: %s",
             label, paste(columns, collapse = ", "), paste(head(duplicated_rows, 20L), collapse = ", "))
  invisible(data)
}

assert_no_nonfinite <- function(x, label) {
  bad <- sum(!is.finite(x[!is.na(x)]))
  if (bad > 0) {
    stop(sprintf("%s contains %d non-finite non-missing values.", label, bad), call. = FALSE)
  }
}

rep_require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  rep_assert(!length(missing), "Missing required R packages: %s", paste(missing, collapse = ", "))
  invisible(packages)
}
