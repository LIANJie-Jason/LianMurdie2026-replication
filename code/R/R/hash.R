rep_command_exists <- function(command) nzchar(Sys.which(command))

rep_sha256_file <- function(path) {
  rep_assert_file(path)
  normalized <- normalizePath(path, mustWork = TRUE)

  if (rep_command_exists("shasum")) {
    output <- system2("shasum", c("-a", "256", shQuote(normalized)), stdout = TRUE, stderr = TRUE)
    status <- attr(output, "status") %||% 0L
    rep_assert(identical(as.integer(status), 0L), "shasum failed for %s: %s", path, paste(output, collapse = " "))
    hash <- strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
  } else if (rep_command_exists("sha256sum")) {
    output <- system2("sha256sum", shQuote(normalized), stdout = TRUE, stderr = TRUE)
    status <- attr(output, "status") %||% 0L
    rep_assert(identical(as.integer(status), 0L), "sha256sum failed for %s: %s", path, paste(output, collapse = " "))
    hash <- strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
  } else {
    rep_abort("No SHA-256 utility found. Install either 'shasum' or 'sha256sum'.")
  }

  hash <- tolower(hash)
  rep_assert(grepl("^[0-9a-f]{64}$", hash), "Invalid SHA-256 output for %s: %s", path, hash)
  unname(hash)
}

rep_file_size <- function(path) {
  rep_assert_file(path)
  unname(file.info(path)$size)
}
