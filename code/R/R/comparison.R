rep_read_csv <- function(path) {
  rep_assert_file(path)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("NA", ""))
}

rep_normalize_scalar <- function(x) {
  ifelse(is.na(x), "<NA>", enc2utf8(as.character(x)))
}

rep_normalize_display_scalar <- function(x) {
  value <- rep_normalize_scalar(x)
  gsub("\\r\\n?", "\\n", value)
}

rep_compare_files_exact <- function(actual, expected) {
  actual_hash <- rep_sha256_file(actual)
  expected_hash <- rep_sha256_file(expected)
  list(
    passed = identical(actual_hash, expected_hash),
    method = "sha256_exact",
    max_abs_diff = NA_real_,
    max_rel_diff = NA_real_,
    detail = if (identical(actual_hash, expected_hash)) "byte-identical" else
      sprintf("SHA-256 differs (actual %s; expected %s)", actual_hash, expected_hash)
  )
}

rep_compare_csv_numeric <- function(actual, expected, abs_tol = 1e-8, rel_tol = 1e-6) {
  lhs <- rep_read_csv(actual)
  rhs <- rep_read_csv(expected)
  if (!identical(dim(lhs), dim(rhs))) {
    return(list(passed = FALSE, method = "csv_numeric", max_abs_diff = Inf,
                max_rel_diff = Inf, detail = sprintf("shape differs: actual %dx%d; expected %dx%d",
                                                     nrow(lhs), ncol(lhs), nrow(rhs), ncol(rhs))))
  }
  if (!identical(names(lhs), names(rhs))) {
    return(list(passed = FALSE, method = "csv_numeric", max_abs_diff = Inf,
                max_rel_diff = Inf, detail = "column names or order differ"))
  }

  max_abs <- 0
  max_rel <- 0
  failures <- character()
  for (column in names(lhs)) {
    lhs_numeric <- suppressWarnings(as.numeric(lhs[[column]]))
    rhs_numeric <- suppressWarnings(as.numeric(rhs[[column]]))
    numeric_like <- all(is.na(lhs[[column]]) | !is.na(lhs_numeric)) &&
      all(is.na(rhs[[column]]) | !is.na(rhs_numeric))

    if (!numeric_like) {
      mismatch <- which(rep_normalize_scalar(lhs[[column]]) != rep_normalize_scalar(rhs[[column]]))
      if (length(mismatch)) failures <- c(failures, sprintf("%s nonnumeric mismatch at rows %s",
                                                             column, paste(head(mismatch, 5L), collapse = ",")))
      next
    }

    same_na <- identical(is.na(lhs_numeric), is.na(rhs_numeric))
    if (!same_na) {
      failures <- c(failures, sprintf("%s NA pattern differs", column))
      next
    }
    usable <- !is.na(lhs_numeric)
    if (!any(usable)) next
    same_infinite <- is.infinite(lhs_numeric) & is.infinite(rhs_numeric) &
      sign(lhs_numeric) == sign(rhs_numeric)
    bad_nonfinite <- usable & (!is.finite(lhs_numeric) | !is.finite(rhs_numeric)) &
      !same_infinite
    if (any(bad_nonfinite)) {
      failures <- c(failures, sprintf("%s non-finite values differ at %d row(s)",
                                      column, sum(bad_nonfinite)))
    }
    finite <- usable & is.finite(lhs_numeric) & is.finite(rhs_numeric)
    if (!any(finite)) next
    absolute <- abs(lhs_numeric[finite] - rhs_numeric[finite])
    scale <- pmax(abs(rhs_numeric[finite]), .Machine$double.eps)
    relative <- absolute / scale
    max_abs <- max(max_abs, absolute)
    max_rel <- max(max_rel, relative)
    bad <- which(absolute > abs_tol + rel_tol * abs(rhs_numeric[finite]))
    if (length(bad)) failures <- c(failures, sprintf("%s exceeds tolerance at %d value(s)", column, length(bad)))
  }

  list(
    passed = !length(failures),
    method = sprintf("csv_numeric(abs=%g,rel=%g)", abs_tol, rel_tol),
    max_abs_diff = max_abs,
    max_rel_diff = max_rel,
    detail = if (!length(failures)) "schema, labels, missingness, and numeric values match" else
      paste(head(failures, 20L), collapse = "; ")
  )
}

rep_compare_csv_exact_strings <- function(actual, expected) {
  rep_assert_file(actual)
  rep_assert_file(expected)
  read_display <- function(path) {
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
             colClasses = "character", na.strings = NULL,
             strip.white = FALSE, blank.lines.skip = FALSE)
  }
  lhs <- read_display(actual)
  rhs <- read_display(expected)
  if (!identical(dim(lhs), dim(rhs))) {
    return(list(passed = FALSE, method = "csv_exact_strings", max_abs_diff = NA_real_,
                max_rel_diff = NA_real_, detail = sprintf("shape differs: actual %dx%d; expected %dx%d",
                                                          nrow(lhs), ncol(lhs), nrow(rhs), ncol(rhs))))
  }
  if (!identical(names(lhs), names(rhs))) {
    return(list(passed = FALSE, method = "csv_exact_strings", max_abs_diff = NA_real_,
                max_rel_diff = NA_real_, detail = "column names or order differ"))
  }
  failures <- character()
  for (column in names(lhs)) {
    mismatch <- which(rep_normalize_display_scalar(lhs[[column]]) !=
                        rep_normalize_display_scalar(rhs[[column]]))
    if (length(mismatch)) {
      failures <- c(failures, sprintf("%s differs at rows %s", column,
                                      paste(head(mismatch, 10L), collapse = ",")))
    }
  }
  list(passed = !length(failures), method = "csv_exact_strings",
       max_abs_diff = NA_real_, max_rel_diff = NA_real_,
       detail = if (!length(failures)) "all displayed strings match" else
         paste(head(failures, 20L), collapse = "; "))
}

rep_compare_display_copy <- function(actual, expected) {
  result <- rep_compare_files_exact(actual, expected)
  result$method <- "accepted_display_copy_identity"
  result$detail <- if (result$passed) {
    "byte-identical accepted display copy; pre-copy builder evidence is in display_materialization.csv"
  } else {
    result$detail
  }
  result
}

rep_compare_output <- function(actual, expected, method = "numeric", abs_tol = 1e-8, rel_tol = 1e-6) {
  if (!file.exists(actual)) return(list(passed = FALSE, method = method, max_abs_diff = Inf,
                                        max_rel_diff = Inf, detail = "actual output is missing"))
  if (!file.exists(expected)) return(list(passed = FALSE, method = method, max_abs_diff = Inf,
                                          max_rel_diff = Inf, detail = "accepted reference is missing"))
  switch(method,
         exact = rep_compare_files_exact(actual, expected),
         csv_exact = rep_compare_csv_exact_strings(actual, expected),
         advisory_csv = rep_compare_csv_exact_strings(actual, expected),
         display_copy = rep_compare_display_copy(actual, expected),
         advisory_numeric = rep_compare_csv_numeric(actual, expected, abs_tol, rel_tol),
         numeric = rep_compare_csv_numeric(actual, expected, abs_tol, rel_tol),
         presence = list(passed = rep_file_size(actual) > 0, method = "presence",
                         max_abs_diff = NA_real_, max_rel_diff = NA_real_,
                         detail = sprintf("output present (%d bytes)", rep_file_size(actual))),
         rep_abort("Unknown comparison method: %s", method))
}
