root <- Sys.getenv("REPLICATION_ROOT", unset = "")
if (!nzchar(root)) {
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
  file_arg <- gsub("~+~", " ", file_arg, fixed = TRUE)
  root <- normalizePath(file.path(dirname(file_arg), "..", ".."), mustWork = TRUE)
}
source(file.path(root, "code", "R", "00_setup.R"))
P <- init_replication(root)
rep_require_packages(c("survival", "logistf", "geepack"))

bundle_files <- c(
  H1 = "H1_models.rds", H21 = "H21_models.rds", H22 = "H22_models.rds",
  H3 = "H3_models.rds", H21_extensions = "H21_extensions_models.rds",
  H3_additional = "H3_additional_moderators_models.rds"
)
bundle_paths <- file.path(P$cache, unname(bundle_files))
rep_assert(all(file.exists(bundle_paths)), "Fit-status audit is missing model bundle(s): %s",
           paste(names(bundle_files)[!file.exists(bundle_paths)], collapse = ", "))
bundles <- setNames(lapply(bundle_paths, readRDS), names(bundle_files))

finite_model_output <- function(model) {
  beta <- tryCatch(coef(model), error = function(error) numeric())
  variance <- tryCatch(vcov(model), error = function(error) NULL)
  length(beta) > 0L && all(is.finite(beta)) && !is.null(variance) &&
    all(is.finite(as.matrix(variance)))
}

audit_one <- function(family, model_id, model) {
  model_class <- class(model)[[1L]]
  finite_output <- finite_model_output(model)
  iterations <- NA_integer_
  iteration_limit <- FALSE
  geese_error <- NA_integer_
  affected_terms <- character()
  reasons <- character()
  converged <- finite_output

  if (inherits(model, "coxph")) {
    iterations <- as.integer(model$iter %||% NA_integer_)
    max_iterations <- survival::coxph.control()$iter.max
    iteration_limit <- is.finite(iterations) && iterations >= max_iterations
    beta <- coef(model)
    monotone_regchange <- "REGCHANGE" %in% names(beta) &&
      is.finite(beta[["REGCHANGE"]]) && abs(beta[["REGCHANGE"]]) > 10
    if (iteration_limit) reasons <- c(reasons, "cox_iteration_limit")
    if (monotone_regchange) {
      reasons <- c(reasons, "cox_monotone_likelihood")
      affected_terms <- c(affected_terms, "REGCHANGE")
    }
    converged <- finite_output && !iteration_limit && !monotone_regchange
  } else if (inherits(model, "survreg")) {
    iterations <- as.integer(model$iter %||% NA_integer_)
    max_iterations <- survival::survreg.control()$iter.max
    iteration_limit <- is.finite(iterations) && iterations >= max_iterations
    fit_failure <- model$fail %||% ""
    table <- tryCatch(summary(model)$table, error = function(error) NULL)
    monotone_regchange <- !is.null(table) && "REGCHANGE" %in% rownames(table) &&
      all(c("Value", "Std. Error") %in% colnames(table)) &&
      is.finite(table["REGCHANGE", "Value"]) &&
      is.finite(table["REGCHANGE", "Std. Error"]) &&
      abs(table["REGCHANGE", "Value"]) > 8 &&
      table["REGCHANGE", "Std. Error"] > 100
    if (iteration_limit) reasons <- c(reasons, "survreg_iteration_limit")
    if (monotone_regchange) {
      reasons <- c(reasons, "survreg_monotone_likelihood")
      affected_terms <- c(affected_terms, "REGCHANGE", "(Intercept)")
    }
    if (nzchar(fit_failure)) reasons <- c(reasons, paste0("survreg_failure:", fit_failure))
    converged <- finite_output && !iteration_limit && !monotone_regchange &&
      !nzchar(fit_failure)
  } else if (inherits(model, "geeglm")) {
    geese_error <- as.integer(model$geese$error %||% NA_integer_)
    if (!identical(geese_error, 0L)) reasons <- c(reasons, "gee_error")
    converged <- finite_output && identical(geese_error, 0L)
  } else if (inherits(model, "logistf")) {
    conv <- abs(model$conv[c("LL change", "max abs score", "beta change")])
    thresholds <- c(model$control$lconv, model$control$gconv, model$control$xconv)
    converged <- finite_output && all(is.finite(conv)) && all(conv <= thresholds)
    if (!converged) reasons <- c(reasons, "logistf_convergence_threshold")
  } else {
    return(NULL)
  }

  if (!finite_output) reasons <- c(reasons, "nonfinite_model_output")
  if (!length(reasons)) reasons <- "none"
  data.frame(
    family = family, model = model_id, estimator_class = model_class,
    iterations = iterations, iteration_limit = iteration_limit,
    geese_error = geese_error,
    affected_terms = paste(unique(affected_terms), collapse = ";"),
    finite_output = finite_output, converged = converged,
    accepted_values_preserved = TRUE,
    reporting_class = if (converged) "accepted_replicated" else
      "accepted_archival_nonconverged",
    reason = paste(unique(reasons), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

rows <- list()
for (family in names(bundles)) {
  bundle <- bundles[[family]]
  for (model_id in names(bundle)) {
    row <- audit_one(family, model_id, bundle[[model_id]])
    if (!is.null(row)) rows[[length(rows) + 1L]] <- row
  }
}
status <- do.call(rbind, rows)
rownames(status) <- NULL
write.csv(status, file.path(P$diagnostics, "accepted_fit_status_all.csv"),
          row.names = FALSE, na = "")
flagged <- status[!status$converged, c("family", "model", "reason", "affected_terms"),
                  drop = FALSE]
cat(sprintf("Fit-status audit: %d supported fits; %d accepted archival/nonconverged.\n",
            nrow(status), nrow(flagged)))
if (nrow(flagged)) print(flagged, row.names = FALSE)
