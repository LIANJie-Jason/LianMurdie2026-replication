# Private helpers for the draft replication table builders.

replication_bootstrap <- function() {
  env_root <- Sys.getenv("REPLICATION_ROOT", unset = "")
  if (nzchar(env_root)) {
    root <- normalizePath(env_root, mustWork = TRUE)
  } else {
    arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(arg) != 1L) {
      stop("Set REPLICATION_ROOT when sourcing a draft script interactively.", call. = FALSE)
    }
    script_path <- sub("^--file=", "", arg)
    script_path <- gsub("~+~", " ", script_path, fixed = TRUE)
    script <- normalizePath(script_path, mustWork = TRUE)
    root <- normalizePath(file.path(dirname(script), "..", ".."), mustWork = TRUE)
  }
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(root)
  sys.source(file.path(root, "code", "R", "00_setup.R"), envir = .GlobalEnv)
  setwd(old_wd)
  init_replication(root)
}

require_file <- function(path, label = basename(path)) {
  if (!file.exists(path)) stop(sprintf("Missing %s: %s", label, path), call. = FALSE)
  path
}

read_models <- function(P, stem) {
  invisible(vapply(c("survival", "logistf", "brglm2", "geepack"), requireNamespace,
                   logical(1), quietly = TRUE))
  path <- require_file(file.path(P$cache, paste0(stem, "_models.rds")), paste0(stem, " model cache"))
  out <- readRDS(path)
  if (!is.list(out) || is.null(names(out))) stop(sprintf("Invalid model bundle: %s", path), call. = FALSE)
  out
}

stars_for_p <- function(p) {
  if (is.na(p)) "" else if (p < .01) "***" else if (p < .05) "**" else if (p < .10) "*" else ""
}

fmt_coef <- function(est, p) sprintf("%.3f%s", est, stars_for_p(p))
fmt_se <- function(se) sprintf("(%.3f)", se)
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < .001, "<0.001", sprintf("%.3f", p)))
fmt_num <- function(x) ifelse(is.na(x), "", ifelse(abs(x) >= 1e4, sprintf("%.3e", x), sprintf("%.3f", x)))
fmt_ci <- function(lo, hi) ifelse(is.na(lo) | is.na(hi), "", sprintf("[%s, %s]", fmt_num(lo), fmt_num(hi)))

extract_term <- function(model, term) {
  coef_names <- names(stats::coef(model))
  if (!term %in% coef_names && grepl(":", term, fixed = TRUE)) {
    swapped <- paste(rev(strsplit(term, ":", fixed = TRUE)[[1]]), collapse = ":")
    if (swapped %in% coef_names) term <- swapped
  }
  if (inherits(model, "survreg")) {
    tab <- summary(model)$table
    if (!term %in% rownames(tab)) return(NULL)
    return(list(estimate = unname(tab[term, "Value"]), se = unname(tab[term, "Std. Error"]), p = unname(tab[term, "p"])))
  }
  if (inherits(model, "logistf")) {
    idx <- match(term, names(model$coefficients))
    if (is.na(idx)) return(NULL)
    return(list(estimate = unname(model$coefficients[idx]), se = unname(sqrt(diag(model$var))[idx]), p = unname(model$prob[idx])))
  }
  if (inherits(model, "coxph")) {
    tab <- summary(model)$coefficients
    if (!term %in% rownames(tab)) return(NULL)
    se_col <- if ("robust se" %in% colnames(tab)) "robust se" else "se(coef)"
    return(list(estimate = unname(tab[term, "coef"]), se = unname(tab[term, se_col]), p = unname(tab[term, "Pr(>|z|)"])))
  }
  if (inherits(model, "bracl")) {
    tab <- summary(model)$coefficients
    if (!term %in% rownames(tab)) return(NULL)
    return(list(estimate = unname(tab[term, "Estimate"]), se = unname(tab[term, "Std. Error"]), p = unname(tab[term, "Pr(>|z|)"])))
  }
  if (inherits(model, "geeglm")) {
    tab <- summary(model)$coefficients
    if (!term %in% rownames(tab)) return(NULL)
    return(list(estimate = unname(tab[term, "Estimate"]), se = unname(tab[term, "Std.err"]), p = unname(tab[term, "Pr(>|W|)"])))
  }
  stop(sprintf("Unsupported model class: %s", paste(class(model), collapse = ",")), call. = FALSE)
}

terms_for_label <- function(label) {
  exact <- list(
    "Proportion of English-Language Signs"=c("eng_c","eng_prob_general","eng_prop_year"),
    "Proportion of English-Language Signs (centered)"="eng_c",
    "Proportion of English-Language Signs Squared"=c("I(eng_c^2)","I(eng_prob_general^2)","I(eng_prop_year^2)"),
    "Proportion of English-Language Signs squared"="I(eng_c^2)",
    "Domestic Autonomy"="domaut_c", "State Authority Over Territory"="stterr_c",
    "Campaign Location FDI Inflows (% of GDP)"="fdiin_c", "International Trade Freedom"="trade_c",
    "Foreign Aid (ln)"="lnaid_c", "Number of IO Donors"="iodonor_c",
    "Campaign Location V-Dem CSO Repression (centered)"="repress_c",
    "Campaign Location V-Dem Civil Liberties (centered)"="civlib_c",
    "Campaign Location V-Dem Political Liberty (centered)"="clpol_c",
    "Campaign Location Polity Score (centered)"="polity_c",
    "V-Dem Political Liberty (centered)"="clpol_c",
    "V-Dem Freedom of Association thick (centered)"="frassoc_c",
    "Eng x CSO Repression"="eng_c:repress_c", "Eng x Civil Liberties"="eng_c:civlib_c",
    "Eng x Political Liberty"="eng_c:clpol_c", "Eng x Polity Score"="eng_c:polity_c",
    "Eng x Freedom of Association"="eng_c:frassoc_c",
    "Campaign Primarily Nonviolent"=c("nonviolent_camp1","nonviolent1"),
    "Campaign Purpose Regime Change"="REGCHANGE", "Campaign Location V-Dem CSO Repression"="v2csreprss",
    "CSO Repression (control)"="v2csreprss", "Polity Score (control)"="p_polity2",
    "Campaign Length, in Years (ln)"="lnlengthofcam", "Time in Campaign"="time_in_campaign",
    "Time in Campaign (years)"="time_in_campaign", "Campaign Length / Time in Campaign"=c("lnlengthofcam","time_in_campaign"),
    "Number of Images on Campaign (ln)"="lnnum_image_sum", "Number of Images in Campaign-Year (ln)"="lnnum_image_sum",
    "Number of Images on Campaign / in Campaign-Year (ln)"="lnnum_image_sum",
    "Campaign Location GDP (ln)"="lngdp", "Campaign Location GDP per capita (ln)"="lngdp", "Campaign-Year GDP per Capita (ln)"="lngdp_yr",
    "Campaign-Year GDP per capita (ln)"="lngdp_yr", "Campaign Location GDP per capita (ln, 1.3)"="lngdp",
    "Campaign-Year GDP per capita (ln, 2.1)"="lngdp_yr",
    "Campaign Location GDP / Campaign-Year GDP per Capita (ln)"=c("lngdp","lngdp_yr"),
    "Campaign Location Population (ln)"="lnpop", "Campaign-Year Population (ln)"="lnpop_yr",
    "English Colonial Legacy"="colonized_english1", "Constant"="(Intercept)"
  )
  if (!is.null(exact[[label]])) return(exact[[label]])
  x <- label
  x <- sub("^Proportion English-Language Signs Squared X ", "I(eng_c^2):", x)
  x <- sub("^Proportion English-Language Signs squared X ", "I(eng_c^2):", x)
  x <- sub("^Proportion English-Language Signs X ", "eng_c:", x)
  replacements <- c("Domestic Autonomy"="domaut_c","State Authority Over Territory"="stterr_c",
    "FDI Inflows (% of GDP)"="fdiin_c","International Trade Freedom"="trade_c","Foreign Aid (ln)"="lnaid_c",
    "Number of IO Donors"="iodonor_c")
  for (pat in names(replacements)) x <- sub(pat,replacements[[pat]],x,fixed=TRUE)
  if (grepl("^(eng_c|I\\(eng_c\\^2\\)):",x)) return(x)
  character()
}

hydrate_scaffold <- function(scaffold_path, models, stats = NULL, stats_ids = names(models)) {
  # The scaffold is a declarative, number-free layout specification. Its data
  # cells may contain only "@" (populate from the fresh model) or blank. The
  # accepted numerical table bodies are never read by a producer.
  body <- read.csv(require_file(scaffold_path,"number-free layout specification"),check.names=FALSE,stringsAsFactors=FALSE,na.strings=NULL)
  cols <- paste0("col",seq_along(models)); if(!all(cols%in%names(body)))stop("Layout/model column mismatch: ",basename(scaffold_path),call.=FALSE)
  bad_layout <- unique(unlist(body[cols], use.names = FALSE))
  bad_layout <- bad_layout[nzchar(bad_layout) & bad_layout != "@"]
  if (length(bad_layout)) stop("Numeric or unsupported content leaked into layout specification: ", basename(scaffold_path), call.=FALSE)
  for(i in seq_len(nrow(body))){label<-body$label[i];if(!nzchar(label))next
    if(label=="Observations"){for(j in seq_along(models)){z<-if(!is.null(stats)&&"N"%in%names(stats))stats[match(stats_ids[j],stats$model),"N"]else NA;body[i,cols[j]]<-if(length(z)&&!is.na(z))as.character(z)else as.character(model_n(models[[j]]))};next}
    if(startsWith(label,"EPV")){if(!is.null(stats)&&"epv"%in%names(stats)){for(j in seq_along(models)){if(!nzchar(body[[cols[j]]][i]))next;z<-stats[match(stats_ids[j],stats$model),"epv"];body[i,cols[j]]<-if(length(z)&&!is.na(z))sprintf("%.1f%s",z,ifelse(z<5,"*",""))else ""}};next}
    if(i==nrow(body)||nzchar(body$label[i+1]))next
    candidates<-terms_for_label(label); if(!length(candidates))stop("No term dictionary entry for label: ",label,call.=FALSE)
    for(j in seq_along(models)){if(!nzchar(body[[cols[j]]][i]))next
      hit<-NULL;for(term in candidates){hit<-extract_term(models[[j]],term);if(!is.null(hit))break}
      if(is.null(hit))stop(sprintf("No candidate term for '%s' in model %s",label,names(models)[j]),call.=FALSE)
      new_coef<-fmt_coef(hit$estimate,hit$p);new_se<-fmt_se(hit$se)
      body[[cols[j]]][i]<-new_coef;body[[cols[j]]][i+1L]<-new_se
    }
  }
  if (any(unlist(body[cols], use.names = FALSE) == "@")) {
    stop("Unpopulated layout marker remains in ", basename(scaffold_path), call. = FALSE)
  }
  body
}

model_n <- function(model) {
  if (inherits(model, "coxph") && !is.null(model$n)) return(as.integer(model$n))
  n <- tryCatch(stats::nobs(model), error = function(e) NA_integer_)
  if (length(n) != 1L || is.na(n)) {
    n <- if (!is.null(model$n)) model$n else if (!is.null(model$model)) nrow(model$model) else NA_integer_
  }
  as.integer(n)
}

model_terms_long <- function(models, exhibit_id) {
  rows <- lapply(seq_along(models), function(index) {
    model <- models[[index]]
    model_id <- names(models)[[index]]
    terms <- names(stats::coef(model))
    values <- lapply(terms, function(term) extract_term(model, term))
    keep <- !vapply(values, is.null, logical(1))
    values <- values[keep]
    terms <- terms[keep]
    data.frame(
      exhibit_id = exhibit_id,
      model_id = model_id,
      estimator = class(model)[[1L]],
      term = terms,
      estimate = vapply(values, `[[`, numeric(1), "estimate"),
      std_error = vapply(values, `[[`, numeric(1), "se"),
      p_value = vapply(values, `[[`, numeric(1), "p"),
      N = model_n(model),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

build_coefficient_body <- function(spec, models, add_epv = NULL) {
  model_ids <- names(models)
  if (!all(model_ids %in% names(spec))) stop("Row specification does not cover every model column.", call. = FALSE)
  rows <- list()
  for (i in seq_len(nrow(spec))) {
    coef_row <- spec$label[i]
    se_row <- ""
    for (id in model_ids) {
      term <- spec[[id]][i]
      if (is.na(term) || !nzchar(term)) {
        coef_row <- c(coef_row, "")
        se_row <- c(se_row, "")
      } else {
        z <- extract_term(models[[id]], term)
        if (is.null(z)) stop(sprintf("Term %s missing from model %s", term, id), call. = FALSE)
        coef_row <- c(coef_row, fmt_coef(z$estimate, z$p))
        se_row <- c(se_row, fmt_se(z$se))
      }
    }
    rows[[length(rows) + 1L]] <- coef_row
    rows[[length(rows) + 1L]] <- se_row
  }
  rows[[length(rows) + 1L]] <- rep("", length(model_ids) + 1L)
  rows[[length(rows) + 1L]] <- c("Observations", vapply(models, function(x) as.character(model_n(x)), character(1)))
  if (!is.null(add_epv)) rows[[length(rows) + 1L]] <- c(add_epv$label, add_epv$values)
  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(out) <- c("label", paste0("col", seq_along(model_ids)))
  out
}

write_body <- function(x, path, expected_rows) {
  if (nrow(x) != expected_rows) stop(sprintf("%s: expected %d rows, found %d", basename(path), expected_rows, nrow(x)), call. = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE, na = "")
  message(sprintf("Wrote %s (%d rows)", path, nrow(x)))
  invisible(path)
}

assert_exact_ids <- function(x, ids, column = "model") {
  got <- as.character(x[[column]])
  if (!identical(got, ids)) stop(sprintf("Unexpected %s order: %s", column, paste(got, collapse = ", ")), call. = FALSE)
  invisible(TRUE)
}
