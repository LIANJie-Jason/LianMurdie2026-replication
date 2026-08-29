###############################################################################
# H3_table4_additional_moderators.R
#
# Standalone extension to Table 4 / H3. Audits the Political Liberty
# substitution and adds one V-Dem domestic-openness companion moderator:
#   v2x_clpol          : Political civil liberties index
#   v2x_frassoc_thick  : Freedom of association thick index
#
# The model specifications mirror the Table 4 main-panel campaign-level models:
#   - Weibull AFT: time to broad campaign success/concession
#   - Firth PML: broad campaign success/concession binary
#   - same controls, mean-centering, and sign convention as H3_code.r
###############################################################################

rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(logistf)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
rep_root <- Sys.getenv("REPLICATION_ROOT", unset = "")
if (!nzchar(rep_root)) {
  script_file <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
  rep_root <- dirname(dirname(dirname(script_file)))
}
source(file.path(rep_root, "code", "R", "00_setup.R"), local = FALSE)
P <- init_replication(rep_root)

navco13_path <- file.path(P$data, "df_final.csv")
moderator_path <- file.path(P$data, "h3_additional_moderators.csv")
out_path     <- file.path(P$estimates, "H3_table4_additional_moderators_results.csv")
h3_results_path <- file.path(P$estimates, "H3_results.csv")
cache_path   <- file.path(P$cache, "H3_additional_moderators_models.rds")

assert_required_cols <- function(dat, required_cols, label) {
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop(
      sprintf("%s is missing required columns: %s", label, paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }
}

assert_unique_keys <- function(dat, keys, label) {
  dupes <- dat[duplicated(dat[keys]) | duplicated(dat[keys], fromLast = TRUE), keys, drop = FALSE]
  if (nrow(dupes) > 0) {
    stop(sprintf("Duplicate %s keys detected.", label), call. = FALSE)
  }
}

assert_no_nonfinite <- function(x, label) {
  bad <- sum(!is.finite(x[!is.na(x)]))
  if (bad > 0) {
    stop(sprintf("%s contains %d non-finite non-missing values.", label, bad), call. = FALSE)
  }
}

stars_for_p <- function(p) {
  if (is.na(p)) {
    ""
  } else if (p < 0.01) {
    "***"
  } else if (p < 0.05) {
    "**"
  } else if (p < 0.10) {
    "*"
  } else {
    ""
  }
}

fmt_coef <- function(est, p) sprintf("%.3f%s", est, stars_for_p(p))
fmt_p <- function(p) ifelse(is.na(p), "NA", sprintf("%.3f", p))
fmt_b <- function(b) ifelse(is.na(b), "NA", sprintf("%.3f", b))

###############################################################################
# Join the frozen campaign-level moderator input. The frozen file is part of
# the replication deposit so this script never depends on a live V-Dem release.
###############################################################################

df <- read.csv(navco13_path)
assert_required_cols(
  df,
  c(
    "CAMPAIGN", "LOCATION", "iso3", "EYEAR", "BYEAR", "success", "limited",
    "nonviolent_camp", "REGCHANGE", "eng_prob_general", "lnnum_image_sum",
    "lnlengthofcam", "v2csreprss", "v2x_civlib", "v2x_clpol", "lngdp"
  ),
  "df_final.csv"
)
assert_unique_keys(df, c("CAMPAIGN", "LOCATION"), "df_final CAMPAIGN-LOCATION")

added_mods <- read.csv(moderator_path)
assert_required_cols(
  added_mods,
  c("CAMPAIGN", "LOCATION", "iso3", "BYEAR", "EYEAR", "v2x_frassoc_thick"),
  "h3_additional_moderators.csv"
)
assert_unique_keys(
  added_mods, c("CAMPAIGN", "LOCATION"),
  "h3_additional_moderators CAMPAIGN-LOCATION"
)
if (nrow(added_mods) != nrow(df)) {
  stop("Frozen H3 moderator input must contain exactly one row per campaign.", call. = FALSE)
}

df <- df %>%
  left_join(
    added_mods,
    by = c("CAMPAIGN", "LOCATION", "iso3", "BYEAR", "EYEAR"),
    relationship = "one-to-one"
  )
if (anyNA(df$v2x_frassoc_thick)) {
  stop("Frozen H3 moderator join produced missing freedom-of-association values.", call. = FALSE)
}

###############################################################################
# Table 4 campaign-level preparation
###############################################################################

df$nonviolent_camp <- as.factor(df$nonviolent_camp)
df$REGCHANGE <- as.numeric(as.character(df$REGCHANGE))

df$duration <- df$EYEAR - df$BYEAR + 1
df$event_concession <- as.integer(
  as.numeric(as.character(df$success)) == 1 |
    as.numeric(as.character(df$limited)) == 1
)
df$concession_bin <- df$event_concession

assert_no_nonfinite(df$eng_prob_general, "df$eng_prob_general")
assert_no_nonfinite(df$lnnum_image_sum, "df$lnnum_image_sum")
assert_no_nonfinite(df$v2x_clpol, "df$v2x_clpol")
assert_no_nonfinite(df$v2x_frassoc_thick, "df$v2x_frassoc_thick")

df$eng_c <- df$eng_prob_general - mean(df$eng_prob_general, na.rm = TRUE)
df$clpol_c <- df$v2x_clpol - mean(df$v2x_clpol, na.rm = TRUE)
df$frassoc_c <- df$v2x_frassoc_thick - mean(df$v2x_frassoc_thick, na.rm = TRUE)

###############################################################################
# Fit helpers: same Table 4 campaign-level specification as H3_code.r.
# Since both added moderators are domestic-openness indices, they use the
# non-repression branch: retain v2csreprss as the domestic control.
###############################################################################

fit_weibull13 <- function(moderator) {
  f <- as.formula(sprintf(
    "Surv(duration, event_concession) ~ eng_c * %s + nonviolent_camp + REGCHANGE + v2csreprss + lnnum_image_sum + lngdp",
    moderator
  ))
  survreg(f, data = df, dist = "weibull")
}

fit_firth13 <- function(moderator) {
  f <- as.formula(sprintf(
    "concession_bin ~ eng_c * %s + nonviolent_camp + REGCHANGE + v2csreprss + lnlengthofcam + lnnum_image_sum + lngdp",
    moderator
  ))
  logistf(
    f,
    data = df,
    control   = logistf.control(maxit = 1000, maxstep = 0.5),
    plcontrol = logistpl.control(maxit = 1000),
    pl        = TRUE
  )
}

cc_stats <- function(dat, vars, event_var) {
  cc_mask <- complete.cases(dat[, vars, drop = FALSE])
  list(
    N = sum(cc_mask),
    n_events = as.integer(sum(dat[[event_var]][cc_mask] == 1, na.rm = TRUE))
  )
}

vars_weibull13 <- function(moderator) {
  c(
    "duration", "event_concession", "eng_c", moderator, "nonviolent_camp",
    "REGCHANGE", "v2csreprss", "lnnum_image_sum", "lngdp"
  )
}

vars_firth13 <- function(moderator) {
  c(
    "concession_bin", "eng_c", moderator, "nonviolent_camp", "REGCHANGE",
    "v2csreprss", "lnlengthofcam", "lnnum_image_sum", "lngdp"
  )
}

extract_term <- function(model, term) {
  if (inherits(model, "survreg")) {
    tab <- summary(model)$table
    ses <- tab[, "Std. Error"]
    return(list(
      estimate = unname(tab[term, "Value"]),
      se = unname(ses[term]),
      p = unname(tab[term, "p"])
    ))
  }

  if (inherits(model, "logistf")) {
    idx <- match(term, names(model$coefficients))
    ses <- sqrt(diag(model$var))
    return(list(
      estimate = unname(model$coefficients[idx]),
      se = unname(ses[idx]),
      p = unname(model$prob[idx])
    ))
  }

  stop(sprintf("Unsupported model class: %s", paste(class(model), collapse = ",")), call. = FALSE)
}

joint_wald_p <- function(model, terms) {
  cf <- coef(model)
  ii <- match(terms, names(cf))
  if (any(is.na(ii))) return(NA_real_)

  V <- tryCatch(vcov(model), error = function(e) NULL)
  if (is.null(V) && inherits(model, "logistf")) V <- model$var
  if (is.null(V)) return(NA_real_)

  b <- cf[ii]
  V2 <- V[ii, ii, drop = FALSE]
  if (!all(is.finite(b)) || !all(is.finite(V2))) return(NA_real_)

  W <- tryCatch(as.numeric(t(b) %*% solve(V2) %*% b), error = function(e) NA_real_)
  if (!is.finite(W)) return(NA_real_)
  pchisq(W, df = length(terms), lower.tail = FALSE)
}

extract_model <- function(tag, estimator, moderator, moderator_label, model, stats) {
  int_term <- paste0("eng_c:", moderator)
  eng <- extract_term(model, "eng_c")
  mod <- extract_term(model, moderator)
  int <- extract_term(model, int_term)
  expected <- if (estimator == "Weibull AFT") "-" else "+"
  sign_consistent <- if (expected == "-") int$estimate < 0 else int$estimate > 0

  data.frame(
    model = tag,
    estimator = estimator,
    data = "NAVCO 1.3",
    response = if (estimator == "Weibull AFT") {
      "Time to broad success/concession"
    } else {
      "Broad success/concession (binary)"
    },
    source_outcome = if (estimator == "Weibull AFT") "event_concession" else "concession_bin",
    moderator = moderator,
    moderator_label = moderator_label,
    beta_eng_c = eng$estimate,
    se_eng_c = eng$se,
    p_eng_c = eng$p,
    beta_mod_c = mod$estimate,
    se_mod_c = mod$se,
    p_mod_c = mod$p,
    beta_int = int$estimate,
    se_int = int$se,
    p_int = int$p,
    joint_wald_p = joint_wald_p(model, c("eng_c", int_term)),
    expected_int_sign = expected,
    sign_consistent = sign_consistent,
    N = stats$N,
    n_events = stats$n_events,
    epv = stats$n_events / length(coef(model)),
    stringsAsFactors = FALSE
  )
}

mods <- list(
  list(name = "clpol_c", label = "Political Liberty (V-Dem v2x_clpol)"),
  list(name = "frassoc_c", label = "Freedom of Association (V-Dem v2x_frassoc_thick)")
)

model_list <- list()
results <- list()

for (i in seq_along(mods)) {
  mod <- mods[[i]]

  weibull <- fit_weibull13(mod$name)
  firth <- fit_firth13(mod$name)
  s_weibull <- cc_stats(df, vars_weibull13(mod$name), "event_concession")
  s_firth <- cc_stats(df, vars_firth13(mod$name), "concession_bin")

  model_list[[paste0("ADD", i, "W")]] <- weibull
  model_list[[paste0("ADD", i, "F")]] <- firth

  results[[length(results) + 1L]] <- extract_model(
    paste0("ADD", i, "W"), "Weibull AFT", mod$name, mod$label, weibull, s_weibull
  )
  results[[length(results) + 1L]] <- extract_model(
    paste0("ADD", i, "F"), "Firth PML", mod$name, mod$label, firth, s_firth
  )
}

results <- bind_rows(results)
results$p_int_holm_added_k4 <- p.adjust(results$p_int, method = "holm", n = nrow(results))

if (!file.exists(h3_results_path)) {
  stop("H3_results.csv is required before fitting the added H3 moderators.", call. = FALSE)
}
h3_results <- read.csv(h3_results_path)
assert_required_cols(h3_results, c("model", "p_int", "main_panel"), "H3_results.csv")
h3_main_p <- h3_results$p_int[h3_results$main_panel %in% TRUE]
if (length(h3_main_p) != 8L || any(!is.finite(h3_main_p))) {
  stop("H3_results.csv must contain exactly eight finite main-panel interaction p-values.", call. = FALSE)
}
augmented_p <- p.adjust(c(h3_main_p, results$p_int), method = "holm", n = 12L)
results$p_int_holm_augmented_k12 <- tail(augmented_p, nrow(results))

results$coef_int_display <- mapply(fmt_coef, results$beta_int, results$p_int)
results$se_int_display <- sprintf("(%.3f)", results$se_int)

write.csv(results, out_path, row.names = FALSE)
stopifnot(identical(names(model_list), results$model))
saveRDS(model_list, cache_path, version = 3)

cat("\n=========== H3 / Table 4 Additional Moderators ===========\n")
cat("Specification: campaign-level Table 4 models; same controls and centering as H3_code.r.\n")
cat("Outcome: broad campaign success/concession (NAVCO success OR limited concession).\n\n")
cat(sprintf(
  "%-7s %-12s %-52s %10s %8s %8s %8s %5s %5s %5s\n",
  "Model", "Estimator", "Moderator", "b_int", "p_int", "p_Holm4", "p_Holm12", "Sign", "N", "Events"
))
cat(paste(rep("-", 135), collapse = ""), "\n")
for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  cat(sprintf(
    "%-7s %-12s %-52s %10s %8s %8s %8s %5s %5d %5d\n",
    r$model,
    r$estimator,
    r$moderator_label,
    fmt_b(r$beta_int),
    fmt_p(r$p_int),
    fmt_p(r$p_int_holm_added_k4),
    fmt_p(r$p_int_holm_augmented_k12),
    ifelse(isTRUE(r$sign_consistent), "cons", "wrong"),
    r$N,
    r$n_events
  ))
}

cat("\nNotes:\n")
cat("  p_Holm4  = Holm-adjusted across the four added Table 4-style models only.\n")
cat("  p_Holm12 = Holm-adjusted after appending these four models to the existing eight Table 4 main-panel models.\n")
cat("  Expected sign: negative for Weibull AFT, positive for Firth PML.\n")
cat(sprintf("\nResults exported to: %s\n", out_path))
cat(sprintf("Model bundle exported to: %s\n", cache_path))
