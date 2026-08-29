# Shared helpers for the power, calibration, diagnostic, and figure scripts.
# These functions intentionally do not read from the accepted-reference tree
# unless a caller explicitly asks for a comparison. Fresh estimates always go
# to output/; accepted evidence remains immutable under reference/accepted/.

resolve_replication_root <- function() {
  env_root <- Sys.getenv("REPLICATION_ROOT", unset = "")
  if (nzchar(env_root)) return(normalizePath(env_root, mustWork = TRUE))
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) != 1L) {
    stop("Set REPLICATION_ROOT when sourcing a replication script interactively.", call. = FALSE)
  }
  script <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
  normalizePath(file.path(dirname(script), "..", ".."), mustWork = TRUE)
}

read_positive_integer <- function(name, default) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) return(as.integer(default))
  value <- suppressWarnings(as.integer(raw))
  if (length(value) != 1L || is.na(value) || value < 1L) {
    stop(sprintf("%s must be a positive integer", name), call. = FALSE)
  }
  value
}

require_files <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop(sprintf("Missing required artifact(s): %s", paste(missing, collapse = ", ")),
         call. = FALSE)
  }
  invisible(paths)
}

rep_model_formula <- function(model) {
  if (inherits(model, "logistf") && inherits(model$formula, "formula")) return(model$formula)
  formula(model)
}

named_seed <- function(family, spec, magnitude = 0, base_seed = 20260524L) {
  token <- paste(base_seed, family, spec,
                 format(magnitude, scientific = FALSE, trim = TRUE), sep = "|")
  # Polynomial hashing is evaluated as a double below 2^53, so every integer
  # operation is exact and platform-stable.  This replaces the old character
  # sum, which produced collisions across distinct family/spec/magnitude cells.
  modulus <- 2147483629
  hash <- 0
  for (byte in utf8ToInt(enc2utf8(token))) hash <- (hash * 131 + byte) %% modulus
  as.integer(hash + 1)
}

power_seed_registry <- function() {
  magnitudes <- c(0, .25, .5, 1, 2)
  expand <- function(family, specs) {
    grid <- expand.grid(spec = specs, magnitude = magnitudes,
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    transform(grid, family = family)[c("family", "spec", "magnitude")]
  }
  power <- rbind(
    expand("H1", c("M1", "M2", "M5")),
    expand("H2.1", c("M1", "M2", "M3", "M4", "M9", "M10")),
    expand("H2.2", paste0("M", 1:8)),
    expand("H3", paste0("M", 1:8))
  )
  calibration <- rbind(
    data.frame(family = "CAL-H1", spec = "M1", magnitude = 0),
    data.frame(family = "CAL-H2.1", spec = c("M1", "M2"), magnitude = 0),
    data.frame(family = "CAL-H2.2", spec = paste0("M", 1:4), magnitude = 0),
    data.frame(family = "CAL-H3", spec = paste0("M", 1:4), magnitude = 0)
  )
  registry <- rbind(power, calibration)
  registry$seed <- mapply(named_seed, registry$family, registry$spec,
                          registry$magnitude, USE.NAMES = FALSE)
  registry$cell_id <- paste(registry$family, registry$spec,
                            format(registry$magnitude, scientific = FALSE, trim = TRUE),
                            sep = "|")
  registry
}

assert_unique_power_seeds <- function() {
  registry <- power_seed_registry()
  duplicated_seed <- duplicated(registry$seed) | duplicated(registry$seed, fromLast = TRUE)
  if (any(duplicated_seed)) {
    collisions <- split(registry$cell_id[duplicated_seed], registry$seed[duplicated_seed])
    stop(sprintf("RNG seed collision(s): %s",
                 paste(vapply(collisions, paste, collapse = " = ", character(1)),
                       collapse = "; ")), call. = FALSE)
  }
  invisible(registry)
}

safe_model_frame <- function(model) {
  mf <- tryCatch(model.frame(rep_model_formula(model)), error = function(e) NULL)
  if (!is.null(mf)) return(mf)
  if (!is.null(model$model) && is.data.frame(model$model)) return(model$model)
  frames <- getOption("replication.analysis.frames")
  data_name <- tryCatch(as.character(model$call$data), error = function(e) "")
  if (is.list(frames) && length(data_name) == 1L && data_name %in% names(frames)) {
    return(model.frame(rep_model_formula(model), data = frames[[data_name]], na.action = stats::na.omit))
  }
  stop(sprintf("Cannot recover model frame for class %s", class(model)[1]), call. = FALSE)
}

load_analysis_frames <- function(cache_dir, family) {
  path <- file.path(cache_dir, paste0(family, "_analysis_frames.rds"))
  require_files(path)
  frames <- readRDS(path)
  if (!is.list(frames) || !all(c("df", "d21") %in% names(frames))) {
    stop(sprintf("%s must contain named df and d21 frames.", basename(path)), call. = FALSE)
  }
  options(replication.analysis.frames = frames)
  invisible(frames)
}

model_simulation_data <- function(model) {
  mf <- safe_model_frame(model)
  model_formula <- rep_model_formula(model)
  response_vars <- all.vars(model_formula[[2]])
  response_col <- attr(terms(model_formula), "response")
  if (is.null(response_col) || !length(response_col)) response_col <- 1L
  predictors <- if (response_col > 0L) mf[, -response_col, drop = FALSE] else mf
  rhs_vars <- all.vars(delete.response(terms(model_formula)))
  absent <- setdiff(rhs_vars, names(predictors))
  for (name in absent) {
    hit <- grep(sprintf("\\(%s\\)$", name), names(predictors), value = TRUE)
    if (length(hit) == 1L) predictors[[name]] <- predictors[[hit]]
  }

  if (inherits(model, c("survreg", "coxph"))) {
    y <- model$y
    if (is.null(y)) y <- model.response(mf)
    if (!(length(response_vars) %in% c(2L, 3L)) || ncol(as.matrix(y)) < 2L) {
      stop("Expected a two- or three-column Surv response.", call. = FALSE)
    }
    values <- if (length(response_vars) == 2L) {
      list(as.numeric(y[, 1]), as.integer(y[, ncol(y)]))
    } else {
      list(as.numeric(y[, 1]), as.numeric(y[, 2]), as.integer(y[, ncol(y)]))
    }
    outcome <- data.frame(setNames(values, response_vars), check.names = FALSE)
  } else {
    if (length(response_vars) != 1L) {
      stop("Expected one response variable for binary power simulation.", call. = FALSE)
    }
    outcome <- data.frame(setNames(list(as.numeric(model.response(mf))), response_vars),
                          check.names = FALSE)
  }
  cbind(outcome, predictors)
}

aligned_model_matrix <- function(model, data) {
  rhs <- delete.response(terms(rep_model_formula(model)))
  mf <- model.frame(rhs, data = data, na.action = na.pass, xlev = model$xlevels)
  X <- model.matrix(rhs, mf, contrasts.arg = model$contrasts)
  beta <- coef(model)
  missing <- setdiff(names(beta), colnames(X))
  if (length(missing)) {
    stop(sprintf("Design matrix misses coefficient(s): %s", paste(missing, collapse = ", ")),
         call. = FALSE)
  }
  X[, names(beta), drop = FALSE]
}

simulate_binary_model <- function(model, focal, magnitude) {
  data <- model_simulation_data(model)
  X <- aligned_model_matrix(model, data)
  beta <- coef(model)
  if (!(focal %in% names(beta))) stop(sprintf("Missing focal term %s", focal), call. = FALSE)
  beta[focal] <- magnitude * beta[focal]
  response <- all.vars(rep_model_formula(model)[[2]])[[1]]
  data[[response]] <- stats::rbinom(nrow(data), 1L, stats::plogis(drop(X %*% beta)))
  data
}

simulate_weibull_model <- function(model, focal, magnitude) {
  data <- model_simulation_data(model)
  X <- aligned_model_matrix(model, data)
  beta <- coef(model)
  if (!(focal %in% names(beta))) stop(sprintf("Missing focal term %s", focal), call. = FALSE)
  beta[focal] <- magnitude * beta[focal]
  response_vars <- all.vars(rep_model_formula(model)[[2]])
  censor_time <- data[[response_vars[[1]]]]
  event_time <- stats::rweibull(
    nrow(data), shape = 1 / model$scale, scale = exp(drop(X %*% beta))
  )
  positive <- censor_time[is.finite(censor_time) & censor_time > 0]
  floor_time <- if (length(positive)) min(positive) / 1000 else 1e-6
  data[[response_vars[[1]]]] <- pmax(pmin(event_time, censor_time), floor_time)
  data[[response_vars[[2]]]] <- as.integer(event_time <= censor_time)
  data
}

assert_refit_safe <- function(model, data, minimum_survival_events = 3L) {
  response_vars <- all.vars(rep_model_formula(model)[[2]])
  if (inherits(model, "survreg")) {
    event <- data[[tail(response_vars, 1L)]]
    duration <- data[[response_vars[[1L]]]]
    if (any(!is.finite(duration)) || any(duration <= 0)) {
      stop("Simulated survival durations must be finite and positive", call. = FALSE)
    }
    event_count <- sum(event == 1L, na.rm = TRUE)
    censor_count <- sum(event == 0L, na.rm = TRUE)
    if (event_count < minimum_survival_events || censor_count < 3L) {
      stop(sprintf("Unsafe survreg refit skipped: events=%d, censored=%d",
                   event_count, censor_count), call. = FALSE)
    }
  } else if (inherits(model, "logistf")) {
    response <- data[[response_vars[[1L]]]]
    if (any(!is.finite(response)) || any(!response %in% c(0, 1))) {
      stop("Simulated binary response must contain only finite 0/1 values",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

refit_focal_p <- function(model, data, focal) {
  assert_refit_safe(model, data)
  fit <- if (inherits(model, "survreg")) {
    survival::survreg(rep_model_formula(model), data = data, dist = model$dist)
  } else if (inherits(model, "logistf")) {
    logistf::logistf(
      rep_model_formula(model), data = data, pl = TRUE,
      control = logistf::logistf.control(maxit = 1000, maxstep = 0.5),
      plcontrol = logistf::logistpl.control(maxit = 1000)
    )
  } else {
    stop(sprintf("Unsupported power estimator: %s", class(model)[1]), call. = FALSE)
  }
  if (inherits(fit, "survreg")) {
    unname(summary(fit)$table[focal, "p"])
  } else {
    unname(fit$prob[focal])
  }
}

wald_joint_p <- function(model, terms) {
  beta <- coef(model)
  if (!all(terms %in% names(beta))) return(NA_real_)
  V <- tryCatch(vcov(model), error = function(e) model$var)
  if (is.null(dimnames(V))) dimnames(V) <- list(names(beta), names(beta))
  block <- V[terms, terms, drop = FALSE]
  b <- beta[terms]
  statistic <- tryCatch(drop(t(b) %*% solve(block, b)), error = function(e) NA_real_)
  if (!is.finite(statistic)) return(NA_real_)
  stats::pchisq(statistic, df = length(terms), lower.tail = FALSE)
}

wald_single_p <- function(model, term) {
  beta <- coef(model)
  if (!(term %in% names(beta))) return(NA_real_)
  V <- tryCatch(vcov(model), error = function(error) model$var)
  if (is.null(dimnames(V))) dimnames(V) <- list(names(beta), names(beta))
  variance <- V[term, term]
  if (!is.finite(beta[[term]]) || !is.finite(variance) || variance <= 0) return(NA_real_)
  2 * stats::pnorm(-abs(beta[[term]] / sqrt(variance)))
}

refit_h1_tests <- function(model, data, linear, quadratic) {
  assert_refit_safe(model, data)
  fit <- if (inherits(model, "survreg")) {
    survival::survreg(rep_model_formula(model), data = data, dist = model$dist)
  } else {
    logistf::logistf(
      rep_model_formula(model), data = data, pl = TRUE,
      control = logistf::logistf.control(maxit = 1000, maxstep = 0.5),
      plcontrol = logistf::logistpl.control(maxit = 1000)
    )
  }
  focal_p <- if (inherits(fit, "survreg")) {
    unname(summary(fit)$table[quadratic, "p"])
  } else {
    wald_single_p(fit, quadratic)
  }
  c(any_eng_effect_eng_eng2 = wald_joint_p(fit, c(linear, quadratic)),
    p_eng2_single_coef = focal_p)
}

run_h1_power <- function(model, spec, estimator, n_reps, k_family = 10L,
                         magnitudes = c(0, .25, .5, 1, 2),
                         alphas = c(.05, .10)) {
  beta <- coef(model)
  quadratic <- grep("\\^2", names(beta), value = TRUE)
  if (length(quadratic) != 1L) stop(sprintf("Could not identify one H1 quadratic term for %s", spec))
  raw_name <- sub("^I\\((.*)\\^2\\)$", "\\1", quadratic)
  linear <- names(beta)[names(beta) == raw_name]
  if (length(linear) != 1L) stop(sprintf("Could not identify H1 linear term for %s", spec))

  rows <- lapply(magnitudes, function(magnitude) {
    set.seed(named_seed("H1", spec, magnitude))
    rep_seeds <- sample.int(.Machine$integer.max, n_reps)
    p <- vapply(rep_seeds, function(seed) {
      set.seed(seed)
      tryCatch({
        data <- model_simulation_data(model)
        X <- aligned_model_matrix(model, data)
        beta_sim <- beta
        beta_sim[c(linear, quadratic)] <- magnitude * beta[c(linear, quadratic)]
        response_vars <- all.vars(rep_model_formula(model)[[2]])
        if (inherits(model, "survreg")) {
          censor <- data[[response_vars[[1]]]]
          event_time <- stats::rweibull(nrow(data), 1 / model$scale,
                                        exp(drop(X %*% beta_sim)))
          data[[response_vars[[1]]]] <- pmax(pmin(event_time, censor), 1e-6)
          data[[response_vars[[2]]]] <- as.integer(event_time <= censor)
        } else {
          data[[response_vars[[1]]]] <- stats::rbinom(nrow(data), 1L,
                                                       stats::plogis(drop(X %*% beta_sim)))
        }
        refit_h1_tests(model, data, linear, quadratic)
      }, error = function(e) c(any_eng_effect_eng_eng2 = NA_real_,
                               p_eng2_single_coef = NA_real_))
    }, numeric(2))
    do.call(rbind, lapply(rownames(p), function(test) {
      values <- p[test, ]
      fit_ok <- mean(is.finite(values))
      do.call(rbind, lapply(alphas, function(alpha) {
        valid <- values[is.finite(values)]
        raw <- sum(values < alpha, na.rm = TRUE) / n_reps
        holm <- sum(values < alpha / k_family, na.rm = TRUE) / n_reps
        data.frame(
          spec = spec, estimator = estimator,
          moderator = "(none - H1 unconditional curvilinear)", test = test,
          beta_eng_hat = unname(beta[linear]), beta_eng2_hat = unname(beta[quadratic]),
          magnitude = magnitude, beta_eng_sim = magnitude * unname(beta[linear]),
          beta_eng2_sim = magnitude * unname(beta[quadratic]), n_reps = n_reps,
          fit_ok_rate = fit_ok, alpha = alpha,
          raw_power_overall = raw, holm_proxy_power_overall = holm,
          raw_power_conditional = if (length(valid)) mean(valid < alpha) else NA_real_,
          holm_proxy_power_conditional = if (length(valid)) mean(valid < alpha / k_family) else NA_real_,
          raw_power = raw, holm_proxy_power = holm, stringsAsFactors = FALSE
        )
      }))
    }))
  })
  do.call(rbind, rows)
}

run_focal_power <- function(model, family, spec, estimator, moderator, focal,
                            k_family, n_reps, magnitudes = c(0, .25, .5, 1, 2),
                            alphas = c(.05, .10)) {
  beta_hat <- unname(coef(model)[focal])
  if (!is.finite(beta_hat)) stop(sprintf("Non-finite fitted focal coefficient for %s %s", family, spec))
  rows <- lapply(magnitudes, function(magnitude) {
    set.seed(named_seed(family, spec, magnitude))
    rep_seeds <- sample.int(.Machine$integer.max, n_reps)
    p <- vapply(rep_seeds, function(seed) {
      set.seed(seed)
      tryCatch({
        sim <- if (inherits(model, "survreg")) {
          simulate_weibull_model(model, focal, magnitude)
        } else {
          simulate_binary_model(model, focal, magnitude)
        }
        refit_focal_p(model, sim, focal)
      }, error = function(e) NA_real_)
    }, numeric(1))
    fit_ok <- mean(is.finite(p))
    do.call(rbind, lapply(alphas, function(alpha) {
      valid <- p[is.finite(p)]
      raw_overall <- sum(p < alpha, na.rm = TRUE) / n_reps
      holm_overall <- sum(p < alpha / k_family, na.rm = TRUE) / n_reps
      data.frame(
        spec = spec, estimator = estimator, moderator = moderator,
        focal_coef = focal, beta_focal_hat = beta_hat, magnitude = magnitude,
        beta_focal_sim = magnitude * beta_hat, n_reps = n_reps,
        fit_ok_rate = fit_ok, alpha = alpha,
        raw_power_overall = raw_overall,
        holm_proxy_power_overall = holm_overall,
        raw_power_conditional = if (length(valid)) mean(valid < alpha) else NA_real_,
        holm_proxy_power_conditional = if (length(valid)) mean(valid < alpha / k_family) else NA_real_,
        raw_power = raw_overall, holm_proxy_power = holm_overall,
        stringsAsFactors = FALSE
      )
    }))
  })
  do.call(rbind, rows)
}

interp_mde <- function(magnitude, power, target = .80) {
  keep <- is.finite(magnitude) & is.finite(power) & magnitude > 0
  magnitude <- magnitude[keep]
  power <- power[keep]
  if (!length(magnitude)) return(NA_real_)
  ord <- order(magnitude)
  magnitude <- magnitude[ord]
  power <- power[ord]
  if (power[[1]] >= target) return(magnitude[[1]])
  crossing <- which(power[-1] >= target & power[-length(power)] < target)
  if (!length(crossing)) return(Inf)
  i <- crossing[[1]]
  magnitude[[i]] + (target - power[[i]]) *
    (magnitude[[i + 1L]] - magnitude[[i]]) / (power[[i + 1L]] - power[[i]])
}

power_summary <- function(power) {
  required <- c("spec", "estimator", "moderator", "beta_focal_hat", "magnitude",
                "alpha", "raw_power", "holm_proxy_power", "fit_ok_rate")
  missing <- setdiff(required, names(power))
  if (length(missing)) {
    stop(sprintf("power_summary() is for focal-interaction families and is missing: %s",
                 paste(missing, collapse = ", ")), call. = FALSE)
  }
  keys <- unique(power[c("spec", "estimator", "moderator")])
  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, ]
    same_moderator <- if (is.na(key$moderator)) is.na(power$moderator) else
      !is.na(power$moderator) & power$moderator == key$moderator
    x <- power[power$spec == key$spec & power$estimator == key$estimator &
                 same_moderator & power$alpha == .05, , drop = FALSE]
    if (!nrow(x)) stop(sprintf("No alpha=.05 power rows for %s", key$spec), call. = FALSE)
    at <- function(m, column = "raw_power") {
      value <- x[x$magnitude == m, column]
      if (length(value)) value[[1]] else NA_real_
    }
    data.frame(
      spec = key$spec, estimator = key$estimator, moderator = key$moderator,
      beta_focal_hat = x$beta_focal_hat[[1]],
      type1_at_null_05 = at(0), raw_power_at_half_05 = at(.5),
      raw_power_at_fitted_05 = at(1), raw_power_at_double_05 = at(2),
      holm_power_at_fitted_05 = at(1, "holm_proxy_power"),
      raw_mde_mult_05 = min(x$magnitude[x$magnitude > 0 & x$raw_power >= .8], Inf),
      holm_mde_mult_05 = min(x$magnitude[x$magnitude > 0 & x$holm_proxy_power >= .8], Inf),
      mde_at_80_raw_interp = interp_mde(x$magnitude, x$raw_power),
      mde_at_80_holm_interp = interp_mde(x$magnitude, x$holm_proxy_power),
      fit_ok_rate_min = min(x$fit_ok_rate), fit_ok_rate_mean = mean(x$fit_ok_rate),
      stringsAsFactors = FALSE
    )
  }))
}

typical_value <- function(x) {
  if (is.factor(x)) return(factor(names(which.max(table(x))), levels = levels(x)))
  values <- unique(stats::na.omit(x))
  if (length(values) <= 2L && all(values %in% c(0, 1))) return(as.numeric(mean(x, na.rm = TRUE) >= .5))
  mean(x, na.rm = TRUE)
}

prediction_grid <- function(model, x_name, moderator = NULL, grid_n = 80L) {
  data <- model_simulation_data(model)
  x <- data[[x_name]]
  x_grid <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = grid_n)
  rhs_vars <- all.vars(delete.response(terms(rep_model_formula(model))))
  base_vars <- setdiff(rhs_vars, c(x_name, moderator))
  base <- as.data.frame(lapply(data[base_vars], typical_value), stringsAsFactors = FALSE)
  levels <- if (is.null(moderator)) NA_real_ else
    as.numeric(stats::quantile(data[[moderator]], c(.10, .50, .90), na.rm = TRUE))
  do.call(rbind, lapply(levels, function(level) {
    nd <- base[rep(1L, grid_n), , drop = FALSE]
    nd[[x_name]] <- x_grid
    if (!is.null(moderator)) nd[[moderator]] <- level
    X <- aligned_model_matrix(model, nd)
    beta <- coef(model)
    if (inherits(model, "survreg")) {
      fit <- exp(drop(X %*% beta)) * log(2)^model$scale
      response_scale <- "Predicted median duration"
    } else if (inherits(model, "logistf")) {
      fit <- stats::plogis(drop(X %*% beta))
      response_scale <- "Predicted probability"
    } else {
      fit <- exp(drop(X %*% beta))
      response_scale <- "Relative hazard"
    }
    data.frame(x = x_grid, moderator_level = level, fit = fit,
               response_scale = response_scale, stringsAsFactors = FALSE)
  }))
}

wald_forest_row <- function(model, term) {
  beta <- coef(model)
  if (!(term %in% names(beta))) return(c(beta = NA, se = NA, lo = NA, hi = NA))
  V <- tryCatch(vcov(model), error = function(e) model$var)
  if (nrow(V) != length(beta)) {
    if (!is.null(rownames(V)) && all(names(beta) %in% rownames(V))) {
      V <- V[names(beta), names(beta), drop = FALSE]
    } else {
      V <- V[seq_along(beta), seq_along(beta), drop = FALSE]
    }
  }
  if (is.null(dimnames(V))) dimnames(V) <- list(names(beta), names(beta))
  se <- sqrt(V[term, term])
  c(beta = unname(beta[term]), se = se,
    lo = unname(beta[term]) - qnorm(.975) * se,
    hi = unname(beta[term]) + qnorm(.975) * se)
}

predict_aft_time_ci <- function(model, newdata, n_sim = 1000L, level = .95) {
  beta <- coef(model)
  V <- vcov(model)
  scale_name <- "Log(scale)"
  if (!(scale_name %in% rownames(V))) stop("survreg vcov lacks Log(scale)")
  mu <- c(beta, setNames(log(model$scale), scale_name))
  X <- aligned_model_matrix(model, newdata)
  draws <- MASS::mvrnorm(n_sim, mu = mu, Sigma = V)
  lp <- X %*% t(draws[, names(beta), drop = FALSE])
  sigma <- exp(draws[, scale_name])
  values <- vapply(seq_len(n_sim), function(i) exp(lp[, i]) * log(2)^sigma[i],
                   numeric(nrow(X)))
  q <- c((1 - level) / 2, 1 - (1 - level) / 2)
  data.frame(fit = rowMeans(values),
             lwr = apply(values, 1, quantile, q[1], names = FALSE),
             upr = apply(values, 1, quantile, q[2], names = FALSE))
}

predict_firth_ci <- function(model, newdata, n_sim = 1000L, level = .95) {
  beta <- coef(model)
  V <- model$var
  if (is.null(dimnames(V))) dimnames(V) <- list(names(beta), names(beta))
  X <- aligned_model_matrix(model, newdata)
  draws <- MASS::mvrnorm(n_sim, mu = beta, Sigma = V)
  values <- plogis(X %*% t(draws))
  q <- c((1 - level) / 2, 1 - (1 - level) / 2)
  data.frame(fit = rowMeans(values),
             lwr = apply(values, 1, quantile, q[1], names = FALSE),
             upr = apply(values, 1, quantile, q[2], names = FALSE))
}

predict_aft_event_probability_ci <- function(model, newdata, reference_time,
                                             n_sim = 1000L, level = .95) {
  beta <- coef(model)
  V <- vcov(model)
  scale_name <- "Log(scale)"
  mu <- c(beta, setNames(log(model$scale), scale_name))
  X <- aligned_model_matrix(model, newdata)
  draws <- MASS::mvrnorm(n_sim, mu = mu, Sigma = V)
  lp <- X %*% t(draws[, names(beta), drop = FALSE])
  sigma <- exp(draws[, scale_name])
  values <- vapply(seq_len(n_sim), function(i) {
    pweibull(reference_time, shape = 1 / sigma[i], scale = exp(lp[, i]))
  }, numeric(nrow(X)))
  q <- c((1 - level) / 2, 1 - (1 - level) / 2)
  data.frame(fit = rowMeans(values),
             lwr = apply(values, 1, quantile, q[1], names = FALSE),
             upr = apply(values, 1, quantile, q[2], names = FALSE))
}

predict_cox_event_probability_ci <- function(model, newdata, reference_time,
                                             level = .95) {
  fit <- survival::survfit(model, newdata = newdata, conf.type = "log-log", conf.int = level)
  at_time <- summary(fit, times = reference_time, extend = TRUE)
  survival <- as.numeric(at_time$surv)
  lower <- as.numeric(at_time$lower)
  upper <- as.numeric(at_time$upper)
  if (length(survival) != nrow(newdata)) stop("Cox prediction count mismatch")
  data.frame(fit = pmin(pmax(1 - survival, 0), 1),
             lwr = pmin(pmax(1 - upper, 0), 1),
             upr = pmin(pmax(1 - lower, 0), 1))
}

turning_point_delta <- function(model, linear, quadratic, level = .95) {
  beta <- coef(model)
  V <- vcov(model)
  b1 <- beta[[linear]]; b2 <- beta[[quadratic]]
  point <- -b1 / (2 * b2)
  gradient <- c(-1 / (2 * b2), b1 / (2 * b2^2))
  variance <- drop(t(gradient) %*% V[c(linear, quadratic), c(linear, quadratic)] %*% gradient)
  se <- sqrt(variance)
  z <- qnorm(1 - (1 - level) / 2)
  list(tp = point, lo = point - z * se, hi = point + z * se)
}

build_conditioned_panel <- function(model, data, raw_moderator, centered_moderator,
                                    controls, estimator, model_id, moderator_label,
                                    grid_n = 80L) {
  required <- unique(c("eng_prob_general", raw_moderator, "eng_c", centered_moderator,
                       controls, all.vars(rep_model_formula(model)[[2]])))
  analytic <- data[complete.cases(data[, required, drop = FALSE]), , drop = FALSE]
  eng_grid <- seq(0, max(analytic$eng_prob_general), length.out = grid_n)
  moderator_values <- unname(quantile(analytic[[raw_moderator]], c(.10, .90)))
  base <- as.data.frame(lapply(analytic[controls], typical_value), stringsAsFactors = FALSE)
  eng_mean <- mean(data$eng_prob_general, na.rm = TRUE)
  moderator_mean <- mean(data[[raw_moderator]], na.rm = TRUE)
  pred <- do.call(rbind, lapply(seq_along(moderator_values), function(i) {
    nd <- base[rep(1L, grid_n), , drop = FALSE]
    nd$eng_c <- eng_grid - eng_mean
    nd[[centered_moderator]] <- moderator_values[i] - moderator_mean
    interval <- if (estimator == "aft") predict_aft_time_ci(model, nd) else predict_firth_ci(model, nd)
    data.frame(model_id = model_id, moderator_label = moderator_label,
               eng = eng_grid, moderator_percentile = factor(c("10th", "90th")[i],
                 levels = c("10th", "90th")), moderator_value = moderator_values[i],
               fit = interval$fit, lwr = interval$lwr, upr = interval$upr)
  }))
  list(pred = pred, rug = data.frame(model_id = model_id, moderator_label = moderator_label,
                                     eng = analytic$eng_prob_general))
}
