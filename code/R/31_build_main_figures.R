root <- Sys.getenv("REPLICATION_ROOT", unset = "")
if (!nzchar(root)) {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  file_arg <- sub("^--file=", "", file_arg[[1L]])
  file_arg <- gsub("~+~", " ", file_arg, fixed = TRUE)
  root <- normalizePath(file.path(dirname(file_arg), "..", ".."), mustWork = TRUE)
}
source(file.path(root, "code", "R", "00_setup.R"))
source(file.path(root, "code", "R", "R", "power_figure_helpers.R"))
P <- init_replication(root)
rep_require_packages(c("survival", "logistf", "MASS", "ggplot2", "patchwork", "png"))
suppressPackageStartupMessages({
  library(survival); library(logistf); library(MASS); library(ggplot2); library(patchwork); library(grid)
})

bundle_paths <- file.path(P$cache, c("H1_models.rds", "H22_models.rds", "H3_models.rds"))
require_files(bundle_paths)
h1 <- readRDS(bundle_paths[1]); h22 <- readRDS(bundle_paths[2]); h3 <- readRDS(bundle_paths[3])

save_figure <- function(plot, stem, plot_data, width, height) {
  write.csv(plot_data, file.path(P$figures, paste0(stem, "_plot_data.csv")), row.names = FALSE, na = "")
  ggsave(file.path(P$figures, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = cairo_pdf)
  ggsave(file.path(P$figures, paste0(stem, ".png")), plot, width = width, height = height,
         units = "in", dpi = 300)
}
base_theme <- theme_minimal(base_size = 12) + theme(
  panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#E7E7E7", linewidth = .35),
  legend.position = "bottom", axis.text = element_text(colour = "#444444"),
  plot.background = element_rect(fill = "white", colour = NA))

# Figure 7: accepted H1 M1/M4 event-probability algorithm.
set.seed(20260506)
h1_frames <- load_analysis_frames(P$cache, "H1")
df_h1 <- h1_frames$df; d21_fc <- h1_frames$d21_first_conc
# survfit.coxph re-evaluates the fitted call's data symbol.
d21_first_conc <- d21_fc
eng_aft <- seq(0, max(df_h1$eng_prob_general, na.rm = TRUE), length.out = 60L)
controls_aft <- c("nonviolent_camp", "REGCHANGE", "v2csreprss", "lnnum_image_sum", "lnpop", "colonized_english")
base_aft <- as.data.frame(lapply(df_h1[controls_aft], typical_value), stringsAsFactors = FALSE)
nd_aft <- cbind(data.frame(eng_prob_general = eng_aft), base_aft[rep(1L, length(eng_aft)), , drop = FALSE])
t_aft <- mean(df_h1$duration[df_h1$event_concession == 1], na.rm = TRUE)
ci_aft <- predict_aft_event_probability_ci(h1$M1, nd_aft, t_aft)
tp_aft <- turning_point_delta(h1$M1, "eng_prob_general", "I(eng_prob_general^2)")
f7_aft <- cbind(panel = "M1: Weibull AFT", eng = eng_aft, ci_aft,
                reference_time = t_aft, turning_point = tp_aft$tp,
                turning_lo = tp_aft$lo, turning_hi = tp_aft$hi)
eng_cox <- seq(0, max(d21_fc$eng_prop_year, na.rm = TRUE), length.out = 60L)
controls_cox <- c("nonviolent", "REGCHANGE", "v2csreprss", "lnnum_image_sum", "lnpop_yr", "colonized_english")
base_cox <- as.data.frame(lapply(d21_fc[controls_cox], typical_value), stringsAsFactors = FALSE)
nd_cox <- cbind(data.frame(eng_prop_year = eng_cox), base_cox[rep(1L, length(eng_cox)), , drop = FALSE], navco21_id = NA)
t_cox <- mean(d21_fc$tstop[d21_fc$event_concession == 1], na.rm = TRUE)
ci_cox <- predict_cox_event_probability_ci(h1$M4, nd_cox, t_cox)
tp_cox <- turning_point_delta(h1$M4, "eng_prop_year", "I(eng_prop_year^2)")
f7_cox <- cbind(panel = "M4: Cox PH", eng = eng_cox, ci_cox,
                reference_time = t_cox, turning_point = tp_cox$tp,
                turning_lo = tp_cox$lo, turning_hi = tp_cox$hi)
f7 <- rbind(f7_aft, f7_cox)
make_h1_panel <- function(data, rug, rug_name, colour) {
  max_x <- max(data$eng); tp <- unique(data[c("turning_point", "turning_lo", "turning_hi")])
  p <- ggplot(data, aes(eng, fit)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = colour, alpha = .20) +
    geom_line(colour = colour, linewidth = .85) +
    geom_rug(data = rug, aes(x = .data[[rug_name]]), inherit.aes = FALSE,
             sides = "b", alpha = .22, length = unit(.010, "npc")) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25), labels = function(x) sprintf("%.2f", x)) +
    scale_x_continuous(limits = c(0, max_x)) +
    labs(x = if (rug_name == "eng_prob_general") "English banner proportion (campaign-level)" else
           "English banner proportion (campaign-year)", y = "Predicted P(success)") + base_theme
  if (is.finite(tp$turning_point) && tp$turning_point >= 0 && tp$turning_point <= max_x &&
      all(is.finite(unlist(tp)))) {
    p <- p + annotate("rect", xmin = max(0, tp$turning_lo), xmax = min(max_x, tp$turning_hi),
                      ymin = -Inf, ymax = Inf, fill = "#F4A261", alpha = .20) +
      geom_vline(xintercept = tp$turning_point, linetype = "longdash", colour = "#E66101", linewidth = .45) +
      annotate("text", x = min(tp$turning_point + max_x * .006, max_x * .82), y = .985,
               label = sprintf("Turning point: %.3f", tp$turning_point), hjust = 0, vjust = 1,
               size = 2.25, colour = "#E66101")
  }
  p
}
p7 <- make_h1_panel(f7_aft, df_h1, "eng_prob_general", "#2C7FB8") +
      make_h1_panel(f7_cox, d21_fc, "eng_prop_year", "#1B9E77") + plot_layout(ncol = 2)
save_figure(p7, "Figure_7", f7, 8.2, 4.2)

# Figure 8: accepted Table 3 M2/M6 algorithm; 10th/90th percentiles only.
set.seed(20260511)
h22_frames <- load_analysis_frames(P$cache, "H22"); df_h22 <- h22_frames$df
trade_specs <- list(
  M2 = list(controls = c("nonviolent_camp", "REGCHANGE", "v2csreprss", "lnnum_image_sum", "lngdp"), estimator = "aft"),
  M6 = list(controls = c("nonviolent_camp", "REGCHANGE", "v2csreprss", "lnlengthofcam", "lnnum_image_sum", "lngdp"), estimator = "firth"))
trade_panels <- lapply(names(trade_specs), function(id) {
  s <- trade_specs[[id]]
  build_conditioned_panel(h22[[id]], df_h22, "fi_ftradeint_pd", "trade_c", s$controls,
                          s$estimator, id, "Trade Freedom", 80L)
})
f8 <- do.call(rbind, lapply(trade_panels, `[[`, "pred")); r8 <- do.call(rbind, lapply(trade_panels, `[[`, "rug"))
f8$panel <- factor(f8$model_id, levels = c("M2", "M6"), labels = c("Weibull AFT", "Firth PML"))
r8$panel <- factor(r8$model_id, levels = c("M2", "M6"), labels = c("Weibull AFT", "Firth PML"))
colours <- c("10th" = "#1F77B4", "90th" = "#D62728")
make_trade <- function(panel, ylabel) {
  pd <- f8[f8$panel == panel, ]; rd <- r8[r8$panel == panel, ]
  ggplot(pd, aes(eng, fit, colour = moderator_percentile, fill = moderator_percentile)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = .13, colour = NA) + geom_line(linewidth = .9) +
    geom_rug(data = rd, aes(x = eng), inherit.aes = FALSE, sides = "b", alpha = .25, length = unit(.012, "npc")) +
    scale_colour_manual(values = colours, name = "Trade Freedom percentile") +
    scale_fill_manual(values = colours, name = "Trade Freedom percentile") +
    labs(title = panel, x = "English banner proportion (campaign-level)", y = ylabel) + base_theme
}
p8 <- (make_trade("Weibull AFT", "Median time to success (years)") +
       make_trade("Firth PML", "Predicted P(success)")) + plot_layout(guides = "collect") +
  plot_annotation(title = "Trade Freedom X English Banner Proportions",
                  theme = theme(plot.title = element_text(hjust = .5, face = "bold", size = 13.5))) &
  theme(legend.position = "bottom")
save_figure(p8, "Figure_8", f8, 11, 4.5)

# Figure 9: accepted Table 4 M1/M5/M3/M7 algorithm; no median condition.
set.seed(20260511)
h3_frames <- load_analysis_frames(P$cache, "H3"); df_h3 <- h3_frames$df
h3_specs <- list(
  M1 = list(raw = "v2csreprss", centered = "repress_c", label = "CSO Repression", estimator = "aft",
            controls = c("nonviolent_camp", "REGCHANGE", "p_polity2", "lnnum_image_sum", "lngdp")),
  M5 = list(raw = "v2csreprss", centered = "repress_c", label = "CSO Repression", estimator = "firth",
            controls = c("nonviolent_camp", "REGCHANGE", "p_polity2", "lnlengthofcam", "lnnum_image_sum", "lngdp")),
  M3 = list(raw = "v2x_clpol", centered = "clpol_c", label = "Political Liberty", estimator = "aft",
            controls = c("nonviolent_camp", "REGCHANGE", "v2csreprss", "lnnum_image_sum", "lngdp")),
  M7 = list(raw = "v2x_clpol", centered = "clpol_c", label = "Political Liberty", estimator = "firth",
            controls = c("nonviolent_camp", "REGCHANGE", "v2csreprss", "lnlengthofcam", "lnnum_image_sum", "lngdp")))
h3_panels <- lapply(names(h3_specs), function(id) {
  s <- h3_specs[[id]]
  build_conditioned_panel(h3[[id]], df_h3, s$raw, s$centered, s$controls, s$estimator,
                          id, s$label, 80L)
})
f9 <- do.call(rbind, lapply(h3_panels, `[[`, "pred")); r9 <- do.call(rbind, lapply(h3_panels, `[[`, "rug"))
make_h3_panel <- function(id, title = NULL, ylabel, show_x = FALSE, y_limits = NULL) {
  pd <- f9[f9$model_id == id, ]; rd <- r9[r9$model_id == id, ]
  p <- ggplot(pd, aes(eng, fit, colour = moderator_percentile, fill = moderator_percentile)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = .13, colour = NA) + geom_line(linewidth = .9) +
    geom_rug(data = rd, aes(x = eng), inherit.aes = FALSE, sides = "b", alpha = .22, length = unit(.010, "npc")) +
    scale_colour_manual(values = colours, name = "Moderator percentile") +
    scale_fill_manual(values = colours, name = "Moderator percentile") +
    labs(title = title, x = if (show_x) "English banner proportion (campaign-level)" else NULL, y = ylabel) + base_theme
  if (!is.null(y_limits)) p <- p + coord_cartesian(ylim = y_limits)
  p
}
row_strip <- function(label) wrap_elements(panel = gTree(children = gList(
  rectGrob(gp = gpar(fill = "#F2F2F2", col = "#BDBDBD", lwd = .7)),
  textGrob(label, rot = -90, gp = gpar(fontsize = 11, col = "#222222")))))
aft_row <- make_h3_panel("M1", "CSO Repression", "Median time to success (years)", y_limits = c(0, 15)) +
           make_h3_panel("M3", "Political Liberty", "Median time to success (years)", y_limits = c(0, 15)) +
           row_strip("AFT model") + plot_layout(widths = c(1, 1, .06))
pml_row <- make_h3_panel("M5", ylabel = "Predicted P(success)", show_x = TRUE) +
           make_h3_panel("M7", ylabel = "Predicted P(success)", show_x = TRUE) +
           row_strip("PML model") + plot_layout(widths = c(1, 1, .06))
p9 <- (aft_row / pml_row) + plot_layout(guides = "collect") +
  plot_annotation(title = "Vulnerability to Domestic Pressure X English Banner Proportions",
                  theme = theme(plot.title = element_text(hjust = .5, face = "bold", size = 13.5))) &
  theme(legend.position = "bottom")
save_figure(p9, "Figure_9", f9, 10.5, 6.2)

# Deterministic raster comparison against the accepted PDFs.
accepted <- c(Figure_7 = "fig_marginal_h1.pdf", Figure_8 = "fig_h22_table3_cols3_4_trade_interaction.pdf",
              Figure_9 = "fig_h3_table4_cols1_2_5_6_interactions.pdf")
raster_mean_abs_threshold <- .05
pdftoppm <- Sys.which("pdftoppm"); if (!nzchar(pdftoppm)) stop("pdftoppm is required for raster audit")
raster_dir <- tempfile("figure-raster-"); dir.create(raster_dir); on.exit(unlink(raster_dir, recursive = TRUE), add = TRUE)
visual <- do.call(rbind, lapply(names(accepted), function(stem) {
  generated_pdf <- file.path(P$figures, paste0(stem, ".pdf")); accepted_pdf <- file.path(P$reference, "figures", accepted[[stem]])
  require_files(c(generated_pdf, accepted_pdf))
  gp <- file.path(raster_dir, paste0(stem, "-generated")); ap <- file.path(raster_dir, paste0(stem, "-accepted"))
  stopifnot(system2(pdftoppm, c("-f", "1", "-singlefile", "-r", "144", "-png",
                                shQuote(generated_pdf), shQuote(gp))) == 0L)
  stopifnot(system2(pdftoppm, c("-f", "1", "-singlefile", "-r", "144", "-png",
                                shQuote(accepted_pdf), shQuote(ap))) == 0L)
  a <- png::readPNG(paste0(gp, ".png")); b <- png::readPNG(paste0(ap, ".png")); same <- identical(dim(a), dim(b))
  mean_difference <- if (same) mean(abs(a - b)) else NA_real_
  data.frame(figure = stem, raster_dpi = 144L, generated_width = dim(a)[2], generated_height = dim(a)[1],
             accepted_width = dim(b)[2], accepted_height = dim(b)[1], same_dimensions = same,
             mean_absolute_pixel_difference = mean_difference,
             max_absolute_pixel_difference = if (same) max(abs(a - b)) else NA_real_,
             predeclared_mean_abs_threshold = raster_mean_abs_threshold,
             verdict = if (same && is.finite(mean_difference) &&
                           mean_difference <= raster_mean_abs_threshold) "PASS" else "FAIL")
}))
write.csv(visual, file.path(P$audit, "figure_visual_comparison.csv"), row.names = FALSE, na = "")
if (any(visual$verdict != "PASS")) {
  warning(sprintf("Advisory figure raster audit exceeded mean absolute pixel threshold %.3f; numerical plot-data validation remains blocking.",
                  raster_mean_abs_threshold), call. = FALSE)
}
cat("Wrote accepted-algorithm Figures 7-9, plot-data sidecars, and raster comparison audit.\n")
