main_prefix <- Sys.getenv("MAIN_SIM_PREFIX", unset = "results/main_n500_rcf")
unbalanced_prefix <- Sys.getenv("IMBALANCED_SIM_PREFIX", unset = "results/unbalanced_aipw_n500_mc1000")
grf_unbalanced_prefix <- Sys.getenv("GRF_IMBALANCED_SIM_PREFIX", unset = "results/grf_unbalanced_focus_n500_mc1000")

sim_summary <- read.csv(paste0(main_prefix, "_summary.csv"), stringsAsFactors = FALSE)
sim_reps <- read.csv(paste0(main_prefix, "_replicates.csv"), stringsAsFactors = FALSE)
unbalanced_summary <- read.csv(paste0(unbalanced_prefix, "_summary.csv"), stringsAsFactors = FALSE)
unbalanced_reps <- read.csv(paste0(unbalanced_prefix, "_replicates.csv"), stringsAsFactors = FALSE)
grf_unbalanced_reps <- read.csv(paste0(grf_unbalanced_prefix, "_replicates.csv"), stringsAsFactors = FALSE)
nsw_summary <- read.csv("results/nsw_table2.csv", stringsAsFactors = FALSE)

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

scenario_map <- c(
  gaussian = "Gaussian",
  heavy_tailed = "Heavy-tailed",
  skewed_mixture = "Skewed mixture",
  misspecified_mean = "Misspecified mean"
)
estimator_levels <- c("Proposed TMLE", "Gaussian OLS", "BART")
estimator_axis_labels <- c("Proposed", "OLS", "BART")
estimator_pch <- c(
  "Proposed TMLE" = 16,
  "Gaussian OLS" = 1,
  "BART" = 17
)
estimator_lty <- c(
  "Proposed TMLE" = 1,
  "Gaussian OLS" = 2,
  "BART" = 3
)
true_ate <- 0.75
unbalanced_estimator_levels <- c("Proposed TMLE", "Gaussian OLS", "AIPW", "BART", "Forest-based TMLE")
unbalanced_axis_labels <- c("Proposed", "OLS", "AIPW", "BART", "Forest TMLE")
unbalanced_pch <- c(
  "Proposed TMLE" = 16,
  "Gaussian OLS" = 1,
  "AIPW" = 2,
  "BART" = 17,
  "Forest-based TMLE" = 15
)
unbalanced_lty <- c(
  "Proposed TMLE" = 1,
  "Gaussian OLS" = 2,
  "AIPW" = 3,
  "BART" = 4,
  "Forest-based TMLE" = 5
)

sim_summary$scenario_label <- unname(scenario_map[sim_summary$scenario])
sim_summary$scenario_label <- factor(
  sim_summary$scenario_label,
  levels = unname(scenario_map)
)
sim_summary$estimator <- factor(sim_summary$estimator, levels = estimator_levels)

sim_reps$scenario_label <- unname(scenario_map[sim_reps$scenario])
sim_reps$scenario_label <- factor(
  sim_reps$scenario_label,
  levels = unname(scenario_map)
)
sim_reps$estimator <- factor(sim_reps$estimator, levels = estimator_levels)

unbalanced_summary$scenario_label <- unname(scenario_map[unbalanced_summary$scenario])
unbalanced_summary$scenario_label <- factor(
  unbalanced_summary$scenario_label,
  levels = unname(scenario_map)
)
unbalanced_summary$estimator <- factor(
  unbalanced_summary$estimator,
  levels = unbalanced_estimator_levels
)

unbalanced_reps$scenario_label <- unname(scenario_map[unbalanced_reps$scenario])
unbalanced_reps$scenario_label <- factor(
  unbalanced_reps$scenario_label,
  levels = unname(scenario_map)
)
unbalanced_reps$estimator <- factor(
  unbalanced_reps$estimator,
  levels = unbalanced_estimator_levels
)

grf_unbalanced_reps$scenario_label <- unname(scenario_map[grf_unbalanced_reps$scenario])
grf_unbalanced_reps$scenario_label <- factor(
  grf_unbalanced_reps$scenario_label,
  levels = unname(scenario_map)
)
grf_unbalanced_reps$estimator <- "Forest-based TMLE"
grf_unbalanced_reps$estimator <- factor(
  grf_unbalanced_reps$estimator,
  levels = unbalanced_estimator_levels
)
grf_unbalanced_reps$treat_rate <- NA_real_

unbalanced_reps <- rbind(
  unbalanced_reps,
  grf_unbalanced_reps[, names(unbalanced_reps)]
)

sim_panel_summary <- aggregate(
  cbind(estimate, ci_lower, ci_upper) ~ scenario + scenario_label + estimator,
  data = sim_reps,
  FUN = function(x) mean(x, na.rm = TRUE)
)

unbalanced_panel_summary <- aggregate(
  cbind(estimate, ci_lower, ci_upper) ~ scenario + scenario_label + estimator,
  data = unbalanced_reps,
  FUN = function(x) mean(x, na.rm = TRUE)
)

sample_remainder <- function(n, scenario) {
  if (scenario == "gaussian") {
    return(rnorm(n))
  }
  if (scenario == "heavy_tailed") {
    return(rt(n, df = 3) / sqrt(3))
  }
  mix <- rbinom(n, 1L, 0.2)
  eps <- rnorm(n, mean = ifelse(mix == 1L, 2, -0.5), sd = ifelse(mix == 1L, 1, 0.5))
  if (scenario == "skewed_mixture") {
    return(eps)
  }
  if (scenario == "misspecified_mean") {
    w2 <- rnorm(n)
    return(0.4 * sin(w2) + eps)
  }
  stop("Unknown scenario: ", scenario)
}

hist_density_ymax <- function(scenarios, n = 4000L, breaks = 24L) {
  ymax <- 0
  for (scenario in scenarios) {
    vals <- sample_remainder(n, scenario)
    h <- hist(vals, breaks = breaks, plot = FALSE)
    ymax <- max(ymax, h$density, na.rm = TRUE)
  }
  1.05 * ymax
}

draw_hist_panel <- function(scenario_key, panel_title, show_y = TRUE) {
  vals <- sample_remainder(4000L, scenario_key)
  hist(
    vals,
    breaks = 24,
    freq = FALSE,
    col = "grey85",
    border = "white",
    main = panel_title,
    xlab = "",
    ylab = if (show_y) "Density" else "",
    yaxt = if (show_y) "s" else "n",
    ylim = c(0, sim_hist_ymax),
    bty = "n",
    cex.main = 1.05,
    cex.lab = 0.95,
    cex.axis = 0.95
  )
  box(bty = "l")
}

draw_ci_panel <- function(scenario_key, show_y = TRUE, show_x = TRUE) {
  dat <- sim_panel_summary[sim_panel_summary$scenario == scenario_key, , drop = FALSE]
  dat <- dat[match(estimator_levels, dat$estimator), , drop = FALSE]
  ypos <- c(3, 2, 1)

  plot(
    NA,
    NA,
    xlim = ci_x_range,
    ylim = c(0.5, 3.5),
    yaxt = "n",
    ylab = "",
    xlab = if (show_x) "ATE" else "",
    xaxt = if (show_x) "s" else "n",
    bty = "n",
    cex.lab = 1.0,
    cex.axis = 0.95
  )
  abline(v = true_ate, lty = 2, col = "grey35")
  abline(h = ypos, col = "grey92", lty = 1)
  if (show_y) {
    axis(2, at = ypos, labels = estimator_axis_labels, las = 1, cex.axis = 0.92)
  } else {
    axis(2, at = ypos, labels = FALSE, tck = 0)
  }
  for (i in seq_along(estimator_levels)) {
    est <- estimator_levels[i]
    row <- dat[dat$estimator == est, , drop = FALSE]
    segments(
      row$ci_lower,
      ypos[i],
      row$ci_upper,
      ypos[i],
      lwd = 1.8,
      lty = estimator_lty[[est]],
      col = "black"
    )
    points(
      row$estimate,
      ypos[i],
      pch = estimator_pch[[est]],
      cex = 1.15,
      col = "black",
      bg = "white"
    )
  }
}

draw_ci_panel_unbalanced <- function(scenario_key, show_y = TRUE, show_x = TRUE, show_right_labels = FALSE) {
  dat <- unbalanced_panel_summary[unbalanced_panel_summary$scenario == scenario_key, , drop = FALSE]
  present_levels <- unbalanced_estimator_levels[unbalanced_estimator_levels %in% as.character(dat$estimator)]
  dat <- dat[match(present_levels, dat$estimator), , drop = FALSE]
  ypos <- rev(seq_len(nrow(dat)))
  axis_labels <- unbalanced_axis_labels[match(present_levels, unbalanced_estimator_levels)]

  plot(
    NA,
    NA,
    xlim = unbalanced_ci_x_range,
    ylim = c(0.5, length(ypos) + 0.5),
    yaxt = "n",
    ylab = "",
    xlab = if (show_x) "ATE" else "",
    xaxt = if (show_x) "s" else "n",
    bty = "n",
    cex.lab = 1.0,
    cex.axis = 0.95
  )
  abline(v = true_ate, lty = 2, col = "grey35")
  abline(h = ypos, col = "grey92", lty = 1)
  if (show_y) {
    axis(2, at = ypos, labels = axis_labels, las = 1, cex.axis = 0.92)
  } else if (show_right_labels) {
    axis(4, at = ypos, labels = axis_labels, las = 1, cex.axis = 0.92)
  } else {
    axis(2, at = ypos, labels = FALSE, tck = 0)
  }
  for (i in seq_along(present_levels)) {
    est <- present_levels[i]
    row <- dat[dat$estimator == est, , drop = FALSE]
    segments(
      row$ci_lower,
      ypos[i],
      row$ci_upper,
      ypos[i],
      lwd = 1.8,
      lty = unbalanced_lty[[est]],
      col = "black"
    )
    points(
      row$estimate,
      ypos[i],
      pch = unbalanced_pch[[est]],
      cex = 1.15,
      col = "black",
      bg = "white"
    )
  }
}

ci_x_range <- range(c(sim_panel_summary$ci_lower, sim_panel_summary$ci_upper, true_ate), finite = TRUE)
ci_pad <- 0.06 * diff(ci_x_range)
if (!is.finite(ci_pad) || ci_pad == 0) {
  ci_pad <- 0.1
}
ci_x_range <- c(ci_x_range[1] - ci_pad, ci_x_range[2] + ci_pad)

unbalanced_ci_x_range <- range(
  c(unbalanced_panel_summary$ci_lower, unbalanced_panel_summary$ci_upper, true_ate),
  finite = TRUE
)
unbalanced_ci_pad <- 0.06 * diff(unbalanced_ci_x_range)
if (!is.finite(unbalanced_ci_pad) || unbalanced_ci_pad == 0) {
  unbalanced_ci_pad <- 0.1
}
unbalanced_ci_x_range <- c(
  unbalanced_ci_x_range[1] - unbalanced_ci_pad,
  unbalanced_ci_x_range[2] + unbalanced_ci_pad
)

set.seed(20260322L)
sim_hist_ymax <- hist_density_ymax(names(scenario_map))
pdf("figures/sim_main_compare.pdf", width = 9.5, height = 10.2)
layout(
  matrix(1:8, nrow = 4, byrow = TRUE),
  heights = c(2.55, 1.4, 2.55, 1.4)
)

par(mar = c(2.2, 4.2, 2.2, 1.1))
draw_hist_panel("gaussian", "Gaussian", show_y = TRUE)
par(mar = c(2.2, 4.2, 2.2, 1.1))
draw_hist_panel("heavy_tailed", "Heavy-tailed", show_y = TRUE)
par(mar = c(3.1, 7.2, 1.1, 1.1))
draw_ci_panel("gaussian", show_y = TRUE, show_x = FALSE)
par(mar = c(3.1, 2.2, 1.1, 1.1))
draw_ci_panel("heavy_tailed", show_y = FALSE, show_x = FALSE)
par(mar = c(2.2, 4.2, 2.2, 1.1))
draw_hist_panel("skewed_mixture", "Skewed mixture", show_y = TRUE)
par(mar = c(2.2, 4.2, 2.2, 1.1))
draw_hist_panel("misspecified_mean", "Misspecified mean", show_y = TRUE)
par(mar = c(3.1, 7.2, 1.1, 1.1))
draw_ci_panel("skewed_mixture", show_y = TRUE, show_x = TRUE)
par(mar = c(3.1, 2.2, 1.1, 1.1))
draw_ci_panel("misspecified_mean", show_y = FALSE, show_x = TRUE)
dev.off()

set.seed(20260322L)
pdf("figures/sim_unbalanced_compare.pdf", width = 9.5, height = 10.8)
layout(
  matrix(1:8, nrow = 4, byrow = TRUE),
  heights = c(2.55, 1.7, 2.55, 1.7)
)

par(mar = c(2.2, 4.2, 2.2, 1.1))
draw_hist_panel("gaussian", "Gaussian", show_y = TRUE)
par(mar = c(2.2, 4.2, 2.2, 1.1))
draw_hist_panel("heavy_tailed", "Heavy-tailed", show_y = TRUE)
par(mar = c(3.5, 8.0, 1.1, 1.1))
draw_ci_panel_unbalanced("gaussian", show_y = TRUE, show_x = FALSE)
par(mar = c(3.5, 2.2, 1.1, 8.0))
draw_ci_panel_unbalanced("heavy_tailed", show_y = FALSE, show_x = FALSE, show_right_labels = TRUE)
par(mar = c(2.2, 4.2, 2.2, 1.1))
draw_hist_panel("skewed_mixture", "Skewed mixture", show_y = TRUE)
par(mar = c(2.2, 4.2, 2.2, 1.1))
draw_hist_panel("misspecified_mean", "Misspecified mean", show_y = TRUE)
par(mar = c(3.5, 8.0, 1.1, 1.1))
draw_ci_panel_unbalanced("skewed_mixture", show_y = TRUE, show_x = TRUE)
par(mar = c(3.5, 2.2, 1.1, 8.0))
draw_ci_panel_unbalanced("misspecified_mean", show_y = FALSE, show_x = TRUE, show_right_labels = TRUE)
dev.off()

nsw_display_order <- c("Proposed TMLE", "AIPW", "BART plug-in", "Forest-based TMLE")
nsw_summary <- nsw_summary[
  match(nsw_display_order, nsw_summary$estimator),
  ,
  drop = FALSE
]
nsw_summary$estimator <- factor(
  nsw_summary$estimator,
  levels = nsw_display_order
)
nsw_ci_lty <- c(
  "Proposed TMLE" = 1,
  "AIPW" = 2,
  "BART plug-in" = 3,
  "Forest-based TMLE" = 4
)
nsw_ci_pch <- c(
  "Proposed TMLE" = 16,
  "AIPW" = 1,
  "BART plug-in" = 17,
  "Forest-based TMLE" = 15
)

pdf("figures/nsw_ci_forest.pdf", width = 10.2, height = 5.2)
op <- par(mar = c(4.5, 7.8, 2.5, 6.2))
ypos <- rev(seq_len(nrow(nsw_summary)))
xlim <- range(c(nsw_summary$ci_lower, nsw_summary$ci_upper))
xpad <- 0.55 * diff(xlim)
if (!is.finite(xpad) || xpad == 0) {
  xpad <- 0.1
}
text_x <- xlim[2] + 0.14 * diff(xlim)

plot(
  NA,
  NA,
  xlim = c(xlim[1] - xpad, xlim[2] + xpad),
  ylim = c(0.5, nrow(nsw_summary) + 0.5),
  yaxt = "n",
  ylab = "",
  xlab = "ATE on asinh(RE78) scale",
  main = "NSW Empirical Comparison",
  bty = "n",
  cex.lab = 1.0,
  cex.axis = 0.95,
  cex.main = 1.05
)
axis(2, at = ypos, labels = as.character(nsw_summary$estimator), las = 1, cex.axis = 0.95)
abline(v = 0, col = "grey80", lty = 2)
abline(h = ypos, col = "grey92", lty = 1)

for (i in seq_len(nrow(nsw_summary))) {
  est_name <- as.character(nsw_summary$estimator[i])
  ci_lty <- nsw_ci_lty[[est_name]]
  ci_pch <- nsw_ci_pch[[est_name]]
  segments(
    nsw_summary$ci_lower[i],
    ypos[i],
    nsw_summary$ci_upper[i],
    ypos[i],
    lwd = 2.5,
    lty = ci_lty,
    col = "black"
  )
  points(nsw_summary$ate[i], ypos[i], pch = ci_pch, cex = 1.3, col = "black", bg = "white")
  text(
    x = text_x,
    y = ypos[i],
    labels = sprintf("w = %.3f", nsw_summary$width[i]),
    adj = c(0, 0.5),
    cex = 0.9,
    col = "grey30",
    xpd = NA
  )
}
par(op)
dev.off()
