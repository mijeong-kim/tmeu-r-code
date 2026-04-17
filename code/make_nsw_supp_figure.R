dat <- read.csv("data/nsw_mixtape.csv", stringsAsFactors = FALSE)
dat$re78_t <- asinh(dat$re78)
zero_share <- mean(dat$re78 == 0)
tbl <- read.csv("results/nsw_table2.csv", stringsAsFactors = FALSE)

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

pdf("figures/nsw_outcome_hist.pdf", width = 8.6, height = 4.2)
op <- par(mfrow = c(1, 2), mar = c(4.5, 4.6, 2.7, 1.2), oma = c(0, 0, 1.2, 0))

hist(
  dat$re78,
  breaks = "FD",
  col = "grey85",
  border = "white",
  main = "Raw RE78",
  xlab = "1978 earnings",
  ylab = "Frequency",
  bty = "n"
)
abline(v = 0, lty = 2, col = "grey35")
box(bty = "l")

hist(
  dat$re78_t,
  breaks = "FD",
  col = "grey85",
  border = "white",
  main = expression(asinh(RE78)),
  xlab = expression(asinh(RE78)),
  ylab = "",
  bty = "n"
)
abline(v = 0, lty = 2, col = "grey35")
box(bty = "l")

mtext(sprintf("Experimental NSW sample; %.1f%% of observed RE78 values are zero.", 100 * zero_share), outer = TRUE, cex = 0.9)
par(op)
dev.off()

tbl$label <- factor(tbl$estimator, levels = rev(tbl$estimator))

pdf("figures/nsw_ci_forest.pdf", width = 8.4, height = 4.6)
op <- par(mar = c(4.5, 8.5, 2.7, 1.2))
ylim <- c(0.5, nrow(tbl) + 0.5)
x_rng <- range(c(tbl$ci_lower, tbl$ci_upper))
plot(
  NA,
  xlim = x_rng,
  ylim = ylim,
  yaxt = "n",
  xlab = expression("ATE on " * asinh * "(RE78) scale"),
  ylab = "",
  main = "NSW empirical comparison",
  bty = "n"
)
axis(2, at = seq_len(nrow(tbl)), labels = rev(tbl$estimator), las = 1)
abline(v = 0, lty = 2, col = "grey45")
segments(tbl$ci_lower, seq_len(nrow(tbl)), tbl$ci_upper, seq_len(nrow(tbl)), lwd = 2, col = "grey35")
points(tbl$ate, seq_len(nrow(tbl)), pch = 19, cex = 1.1)
box(bty = "l")
par(op)
dev.off()
