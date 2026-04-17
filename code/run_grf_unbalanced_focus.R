args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(x) {
  out <- list(
    n = 500L,
    reps = 1000L,
    seed = 20260415L,
    grf_trees = 400L,
    out_prefix = "results/grf_unbalanced_focus_n500_mc1000",
    scenarios = c("heavy_tailed", "skewed_mixture")
  )
  if (length(x) == 0L) {
    return(out)
  }
  for (arg in x) {
    parts <- strsplit(arg, "=", fixed = TRUE)[[1L]]
    if (length(parts) != 2L) next
    key <- parts[1L]
    val <- parts[2L]
    if (key %in% c("n", "reps", "seed", "grf_trees")) {
      out[[key]] <- as.integer(val)
    } else if (key == "scenarios") {
      out[[key]] <- strsplit(val, ",", fixed = TRUE)[[1L]]
    } else if (key == "out_prefix") {
      out[[key]] <- val
    }
  }
  out
}

cfg <- parse_args(args)

if (!requireNamespace("grf", quietly = TRUE)) {
  stop("'grf' is required for this script.")
}

dir.create(dirname(cfg$out_prefix), showWarnings = FALSE, recursive = TRUE)

true_ate <- 0.75
scenario_levels <- c("heavy_tailed", "skewed_mixture")
cfg$scenarios <- intersect(cfg$scenarios, scenario_levels)
if (length(cfg$scenarios) == 0L) {
  stop("No valid scenarios requested.")
}

generate_errors <- function(n, scenario) {
  if (scenario == "heavy_tailed") {
    return(rt(n, df = 3) / sqrt(3))
  }
  if (scenario == "skewed_mixture") {
    mix <- rbinom(n, 1L, 0.2)
    err <- rnorm(n, mean = ifelse(mix == 1L, 2, -0.5), sd = ifelse(mix == 1L, 1, 0.5))
    return(err - mean(err))
  }
  stop("Unknown scenario: ", scenario)
}

generate_data <- function(n, scenario) {
  W1 <- rnorm(n)
  W2 <- rnorm(n)
  W3 <- rbinom(n, 1L, 0.5)
  W4 <- runif(n, min = -1, max = 1)
  p_raw <- plogis(-1.0 + 0.9 * W1 - 0.8 * W2 + 0.8 * W3 - 0.6 * W4)
  p <- pmin(pmax(p_raw, 0.08), 0.92)
  A <- rbinom(n, 1L, p)
  eps <- generate_errors(n, scenario)
  mean_part <- 1 * A + 0.8 * W1 - 0.6 * W2 + 0.5 * W3 + 0.4 * W4 + 0.7 * A * W1 - 0.5 * A * W3
  Y <- mean_part + eps
  data.frame(Y = Y, A = A, W1 = W1, W2 = W2, W3 = W3, W4 = W4)
}

fit_grf_tmle <- function(dat, num.trees) {
  x <- as.matrix(dat[c("W1", "W2", "W3", "W4")])
  fit <- grf::causal_forest(
    X = x,
    Y = dat$Y,
    W = dat$A,
    num.trees = num.trees,
    num.threads = 1L,
    seed = sample.int(.Machine$integer.max, 1L)
  )
  ate <- grf::average_treatment_effect(fit, target.sample = "all", method = "TMLE")
  est <- unname(ate[["estimate"]])
  se <- unname(ate[["std.err"]])
  c(
    estimate = est,
    se = se,
    ci_lower = est - 1.96 * se,
    ci_upper = est + 1.96 * se
  )
}

summarize_results <- function(results, truth) {
  split_key <- interaction(results$scenario, results$estimator, drop = TRUE)
  pieces <- split(results, split_key)
  do.call(rbind, lapply(pieces, function(df) {
    est <- df$estimate
    width <- df$ci_upper - df$ci_lower
    data.frame(
      scenario = df$scenario[1L],
      estimator = df$estimator[1L],
      bias = mean(est - truth, na.rm = TRUE),
      esd = stats::sd(est, na.rm = TRUE),
      rmse = sqrt(mean((est - truth)^2, na.rm = TRUE)),
      cover = mean(df$ci_lower <= truth & df$ci_upper >= truth, na.rm = TRUE),
      width = mean(width, na.rm = TRUE),
      failures = sum(!is.finite(est)),
      stringsAsFactors = FALSE
    )
  }))
}

set.seed(cfg$seed)

rows <- vector("list", length(cfg$scenarios) * cfg$reps)
counter <- 1L
for (scenario in cfg$scenarios) {
  for (rep_id in seq_len(cfg$reps)) {
    dat <- generate_data(cfg$n, scenario)
    vals <- tryCatch(
      fit_grf_tmle(dat, cfg$grf_trees),
      error = function(e) c(estimate = NA_real_, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_)
    )
    rows[[counter]] <- data.frame(
      scenario = scenario,
      rep = rep_id,
      estimator = "GRF TMLE",
      estimate = vals[["estimate"]],
      se = vals[["se"]],
      ci_lower = vals[["ci_lower"]],
      ci_upper = vals[["ci_upper"]],
      stringsAsFactors = FALSE
    )
    counter <- counter + 1L
    if (rep_id %% 10L == 0L) {
      message(sprintf("[%s] completed %d/%d", scenario, rep_id, cfg$reps))
    }
  }
}

results <- do.call(rbind, rows)
summary_df <- summarize_results(results, true_ate)
summary_df$scenario <- factor(summary_df$scenario, levels = scenario_levels)
summary_df <- summary_df[order(summary_df$scenario), ]

write.csv(results, paste0(cfg$out_prefix, "_replicates.csv"), row.names = FALSE)
write.csv(summary_df, paste0(cfg$out_prefix, "_summary.csv"), row.names = FALSE)
print(summary_df)
