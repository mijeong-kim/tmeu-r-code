args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(x) {
  out <- list(
    n = 500L,
    reps = 1000L,
    seed = 20260325L,
    split_reps = 20L,
    folds = 5L,
    maxit = 6L,
    out_prefix = "results/main_n500_rcf",
    scenarios = c("gaussian", "heavy_tailed", "skewed_mixture", "misspecified_mean")
  )
  if (length(x) == 0L) {
    return(out)
  }
  for (arg in x) {
    parts <- strsplit(arg, "=", fixed = TRUE)[[1L]]
    if (length(parts) != 2L) {
      next
    }
    key <- parts[1L]
    val <- parts[2L]
    if (key %in% c("n", "reps", "seed", "split_reps", "folds", "maxit")) {
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

suppressPackageStartupMessages({
  library(MASS)
})

dir.create(dirname(cfg$out_prefix), showWarnings = FALSE, recursive = TRUE)

true_ate <- 0.75
scenario_levels <- c("gaussian", "heavy_tailed", "skewed_mixture", "misspecified_mean")
cfg$scenarios <- intersect(cfg$scenarios, scenario_levels)
if (length(cfg$scenarios) == 0L) {
  stop("No valid scenarios requested.")
}
rhs_formula <- stats::delete.response(stats::terms(Y ~ A + W1 + W2 + W3 + W4 + A:W1 + A:W3))
fit_formula <- Y ~ A + W1 + W2 + W3 + W4 + A:W1 + A:W3

generate_errors <- function(n, scenario) {
  if (scenario == "gaussian") {
    return(rnorm(n))
  }
  if (scenario == "heavy_tailed") {
    return(rt(n, df = 3) / sqrt(3))
  }
  if (scenario %in% c("skewed_mixture", "misspecified_mean")) {
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
  p <- plogis(-0.2 + 0.5 * W1 - 0.4 * W2 + 0.6 * W3 - 0.3 * W4)
  A <- rbinom(n, 1L, p)
  eps <- generate_errors(n, scenario)
  mean_part <- 1 * A + 0.8 * W1 - 0.6 * W2 + 0.5 * W3 + 0.4 * W4 + 0.7 * A * W1 - 0.5 * A * W3
  if (scenario == "misspecified_mean") {
    mean_part <- mean_part + 0.4 * sin(W2)
  }
  Y <- mean_part + eps
  data.frame(Y = Y, A = A, W1 = W1, W2 = W2, W3 = W3, W4 = W4)
}

kernel_scores <- function(e, r, h) {
  u <- outer(as.numeric(e), as.numeric(r), "-") / h
  K <- dnorm(u)
  f <- rowMeans(K) / h
  f <- pmax(f, 1e-8)
  fp <- -rowMeans(u * K) / (h^2)
  fpp <- rowMeans((u^2 - 1) * K) / (h^3)
  l1 <- fp / f
  l2 <- fpp / f - l1^2
  list(l1 = l1, l2 = l2)
}

ate_gradient <- function(dat, coef_names) {
  g <- rep(0, length(coef_names))
  names(g) <- coef_names
  g["A"] <- 1
  g["A:W1"] <- mean(dat$W1)
  g["A:W3"] <- mean(dat$W3)
  g
}

delta_from_beta <- function(beta, dat) {
  dat1 <- dat
  dat0 <- dat
  dat1$A <- 1
  dat0$A <- 0
  x1 <- model.matrix(rhs_formula, data = dat1)
  x0 <- model.matrix(rhs_formula, data = dat0)
  as.numeric((x1 - x0) %*% beta)
}

fit_proposed_core <- function(dat, maxit = 6L) {
  fit0 <- lm(fit_formula, data = dat)
  X <- model.matrix(fit0)
  y <- dat$Y
  beta <- coef(fit0)
  p <- ncol(X)

  for (it in seq_len(maxit)) {
    r <- as.numeric(y - X %*% beta)
    h <- max(1.06 * stats::sd(r) * length(r)^(-1 / 5), 0.15)
    ks <- kernel_scores(r, r, h)
    U <- crossprod(X, ks$l1)
    J <- -crossprod(X * as.numeric(ks$l2), X) + diag(1e-6, p)
    step <- tryCatch(solve(J, U), error = function(e) MASS::ginv(J) %*% U)
    beta_new <- as.numeric(beta - as.numeric(step))
    names(beta_new) <- colnames(X)
    if (!all(is.finite(beta_new))) {
      break
    }
    if (max(abs(beta_new - beta)) < 1e-6) {
      beta <- beta_new
      break
    }
    beta <- beta_new
  }

  r <- as.numeric(y - X %*% beta)
  h <- max(1.06 * stats::sd(r) * length(r)^(-1 / 5), 0.15)
  ks <- kernel_scores(r, r, h)
  Ahat <- crossprod(X * as.numeric(ks$l2), X) / nrow(X)
  Ahat <- Ahat + diag(1e-6, ncol(Ahat))
  Ainv <- tryCatch(solve(Ahat), error = function(e) MASS::ginv(Ahat))

  list(
    beta = beta,
    Ainv = Ainv,
    train_resid = r,
    h = h,
    coef_names = colnames(X)
  )
}

score_l1_new <- function(e_new, ref_resid, h) {
  u <- outer(as.numeric(e_new), as.numeric(ref_resid), "-") / h
  K <- dnorm(u)
  f <- rowMeans(K) / h
  f <- pmax(f, 1e-8)
  fp <- -rowMeans(u * K) / (h^2)
  fp / f
}

fit_repeated_cf <- function(dat, split_reps = 20L, folds = 5L, maxit = 6L) {
  n <- nrow(dat)
  g_full <- ate_gradient(dat, colnames(model.matrix(fit_formula, data = dat)))
  split_est <- numeric(split_reps)
  split_var <- numeric(split_reps)

  for (b in seq_len(split_reps)) {
    fold_id <- sample(rep(seq_len(folds), length.out = n))
    delta_all <- numeric(n)
    infl_all <- matrix(NA_real_, nrow = n, ncol = length(g_full))
    colnames(infl_all) <- names(g_full)

    for (k in seq_len(folds)) {
      test_idx <- which(fold_id == k)
      train_idx <- which(fold_id != k)
      train <- dat[train_idx, , drop = FALSE]
      test <- dat[test_idx, , drop = FALSE]
      fit <- fit_proposed_core(train, maxit = maxit)

      delta_all[test_idx] <- delta_from_beta(fit$beta, test)

      X_test <- model.matrix(rhs_formula, data = test)
      e_test <- as.numeric(test$Y - X_test %*% fit$beta)
      l1_test <- score_l1_new(e_test, fit$train_resid, fit$h)
      infl_all[test_idx, ] <- (X_test * as.numeric(l1_test)) %*% t(fit$Ainv)
    }

    split_est[b] <- mean(delta_all)
    eif <- delta_all - split_est[b] + as.numeric(infl_all %*% g_full)
    split_var[b] <- mean(eif^2) / n
  }

  est <- mean(split_est)
  se <- sqrt(mean(split_var) + stats::var(split_est))
  c(
    estimate = est,
    se = se,
    ci_lower = est - 1.96 * se,
    ci_upper = est + 1.96 * se
  )
}

run_one <- function(scenario, rep_id, cfg) {
  dat <- generate_data(cfg$n, scenario)
  vals <- tryCatch(
    fit_repeated_cf(dat, split_reps = cfg$split_reps, folds = cfg$folds, maxit = cfg$maxit),
    error = function(e) c(estimate = NA_real_, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_)
  )
  data.frame(
    scenario = scenario,
    rep = rep_id,
    estimator = "Proposed TMLE",
    estimate = vals[["estimate"]],
    se = vals[["se"]],
    ci_lower = vals[["ci_lower"]],
    ci_upper = vals[["ci_upper"]],
    stringsAsFactors = FALSE
  )
}

summarize_results <- function(results, truth) {
  pieces <- split(results, results$scenario)
  out <- lapply(pieces, function(df) {
    est <- df$estimate
    width <- df$ci_upper - df$ci_lower
    data.frame(
      scenario = df$scenario[1L],
      estimator = "Proposed TMLE",
      bias = mean(est - truth, na.rm = TRUE),
      esd = stats::sd(est, na.rm = TRUE),
      rmse = sqrt(mean((est - truth)^2, na.rm = TRUE)),
      cover = mean(df$ci_lower <= truth & df$ci_upper >= truth, na.rm = TRUE),
      width = mean(width, na.rm = TRUE),
      failures = sum(!is.finite(est) | !is.finite(df$ci_lower) | !is.finite(df$ci_upper)),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  out$scenario <- factor(out$scenario, levels = scenario_levels)
  out[order(out$scenario), ]
}

set.seed(cfg$seed)
results_list <- vector("list", length = cfg$reps * length(cfg$scenarios))
counter <- 1L

for (scenario in cfg$scenarios) {
  for (rep_id in seq_len(cfg$reps)) {
    results_list[[counter]] <- run_one(scenario, rep_id, cfg)
    counter <- counter + 1L
    if (rep_id %% 25L == 0L) {
      message(sprintf("[%s] completed %d/%d", scenario, rep_id, cfg$reps))
    }
  }
}

results <- do.call(rbind, results_list)
summary_df <- summarize_results(results, truth = true_ate)

write.csv(results, paste0(cfg$out_prefix, "_replicates.csv"), row.names = FALSE)
write.csv(summary_df, paste0(cfg$out_prefix, "_summary.csv"), row.names = FALSE)

print(summary_df)
