args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(x) {
  out <- list(
    n = 500L,
    reps = 1000L,
    seed = 20260325L,
    split_reps = 20L,
    folds = 5L,
    bart_samples = 60L,
    bart_burn = 20L,
    bart_chains = 1L,
    grf_trees = 800L,
    maxit = 6L,
    include_bart = TRUE,
    include_grf = TRUE,
    out_prefix = "results/unbalanced_aipw_n500_mc1000",
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
    if (key %in% c("n", "reps", "seed", "split_reps", "folds", "bart_samples", "bart_burn", "bart_chains", "grf_trees", "maxit")) {
      out[[key]] <- as.integer(val)
    } else if (key %in% c("include_bart", "include_grf")) {
      out[[key]] <- tolower(val) %in% c("1", "true", "t", "yes", "y")
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

grf_available <- FALSE
if (isTRUE(cfg$include_grf)) {
  grf_available <- requireNamespace("grf", quietly = TRUE)
  if (!grf_available) {
    stop("'grf' is not installed but include_grf=TRUE was requested.")
  }
}

bart_backend <- NULL
if (isTRUE(cfg$include_bart)) {
  bart_backend <- if (requireNamespace("bartCause", quietly = TRUE)) {
    "bartCause"
  } else if (requireNamespace("bartMachine", quietly = TRUE)) {
    bartMachine::set_bart_machine_num_cores(1L)
    "bartMachine"
  } else {
    stop("Neither 'bartCause' nor 'bartMachine' is installed.")
  }
}

if (is.null(bart_backend)) {
  message("BART benchmark skipped for this run.")
} else {
  message(sprintf("Using %s for the BART benchmark.", bart_backend))
}
if (grf_available) {
  message(sprintf("Using grf with %d trees for the forest-based TMLE benchmark.", cfg$grf_trees))
} else {
  message("GRF benchmark skipped for this run.")
}
dir.create(dirname(cfg$out_prefix), showWarnings = FALSE, recursive = TRUE)

true_ate <- 0.75
scenario_levels <- c("gaussian", "heavy_tailed", "skewed_mixture", "misspecified_mean")
cfg$scenarios <- intersect(cfg$scenarios, scenario_levels)
if (length(cfg$scenarios) == 0L) {
  stop("No valid scenarios requested.")
}
fit_formula <- Y ~ A + W1 + W2 + W3 + W4 + A:W1 + A:W3
rhs_formula <- stats::delete.response(stats::terms(fit_formula))
wm_formula <- Y ~ W1 + W2 + W3 + W4
ps_formula <- A ~ W1 + W2 + W3 + W4

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
  p_raw <- plogis(-1.0 + 0.9 * W1 - 0.8 * W2 + 0.8 * W3 - 0.6 * W4)
  p <- pmin(pmax(p_raw, 0.08), 0.92)
  A <- rbinom(n, 1L, p)
  eps <- generate_errors(n, scenario)
  mean_part <- 1 * A + 0.8 * W1 - 0.6 * W2 + 0.5 * W3 + 0.4 * W4 + 0.7 * A * W1 - 0.5 * A * W3
  if (scenario == "misspecified_mean") {
    mean_part <- mean_part + 0.4 * sin(W2)
  }
  Y <- mean_part + eps
  data.frame(Y = Y, A = A, W1 = W1, W2 = W2, W3 = W3, W4 = W4, p_true = p)
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

fit_repeated_cf_proposed <- function(dat, split_reps = 10L, folds = 5L, maxit = 6L) {
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

fit_aipw_repeated_cf <- function(dat, split_reps = 10L, folds = 5L) {
  n <- nrow(dat)
  split_est <- numeric(split_reps)
  split_var <- numeric(split_reps)

  for (b in seq_len(split_reps)) {
    fold_id <- sample(rep(seq_len(folds), length.out = n))
    pseudo <- numeric(n)

    for (k in seq_len(folds)) {
      test_idx <- which(fold_id == k)
      train_idx <- which(fold_id != k)
      train <- dat[train_idx, , drop = FALSE]
      test <- dat[test_idx, , drop = FALSE]

      if (sum(train$A == 1) < 10L || sum(train$A == 0) < 10L) {
        return(c(estimate = NA_real_, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_))
      }

      ps_fit <- glm(ps_formula, data = train, family = binomial())
      e_hat <- predict(ps_fit, newdata = test, type = "response")
      e_hat <- pmin(pmax(e_hat, 0.02), 0.98)

      mu1_fit <- lm(wm_formula, data = train[train$A == 1, , drop = FALSE])
      mu0_fit <- lm(wm_formula, data = train[train$A == 0, , drop = FALSE])
      mu1_hat <- predict(mu1_fit, newdata = test)
      mu0_hat <- predict(mu0_fit, newdata = test)

      pseudo[test_idx] <- mu1_hat - mu0_hat +
        test$A * (test$Y - mu1_hat) / e_hat -
        (1 - test$A) * (test$Y - mu0_hat) / (1 - e_hat)
    }

    split_est[b] <- mean(pseudo)
    split_var[b] <- stats::var(pseudo) / n
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

fit_ols <- function(dat) {
  fit <- lm(fit_formula, data = dat)
  beta <- coef(fit)
  g <- ate_gradient(dat, names(beta))
  ate <- as.numeric(sum(g * beta))
  se <- sqrt(drop(t(g) %*% vcov(fit) %*% g))
  c(
    estimate = ate,
    se = se,
    ci_lower = ate - 1.96 * se,
    ci_upper = ate + 1.96 * se
  )
}

fit_bart <- function(dat, bart_samples, bart_burn, bart_chains) {
  if (bart_backend == "bartCause") {
    x <- dat[c("W1", "W2", "W3", "W4")]
    fit <- bartCause::bartc(
      response = dat$Y,
      treatment = dat$A,
      confounders = x,
      method.rsp = "bart",
      method.trt = "glm",
      estimand = "ate",
      verbose = FALSE,
      n.samples = bart_samples,
      n.burn = bart_burn,
      n.chains = bart_chains,
      n.threads = 1L,
      seed = sample.int(.Machine$integer.max, 1L)
    )
    est <- summary(fit)$estimates[1L, ]
    return(
      c(
        estimate = as.numeric(est[["estimate"]]),
        se = as.numeric(est[["sd"]]),
        ci_lower = as.numeric(est[["ci.lower"]]),
        ci_upper = as.numeric(est[["ci.upper"]])
      )
    )
  }

  x <- dat[c("A", "W1", "W2", "W3", "W4")]
  fit <- bartMachine::bartMachine(
    X = x,
    y = dat$Y,
    num_trees = 50L,
    num_burn_in = bart_burn,
    num_iterations_after_burn_in = bart_samples,
    seed = sample.int(.Machine$integer.max, 1L),
    verbose = FALSE,
    mem_cache_for_speed = FALSE,
    flush_indices_to_save_RAM = TRUE
  )
  x1 <- x
  x0 <- x
  x1$A <- 1
  x0$A <- 0
  post1 <- bartMachine::bart_machine_get_posterior(fit, x1)$y_hat_posterior_samples
  post0 <- bartMachine::bart_machine_get_posterior(fit, x0)$y_hat_posterior_samples
  ate_draws <- colMeans(post1 - post0)
  ci <- stats::quantile(ate_draws, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  out <- c(
    estimate = mean(ate_draws, na.rm = TRUE),
    se = stats::sd(ate_draws, na.rm = TRUE),
    ci_lower = ci[1L],
    ci_upper = ci[2L]
  )
  rm(fit, post1, post0, ate_draws)
  gc()
  out
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

run_one <- function(scenario, rep_id, cfg) {
  dat <- generate_data(cfg$n, scenario)
  out <- list()

  out[["Proposed TMLE"]] <- tryCatch(
    fit_repeated_cf_proposed(dat, split_reps = cfg$split_reps, folds = cfg$folds, maxit = cfg$maxit),
    error = function(e) c(estimate = NA_real_, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_)
  )
  out[["AIPW"]] <- tryCatch(
    fit_aipw_repeated_cf(dat, split_reps = cfg$split_reps, folds = cfg$folds),
    error = function(e) c(estimate = NA_real_, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_)
  )
  out[["Gaussian OLS"]] <- tryCatch(
    fit_ols(dat),
    error = function(e) c(estimate = NA_real_, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_)
  )
  if (grf_available) {
    out[["GRF TMLE"]] <- tryCatch(
      fit_grf_tmle(dat, cfg$grf_trees),
      error = function(e) c(estimate = NA_real_, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_)
    )
  }
  if (!is.null(bart_backend)) {
    out[["BART"]] <- tryCatch(
      fit_bart(dat, cfg$bart_samples, cfg$bart_burn, cfg$bart_chains),
      error = function(e) c(estimate = NA_real_, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_)
    )
  }

  rows <- lapply(names(out), function(estimator) {
    vals <- out[[estimator]]
    data.frame(
      scenario = scenario,
      rep = rep_id,
      estimator = estimator,
      estimate = vals[["estimate"]],
      se = vals[["se"]],
      ci_lower = vals[["ci_lower"]],
      ci_upper = vals[["ci_upper"]],
      treat_rate = mean(dat$A),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

summarize_results <- function(results, truth) {
  estimator_levels <- c("Proposed TMLE", "Gaussian OLS", "AIPW")
  if (grf_available) {
    estimator_levels <- c(estimator_levels, "GRF TMLE")
  }
  if (!is.null(bart_backend)) {
    estimator_levels <- c(estimator_levels, "BART")
  }
  split_key <- interaction(results$scenario, results$estimator, drop = TRUE)
  pieces <- split(results, split_key)
  out <- lapply(pieces, function(df) {
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
      treat_rate = mean(df$treat_rate, na.rm = TRUE),
      failures = sum(!is.finite(est)),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  out$scenario <- factor(out$scenario, levels = scenario_levels)
  out$estimator <- factor(out$estimator, levels = estimator_levels)
  out[order(out$scenario, out$estimator), ]
}

format_num <- function(x, digits = 3) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
}

bold_cell <- function(x, flag, digits = 3) {
  txt <- format_num(x, digits = digits)
  ifelse(flag, paste0("\\textbf{", txt, "}"), txt)
}

best_mask <- function(x, score = identity) {
  s <- score(x)
  ok <- is.finite(s)
  out <- rep(FALSE, length(x))
  if (any(ok)) {
    out[ok] <- s[ok] == min(s[ok], na.rm = TRUE)
  }
  out
}

write_latex_table <- function(summary_df, path, n, reps) {
  scenario_labels <- c(
    gaussian = "Gaussian",
    heavy_tailed = "Heavy-tailed",
    skewed_mixture = "Skewed mixture",
    misspecified_mean = "Misspecified mean"
  )
  note_parts <- c("The proposed estimator and AIPW are reported with repeated cross-fitting.")
  if (any(summary_df$estimator == "GRF TMLE", na.rm = TRUE)) {
    note_parts <- c(note_parts, "GRF TMLE uses the forest-based targeted estimator from \\texttt{grf::average\\_treatment\\_effect(..., method = \"TMLE\")}.")
  }
  if (any(summary_df$estimator == "BART", na.rm = TRUE)) {
    note_parts <- c(note_parts, "Gaussian OLS and BART use their standard interval constructions.")
  } else {
    note_parts <- c(note_parts, "Gaussian OLS uses its standard interval construction.")
  }
  estimator_note <- paste(note_parts, collapse = " ")
  lines <- c(
    sprintf("%% Generated on %s", Sys.time()),
    "\\begin{table}[t]",
    "\\centering",
    sprintf("\\caption{Simulation results for $n=%d$ under unbalanced treatment assignment over %d Monte Carlo replications.}", n, reps),
    "\\label{tab:sim_unbalanced_aipw}",
    "\\small",
    "\\begin{tabular}{llccccc}",
    "\\toprule",
    "Scenario & Estimator & Bias & ESD & RMSE & Cover. & Width \\\\",
    "\\midrule"
  )

  for (sc in scenario_levels) {
    block <- summary_df[summary_df$scenario == sc, ]
    best_bias <- best_mask(block$bias, score = abs)
    best_esd <- best_mask(block$esd)
    best_rmse <- best_mask(block$rmse)
    best_cover <- best_mask(block$cover, score = function(x) abs(x - 0.95))
    best_width <- best_mask(block$width)
    for (i in seq_len(nrow(block))) {
      row <- block[i, ]
      lines <- c(
        lines,
        sprintf(
          "%s & %s & %s & %s & %s & %s & %s \\\\",
          scenario_labels[[sc]],
          as.character(row$estimator),
          bold_cell(row$bias, best_bias[i]),
          bold_cell(row$esd, best_esd[i]),
          bold_cell(row$rmse, best_rmse[i]),
          bold_cell(row$cover, best_cover[i]),
          bold_cell(row$width, best_width[i])
        )
      )
    }
    if (sc != tail(scenario_levels, 1L)) {
      lines <- c(lines, "\\addlinespace")
    }
  }

  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\par\\medskip",
    sprintf("\\footnotesize Treatment assignment is generated from a covariate-dependent propensity score with average treated fraction %.3f across replications. ESD, empirical standard deviation; Cover., empirical coverage of the nominal $95\\%%$ confidence interval; Width, average confidence interval width. %s Within each scenario, boldface marks the best displayed value in each column, using smallest absolute bias, smallest ESD, smallest RMSE, coverage closest to $0.95$, and smallest interval width; ties are both boldfaced.", mean(unique(summary_df$treat_rate), na.rm = TRUE), estimator_note),
    "\\end{table}"
  )
  writeLines(lines, path)
}

set.seed(cfg$seed)

results_list <- vector("list", length = length(cfg$scenarios) * cfg$reps)
counter <- 1L
for (scenario in cfg$scenarios) {
  for (rep_id in seq_len(cfg$reps)) {
    results_list[[counter]] <- run_one(scenario, rep_id, cfg)
    counter <- counter + 1L
    if (rep_id %% 10L == 0L) {
      message(sprintf("[%s] completed %d/%d", scenario, rep_id, cfg$reps))
    }
  }
}

results <- do.call(rbind, results_list)
summary_df <- summarize_results(results, true_ate)

write.csv(results, paste0(cfg$out_prefix, "_replicates.csv"), row.names = FALSE)
write.csv(summary_df, paste0(cfg$out_prefix, "_summary.csv"), row.names = FALSE)
write_latex_table(summary_df, paste0(cfg$out_prefix, "_table.tex"), cfg$n, cfg$reps)

print(summary_df)
