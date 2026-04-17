args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(x) {
  out <- list(
    n = 500L,
    reps = 1000L,
    seed = 20260415L,
    mc_cores = 4L,
    split_reps = 10L,
    folds = 5L,
    bart_samples = 40L,
    bart_burn = 15L,
    bart_chains = 1L,
    maxit = 6L,
    out_prefix = "results/saturation_spline_n500_mc1000",
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
    if (key %in% c("n", "reps", "seed", "mc_cores", "split_reps", "folds", "bart_samples", "bart_burn", "bart_chains", "maxit")) {
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

score_clip <- 8
info_ridge <- 1e-3
step_clip <- 0.5
min_bandwidth <- 0.25
sanity_bound <- 5

suppressPackageStartupMessages({
  library(MASS)
  library(parallel)
  library(splines)
})

bart_backend <- if (requireNamespace("bartCause", quietly = TRUE)) {
  "bartCause"
} else if (requireNamespace("bartMachine", quietly = TRUE)) {
  bartMachine::set_bart_machine_num_cores(1L)
  "bartMachine"
} else {
  stop("Neither 'bartCause' nor 'bartMachine' is installed.")
}

dir.create(dirname(cfg$out_prefix), showWarnings = FALSE, recursive = TRUE)
message(sprintf("Using %s for the BART benchmark.", bart_backend))

scenario_levels <- c("heavy_tailed", "skewed_mixture")
cfg$scenarios <- intersect(cfg$scenarios, scenario_levels)
if (length(cfg$scenarios) == 0L) {
  stop("No valid scenarios requested.")
}

scenario_labels <- c(
  heavy_tailed = "Heavy-tailed",
  skewed_mixture = "Skewed mixture"
)

sat_mean <- 1 - log(5) / 4
true_ate <- 0.875
tau_intercept <- true_ate - 0.35 * sat_mean + 0.30 * 0.5

make_spline_df <- function(u) {
  out <- as.data.frame(
    splines::ns(
      u,
      knots = c(1, 2, 3),
      Boundary.knots = c(0, 4),
      intercept = FALSE
    )
  )
  names(out) <- paste0("S", seq_len(ncol(out)))
  out
}

fit_formula <- Y ~ A + S1 + S2 + S3 + S4 + W2 + W3 + W4 + A:S1 + A:S2 + A:S3 + A:S4 + A:W3
rhs_formula <- stats::delete.response(stats::terms(fit_formula))
wm_formula <- Y ~ S1 + S2 + S3 + S4 + W2 + W3 + W4
ps_formula <- A ~ S1 + S2 + S3 + S4 + W2 + W3 + W4

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
  W1 <- runif(n, min = 0, max = 4)
  W2 <- rnorm(n)
  W3 <- rbinom(n, 1L, 0.5)
  W4 <- runif(n, min = -1, max = 1)
  sat_u <- W1 / (1 + W1)
  spline_df <- make_spline_df(W1)

  p_raw <- plogis(-0.20 + 1.00 * sat_u - 0.55 * W2 + 0.60 * W3 - 0.45 * W4)
  p <- pmin(pmax(p_raw, 0.05), 0.95)
  A <- rbinom(n, 1L, p)
  eps <- generate_errors(n, scenario)

  mu0 <- 0.25 + 0.70 * sat_u - 0.40 * W2 + 0.35 * W3 + 0.25 * W4
  tau <- tau_intercept + 0.35 * sat_u - 0.30 * W3
  Y <- mu0 + A * tau + eps

  cbind(
    data.frame(Y = Y, A = A, W1 = W1, W2 = W2, W3 = W3, W4 = W4, sat_u = sat_u),
    spline_df
  )
}

kernel_scores <- function(e, r, h) {
  u <- outer(as.numeric(e), as.numeric(r), "-") / h
  K <- dnorm(u)
  f <- rowMeans(K) / h
  f <- pmax(f, 1e-8)
  fp <- -rowMeans(u * K) / (h^2)
  fpp <- rowMeans((u^2 - 1) * K) / (h^3)
  l1 <- fp / f
  l1 <- pmax(pmin(l1, score_clip), -score_clip)
  l2 <- fpp / f - l1^2
  list(l1 = l1, l2 = l2)
}

ate_gradient <- function(dat, coef_names) {
  g <- rep(0, length(coef_names))
  names(g) <- coef_names
  g["A"] <- 1
  g["A:S1"] <- mean(dat$S1)
  g["A:S2"] <- mean(dat$S2)
  g["A:S3"] <- mean(dat$S3)
  g["A:S4"] <- mean(dat$S4)
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

fit_proposed_core <- function(dat, maxit = cfg$maxit) {
  fit0 <- lm(fit_formula, data = dat)
  X <- model.matrix(fit0)
  y <- dat$Y
  beta <- coef(fit0)
  beta_init <- beta
  p <- ncol(X)

  for (it in seq_len(maxit)) {
    r <- as.numeric(y - X %*% beta)
    h <- max(1.06 * stats::sd(r) * length(r)^(-1 / 5), min_bandwidth)
    ks <- kernel_scores(r, r, h)
    U <- crossprod(X, ks$l1)
    J <- -crossprod(X * as.numeric(ks$l2), X) + diag(info_ridge, p)
    step <- tryCatch(solve(J, U), error = function(e) MASS::ginv(J) %*% U)
    step <- pmax(pmin(as.numeric(step), step_clip), -step_clip)
    beta_new <- as.numeric(beta - as.numeric(step))
    names(beta_new) <- colnames(X)
    if (!all(is.finite(beta_new)) || max(abs(beta_new)) > 20) {
      beta <- beta_init
      break
    }
    if (max(abs(beta_new - beta)) < 1e-6) {
      beta <- beta_new
      break
    }
    beta <- beta_new
  }

  r <- as.numeric(y - X %*% beta)
  h <- max(1.06 * stats::sd(r) * length(r)^(-1 / 5), min_bandwidth)
  ks <- kernel_scores(r, r, h)
  Ahat <- crossprod(X * as.numeric(ks$l2), X) / nrow(X)
  Ahat <- Ahat + diag(info_ridge, ncol(Ahat))
  Ainv <- tryCatch(solve(Ahat), error = function(e) MASS::ginv(Ahat))

  list(
    beta = beta,
    Ainv = Ainv,
    train_resid = r,
    h = h
  )
}

score_l1_new <- function(e_new, ref_resid, h) {
  u <- outer(as.numeric(e_new), as.numeric(ref_resid), "-") / h
  K <- dnorm(u)
  f <- rowMeans(K) / h
  f <- pmax(f, 1e-8)
  fp <- -rowMeans(u * K) / (h^2)
  l1 <- fp / f
  pmax(pmin(l1, score_clip), -score_clip)
}

fit_repeated_cf_proposed <- function(dat, split_reps = cfg$split_reps, folds = cfg$folds, maxit = cfg$maxit) {
  n <- nrow(dat)
  g_full <- ate_gradient(dat, colnames(model.matrix(fit_formula, data = dat)))
  split_est <- numeric(split_reps)
  split_var <- numeric(split_reps)

  for (b in seq_len(split_reps)) {
    fold_id <- sample(rep(seq_len(folds), length.out = n))
    delta_all <- numeric(n)
    infl_all <- matrix(NA_real_, nrow = n, ncol = length(g_full))

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
  if (!is.finite(est) || !is.finite(se) || abs(est) > sanity_bound || se > sanity_bound) {
    return(c(estimate = NA_real_, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_))
  }
  c(
    estimate = est,
    se = se,
    ci_lower = est - 1.96 * se,
    ci_upper = est + 1.96 * se
  )
}

fit_aipw_repeated_cf <- function(dat, split_reps = cfg$split_reps, folds = cfg$folds) {
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

run_one <- function(scenario, rep_id, cfg) {
  dat <- generate_data(cfg$n, scenario)
  out <- list()

  out[["Proposed TMLE"]] <- tryCatch(
    fit_repeated_cf_proposed(dat, split_reps = cfg$split_reps, folds = cfg$folds, maxit = cfg$maxit),
    error = function(e) c(estimate = NA_real_, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_)
  )
  out[["Gaussian OLS"]] <- tryCatch(
    fit_ols(dat),
    error = function(e) c(estimate = NA_real_, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_)
  )
  out[["AIPW"]] <- tryCatch(
    fit_aipw_repeated_cf(dat, split_reps = cfg$split_reps, folds = cfg$folds),
    error = function(e) c(estimate = NA_real_, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_)
  )
  out[["BART"]] <- tryCatch(
    fit_bart(dat, cfg$bart_samples, cfg$bart_burn, cfg$bart_chains),
    error = function(e) c(estimate = NA_real_, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_)
  )

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
  estimator_levels <- c("Proposed TMLE", "Gaussian OLS", "AIPW", "BART")
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
  lines <- c(
    sprintf("%% Generated on %s", Sys.time()),
    "\\begin{table}[t]",
    "\\centering",
    sprintf("\\caption{Saturation-model simulation with a common natural cubic spline working basis at $n=%d$ over %d Monte Carlo replications. The true conditional mean uses the nonlinear saturation term $W_1/(1+W_1)$, whereas the proposed estimator, Gaussian OLS, and AIPW all receive only the same four-dimensional natural spline basis in $W_1$.}", n, reps),
    "\\label{tab:saturation_spline}",
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
    sprintf("\\footnotesize The saturation covariate is generated as $W_1\\sim\\mathrm{Unif}(0,4)$, and the treatment effect is averaged to the true ATE $%.3f$. The proposed estimator, Gaussian OLS, and AIPW all use the same fixed-knot natural cubic spline basis in $W_1$ with internal knots at $1$, $2$, and $3$; BART receives the raw covariates $(W_1,W_2,W_3,W_4)$. ESD, empirical standard deviation; Cover., empirical coverage of the nominal $95\\%%$ confidence interval; Width, average confidence interval width. Within each scenario, boldface marks the best displayed value in each column, using smallest absolute bias, smallest ESD, smallest RMSE, coverage closest to $0.95$, and smallest interval width; ties are both boldfaced.", true_ate),
    "\\end{table}"
  )
  writeLines(lines, path)
}

set.seed(cfg$seed)
tasks <- expand.grid(
  scenario = cfg$scenarios,
  rep_id = seq_len(cfg$reps),
  stringsAsFactors = FALSE
)
task_seeds <- sample.int(.Machine$integer.max, nrow(tasks))

worker <- function(i) {
  set.seed(task_seeds[i])
  run_one(tasks$scenario[i], tasks$rep_id[i], cfg)
}

if (.Platform$OS.type == "unix" && cfg$mc_cores > 1L) {
  message(sprintf("Running %d saturation-spline tasks with %d cores.", nrow(tasks), cfg$mc_cores))
  results_list <- parallel::mclapply(seq_len(nrow(tasks)), worker, mc.cores = cfg$mc_cores, mc.set.seed = FALSE)
} else {
  results_list <- vector("list", nrow(tasks))
  for (i in seq_len(nrow(tasks))) {
    results_list[[i]] <- worker(i)
    if (tasks$rep_id[i] %% 10L == 0L) {
      message(sprintf("[%s] completed %d/%d", tasks$scenario[i], tasks$rep_id[i], cfg$reps))
    }
  }
}

results <- do.call(rbind, results_list)
summary_df <- summarize_results(results, true_ate)

write.csv(results, paste0(cfg$out_prefix, "_replicates.csv"), row.names = FALSE)
write.csv(summary_df, paste0(cfg$out_prefix, "_summary.csv"), row.names = FALSE)
write_latex_table(summary_df, paste0(cfg$out_prefix, "_table.tex"), cfg$n, cfg$reps)

print(summary_df)
