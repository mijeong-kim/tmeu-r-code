suppressPackageStartupMessages({
  library(MASS)
})

prepare_data <- function() {
  if (requireNamespace("causaldata", quietly = TRUE)) {
    dat <- as.data.frame(causaldata::nsw_mixtape)
  } else {
    dat <- read.csv("data/nsw_mixtape.csv", stringsAsFactors = FALSE)
  }
  dat$Y <- asinh(dat$re78)
  dat$A <- dat$treat
  dat$re74_t <- asinh(dat$re74)
  dat$re75_t <- asinh(dat$re75)
  dat
}

kernel_scores <- function(e, r, h) {
  u <- outer(as.numeric(e), as.numeric(r), "-") / h
  K <- dnorm(u)
  f <- pmax(rowMeans(K) / h, 1e-8)
  fp <- -rowMeans(u * K) / (h^2)
  fpp <- rowMeans((u^2 - 1) * K) / (h^3)
  l1 <- fp / f
  l2 <- fpp / f - l1^2
  list(l1 = l1, l2 = l2)
}

fit_formula <- Y ~ A + age + educ + black + hisp + marr + nodegree + re74_t + re75_t + A:re74_t + A:re75_t
wm_formula <- Y ~ age + educ + black + hisp + marr + nodegree + re74_t + re75_t
ps_formula <- A ~ age + educ + black + hisp + marr + nodegree + re74_t + re75_t

ate_gradient <- function(dat, coef_names) {
  g <- rep(0, length(coef_names))
  names(g) <- coef_names
  g["A"] <- 1
  g["A:re74_t"] <- mean(dat$re74_t)
  g["A:re75_t"] <- mean(dat$re75_t)
  g
}

delta_from_beta <- function(beta, dat) {
  beta <- beta[names(beta)]
  beta["A"] + beta["A:re74_t"] * dat$re74_t + beta["A:re75_t"] * dat$re75_t
}

score_l1_new <- function(e_new, ref_resid, h) {
  u <- outer(as.numeric(e_new), as.numeric(ref_resid), "-") / h
  K <- dnorm(u)
  f <- pmax(rowMeans(K) / h, 1e-8)
  fp <- -rowMeans(u * K) / (h^2)
  fp / f
}

fit_proposed_stable <- function(dat, lambda = 0.01) {
  fit0 <- lm(fit_formula, data = dat)
  X <- model.matrix(fit0)
  y <- dat$Y
  beta0 <- coef(fit0)

  r <- as.numeric(y - X %*% beta0)
  h <- max(1.06 * stats::sd(r) * length(r)^(-1 / 5), 0.15)
  ks <- kernel_scores(r, r, h)
  U <- crossprod(X, ks$l1)
  J <- -crossprod(X * as.numeric(ks$l2), X) + diag(1e-6, ncol(X))
  step <- tryCatch(solve(J, U), error = function(e) MASS::ginv(J) %*% U)
  beta <- as.numeric(beta0 - lambda * as.numeric(step))
  names(beta) <- colnames(X)

  r2 <- as.numeric(y - X %*% beta)
  h2 <- max(1.06 * stats::sd(r2) * length(r2)^(-1 / 5), 0.15)
  ks2 <- kernel_scores(r2, r2, h2)
  Ahat <- crossprod(X * as.numeric(ks2$l2), X) / nrow(X) + diag(1e-6, ncol(X))
  Bhat <- crossprod(X * as.numeric(ks2$l1), X * as.numeric(ks2$l1)) / nrow(X)
  Ainv <- tryCatch(solve(Ahat), error = function(e) MASS::ginv(Ahat))
  vcov_beta <- Ainv %*% Bhat %*% t(Ainv) / nrow(X)

  g <- ate_gradient(dat, colnames(X))
  ate <- as.numeric(sum(g * beta))
  se <- sqrt(drop(t(g) %*% vcov_beta %*% g))

  list(
    estimate = ate,
    se = se,
    ci_lower = ate - 1.96 * se,
    ci_upper = ate + 1.96 * se,
    beta = beta,
    Ainv = Ainv,
    train_resid = r2,
    h = h2
  )
}

fit_proposed_rcf <- function(dat, lambda = 0.01, K = 5L, split_reps = 20L, seed = 1L) {
  set.seed(seed)
  n <- nrow(dat)
  coef_names <- colnames(model.matrix(fit_formula, data = dat))
  g_full <- ate_gradient(dat, coef_names)
  split_est <- numeric(split_reps)
  split_var <- numeric(split_reps)

  for (b in seq_len(split_reps)) {
    fold <- sample(rep(seq_len(K), length.out = n))
    delta_all <- numeric(n)
    infl_all <- matrix(NA_real_, nrow = n, ncol = length(g_full))
    colnames(infl_all) <- names(g_full)

    for (k in seq_len(K)) {
      train <- dat[fold != k, , drop = FALSE]
      test_idx <- which(fold == k)
      test <- dat[test_idx, , drop = FALSE]
      fit <- fit_proposed_stable(train, lambda = lambda)

      delta_all[test_idx] <- delta_from_beta(fit$beta, test)

      X_test <- model.matrix(fit_formula, data = test)
      e_test <- as.numeric(test$Y - X_test %*% fit$beta)
      l1_test <- score_l1_new(e_test, fit$train_resid, fit$h)
      infl_all[test_idx, ] <- (X_test * as.numeric(l1_test)) %*% t(fit$Ainv)
    }

    split_est[b] <- mean(delta_all)
    eif <- delta_all - split_est[b] + as.numeric(infl_all %*% g_full)
    split_var[b] <- mean(eif^2) / n
  }

  estimate <- mean(split_est)
  se <- sqrt(mean(split_var) + stats::var(split_est))
  list(
    estimate = estimate,
    se = se,
    ci_lower = estimate - 1.96 * se,
    ci_upper = estimate + 1.96 * se,
    split_sd = stats::sd(split_est)
  )
}

fit_aipw_rcf <- function(dat, K = 5L, split_reps = 20L, seed = 1L) {
  set.seed(seed)
  n <- nrow(dat)
  split_est <- numeric(split_reps)
  split_var <- numeric(split_reps)

  for (b in seq_len(split_reps)) {
    fold <- sample(rep(seq_len(K), length.out = n))
    pseudo <- numeric(n)

    for (k in seq_len(K)) {
      train <- dat[fold != k, , drop = FALSE]
      test_idx <- which(fold == k)
      test <- dat[test_idx, , drop = FALSE]

      if (sum(train$A == 1) < 10L || sum(train$A == 0) < 10L) {
        return(list(
          estimate = NA_real_,
          se = NA_real_,
          ci_lower = NA_real_,
          ci_upper = NA_real_,
          split_sd = NA_real_
        ))
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

  estimate <- mean(split_est)
  se <- sqrt(mean(split_var) + stats::var(split_est))
  list(
    estimate = estimate,
    se = se,
    ci_lower = estimate - 1.96 * se,
    ci_upper = estimate + 1.96 * se,
    split_sd = stats::sd(split_est)
  )
}

fit_bart_once <- function(dat, seed = 1L, n.samples = 100L, n.burn = 30L) {
  x <- dat[c("age", "educ", "black", "hisp", "marr", "nodegree", "re74_t", "re75_t")]
  fit <- bartCause::bartc(
    response = dat$Y,
    treatment = dat$A,
    confounders = x,
    method.rsp = "bart",
    method.trt = "none",
    p.scoreAsCovariate = FALSE,
    estimand = "ate",
    verbose = FALSE,
    n.samples = n.samples,
    n.burn = n.burn,
    n.chains = 1L,
    n.threads = 1L,
    seed = seed
  )
  est <- summary(fit)$estimates[1L, ]
  c(
    estimate = as.numeric(est[["estimate"]]),
    se = as.numeric(est[["sd"]]),
    ci_lower = as.numeric(est[["ci.lower"]]),
    ci_upper = as.numeric(est[["ci.upper"]])
  )
}

fit_grf_once <- function(dat, seed = 1L, num.trees = 800L) {
  if (!requireNamespace("grf", quietly = TRUE)) {
    stop("Package 'grf' is required for the forest-based TMLE benchmark.")
  }
  x <- as.matrix(dat[c("age", "educ", "black", "hisp", "marr", "nodegree", "re74_t", "re75_t")])
  fit <- grf::causal_forest(
    X = x,
    Y = dat$Y,
    W = dat$A,
    num.trees = num.trees,
    num.threads = 1L,
    seed = seed
  )
  ate <- grf::average_treatment_effect(fit, target.sample = "all", method = "TMLE")
  c(
    estimate = as.numeric(ate[["estimate"]]),
    se = as.numeric(ate[["std.err"]]),
    ci_lower = as.numeric(ate[["estimate"]]) - 1.96 * as.numeric(ate[["std.err"]]),
    ci_upper = as.numeric(ate[["estimate"]]) + 1.96 * as.numeric(ate[["std.err"]])
  )
}

dat <- prepare_data()

prop <- fit_proposed_rcf(dat, lambda = 0.01, K = 5L, split_reps = 20L, seed = 20260324L)
set.seed(20260326L)
prop_split_sd <- stats::sd(as.numeric(replicate(
  20,
  fit_proposed_rcf(dat, lambda = 0.01, K = 5L, split_reps = 20L, seed = sample.int(1e6, 1L))[["estimate"]]
)))

aipw <- fit_aipw_rcf(dat, K = 5L, split_reps = 20L, seed = 20260327L)
set.seed(20260328L)
aipw_split_sd <- stats::sd(as.numeric(replicate(
  20,
  fit_aipw_rcf(dat, K = 5L, split_reps = 20L, seed = sample.int(1e6, 1L))[["estimate"]]
)))

if (requireNamespace("bartCause", quietly = TRUE)) {
  bart_main <- fit_bart_once(dat, seed = 20260322L, n.samples = 100L, n.burn = 30L)
  set.seed(20260325L)
  bart_split_sd <- stats::sd(as.numeric(replicate(20, fit_bart_once(dat, seed = sample.int(1e6, 1L), n.samples = 80L, n.burn = 20L)["estimate"])))
} else {
  message("Package 'bartCause' not available; reusing the stored BART benchmark from results/nsw_table2.csv.")
  prev <- read.csv("results/nsw_table2.csv", stringsAsFactors = FALSE)
  prev_bart <- prev[prev$estimator == "BART plug-in", , drop = FALSE]
  if (nrow(prev_bart) != 1L) {
    stop("BART benchmark unavailable: install 'bartCause' or provide a stored BART row in results/nsw_table2.csv.")
  }
  bart_main <- c(
    estimate = prev_bart$ate[1],
    se = prev_bart$se[1],
    ci_lower = prev_bart$ci_lower[1],
    ci_upper = prev_bart$ci_upper[1]
  )
  bart_split_sd <- prev_bart$split_sd[1]
}

if (requireNamespace("grf", quietly = TRUE)) {
  grf_main <- fit_grf_once(dat, seed = 20260416L, num.trees = 800L)
  set.seed(20260417L)
  grf_split_sd <- stats::sd(as.numeric(replicate(
    20,
    fit_grf_once(dat, seed = sample.int(1e6, 1L), num.trees = 600L)[["estimate"]]
  )))
} else {
  grf_main <- c(estimate = NA_real_, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_)
  grf_split_sd <- NA_real_
}

res <- data.frame(
  estimator = c("Proposed TMLE", "AIPW", "BART plug-in", "Forest-based TMLE"),
  ate = c(prop$estimate, aipw$estimate, bart_main[["estimate"]], grf_main[["estimate"]]),
  se = c(prop$se, aipw$se, bart_main[["se"]], grf_main[["se"]]),
  ci_lower = c(prop$ci_lower, aipw$ci_lower, bart_main[["ci_lower"]], grf_main[["ci_lower"]]),
  ci_upper = c(prop$ci_upper, aipw$ci_upper, bart_main[["ci_upper"]], grf_main[["ci_upper"]]),
  width = c(
    prop$ci_upper - prop$ci_lower,
    aipw$ci_upper - aipw$ci_lower,
    bart_main[["ci_upper"]] - bart_main[["ci_lower"]],
    grf_main[["ci_upper"]] - grf_main[["ci_lower"]]
  ),
  split_sd = c(prop_split_sd, aipw_split_sd, bart_split_sd, grf_split_sd),
  stringsAsFactors = FALSE
)

dir.create("results", showWarnings = FALSE, recursive = TRUE)
write.csv(res, "results/nsw_table2.csv", row.names = FALSE)
print(res)
