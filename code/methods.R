
# Methods for galuppi.qmd

# Convention of the ridge inference below: X carries no intercept column and
# has unit-sd columns, and lambda is in "theory" scale (X'X + 2 * lambda * I),
# so callers pass (n - 1) * lambda_glmnet / 2, which is exact: glmnet
# standardizes internally by the 1 / n variance while the recipe delivers
# unit (n - 1)-sd columns, rescaling the penalty by (n - 1) / n. Reading
# eq. (2) without that factor changes the standard errors by at most 0.3%.

# Notation map to Section 3.3 and Appendix A of the paper:
# - H(beta) of eq. (hess) is the summed crossprod XtWX below.
# - The penalty matrix P_lambda is `penalty`/`dpen`.
# - Sigma_lambda of eq. (Sigmalambda) is what the *_cov() functions return.
# - Coefficients are class-major, as.vector(t(beta_mat)).
# - Lambda chain: the eq-(2) lambda is lambda_glmnet / 2 for the multinomial
#   fits and lambda_glmnet for the binomial one, and the theory-scale lambda
#   here is (n - 1) * lambda_glmnet / 2, both conversions pinned in
#   test_methods.R.

## Logistic ridge inference
{

# Le Cessie & Van Houwelingen (1992, Section 2.1, p. 194, an unnumbered
# equation after (2.3)) sandwich for the logistic ridge,
# (X'WX + P)^-1 X'WX (X'WX + P)^-1, with P = diag(0, 2 * lambda, ..., 2 *
# lambda). Their design carries no constant term (p. 193); the leading 0 of P
# adds an unpenalized intercept to their formula.
logis_ridge_cov <- function(X, beta_hat, lambda) {

  # Needs standardized data (no intercept column)
  X <- as.matrix(X)
  p <- ncol(X)
  stopifnot(max(abs(apply(X, 2, sd) - 1)) < 1e-6)

  # Binomial weights of the Fisher information
  eta_hat <- drop(beta_hat[1] + X %*% beta_hat[-1])
  p_hat <- 1 / (1 + exp(-eta_hat))
  w_hat <- p_hat * (1 - p_hat)

  # Augmented design with intercept, intercept unpenalized
  Xt <- cbind("(Intercept)" = 1, X)
  H <- crossprod(Xt * sqrt(w_hat))
  penalty <- diag(c(0, rep(2 * lambda, p)))
  H_lambda_inv <- solve(H + penalty)

  # Asymptotic covariance, (p + 1) x (p + 1) including the intercept
  H_lambda_inv %*% H %*% H_lambda_inv

}

# Wald CIs for the slopes of the logistic ridge, Bonferroni over the family
# size m when is given (typically the number of slopes) and marginal when m
# is NULL. z_value and p_value are always marginal; only p_value_adj and the
# interval widths carry m.
logis_ridge_ci <- function(X, beta_hat, lambda, level = 0.95, m = NULL) {

  # Needs standardized data (no intercept column)
  X <- as.matrix(X)
  p <- ncol(X)
  stopifnot(max(abs(apply(X, 2, sd) - 1)) < 1e-6)

  # Asymptotic covariance of slopes and intercept
  Sigma_beta <- logis_ridge_cov(X = X, beta_hat = beta_hat, lambda = lambda)

  # Standard errors and z-based CIs for the slopes
  se <- sqrt(diag(Sigma_beta))[-1]
  alpha_slopes <- if (is.null(m)) 1 - level else (1 - level) / m
  z_crit <- qnorm(1 - alpha_slopes / 2)
  ci_low <- beta_hat[-1] - z_crit * se
  ci_up <- beta_hat[-1] + z_crit * se

  # Adjusted p-values are Bonferroni on the slopes
  tibble::tibble(
    term = colnames(X) %||% paste0("x", seq_len(p)),
    estimate = beta_hat[-1],
    std_error = se,
    ci_low = ci_low,
    ci_up = ci_up,
    z_value = estimate / std_error,
    p_value = 2 * pnorm(-abs(z_value)),
    p_value_adj = pmin(1, (if (is.null(m)) 1 else m) * p_value)
  )

}

# Prediction CIs for the logistic ridge: a Wald interval for
# eta_0 = beta_0 + x_0'beta, mapped through the logistic link.
logis_ridge_pred_ci <- function(X, beta_hat, lambda, new_X, level = 0.95) {

  # Needs standardized data (no intercept column)
  X <- as.matrix(X)
  stopifnot(max(abs(apply(X, 2, sd) - 1)) < 1e-6)
  new_X <- as.matrix(new_X)

  # Asymptotic covariance of slopes and intercept
  Sigma_beta <- logis_ridge_cov(X = X, beta_hat = beta_hat, lambda = lambda)

  # Predicted linear predictor and probability at the new points
  Xt_new <- cbind(1, new_X)
  eta_new <- drop(Xt_new %*% beta_hat)
  p_new <- 1 / (1 + exp(-eta_new))

  # Var(eta_new_i) = x_i' Sigma_beta x_i for each row i of Xt_new
  var_eta_new <- rowSums((Xt_new %*% Sigma_beta) * Xt_new)
  se_eta_new <- sqrt(var_eta_new)

  # Wald interval on the linear predictor
  z_crit <- qnorm(1 - (1 - level) / 2)
  eta_low <- eta_new - z_crit * se_eta_new
  eta_up <- eta_new + z_crit * se_eta_new

  # Wald interval mapped to the probability
  ci_low <- 1 / (1 + exp(-eta_low))
  ci_up <- 1 / (1 + exp(-eta_up))
  tibble::tibble(
    .row = seq_len(nrow(new_X)),
    eta = eta_new,
    eta_se = se_eta_new,
    p_hat = p_new,
    p_ci_low = ci_low,
    p_ci_up = ci_up
  )

}

}

## Multinomial logistic ridge inference
{

# Asymptotic covariance Sigma = (H + P)^+ H (H + P)^+ of the multinomial
# ridge in glmnet's symmetric parametrization, returned in class-major order,
# i.e. for as.vector(t(beta_mat)) (Zahid & Tutz, 2013, Section 2.1, p. 1025).
# Zahid & Tutz use a plain inverse because they work in the K - 1 free
# sum-to-zero coordinates (logits against the geometric mean, not
# reference-category logits); the whole covariance is equal across the two
# parametrizations, E Sigma_free E' with E the embedding, see the Zahid & Tutz
# test in test_methods.R.
# K = 2 is the binary case at half the penalty, see test_methods.R.
# Conservative by construction, see coverage_study.R.
mult_ridge_cov <- function(X, beta_mat, lambda) {

  # Needs standardized data (no intercept column)
  X <- as.matrix(X)
  stopifnot(max(abs(apply(X, 2, sd) - 1)) < 1e-6)

  # beta_mat is K x (p + 1), row k = (intercept_k, slopes_k)
  n <- nrow(X)
  p <- ncol(X)
  K <- nrow(beta_mat)

  # Add column of ones
  q <- p + 1
  Xt <- cbind(1, X) # n x q

  # Fitted class probabilities as softmax of the linear predictors, n x K
  eta <- Xt %*% t(beta_mat)
  eta <- eta - apply(eta, 1, max) # For numerical stabilization
  ex <- exp(eta)
  P_hat <- ex / rowSums(ex)

  # Multinomial Fisher information H by blocks, kth block is
  # l = Xt' diag(P_k (delta_kl - P_l)) Xt)
  H <- matrix(0, K * q, K * q)
  for (k in 1:K) {
    for (l in 1:K) {
      w <- P_hat[, k] * ((k == l) - P_hat[, l]) # length n weights
      block <- crossprod(Xt, Xt * w) # q x q = Xt' diag(w) Xt
      H[((k - 1) * q + 1):(k * q), ((l - 1) * q + 1):(l * q)] <- block
    }
  }

  # Ridge penalty Hessian, 2 * lambda per class on the slopes
  penalty <- diag(rep(c(0, rep(2 * lambda, p)), times = K))

  # Pseudo-inverse: (H + P) has a 1-dim null space (the common shift of all
  # class intercepts), harmless as every reported quantity is shift-invariant.
  H_lambda_inv <- MASS::ginv(H + penalty)
  H_lambda_inv %*% H %*% H_lambda_inv

}

# CIs for the slopes of the multinomial ridge. Mirrors logis_ridge_ci:
# intercepts are excluded from the rows and from the Bonferroni family m,
# which is K * p when all slopes are reported at once, and z_value and p_value
# stay marginal.
mult_ridge_ci <- function(X, beta_mat, lambda, level = 0.95, m = NULL,
                          class_names = NULL) {

  # Needs standardized data (no intercept column)
  X <- as.matrix(X)
  stopifnot(max(abs(apply(X, 2, sd) - 1)) < 1e-6)

  # Asymptotic covariance of slopes and intercepts, in class-major order
  p <- ncol(X)
  K <- nrow(beta_mat)
  q <- p + 1
  Sigma_beta <- mult_ridge_cov(X = X, beta_mat = beta_mat, lambda = lambda)

  # Extract coefficients and standard errors of the slopes in class-major order,
  # intercept entries dropped
  keep <- rep(c(FALSE, rep(TRUE, p)), times = K)
  se <- sqrt(pmax(diag(Sigma_beta), 0))[keep]
  beta_vec <- as.vector(t(beta_mat))[keep]

  # Standard errors and z-based CIs for the slopes
  alpha_slopes <- if (is.null(m)) 1 - level else (1 - level) / m
  z_crit <- qnorm(1 - alpha_slopes / 2)
  ci_low <- beta_vec - z_crit * se
  ci_up <- beta_vec + z_crit * se
  terms <- colnames(X) %||% paste0("x", seq_len(p))
  class_names <- class_names %||% paste0("class", seq_len(K))

  # Adjusted p-values are Bonferroni on the slopes
  tibble::tibble(
    class = rep(class_names, each = p),
    term = rep(terms, times = K),
    estimate = beta_vec,
    std_error = se,
    ci_low = ci_low,
    ci_up = ci_up,
    z_value = estimate / std_error,
    p_value = 2 * pnorm(-abs(z_value)),
    p_value_adj = pmin(1, (if (is.null(m)) 1 else m) * p_value)
  )

}

# Linear predictors, class probabilities, and Var(eta_0) at one new point x0
# (with its leading 1).
mult_ridge_eta <- function(x0, beta_mat, Sigma_beta) {

  # Dimensions
  K <- nrow(beta_mat)
  q <- length(x0)

  # Fitted class probabilities as softmax of the linear predictors
  eta0 <- as.numeric(beta_mat %*% x0)
  eta0 <- eta0 - max(eta0) # For numerical stabilization
  ex <- exp(eta0)
  p0 <- ex / sum(ex)

  # Var(eta0) is K x K: entry (k, l) = x0' Sigma_{kl} x0
  var_eta <- matrix(0, K, K)
  for (k in 1:K) {
    for (l in 1:K) {
      block <- Sigma_beta[((k - 1) * q + 1):(k * q), ((l - 1) * q + 1):(l * q)]
      var_eta[k, l] <- as.numeric(t(x0) %*% block %*% x0)
    }
  }

  # Result
  list(eta0 = eta0, p0 = p0, var_eta = var_eta)

}

# Prediction CIs for the multinomial ridge: the softmax probabilities are
# propagated by the delta method through J = diag(p_0) - p_0 p_0', and the
# Wald interval is built on logit(p_k) and mapped back.
mult_ridge_pred_ci <- function(X, beta_mat, lambda, new_X, level = 0.95,
                               class_names = NULL) {

  # Prepare needed objects
  X <- as.matrix(X)
  new_X <- as.matrix(new_X)
  K <- nrow(beta_mat)
  Sigma_beta <- mult_ridge_cov(X = X, beta_mat = beta_mat, lambda = lambda)
  Xt_new <- cbind(1, new_X)
  n_new <- nrow(Xt_new)
  z_crit <- qnorm(1 - (1 - level) / 2)
  class_names <- class_names %||% paste0("class", seq_len(K))

  # Compute the prediction intervals one row at a time
  out <- vector("list", n_new)
  for (i in 1:n_new) {

    # Predicted linear predictor, probability, and Var(eta_0) at the new point
    pred0 <- mult_ridge_eta(x0 = Xt_new[i, ], beta_mat = beta_mat,
                            Sigma_beta = Sigma_beta)
    p0 <- pred0$p0
    var_eta <- pred0$var_eta

    # Delta method through the softmax Jacobian
    jac <- diag(p0) - tcrossprod(p0)
    var_p <- jac %*% var_eta %*% t(jac)
    se_p <- sqrt(pmax(diag(var_p), 0))

    # Second delta step to the logit scale, where the Wald interval lives
    se_logit <- se_p / (p0 * (1 - p0))
    logit_p0 <- log(p0 / (1 - p0))
    out[[i]] <- tibble::tibble(
      .row = i,
      class = class_names,
      p_hat = p0,
      p_ci_low = 1 / (1 + exp(-(logit_p0 - z_crit * se_logit))),
      p_ci_up = 1 / (1 + exp(-(logit_p0 + z_crit * se_logit)))
    )

  }

  # Results
  dplyr::bind_rows(out)

}

# CI for the contrast between two class probabilities at a new x_0, since
# overlapping marginal intervals are not a test of pi_j - pi_k. The test lives
# on the log-odds contrast eta_j - eta_k, which is linear and identified, so
# its Wald interval needs no delta step and no truncation; the probability
# difference is returned as the effect size.
mult_ridge_contrast_ci <- function(X, beta_mat, lambda, new_X,
                                   class_j, class_k, level = 0.95,
                                   class_names = NULL) {

  # Prepare needed objects
  X <- as.matrix(X)
  new_X <- as.matrix(new_X)
  K <- nrow(beta_mat)
  Sigma_beta <- mult_ridge_cov(X = X, beta_mat = beta_mat, lambda = lambda)
  Xt_new <- cbind(1, new_X)
  n_new <- nrow(Xt_new)
  z_crit <- qnorm(1 - (1 - level) / 2)
  class_names <- class_names %||% paste0("class", seq_len(K))

  # class_j, class_k: indices or names of the two classes to contrast
  j <- if (is.character(class_j)) match(class_j, class_names) else class_j
  k <- if (is.character(class_k)) match(class_k, class_names) else class_k
  stopifnot(!is.na(j), !is.na(k), j != k)

  # Compute the prediction intervals one row at a time
  out <- vector("list", n_new)
  for (i in seq_len(n_new)) {

    # Predicted linear predictor, probability, and Var(eta_0) at the new point
    pred0 <- mult_ridge_eta(x0 = Xt_new[i, ], beta_mat = beta_mat,
                            Sigma_beta = Sigma_beta)
    p0 <- pred0$p0
    eta0 <- pred0$eta0
    var_eta <- pred0$var_eta

    # Probability difference, the effect size
    diff_jk <- p0[j] - p0[k]

    # Identified log-odds contrast, invariant to the common shift of eta0
    logodds_jk <- eta0[j] - eta0[k]
    se_logodds_jk <- sqrt(max(var_eta[j, j] + var_eta[k, k] -
                                2 * var_eta[j, k], 0))
    out[[i]] <- tibble::tibble(
      .row = i,
      class_j = class_names[j],
      class_k = class_names[k],
      diff = diff_jk,
      logodds = logodds_jk,
      logodds_se = se_logodds_jk,
      logodds_ci_low = logodds_jk - z_crit * se_logodds_jk,
      logodds_ci_up = logodds_jk + z_crit * se_logodds_jk
    )

  }

  # Results
  dplyr::bind_rows(out)

}

}

## FDR inference
{

# Benjamini-Hochberg selection at FDR alpha over the m slope p-values of a
# *_ridge_ci() table, with the Benjamini & Yekutieli (2005) false-coverage-rate
# intervals of the R selected coefficients: marginal Wald CIs at level
# 1 - R * alpha / m, i.e. critical value z_{1 - R * alpha / (2 m)}, so that
# the expected proportion of selected intervals missing their target is at
# most alpha. R = m recovers the marginal intervals; R = 1 the Bonferroni one.
fcr_ci <- function(ci_df, m = nrow(ci_df), alpha = 0.05) {

  # BH selection on the marginal p-values, at family size m (>= nrow when the
  # table is a subset of the family)
  p_bh <- pmin(1, p.adjust(ci_df$p_value, method = "BH", n = m))
  R <- sum(p_bh < alpha)
  z_fcr <- qnorm(1 - max(R, 1) * alpha / (2 * m))

  # FCR intervals for every row; selection flags the BH discoveries
  out <- ci_df
  out$p_bh <- p_bh
  out$ci_low <- out$estimate - z_fcr * out$std_error
  out$ci_up <- out$estimate + z_fcr * out$std_error
  out$selected <- p_bh < alpha

  # Result
  list(ci = out, R = R, z_fcr = z_fcr, alpha = alpha, m = m)

}

}

## Recipes
{

# Preprocessing recipe shared by the binary and multinomial problems
preprocess_recipe <- function(data) {

  recipe(Composer ~ ., data = data) |>
    # Remove ID variable
    step_rm(AriaId) |>
    # Preprocess all predictors BEFORE dummying
    step_nzv(all_predictors()) |>
    step_filter_missing(all_predictors(), threshold = 0.10) |>
    step_impute_knn(all_predictors(), neighbors = 10) |>
    # Collapse rare factor levels (BEFORE dummy creation)
    step_other(all_nominal_predictors(), threshold = 0.10) |>
    # Preprocess numeric predictors only
    step_corr(all_numeric_predictors(), threshold = 0.80) |>
    step_normalize(all_numeric_predictors()) |>
    # Dummies AFTER normalization, one-hot so that the "other" level can be
    # dropped below and become the reference.
    step_dummy(all_nominal_predictors(), one_hot = TRUE,
               naming = function(...) paste0("dum_", dummy_names(...))) |>
    # Drop the pooled "other" dummies
    step_rm(starts_with("dum_") & ends_with("_other")) |>
    # Post-dummy cleanup: the correlation filter now also deduplicates the
    # dummies. No step_lincomb() because that truncated the design to n - 1;
    # the exact structural dependencies that survive are handled with
    # regularization.
    step_nzv(all_predictors()) |>
    step_corr(all_predictors(), threshold = 0.80) |>
    # Standardize last so the dummies enter glmnet on unit-sd scale too, as the
    # uniform 2 * lambda penalty of the sandwich covariance requires.
    step_normalize(all_predictors())

}


# Top correlates of each selected marker among the features the preprocessing
# filtered out (setdiff of the pre-correlation pool and the final design):
# those companions never had the chance to be selected, and the post-encoding
# 0.8 filter guarantees no kept pair exceeds that threshold, so the removed
# pool is complete at the |rho| >= 0.8 level.
marker_correlates <- function(X_pool, final_cols, markers, top_k = 10) {

  X_pool <- as.matrix(X_pool)
  removed <- setdiff(colnames(X_pool), final_cols)
  missing_markers <- setdiff(markers, colnames(X_pool))
  if (length(missing_markers) > 0) {

    warning("Markers absent from the pool: ",
            paste(missing_markers, collapse = ", "))
    markers <- setdiff(markers, missing_markers)

  }

  # Marker-by-removed correlations, top_k per marker by absolute value
  rho <- cor(X_pool[, markers, drop = FALSE],
             X_pool[, removed, drop = FALSE])
  out <- tibble::tibble(
    marker = rep(rownames(rho), times = ncol(rho)),
    correlate = rep(colnames(rho), each = nrow(rho)),
    rho = as.vector(rho)
  ) |>
    dplyr::arrange(marker, dplyr::desc(abs(rho))) |>
    dplyr::group_by(marker) |>
    dplyr::slice_head(n = top_k) |>
    dplyr::ungroup()

  # Result
  return(out)

}

}

## Resampling
{

# Fold ids from an rsample object, so cv.glmnet and h2o share the same
# folds determined by tidymodels.
make_foldid <- function(folds, n) {

  foldid <- integer(n)
  for (i in seq_along(folds$splits)) {
    idx <- complement(folds$splits[[i]]) # Indexes in the split
    foldid[idx] <- i
  }
  stopifnot(all(foldid %in% seq_along(folds$splits)))
  foldid

}

# Shared resampling indices, so every model comparison uses the same B
# bootstrap samples of the test set.
boot_idx <- function(n, B = 1e5, seed = 987204452) {

  # Do not pollute the caller's RNG state with the internal seed
  withr::local_seed(seed)
  matrix(sample(n, size = n * B, replace = TRUE), nrow = n, ncol = B)

}

}

## Model fitting
{

# Fit and test-predict a penalized logistic regression, multinomial for K > 2,
# at the one-standard-error penalty of glmnet's own cross-validated deviance.
# We use cv.glmnet and not tune_grid because the tidymodels one-SE rule
# measures its standard error across folds rather than across observations.
fit_glmnet_model <- function(train, test, folds, mixture, penalty_range) {

  # Call to cv.glmnet
  K <- nlevels(train$Composer)
  family <- switch((K > 2) + 1, "binomial", "multinomial")
  reg <- switch((K > 2) + 1, logistic_reg, multinom_reg)
  lambda_seq <- exp(seq(max(penalty_range), min(penalty_range), l = 300))
  cv <- cv.glmnet(x = train |>
                    select(-Composer) |>
                    as.matrix(),
                  y = train$Composer,
                  family = family,
                  type.measure = "deviance",
                  alpha = mixture,
                  foldid = make_foldid(folds, nrow(train)),
                  parallel = TRUE,
                  lambda = lambda_seq)

  # Both the deviance minimizer and the one-standard-error penalty must be
  # strictly inside the grid, error if not and force changing penalty_range:
  # a railed lambda.min means the grid does not bracket the minimum, so the
  # one-standard-error anchor is not trustworthy either.
  if (cv$lambda.min <= min(lambda_seq) || cv$lambda.min >= max(lambda_seq) ||
      cv$lambda.1se <= min(lambda_seq) || cv$lambda.1se >= max(lambda_seq)) {
    stop("lambda.min or lambda.1se is at the edge of the grid, ",
         "change penalty_range = c(", signif(min(lambda_seq), 4), ", ",
         signif(max(lambda_seq), 4), ")")
  }
  best <- tibble(penalty = cv$lambda.1se)

  # path_values pins the workflow's glmnet path to the tuning grid, so the
  # fit evaluated at lambda.1se is exact rather than path-interpolated.
  spec <- reg(penalty = tune(), mixture = mixture) |>
    set_mode("classification") |>
    set_engine("glmnet", path_values = lambda_seq)
  wf <- workflow() |>
    add_model(spec) |>
    add_formula(Composer ~ .)
  fit <- finalize_workflow(wf, best) |>
    fit(train)
  list(fit = fit,
       params = paste0("penalty = ", signif(best$penalty, 4)),
       penalty = best$penalty,
       class = predict(fit, test)$.pred_class,
       prob = predict(fit, test, type = "prob")[[2]])

}

# Tune, fit, and test-predict a random forest with 5000 trees
fit_rf_model <- function(train, test, folds) {

  p <- ncol(train) - 1
  spec <- rand_forest(mtry = tune(), trees = 5000, min_n = tune()) |>
    set_mode("classification") |>
    set_engine("ranger", seed = 42)
  grid <- grid_regular(mtry(range = c(1, p)),
                       min_n(range = c(2, 10)),
                       levels = c(10, 3))
  wf <- workflow() |>
    add_model(spec) |>
    add_formula(Composer ~ .)

  # Seeded locally, so the tuning does not advance the caller's RNG state
  withr::local_seed(42)
  tuned <- tune_grid(wf, resamples = folds, grid = grid)
  best <- select_best(tuned, metric = "accuracy")
  fit <- finalize_workflow(wf, best) |>
    fit(train)
  list(fit = fit,
       params = paste0("mtry = ", best$mtry, ", min_n = ", best$min_n),
       class = predict(fit, test)$.pred_class,
       prob = predict(fit, test, type = "prob")[[2]])

}

# Fit and test-predict an h2o AutoML search over the full algorithm set, at a
# fixed model budget and seed. The search runs on all the threads of the h2o
# cluster and is not run-reproducible; cache_file makes the reported run
# transferable: when the file exists the search is skipped and its stored
# results are returned.
fit_automl_model <- function(train, test, unk, foldid, cache_file = NULL,
                             max_models = 100) {

  # Reuse the stored search if present
  if (!is.null(cache_file) && file.exists(cache_file)) {

    return(readRDS(cache_file))

  }

  # Add foldid columns to the train/test sets, so h2o can use the same folds as
  # the other models.
  train_foldid <- train |>
    mutate(foldid = foldid)
  test_foldid <- test |>
    mutate(foldid = 0L)

  # AUC is a binomial-only leaderboard metric in h2o
  lvls <- levels(train$Composer)
  metric <- switch((length(lvls) > 2) + 1, "AUC", "mean_per_class_error")

  # Call to h2o::h2o.automl
  spec <- auto_ml() |>
    set_engine("h2o", max_models = max_models, seed = 42,
               stopping_metric = metric, sort_metric = metric,
               fold_column = "foldid") |>
    set_mode("classification")
  fit <- workflow() |>
    add_model(spec) |>
    add_formula(Composer ~ .) |>
    fit(train_foldid)

  # Extract the leader model and predict on the test/unknown sets
  leader <- as.data.frame(fit$fit$fit$fit@leaderboard)$model_id[1]
  prob_all <- predict(fit, test_foldid, type = "prob")
  prob <- prob_all[[2]]
  unk_prob_all <- predict(fit, unk |> mutate(foldid = 0L), type = "prob")
  unk_prob <- unk_prob_all[[2]]

  # Results, with the h2o version so the stored run keeps its provenance
  res <- list(params = leader,
              class = switch((length(lvls) > 2) + 1,
                             factor(lvls[1 + (prob > 0.5)], levels = lvls),
                             factor(lvls[max.col(as.matrix(prob_all),
                                                 ties.method = "first")],
                                    levels = lvls)),
              prob = prob,
              unk_prob = unk_prob,
              h2o_version = as.character(packageVersion("h2o")))
  if (!is.null(cache_file)) saveRDS(res, file = cache_file)
  res

}

# Fit the three deterministic candidates on one preprocessed dataset, so the
# binary problem and the two screenings are compared on the same model set.
# AutoML is fitted apart by fit_automl_model(), in its own report chunk, so a
# failed search does not invalidate the cache of these fits.
fit_base_models <- function(train, test, folds) {

  list(
    LL = fit_glmnet_model(train = train, test = test, folds = folds,
                          mixture = 1, penalty_range = c(-8, 1)),
    RL = fit_glmnet_model(train = train, test = test, folds = folds,
                          mixture = 0, penalty_range = c(-6, 3)),
    RF = fit_rf_model(train = train, test = test, folds = folds)
  )

}

# Fit the deterministic candidates and the majority-class baseline of a
# screening, with the confusion matrix and per-class recalls of its ridge. The
# folds are returned so the report can show how rsample balanced the rare
# composers, and their foldid representation for the separate AutoML search.
fit_screen_models <- function(train, test, seed = 42) {

  set.seed(seed)
  folds <- vfold_cv(train, v = 10, strata = Composer)
  models <- fit_base_models(train = train, test = test, folds = folds)
  cm <- conf_mat(data = tibble(Composer = test$Composer,
                               est = models$RL$class),
                 truth = Composer, estimate = est)

  # Constant prediction at the majority class of the training set
  mc_class <- factor(rep(names(which.max(table(train$Composer))), nrow(test)),
                     levels = levels(train$Composer))

  list(models = models, folds = folds,
       foldid = make_foldid(folds, nrow(train)), cm = cm,
       recalls = diag(cm$table) / colSums(cm$table), mc_class = mc_class)

}

}

## Model comparison
{

# Macro-F1, with zero-filling for classes that are never predicted. Using the
# explicit lvls keeps the class set fixed even if a bootstrap resample omits a
# rare class. y_hat may be an n x B matrix of bootstrap replicates, scored one
# column at a time and returned as a vector of length B; y is then either the
# matching n x B matrix of truths or a single vector recycled down the columns.
f1 <- function(y, y_hat, lvls = levels(factor(y))) {

  # Class indices, NA outside lvls, exactly as factor(., levels = lvls) gives
  # them, and a plain vector of predictions as the one-column case
  y <- match(y, lvls)
  y_hat <- matrix(match(y_hat, lvls), nrow = NROW(y_hat), ncol = NCOL(y_hat))

  # Class-wise counts, dropping the pairs with an NA as table() does
  f1_sum <- 0
  for (k in seq_along(lvls)) {

    tp <- colSums(y_hat == k & y == k, na.rm = TRUE)
    fp <- colSums(y_hat == k & y != k, na.rm = TRUE)
    fn <- colSums(y_hat != k & y == k, na.rm = TRUE)

    # A class never predicted or never present gives 0 / 0: zero-filled
    prec <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
    rec <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
    f1_sum <- f1_sum + ifelse(prec + rec == 0, 0,
                              2 * prec * rec / (prec + rec))

  }

  drop(f1_sum / length(lvls))

}

# Percentile bootstrap for the differences in accuracy and macro-F1 of two
# models, resampling the test-set predictions without retraining. p-values
# test H1: M1 > M2.
boot_diff <- function(truth, pred1, pred2, idx, alpha = 0.05) {

  # Bootstrap replicates of the test-set predictions; resamples at rows
  B <- ncol(idx)
  lvls <- levels(truth)
  tc <- matrix(as.integer(truth)[idx], nrow = nrow(idx))
  p1 <- matrix(as.integer(factor(pred1, levels = lvls))[idx], nrow = nrow(idx))
  p2 <- matrix(as.integer(factor(pred2, levels = lvls))[idx], nrow = nrow(idx))

  # Bootstrap differences of accuracies and zero-filled macro-F1 scores
  d_acc <- colMeans(p1 == tc) - colMeans(p2 == tc)
  d_f1 <- f1(y = tc, y_hat = p1, lvls = seq_along(lvls)) -
    f1(y = tc, y_hat = p2, lvls = seq_along(lvls))

  # Order-statistic percentile CI and one-sided p-value
  ci <- function(d) sort(d)[c(ceiling((B + 1) * alpha / 2),
                              floor((B + 1) * (1 - alpha / 2)))]
  pv <- function(d) (1 + sum(d <= 0)) / (B + 1)
  ci_acc <- ci(d_acc)
  ci_f1 <- ci(d_f1)
  list(acc_ci_low = ci_acc[1], acc_ci_up = ci_acc[2], acc_p = pv(d_acc),
       f1_ci_low = ci_f1[1], f1_ci_up = ci_f1[2], f1_p = pv(d_f1))

}

# Test metrics and bootstrap comparisons of each model against the incumbent
# and the baseline, on shared resamples. AUC is filled only when probabilities
# are given.
compare_models <- function(truth, all_class, incumbent, baseline, idx,
                           all_prob = NULL) {

  # Metrics table: accuracy, macro-F1, and AUC if probabilities are given
  metrics_tab <- tibble(
    combo = names(all_class),
    Acc = map_dbl(all_class, ~ mean(.x == truth)),
    F1 = map_dbl(all_class, ~ f1(truth, .x))
  )
  if (!is.null(all_prob)) {

    # AUC with the second factor level as the event
    metrics_tab$AUC <- map_dbl(all_prob, function(p)
      yardstick::roc_auc_vec(truth = truth, estimate = drop(p),
                             event_level = "second"))
    metrics_tab$AUC[metrics_tab$combo == baseline] <- 0.5
    metrics_tab <- metrics_tab |>
      select(combo, Acc, AUC, F1)

  }

  # Bootstrap comparisons of each model against the incumbent and the baseline
  boot_tab <- map_dfr(names(all_class), function(m) {

    vs_rl <- if (m == incumbent) NULL else
      boot_diff(truth, all_class[[m]], all_class[[incumbent]], idx = idx)
    vs_mc <- if (m == baseline) NULL else
      boot_diff(truth, all_class[[m]], all_class[[baseline]], idx = idx)
    fmt_ci <- function(x, lo, up)
      if (is.null(x)) NA_character_ else
        sprintf("(%.4f, %.4f)", x[[lo]], x[[up]])
    fmt_p <- function(x, p) if (is.null(x)) NA_real_ else x[[p]]
    tibble(
      combo = m,
      ci_rl = fmt_ci(vs_rl, "acc_ci_low", "acc_ci_up"),
      p_rl = fmt_p(vs_rl, "acc_p"),
      ci_f1_rl = fmt_ci(vs_rl, "f1_ci_low", "f1_ci_up"),
      p_f1_rl = fmt_p(vs_rl, "f1_p"),
      ci_mc = fmt_ci(vs_mc, "acc_ci_low", "acc_ci_up"),
      p_mc = fmt_p(vs_mc, "acc_p"),
      ci_f1_mc = fmt_ci(vs_mc, "f1_ci_low", "f1_ci_up"),
      p_f1_mc = fmt_p(vs_mc, "f1_p")
    )

  })

  left_join(metrics_tab, boot_tab, by = "combo")

}

}

## Final screening fit
{

# Final multinomial ridge of a screening on all its known arias, with the
# class probabilities of its doubtful arias, their prediction CIs, and one
# contrast of class_j against each entry of contrasts. Exports the design,
# coefficients, and penalty for the coverage study of coverage_study.R.
final_screen <- function(all_clean, unk_clean, unk_ids, aria_labels,
                         contrasts, rds_file, class_j = "Baldassare Galuppi",
                         level = 0.95) {

  # Fold ids on all known arias for a reproducible cv.glmnet, seeded locally
  withr::local_seed(42)
  folds <- vfold_cv(all_clean, v = 10, strata = Composer)
  foldid <- make_foldid(folds, nrow(all_clean))
  X_all <- all_clean |>
    select(-Composer) |>
    as.matrix()
  cv <- cv.glmnet(x = X_all,
                  y = all_clean$Composer,
                  family = "multinomial",
                  type.measure = "deviance",
                  alpha = 0,
                  foldid = foldid,
                  parallel = TRUE,
                  lambda = exp(seq(3, -6, l = 300)))

  # Both the deviance minimizer and the selected penalty must be strictly
  # inside the grid
  stopifnot(cv$lambda.min > min(cv$lambda), cv$lambda.min < max(cv$lambda),
            cv$lambda.1se > min(cv$lambda), cv$lambda.1se < max(cv$lambda))
  lambda_th <- (nrow(X_all) - 1) * cv$lambda.1se / 2

  # rbind keeps the single-aria case a named 1 x K matrix, and is a no-op
  # otherwise.
  X_unk <- unk_clean |>
    select(-Composer) |>
    as.matrix()
  probs <- rbind(predict(cv, newx = X_unk, s = cv$lambda.1se,
                         type = "response")[, , 1])
  rownames(probs) <- aria_labels[as.character(unk_ids)]

  # Order the beta rows to match the columns of probs
  class_order <- colnames(probs)
  beta_list <- coef(cv, s = cv$lambda.1se)
  beta_mat <- t(sapply(class_order, function(cl) as.numeric(beta_list[[cl]])))
  ci <- mult_ridge_pred_ci(X = X_all, beta_mat = beta_mat,
                           lambda = lambda_th, new_X = X_unk,
                           level = level,
                           class_names = class_order) |>
    mutate(aria_id = as.character(unk_ids[.row]),
           aria = aria_labels[aria_id])

  # Sanity: machinery point probabilities match glmnet's
  stopifnot(max(abs(ci$p_hat - as.vector(t(probs)))) < 1e-6)

  diffs <- lapply(contrasts, function(class_k)
    mult_ridge_contrast_ci(X = X_all, beta_mat = beta_mat,
                           lambda = lambda_th, new_X = X_unk,
                           class_j = class_j, class_k = class_k,
                           level = level,
                           class_names = class_order) |>
      mutate(aria_id = as.character(unk_ids[.row])))

  # Gitignored, read back by coverage_study.R
  saveRDS(list(X = X_all, beta_mat = beta_mat,
               lambda_glmnet = cv$lambda.1se), file = rds_file)

  # Result
  list(cv = cv, X_all = X_all, probs = probs, class_order = class_order,
       beta_mat = beta_mat, ci = ci, diffs = diffs)

}

}

## Plots
{

# Print the logistic ridge coefficients as percent changes in odds, adding
# horizontal CIs.
logis_ridge_ci_plot <- function(ci_df) {

  ci_df <- ci_df |>
    mutate(
      pct_change = (exp(estimate) - 1) * 100,
      pct_low = (exp(ci_low) - 1) * 100,
      pct_up = (exp(ci_up) - 1) * 100
    ) |>
    mutate(
      term = factor(term,
                    levels = term[order(pct_change, decreasing = FALSE)])
    )
  ci_df |>
    select(term, pct_change, pct_low, pct_up, p_bh) |>
    arrange(desc(pct_change)) |>
    knitr::kable(digits = c(2, 2, 2, 2, 4)) |>
    print()

  ggplot(ci_df, aes(x = term, y = pct_change, ymin = pct_low, ymax = pct_up)) +
    geom_pointrange() +
    geom_hline(yintercept = 0, linetype = 2) +
    coord_flip() +
    labs(x = NULL,
         y = paste0("Percent change in odds of Perez over Galuppi\n",
                    "(per 1-SD increase in each standardized predictor)")
    ) +
    theme_bw()

}

# Print the multinomial ridge coefficients as percent changes in the
# probability ratio to the geometric-mean baseline (the sum-to-zero
# parametrization), adding horizontal CIs, one facet per class.
mult_ridge_ci_plot <- function(ci_df, ncol = 2) {

  # Order the coefficients within each class facet by pasting the class into
  # the level names, then stripping it from the axis labels
  ci_df <- ci_df |>
    mutate(
      pct_change = (exp(estimate) - 1) * 100,
      pct_low = (exp(ci_low) - 1) * 100,
      pct_up = (exp(ci_up) - 1) * 100,
      item = paste(class, term, sep = "___")
    ) |>
    arrange(pct_change) |>
    mutate(item = factor(item, levels = unique(item)))
  ci_df |>
    select(class, term, pct_change, pct_low, pct_up, p_bh) |>
    arrange(class, desc(pct_change)) |>
    knitr::kable(digits = c(2, 2, 2, 2, 2, 4)) |>
    print()

  ggplot(ci_df, aes(x = item, y = pct_change, ymin = pct_low,
                    ymax = pct_up)) +
    geom_pointrange() +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_x_discrete(labels = function(x) sub("^.*___", "", x)) +
    coord_flip() +
    facet_wrap(~class, scales = "free_y", ncol = ncol) +
    labs(x = NULL,
         y = paste0("Percent change in the probability ratio to the ",
                    "geometric-mean baseline\n(per 1-SD increase in each ",
                    "standardized predictor)")) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 7))

}

# Class-wise kernel density estimate of probabilities on [0, 1], each class
# scaled by its sample share so the curves add up to the mixture density.
kde_probs <- function(p, Composer, gridsize = 401) {

  # The probit transformation needs probabilities strictly inside (0, 1)
  stopifnot(all(p > 0 & p < 1))

  n <- length(p)
  bind_rows(lapply(split(p, Composer), function(pj) {

    # Throw a warning if not enough points to compute the kde
    if (length(unique(pj)) < 2) {

      warning("kde_probs: class with < 2 distinct probabilities dropped")
      return(NULL)

    }

    # ks::kde(unit.interval = TRUE) estimates the density of qnorm(p) and
    # back-transforms, with the plug-in bandwidth returned in h.
    fh <- ks::kde(x = pj, gridsize = gridsize, unit.interval = TRUE)
    tibble(p = fh$eval.points, h = fh$h,
           dens = pmax(fh$estimate, 0) * length(pj) / n)

  }), .id = "Composer")

}

# Class-wise density figure of the predicted Galuppi probabilities
# Faceted signed-bar chart of the filtered-out correlates of each selected
# marker: one facet per marker, bars of |rho| for its top correlates, fill by
# the sign of the correlation, and the 0.8 filter threshold marked. The blue
# and orange hues are colorblind-safe.
plot_marker_correlates <- function(corr_df, threshold = 0.8, ncol = 3) {

  # Order the bars within each facet by pasting the marker into the level
  # names, then stripping it from the axis labels
  corr_df <- corr_df |>
    mutate(sign = factor(ifelse(rho >= 0, "Positive", "Negative"),
                         levels = c("Positive", "Negative")),
           item = paste(marker, correlate, sep = "___")) |>
    arrange(abs(rho)) |>
    mutate(item = factor(item, levels = unique(item)))
  ggplot(corr_df, aes(x = item, y = abs(rho), fill = sign)) +
    geom_col() +
    geom_hline(yintercept = threshold, linetype = 2) +
    scale_x_discrete(labels = function(x) sub("^.*___", "", x)) +
    scale_fill_manual(values = c("Positive" = "#0072B2",
                                 "Negative" = "#E69F00"),
                      name = "Correlation sign") +
    coord_flip() +
    facet_wrap(~marker, scales = "free_y", ncol = ncol) +
    labs(x = NULL, y = "Absolute correlation with the selected marker") +
    theme_bw() +
    theme(legend.position = "top",
          axis.text.y = element_text(size = 7))

}


plot_probs_kde <- function(kde_df, rug_df, arias_df, xlab, facet = NULL,
                           labeller = "label_value",
                           legend_box = "vertical") {

  # Composer legend first and aria legend second, side by side
  # (legend_box = "horizontal") or stacked (legend_box = "vertical")
  pl <- ggplot(kde_df, aes(x = p, y = dens, fill = Composer,
                           colour = Composer)) +
    geom_area(alpha = 0.3, position = "identity") +
    geom_line() +
    geom_rug(data = rug_df, aes(x = p_galuppi, colour = Composer),
             inherit.aes = FALSE, sides = "b", length = unit(0.03, "npc"),
             alpha = 0.7) +
    geom_vline(data = arias_df, aes(xintercept = p_galuppi, linetype = Aria),
               inherit.aes = FALSE) +
    labs(x = xlab, y = "Density (test arias)") +
    guides(fill = guide_legend(order = 1), colour = guide_legend(order = 1),
           linetype = guide_legend(order = 2)) +
    theme_bw() +
    theme(legend.position = "bottom", legend.box = legend_box)

  # The screening figure gets one panel per problem, each on its own y scale
  if (!is.null(facet)) {

    pl <- pl + facet_wrap(facet, ncol = 1, scales = "free_y",
                          labeller = labeller)

  }

  pl

}

}

## Paper numbers and macros
{

# Metric and bootstrap-comparison rows of paper_numbers for one problem, named
# "<prefix><stem><combo>". cols selects and orders the columns of a
# compare_models() table; their values pass through unconverted, so one call
# must not mix the numeric metrics with the character CI strings.
comparison_numbers <- function(tab, cols, prefix = "") {

  stems <- c(Acc = "acc_", AUC = "auc_", F1 = "f1_",
             p_rl = "p_rl_", p_f1_rl = "p_f1_rl_",
             p_mc = "p_mc_", p_f1_mc = "p_f1_mc_",
             ci_rl = "ci_rl_", ci_f1_rl = "ci_f1_rl_",
             ci_mc = "ci_mc_", ci_f1_mc = "ci_f1_mc_")
  combo <- gsub("-", "_", tab$combo)
  bind_rows(lapply(cols, function(cl)
    tibble(name = paste0(prefix, stems[[cl]], combo), value = tab[[cl]])))

}

# Class-probability and prediction-CI rows of paper_numbers for one screening,
# named "<prefix>_p_", "<prefix>_cilow_", "<prefix>_ciup_" plus the class and
# the aria id.
screen_numbers <- function(probs, ci, unk_ids, prefix) {

  alnum <- function(x) gsub("[^A-Za-z0-9]", "", x)
  bind_rows(
    as.data.frame(probs) |>
      mutate(aria_id = as.character(unk_ids)) |>
      pivot_longer(-aria_id, names_to = "composer", values_to = "prob") |>
      transmute(name = paste0(prefix, "_p_", alnum(composer), "_", aria_id),
                value = prob),
    ci |>
      transmute(name = paste0(prefix, "_cilow_", alnum(class), "_", aria_id),
                value = p_ci_low),
    ci |>
      transmute(name = paste0(prefix, "_ciup_", alnum(class), "_", aria_id),
                value = p_ci_up)
  )

}

# One \gp<Name> macro per quantity for \input{paper_numbers.tex}, so the
# manuscript numbers update with every render. Macro names must be purely
# alphabetic, hence the spelled-out aria ids and digits.
tex_name <- function(x) {

  x <- gsub("1239", "_nonsarei", x)
  x <- gsub("1241", "_sperai", x)
  x <- gsub("1242", "_tusai", x)
  x <- gsub("1se", "onese", x)
  x <- gsub("f1", "fone", x)
  # Marker names carry interval and degree numbers
  digit_words <- c("0" = "Zero", "1" = "One", "2" = "Two", "3" = "Three",
                   "4" = "Four", "5" = "Five", "6" = "Six", "7" = "Seven",
                   "8" = "Eight", "9" = "Nine")
  for (d in names(digit_words)) x <- gsub(d, digit_words[[d]], x, fixed = TRUE)
  parts <- strsplit(x, "_", fixed = TRUE)[[1]]
  parts <- parts[parts != ""]
  out <- paste0("gp", paste0(toupper(substring(parts, 1, 1)),
                             substring(parts, 2), collapse = ""))
  stopifnot(grepl("^[A-Za-z]+$", out))
  out

}

# Formats as reported in the paper: counts with the math-safe {,} separator,
# metrics and p-values to four decimals, with p-values below 1e-4 (possible
# because of the 1 / (B + 1) floor) reported as {<}0.0001, recalls to two,
# the rest to three; CI strings pass through unchanged.
tex_value <- function(name, value) {

  num <- suppressWarnings(as.numeric(value))
  if (is.na(num)) return(value)
  if (grepl("^pctl_", name)) {
    k <- round(num)
    suf <- if (k %% 100 %in% 11:13) "th" else
      switch(as.character(k %% 10), "1" = "st", "2" = "nd", "3" = "rd", "th")
    return(paste0("$", k, "$", suf))
  }
  if (grepl(paste0("^(n_|multi_n_|ale_n_|p_train$|p_all$|",
                   "multi_p_train$|multi_p_all$|ale_p_train$|ale_p_all$|",
                   "cov_reps$)"), name)) {
    return(format(round(num), big.mark = "{,}", scientific = FALSE))
  }
  if (grepl("^(acc_|auc_|f1_|multi_acc|multi_f1|ale_acc|ale_f1)", name)) {
    return(sprintf("%.4f", num))
  }
  if (grepl("(^|_)p_(rl|mc|f1)", name)) {
    return(switch((num < 1e-4) + 1, sprintf("%.4f", num), "{<}0.0001"))
  }
  if (grepl("^(multi|ale)_recall_", name)) return(sprintf("%.2f", num))
  if (name %in% c("max_pct_change", "null_expected") || grepl("_pct$", name) ||
      grepl("^(cov_nsel_|ex_dyn_)", name))
    return(sprintf("%.1f", num))
  sprintf("%.3f", num)

}

# One \gp<Name> macro per row of paper_numbers, written to file for
# \input{paper_numbers.tex} and returned for inspection.
write_paper_macros <- function(paper_numbers, file) {

  macros <- paper_numbers |>
    rowwise() |>
    mutate(macro = sprintf("\\newcommand{\\%s}{%s}", tex_name(name),
                           tex_value(name, value))) |>
    pull(macro)
  writeLines(macros, con = file)
  macros

}

}
