
# Coverage study of the ridge inference in methods.R, at the settings of the
# three reported problems (binary, demofoonte with K = 5, alessandro with
# K = 6). Each problem is run in three phases:
#
#  1. The population ridge coefficient beta_lambda is computed exactly, by
#     fitting glmnet to the true class probabilities. This, and not the
#     data-generating beta, is the estimand the intervals are scored against.
#  2. The design is held fixed and only the responses are resampled from the
#     true class probabilities, R_cov times; each replicate is refitted and
#     handed to the sandwich covariance of methods.R.
#  3. The marginal Wald intervals of all the slopes are scored against
#     beta_lambda at the 95%, 90%, and 80% nominal levels, and the standard
#     errors are compared with the Monte Carlo spread of the estimates.
#
# It runs on the designs, coefficients, and penalties that galuppi.qmd exports
# to coverage_{binary,multi,ale}.rds. Its output is the appendix table, read
# by galuppi.qmd into the \gpCov macros.

# Load required methods
library(glmnet)
library(foreach)
library(future)
library(doFuture)
library(progressr)
source("methods.R")

## Setup
{

# Monte Carlo replicates per problem, and cores used
R_cov <- 1000
n_cores <- 10

# Nominal confidence levels at which the marginal coverage is scored, on the
# same simulations, in the order of the appendix table
conf <- c(0.95, 0.90, 0.80)

# Reproducible parallel replicates
set.seed(5)
plan(multisession, workers = n_cores)

# Progress bar
options(progressr.enable = TRUE)
handlers(handler_progress(
  format = ":spin [:bar] :percent Total: :elapsedfull End \u2248 :eta",
  clear = FALSE))

# Row-wise softmax, shifted by the row maximum for numerical stability
softmax <- function(E) {

  E <- E - apply(E, 1, max)
  ex <- exp(E)
  ex / rowSums(ex)

}

# The three reported problems and the exports they are run on; n, p, and the
# penalty are read off each design below, not declared here.
cfg <- data.frame(
  problem = c("binary", "demofoonte", "alessandro"),
  K = c(2, 5, 6),
  rds = c("coverage_binary.rds", "coverage_multi.rds", "coverage_ale.rds")
)

# Check the designs are present
if (!all(file.exists(cfg$rds))) {

  stop("Missing coverage designs (",
       paste(cfg$rds[!file.exists(cfg$rds)], collapse = ", "),
       "), render galuppi.qmd first.")

}

}

## Coverage study
{

cov_rows <- with_progress({

# One bar over all refits
prog <- progressor(steps = nrow(cfg) * R_cov)

lapply(seq_len(nrow(cfg)), function(i) {

  ## Design and data-generating truth

  # Design, penalty, and coefficients of the reported fit, exported by the qmd;
  # the coefficients are the truth the responses are resampled from.
  K <- cfg$K[i]
  obj <- readRDS(cfg$rds[i])
  X <- as.matrix(obj$X)
  lambda_glmnet <- obj$lambda_glmnet
  if (K == 2) beta_true <- as.numeric(obj$beta) else
    beta_mat_true <- obj$beta_mat

  # Design with intercept, and the theory-scale penalty of methods.R
  n <- nrow(X)
  p <- ncol(X)
  Xt <- cbind(1, X)
  lambda_th <- (n - 1) * lambda_glmnet / 2

  # Marginal critical values at each conf level
  z_marg <- qnorm(1 - (1 - conf) / 2)
  message(sprintf("%s: K = %d, n = %d, p = %d, lambda = %.4g",
                  cfg$problem[i], K, n, p, lambda_glmnet))

  if (K == 2) {

    # Binary problem, with the binomial machinery deployed in the analysis
    pi_true <- drop(1 / (1 + exp(-Xt %*% beta_true)))

    ## Population ridge target beta_lambda

    # The estimand is not beta_true but the population ridge coefficient
    #
    #   beta_lambda = argmin_b { E[-loglik(b; Y, X)] / n + lambda * ||b||^2 }.
    #
    # It is obtained exactly, with no Monte Carlo error, because the expected
    # deviance equals the deviance evaluated at the true class probabilities
    # used as fractional responses: handing cbind(1 - pi_true, pi_true) to
    # glmnet as a two-column response minimizes precisely that criterion.
    fit_t <- glmnet(X, cbind(1 - pi_true, pi_true), family = "binomial",
                    alpha = 0, lambda = lambda_glmnet)
    target <- as.numeric(coef(fit_t, s = lambda_glmnet))[-1]

    ## Monte Carlo replicates

    # Fixed design, responses resampled from pi_true. One future per replicate
    # (chunk.size = 1), as progress is only relayed when a future resolves.
    reps <- foreach(r = seq_len(R_cov),
                    .options.future = list(seed = TRUE,
                                           chunk.size = 1)) %dofuture% {

      prog()
      y <- factor(rbinom(n, 1, pi_true), levels = 0:1)

      # Degenerate resamples have no fit and are dropped from R_eff below
      if (nlevels(droplevels(y)) < 2) NULL else {

        beta_hat <- as.numeric(coef(glmnet(X, y, family = "binomial",
                                           alpha = 0, lambda = lambda_glmnet),
                                    s = lambda_glmnet))
        ci_r <- logis_ridge_ci(X, beta_hat = beta_hat, lambda = lambda_th)
        list(est = ci_r$estimate, se = ci_r$std_error)

      }

    }

  } else {

    # Screenings, with the multinomial machinery deployed in the analysis
    pi_true <- softmax(Xt %*% t(beta_mat_true))
    class_names <- paste0("c", seq_len(K))

    ## Population ridge target beta_lambda

    # As above, with pi_true as the fractional response matrix. Coefficients are
    # class-major, as.vector(t(beta_mat)); keep drops the K intercepts, which
    # are not reported and are identified only up to a common shift.
    fit_t <- glmnet(X, pi_true, family = "multinomial", alpha = 0,
                    lambda = lambda_glmnet, type.multinomial = "ungrouped")
    beta_list_true <- coef(fit_t, s = lambda_glmnet)
    beta_mat_target <-
      t(sapply(seq_len(K), function(k) as.numeric(beta_list_true[[k]])))
    keep <- rep(c(FALSE, rep(TRUE, p)), times = K)
    target <- as.vector(t(beta_mat_target))[keep]

    ## Monte Carlo replicates

    # Same, with each label drawn from its own row of pi_true;
    # mult_ridge_ci() returns the slopes in class-major order, matching target.
    reps <- foreach(r = seq_len(R_cov),
                    .options.future = list(seed = TRUE,
                                           chunk.size = 1)) %dofuture% {

      prog()
      y <- factor(class_names[apply(pi_true, 1,
                                    function(q) sample(K, 1, prob = q))],
                  levels = class_names)

      # Degenerate resamples, and rare-class refits on which glmnet fails
      # to converge at the single lambda, have no fit and are dropped from
      # R_eff below
      if (nlevels(droplevels(y)) < K) NULL else tryCatch({

        fit <- glmnet(X, y, family = "multinomial", alpha = 0,
                      lambda = lambda_glmnet, type.multinomial = "ungrouped")
        beta_list <- coef(fit, s = lambda_glmnet)
        beta_mat_hat <-
          t(sapply(seq_len(K), function(k) as.numeric(beta_list[[k]])))
        ci_r <- mult_ridge_ci(X, beta_mat = beta_mat_hat, lambda = lambda_th)
        list(est = ci_r$estimate, se = ci_r$std_error)

      }, error = function(e) NULL)

    }

  }

  ## Coverage summaries

  # Drop degenerate resamples and stack; R_eff, not R_cov, is the reported R
  reps <- Filter(Negate(is.null), reps)
  R_eff <- length(reps)
  est <- do.call(rbind, lapply(reps, `[[`, "est"))
  se <- do.call(rbind, lapply(reps, `[[`, "se"))

  # cov_marginal_<level> covers entrywise at each nominal level, reusing the
  # same replicates; se_ratio is the median over the slopes of the mean
  # standard error over the Monte Carlo standard deviation, above one meaning
  # conservative.
  tgt <- rep(target, each = R_eff)
  cov_marg <- sapply(z_marg, function(z) mean(abs(est - tgt) <= z * se))
  names(cov_marg) <- paste0("cov_marginal_", round(100 * conf))
  se_ratio <- median(colMeans(se) / apply(est, 2, sd))
  data.frame(problem = cfg$problem[i], K = K, n = n, p = p, R = R_eff,
             t(round(cov_marg, 4)), se_ratio = round(se_ratio, 3))

})

})

# Release the workers
plan(sequential)

coverage_tab <- do.call(rbind, cov_rows)
print(coverage_tab)
write.csv(coverage_tab, file = "coverage.csv", row.names = FALSE)

}
