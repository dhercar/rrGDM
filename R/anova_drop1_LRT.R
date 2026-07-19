#' Wald Tests for Environmental Loadings and Dispersion Effects
#'
#' @description
# 'Wald tests for sets of predictor effects in a fitted \code{rrGDM} model:
#' \itemize{
#'   \item \strong{Environmental loadings} (\code{lambda}): for each
#'     predictor, a multivariate (multi-df) Wald test of the null that all
#'     free loadings of that predictor (one per latent dimension it loads
#'     on) are jointly zero. This accounts for the covariance between
#'     loadings within a predictor, which can be substantial due to
#'     rotational invariance.
#'   \item \strong{Dispersion effects} (\code{beta}): for each predictor in
#'     the dispersion model, a univariate Wald z-test, since each
#'     coefficient is a scalar effect on \eqn{\log(\tau_i)}.
#' }
#'
#' @param object A fitted model of class \code{rrGDM}, fitted with
#'   \code{se = TRUE}.
#' @param include_intercept Logical. Should the intercept term of
#'   \code{formula_sd} be included in the \code{beta} output? Default is
#'   \code{FALSE}.
#'
#' @details
#' This test is preferred over testing each latent dimension separately
#' (e.g., with per-dimension p-values corrected for multiple comparisons),
#' because it explicitly accounts for the covariance between loadings
#' within a predictor rather than discarding it. If
#' \eqn{\boldsymbol{\Sigma}_i} is computationally singular (which can occur
#' under strong collinearity between latent dimensions for a given
#' predictor), a Moore-Penrose pseudoinverse is used instead, with a warning.
#'
#' The degrees of freedom of the \code{lambda} test depend on
#' \code{object$constraint}:
#' \itemize{
#'   \item Under \code{constraint == 1}, the loadings matrix is constrained
#'     to be lower triangular: predictor (row) \code{i} only has free
#'     loadings in latent dimensions \code{1:min(i, k)}, the remainder being
#'     structural zeros that are fixed (not estimated) when the model is
#'     fit. The Wald test for predictor \code{i} therefore uses
#'     \code{Df = min(i, k)}, and only the free loadings (and their
#'     covariance sub-block) enter the test statistic.
#'   \item Under \code{constraint \%in\% c(0, 2)}, \code{lambda} is
#'     unconstrained, so every predictor is tested with \code{Df = k}.
#'     Note that under \code{constraint == 0} the model is not identified
#'     under rotation of the latent space, so individual loadings (and therefore this test)
#'     should be interpreted with caution.
#' }
#'
#' @return
#' A list of class \code{wald_rrGDM} with two elements:
#' \itemize{
#'   \item \code{lambda}: a \code{data.frame} with one row per predictor in
#'     \code{X_design}, with columns \code{term}, \code{Wald}, \code{Df}
#'     (see Details -- \code{k} under \code{constraint \%in\% c(0, 2)},
#'     \code{min(i, k)} for predictor \code{i} under \code{constraint == 1}),
#'     \code{Pr_Chisq}, \code{Sig}.
#'   \item \code{beta}: a \code{data.frame} with one row per predictor in
#'     \code{W_design} (excluding the intercept unless
#'     \code{include_intercept = TRUE}), with columns \code{term},
#'     \code{Estimate}, \code{SE}, \code{z}, \code{Pr_z}, \code{Sig}.
#' }
#'
#' @export
#' @importFrom stats pchisq pnorm
wald_test <- function(object, include_intercept = FALSE) {

  if (!inherits(object, "rrGDM")) {
    cli::cli_abort("{.arg object} must be a fitted {.cls rrGDM} model.")
  }

  if (is.null(object$sd_rep)) {
    cli::cli_abort(c(
      "x" = "{.arg object} has no {.field sd_rep}.",
      "i" = "Refit the model with {.code se = TRUE}, or run {.fn include_se} on the fitted object."
    ))
  }

  if (is.null(object$constraint)) {
    cli::cli_abort(c(
      "x" = "{.arg object} has no {.field constraint} recorded.",
      "i" = "Refit the model with the current version of {.fn rrGDM}, which stores {.field constraint} on the fitted object.",
      "i" = "The degrees of freedom for the {.field lambda} Wald test depend on which identifiability constraint was used."
    ))
  }

  constraint <- object$constraint

  sd_rep <- object$sd_rep
  sdr_vals <- sd_rep$value
  sdr_cov  <- sd_rep$cov

  # LOADINGS (lambda)

  n_env <- ncol(object$X_design)

  if (is.null(n_env) || n_env == 0) {
    cli::cli_abort(c(
      "x" = "The model contains no predictors in {.field X_design}.",
      "i" = "{.fn wald_test} requires at least one predictor with an estimated loading."
    ))
  }

  k  <- object$k
  env_nms <- colnames(object$X_design)

  lam_idx <- which(names(sdr_vals) == "lambda")

  if (length(lam_idx) != n_env * k) {
    cli::cli_abort(c(
      "x" = "Could not locate the expected number of {.field lambda} entries in {.field sd_rep}.",
      "i" = "Expected {n_env * k}, found {length(lam_idx)}.",
      "!" = "This may indicate a mismatch between the fitted model and {.arg object}."
    ))
  }

  lambda_hat <- matrix(sdr_vals[lam_idx], nrow = n_env, ncol = k)
  Sigma_full <- sdr_cov[lam_idx, lam_idx, drop = FALSE]

  report_lambda <- tryCatch(object$obj$report(object$opt$par)$lambda, error = function(e) NULL)
  if (!is.null(report_lambda) &&
      all(dim(report_lambda) == dim(lambda_hat)) &&
      !isTRUE(all.equal(unname(lambda_hat), unname(report_lambda), tolerance = 1e-6))) {
    cli::cli_warn(c(
      "!" = "ADREPORT-derived {.field lambda} does not match {.code object$obj$report()$lambda}.",
      "i" = "Results may be misaligned; proceed with caution."
    ))
  }

  if (constraint == 0) {
    cli::cli_warn(c(
      "!" = "Model was fit with {.code constraint = 0} (no identifiability constraint).",
      "i" = "Individual loadings in {.field lambda} are not identified under rotation of the latent space; interpret per-predictor Wald tests with caution."
    ))
  }

  # excluded from both the test statistic and its degrees of freedom.
  get_free_dims <- function(i) {
    if (constraint == 1) seq_len(min(i, k)) else seq_len(k)
  }

  singular_terms <- character(0)

  lambda_res <- lapply(seq_len(n_env), function(i) {
    fd <- get_free_dims(i)
    zd <- setdiff(seq_len(k), fd)  # dimensions expected to be structural zeros for this predictor

    if (length(zd) > 0 &&
        !isTRUE(all.equal(unname(lambda_hat[i, zd]), rep(0, length(zd)), tolerance = 1e-8))) {
      cli::cli_abort(c(
        "x" = "Unexpected non-zero loading found where a structural zero was expected for {.val {env_nms[i]}}.",
        "i" = "This suggests {.field object$constraint} does not match how the model was actually fitted."
      ))
    }

    df_i  <- length(fd)
    idx_i  <- i + (fd - 1) * n_env
    lam_i  <- lambda_hat[i, fd]
    Sigma_i <- Sigma_full[idx_i, idx_i, drop = FALSE]

    Sigma_inv <- tryCatch(
      solve(Sigma_i),
      error = function(e) {
        singular_terms <<- c(singular_terms, env_nms[i])
        MASS::ginv(Sigma_i)
      }
    )

    W  <- as.numeric(t(lam_i) %*% Sigma_inv %*% lam_i)
    pv <- stats::pchisq(W, df = df_i, lower.tail = FALSE)

    data.frame(term = env_nms[i], Wald = W, Df = df_i, Pr_Chisq = pv,
               stringsAsFactors = FALSE)
  })

  lambda_df <- do.call(rbind, lambda_res)

  if (length(singular_terms) > 0) {
    cli::cli_warn(c(
      "!" = "Covariance matrix was computationally singular for loading term{?s}: {.val {singular_terms}}.",
      "i" = "A Moore-Penrose pseudoinverse was used instead; interpret {.field Wald}/{.field Pr_Chisq} for {?this/these} term{?s} with caution."
    ))
  }

  lambda_df$Sig <- cut(
    lambda_df$Pr_Chisq,
    breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
    labels = c("***", "**", "*", ".", " ")
  )
  rownames(lambda_df) <- NULL

  #  DISPERSION (beta)

  w_nms <- colnames(object$W_design)

  if (is.null(w_nms) || length(w_nms) == 0) {
    cli::cli_abort("The model contains no terms in {.field W_design}.")
  }

  beta_idx <- which(names(sdr_vals) == "beta")

  if (length(beta_idx) != length(w_nms)) {
    cli::cli_abort(c(
      "x" = "Could not locate the expected number of {.field beta} entries in {.field sd_rep}.",
      "i" = "Expected {length(w_nms)}, found {length(beta_idx)}.",
      "!" = "This may indicate a mismatch between the fitted model and {.arg object}."
    ))
  }

  beta_hat <- sdr_vals[beta_idx]
  se_hat   <- sqrt(diag(sdr_cov[beta_idx, beta_idx, drop = FALSE]))

  # Sanity check alignment against obj$report()$beta, evaluated explicitly
  # at the fitted parameters (object$opt$par) -- see note above re: obj's
  # mutable internal state.
  report_beta <- tryCatch(object$obj$report(object$opt$par)$beta, error = function(e) NULL)
  if (!is.null(report_beta) && length(report_beta) == length(beta_hat) &&
      !isTRUE(all.equal(unname(beta_hat), unname(report_beta), tolerance = 1e-6))) {
    cli::cli_warn(c(
      "!" = "ADREPORT-derived {.field beta} does not match {.code object$obj$report()$beta}.",
      "i" = "Results may be misaligned."
    ))
  }

  z    <- beta_hat / se_hat
  pval <- 2 * stats::pnorm(-abs(z))

  beta_df <- data.frame(
    term     = w_nms,
    Estimate = beta_hat,
    SE       = se_hat,
    z        = z,
    Pr_z     = pval,
    stringsAsFactors = FALSE
  )

  if (!include_intercept) {
    beta_df <- beta_df[w_nms != "(Intercept)", , drop = FALSE]
  }

  beta_df$Sig <- cut(
    beta_df$Pr_z,
    breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
    labels = c("***", "**", "*", ".", " ")
  )
  rownames(beta_df) <- NULL

  #  RESULT

  res <- list(lambda = lambda_df, beta = beta_df)
  class(res) <- "wald_rrGDM"
  res
}

#' Print method for wald_rrGDM objects
#'
#' @param x An object of class \code{wald_rrGDM}, as returned by
#'   \code{\link{wald_test}}.
#' @param ... Currently unused.
#' @export
#' @method print wald_rrGDM
print.wald_rrGDM <- function(x, ...) {
  cli::cli_h2("Wald tests: environmental loadings")
  print(x$lambda, row.names = FALSE)
  cli::cli_h2("Wald tests: dispersion effects")
  print(x$beta, row.names = FALSE)
  invisible(x)
}
