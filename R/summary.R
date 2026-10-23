#' Summary method for rrGDM models
#'
#' @description
#' Summarises a fitted \code{rrGDM} object. The summary reports:
#' \itemize{
#'   \item the model call and the (approximate marginal) log-likelihood;
#'   \item for each predictor in the mean model, the effect size
#'     (\code{Size}, the modulus of its loadings across the latent
#'     dimensions) and a multivariate Wald test (see \code{\link{wald_test}});
#'   \item the intercept and slopes of the dispersion model (effects on
#'     \eqn{\log(\tau)}, where \eqn{\tau} is the sample-level SD), with
#'     univariate Wald z-tests;
#'   \item the standard deviation of additional random effects, if any;
#'   \item the observation-level (pairwise) parameter of the response
#'     distribution, if the family has one (\code{phi} for beta and
#'     beta-binomial, \code{SD} for Gaussian)
#'   \item the conditional stress and corR2 of the model.
#' }
#'
#' @param object A fitted model of class \code{rrGDM}.
#' @param ... Currently unused.
#'
#' @details
#' The log-likelihood is the approximate marginal log-likelihood obtained
#' with the Laplace approximation, and \code{df} is the number of
#' (non-mapped) fixed parameters estimated.
#'
#' Two stress values are reported, both computed from the conditional latent
#' coordinates (see \code{\link{stress}}): the metric stress on the
#' dissimilarity (response) scale, and Kruskal's non-metric stress on the
#' latent distances (link scale), which is the one comparable with NMDS software.
#' The corR2 (see \code{\link{cor_r2}}) is the conditional squared correlation
#' between observed and fitted dissimilarities. Both indicate how well the
#' latent space captures the dissimilarities, but not the explanatory power of
#' the predictors.
#'
#' If the model was fitted without standard errors (\code{se = FALSE}), the
#' Hessian at the optimum is not positive definite, or the covariance matrix of
#' the estimates is missing, a warning is issued and only point estimates are
#' reported.
#'
#' If the stress or corR2 cannot be computed, a warning is given and the value
#' is \code{NA}.
#'
#' @return
#' A list of class \code{summary.rrGDM}, with elements including
#' \code{call}, \code{loglik}, \code{df}, \code{aic}, \code{lambda} (see
#' \code{\link{wald_test}}), \code{beta}, \code{tau0} (sample-level SD when all
#' dispersion predictors are at their reference value), \code{sigma_re},
#' \code{scale_par}, \code{stress_metric}, \code{stress_nonmetric} and
#' \code{cor_r2}.
#'
#' @seealso \code{\link{wald_test}}, \code{\link{stress}}, \code{\link{cor_r2}}
#'
#' @method summary rrGDM
#' @export
summary.rrGDM <- function(object, ...) {

  if (!inherits(object, "rrGDM")) {
    cli::cli_abort("{.arg object} must be a fitted {.cls rrGDM} model.")
  }

  report <- tryCatch(get_report(object), error = function(e) {
    cli::cli_abort(c(
      "x" = "Could not extract model report from {.arg object$obj}.",
      "i" = "Original error: {conditionMessage(e)}"
    ))
  })

  sd_rep <- object$sd_rep
  has_sd_rep <- !is.null(sd_rep)
  pd_ok <- has_sd_rep && isTRUE(hess_ok(object))
  has_cov <- has_sd_rep && is.matrix(sd_rep$cov)
  # no SEs / Wald tests if there is no sd_rep, the Hessian is not usable or the covariance is missing
  has_se <- pd_ok && has_cov
  has_sample_re <- !is.null(report$u_sample)

  n_env <- ncol(object$X_design)
  if (is.null(n_env)) n_env <- 0L

  if (!has_sd_rep) {
    cli::cli_warn(c(
      "!" = "{.arg object} has no {.field sd_rep}: standard errors and Wald tests are not available.",
      "i" = "Refit the model with {.code se = TRUE}, or run {.fn include_se} on the fitted object."
    ))
  } else if (!pd_ok) {
    cli::cli_warn(c(
      "!" = "The Hessian at the optimum is not positive definite (or is nearly singular): standard errors and Wald tests are omitted.",
      "i" = "Only point estimates are reported."
    ))
  } else if (!has_cov) {
    cli::cli_warn(c(
      "!" = "The covariance matrix of the {.field sd_rep} is missing: standard errors and Wald tests are omitted.",
      "i" = "Run {.fn include_se} on the fitted object without {.code getReportCovariance = FALSE}."
    ))
  }

  # LOCATION AND DISPERSION EFFECTS
  if (has_se) {
    wald <- wald_test(object, include_intercept = TRUE)
    lambda_df <- wald$lambda
    beta_df <- wald$beta

  } else {
    lambda_df <- data.frame(
      term = if (n_env > 0) colnames(object$X_design) else character(0),
      Size = as.numeric(report$lambda_mod),
      Wald = rep(NA_real_, n_env),
      Df = rep(NA_integer_, n_env),
      Pr_Chisq = rep(NA_real_, n_env),
      Sig = rep(NA_character_, n_env),
      stringsAsFactors = FALSE)

    beta_df <- data.frame(
      term = colnames(object$W_design),
      Estimate = as.numeric(report$beta),
      SE = NA_real_,
      z = NA_real_,
      Pr_z = NA_real_,
      Sig = NA_character_,
      stringsAsFactors = FALSE)
  }

  # No dispersion model if sample-level errors are not included
  if (!has_sample_re) beta_df <- NULL

  # Sample-level SD at baseline (all dispersion predictors at 0 / reference level)
  tau0 <- NA_real_
  if (!is.null(beta_df)) {
    b0 <- beta_df$Estimate[beta_df$term == "(Intercept)"]
    if (length(b0) == 1) tau0 <- exp(b0)
  }

  # RANDOM EFFECTS (SD)
  re_df <- NULL
  if (length(object$re_vars) > 0) {
    sigma_re <- as.numeric(report$sigma_re)

    re_df <- data.frame(
      term = object$re_vars,
      SD = sigma_re,
      stringsAsFactors = FALSE)
  }

  # OBSERVATION-LEVEL (PAIRWISE) PARAMETER
  scale_name <- switch(object$family,
                       "gaussian" = "SD",
                       "beta" = "phi",
                       "beta.binomial" = "phi",
                       NULL)
  scale_df <- NULL

  if (!is.null(scale_name)) {
    log_scale <- object$opt$par[names(object$opt$par) == "log_scale"]
    scale_df <- data.frame(
      term = scale_name,
      Estimate = if (length(log_scale) == 1) exp(as.numeric(log_scale)) else NA_real_,
      stringsAsFactors = FALSE)
  }

  # Likelihood
  loglik <- -object$opt$objective
  df <- length(object$opt$par)
  aic <- 2 * df - 2 * loglik

  #  FIT OF THE LATENT SPACE (i.e., conditional fit)
  safe <- function(expr, what) {
    tryCatch(expr, error = function(e) {
      cli::cli_warn("Could not compute {what}: {conditionMessage(e)}")
      NA_real_
    })
  }
  stress_metric <- safe(stress(object, type = "metric", scale = "response"), "the metric stress")
  stress_nonmetric <- safe(stress(object, type = "nonmetric", scale = "link"), "the non-metric stress")
  r2 <- safe(cor_r2(object, message = FALSE), "corR2")

  res <- list(
    call = object$call,
    family = object$family,
    link = object$link,
    k = object$k,
    n_samples = nrow(object$X_design),
    n_pairs = length(object$Y_diss),
    loglik = loglik,
    df = df,
    aic = aic,
    lambda = lambda_df,
    beta = beta_df,
    tau0 = tau0,
    sigma_re = re_df,
    scale_par = scale_df,
    stress_metric = stress_metric,
    stress_nonmetric = stress_nonmetric,
    cor_r2 = r2,
    has_se = has_se)

  class(res) <- "summary.rrGDM"
  res
}

#' Print method for summary.rrGDM objects
#'
#' @param x An object of class \code{summary.rrGDM}, as returned by
#'   \code{\link{summary.rrGDM}}.
#' @param digits Number of significant digits for test statistics,
#'   estimates and effect sizes. Default is 4.
#' @param ... Currently unused.
#'
#' @return Invisibly returns \code{x}.
#' @export
#' @method print summary.rrGDM
print.summary.rrGDM <- function(x, digits = 4, ...) {
  fmt_num <- function(v) formatC(v, digits = digits, format = "fg", flag = "#")

  cat("Call:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")

  cat(sprintf("Family: %s (link: %s); %d latent dimension%s; %d samples (%d pairs)\n",
              x$family, x$link, x$k, if (x$k == 1) "" else "s", x$n_samples, x$n_pairs))
  cat(sprintf("Log-likelihood: %s (df = %d); AIC: %s\n",
              formatC(x$loglik, digits = 2, format = "f"), x$df,
              formatC(x$aic, digits = 2, format = "f")))

  # Location effects
  cli::cli_h2("Effects on composition (location)")
  if (nrow(x$lambda) > 0) {
    d <- x$lambda
    d$Size <- fmt_num(d$Size)
    if (x$has_se) {
      d$Wald <- round(d$Wald, 3)
      d$Pr_Chisq <- format_pval(d$Pr_Chisq)
      d$Sig[is.na(d$Sig)] <- " "   # NA is printed as blank
    } else {
      d <- d[, c("term", "Size"), drop = FALSE]
    }
    print(d, row.names = FALSE, right = TRUE)
  } else {
    cli::cli_alert_info("No predictors in the mean model.")
  }

  #  Dispersion effects
  cli::cli_h2("Effects on dispersion (log of sample-level SD)")
  if (!is.null(x$beta) && nrow(x$beta) > 0) {
    d <- x$beta
    d$Estimate <- fmt_num(d$Estimate)
    if (x$has_se) {
      d$SE <- fmt_num(d$SE)
      d$z <- fmt_num(d$z)
      d$Pr_z <- format_pval(d$Pr_z)
      d$Sig[is.na(d$Sig)] <- " "
    } else {
      d <- d[, c("term", "Estimate"), drop = FALSE]
    }
    print(d, row.names = FALSE, right = TRUE)
    cat(sprintf("\nSample-level SD at baseline: tau = exp(Intercept) = %s\n", fmt_num(x$tau0)))
  } else {
    cli::cli_alert_info("Sample-level errors not included in the model.")
  }

  #  Random effects
  if (!is.null(x$sigma_re)) {
    cli::cli_h2("Random effects (SD)")
    d <- x$sigma_re
    d$SD <- fmt_num(d$SD)
    print(d, row.names = FALSE, right = TRUE)
  }

  #  Observation-level error
  if (!is.null(x$scale_par)) {
    cli::cli_h2("Pairwise (observation-level) error")
    d <- x$scale_par
    d$Estimate <- fmt_num(d$Estimate)
    print(d, row.names = FALSE, right = TRUE)
  }

  # Fit of latent space
  cli::cli_h2("Fit of the latent space (conditional)")
  cat(sprintf("  Stress (metric, response scale):    %s\n", fmt_num(x$stress_metric)))
  cat(sprintf("  Stress (non-metric, link scale):    %s\n", fmt_num(x$stress_nonmetric)))
  cat(sprintf("  corR2:                              %s\n", fmt_num(x$cor_r2)))
  cat("  Stress and corR2 show how well the latent space captures dissimilarities,\n")
  cat("  not the explanatory power of the predictors.\n")

  cat("\n---\n")
  if (x$has_se) {
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  } else {
    cat("Standard errors not available: Wald tests omitted.\n")
  }
  invisible(x)
}
