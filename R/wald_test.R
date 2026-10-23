#' Wald Tests for Environmental Loadings and Dispersion Effects
#'
#' @description
#' Wald tests for sets of predictor effects in a fitted \code{rrGDM} model:
#' \itemize{
#'   \item \strong{Environmental loadings} (\code{lambda}): for each
#'     predictor, a multivariate (multi-df) Wald test of the null that all
#'     free loadings of that predictor (one per latent dimension) are jointly zero.
#'     This accounts for the potential covariance between
#'     loadings within a predictor. The \code{Size} column gives the
#'     modulus (Euclidean norm) of the loading vector across the \code{k}
#'     latent dimensions, a rotation-invariant measure of the predictor's
#'     total effect on community composition.
#'   \item \strong{Dispersion effects} (\code{beta}): for each predictor in
#'     the dispersion model, a univariate Wald z-test, since each
#'     coefficient is a scalar effect on \eqn{\log(\tau)}.
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
#' because it explicitly accounts for the covariance between loadings.
#' If the covariance matrix \eqn{\boldsymbol{\Sigma}_p} of a predictor's
#' loadings is singular or ill-conditioned (judged on the correlation scale, so
#' independently of the units of the coefficients; reciprocal condition number
#' below \code{1e-8}), the Wald statistic and p-value for that predictor are set
#' to \code{NA} with a warning. This typically indicates unidentified or nearly
#' collinear coefficients, for which the test is not meaningful. If the Hessian
#' at the optimum was not positive definite, the function stops with an error.
#'
#' The degrees of freedom of the \code{lambda} test depend on
#' \code{object$constraint}:
#' \itemize{
#'   \item Under \code{constraint == 1}, the loadings matrix is constrained
#'     to be lower triangular: predictor (row) \code{p} only has free
#'     loadings in latent dimensions \code{1:min(p, k)}, the remainder being
#'     structural zeros that are fixed (not estimated) when the model is
#'     fit. The Wald test for predictor \code{p} therefore uses
#'     \code{Df = min(p, k)}, and only the free loadings (and their
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
#'     \code{X_design} (zero rows if the mean model has no predictors), with
#'     columns \code{term}, \code{Size}, \code{Wald},
#'     \code{Df} (see Details), \code{Pr_Chisq}, \code{Sig}.
#'   \item \code{beta}: a \code{data.frame} with one row per predictor in
#'     \code{W_design} (excluding the intercept unless
#'     \code{include_intercept = TRUE}), with columns \code{term},
#'     \code{Estimate}, \code{SE}, \code{z}, \code{Pr_z}, \code{Sig}.
#' }
#'
#' @seealso \code{\link{anova.rrGDM}} for tests by variable,
#'   \code{\link{summary.rrGDM}}
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

  check_pd_hessian(object)

  if (is.null(object$constraint)) {
    cli::cli_abort("{.arg object} has no {.field constraint} recorded; try refitting the model with {.fn rrGDM}. Please report if the issue persists.")
  }

  constraint <- object$constraint
  sd_rep <- object$sd_rep
  sdr_vals <- sd_rep$value
  sdr_cov  <- sd_rep$cov

  n_env <- ncol(object$X_design)
  if (is.null(n_env)) n_env <- 0L

  k  <- object$k
  env_nms <- colnames(object$X_design)
  lam_idx <- which(names(sdr_vals) == "lambda")

  if (length(lam_idx) != n_env * k) {
    cli::cli_abort(c(
      "x" = "Could not locate the expected number of {.field lambda} entries in {.field sd_rep}."
    ))
  }

  lambda_hat <- matrix(sdr_vals[lam_idx], nrow = n_env, ncol = k)
  Sigma_full <- sdr_cov[lam_idx, lam_idx, drop = FALSE]

  # Check values match with sdreport
  report_lambda <- tryCatch(get_report(object)$lambda, error = function(e) NULL)
  if (!is.null(report_lambda) &&
      all(dim(report_lambda) == dim(lambda_hat)) &&
      !isTRUE(all.equal(unname(lambda_hat), unname(report_lambda), tolerance = 1e-6))) {
    cli::cli_abort(c(
      "x" = "ADREPORT-derived {.field lambda} does not match {.code get_report(object)$lambda}.",
      "i" = "If the problem persists, please report the issue."
    ))
  }

  if (constraint == 0 && n_env > 0) {
    cli::cli_warn(c(
      "!" = "Model was fit with {.code constraint = 0} (no identifiability constraint).",
      "i" = "Individual loadings in {.field lambda} are not identified under rotation of the latent space; interpret per-predictor Wald tests with caution."
    ))
  }

  # Structural zeros are excluded from both the test statistic and its df.
  get_free_dims <- function(i) {
    if (constraint == 1) seq_len(min(i, k)) else seq_len(k)
  }

  singular_terms <- character(0)

  # Wald test for each predictor
  lambda_res <- lapply(seq_len(n_env), function(i) {
    fd <- get_free_dims(i)
    zd <- setdiff(seq_len(k), fd)

    if (length(zd) > 0 &&
        !isTRUE(all.equal(unname(lambda_hat[i, zd]), rep(0, length(zd)), tolerance = 1e-8))) {
      cli::cli_abort(c(
        "x" = "Unexpected non-zero loading found where a structural zero was expected for {.val {env_nms[i]}}."
      ))}

    df_i <- length(fd)
    idx_i <- i + (fd - 1) * n_env
    lam_i <- lambda_hat[i, fd]
    Sigma_i <- Sigma_full[idx_i, idx_i, drop = FALSE]

    # Effect size
    size_i <- sqrt(sum(lam_i^2))

    # Wald statistic (NA if Sigma_i is singular or ill-conditioned) and pval from chi dist.
    wj <- wald_joint(lam_i, Sigma_i)
    if (!wj$ok) singular_terms <<- c(singular_terms, env_nms[i])
    W <- wj$W
    pv <- if (is.na(W)) NA_real_ else stats::pchisq(W, df = df_i, lower.tail = FALSE)

    data.frame(term = env_nms[i], Size = size_i, Wald = W, Df = df_i,
               Pr_Chisq = pv, stringsAsFactors = FALSE)
  })

  lambda_df <- do.call(rbind, lambda_res)

  # Empty loadings table when there are no predictors in the mean model
  if (n_env == 0) {
    lambda_df <- data.frame(term = character(0), Size = numeric(0), Wald = numeric(0),
                            Df = integer(0), Pr_Chisq = numeric(0),
                            stringsAsFactors = FALSE)
  }

  if (length(singular_terms) > 0) {
    cli::cli_warn(c(
      "!" = "Covariance matrix was singular or ill-conditioned for term{?s}: {.val {singular_terms}}.",
      "i" = "{.field Wald} and {.field Pr_Chisq} were set to {.code NA} for {cli::qty(length(singular_terms))}{?this/these} term{?s}.",
      "i" = "This usually indicates unidentified or nearly collinear coefficients."
    ))
  }

  lambda_df$Sig <- cut(
    lambda_df$Pr_Chisq,
    breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
    labels = c("***", "**", "*", ".", " "))

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
      "i" = "Expected {length(w_nms)}, found {length(beta_idx)}."
    ))
  }

  beta_hat <- sdr_vals[beta_idx]
  se_hat <- sqrt(diag(sdr_cov[beta_idx, beta_idx, drop = FALSE]))

  # Sanity check that values match
  report_beta <- tryCatch(get_report(object)$beta, error = function(e) NULL)
  if (!is.null(report_beta) && length(report_beta) == length(beta_hat) &&
      !isTRUE(all.equal(unname(beta_hat), unname(report_beta), tolerance = 1e-6))) {
    cli::cli_abort(c(
      "x" = "ADREPORT-derived {.field beta} does not match {.code get_report(object)$beta}.",
      "i" = "This is likely a bug, please report if the problem persists."
    ))
  }

  z <- beta_hat / se_hat
  pval <- 2 * stats::pnorm(-abs(z))

  beta_df <- data.frame(
    term = w_nms,
    Estimate = beta_hat,
    SE = se_hat,
    z = z,
    Pr_z = pval,
    stringsAsFactors = FALSE)

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

#' Format p-values for printing (internal)
#'
#' Fixed notation for moderate/large values (3 decimals), scientific notation
#' for small values (<0.001), and a floor of "<1e-06" for very small values.
#' @param p Numeric vector of p-values.
#' @param floor Values below this are printed as "<floor".
#' @noRd
format_pval <- function(p, floor = 1e-6) {
  out <- character(length(p))
  na <- is.na(p)
  lo <- !na & p < floor
  sci <- !na & p >= floor & p < 1e-3
  fix <- !na & p >= 1e-3
  out[na] <- "NA"
  out[lo] <- paste0("<", format(floor, scientific = TRUE))
  out[sci] <- formatC(p[sci], format = "e", digits = 1)
  out[fix] <- formatC(p[fix], format = "f", digits = 3)
  out
}

#' Print method for wald_rrGDM objects
#'
#' @param x An object of class \code{wald_rrGDM}, as returned by
#'   \code{\link{wald_test}}.
#' @param digits Number of significant digits for test statistics and
#'   effect sizes. Default is 4.
#' @param ... Currently unused.
#' @return Invisibly returns \code{x}.
#' @export
#' @method print wald_rrGDM
print.wald_rrGDM <- function(x, digits = 4, ...) {
  fmt_num <- function(v) formatC(v, digits = digits, format = "fg", flag = "#")

  cli::cli_h2("Wald tests: environmental loadings")
  if (nrow(x$lambda) > 0) {
    d <- x$lambda
    d$Size <- fmt_num(d$Size)
    d$Wald <- round(d$Wald, 3)
    d$Pr_Chisq <- format_pval(d$Pr_Chisq)
    d$Sig[is.na(d$Sig)] <- " "   # NA is printed as blank
    print(d, row.names = FALSE, right = TRUE)
  } else {
    cli::cli_alert_info("No predictors to test.")
  }

  cli::cli_h2("Wald tests: dispersion effects")
  if (nrow(x$beta) > 0) {
    d <- x$beta
    d$Estimate <- fmt_num(d$Estimate)
    d$SE <- fmt_num(d$SE)
    d$z <- fmt_num(d$z)
    d$Pr_z <- format_pval(d$Pr_z)
    d$Sig[is.na(d$Sig)] <- " "
    print(d, row.names = FALSE, right = TRUE)
  } else {
    cli::cli_alert_info("No dispersion terms.")
  }

  cat("\n---\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  cat("Wald tests on loadings can be anti-conservative with few samples.\n")
  invisible(x)
}




#' Wald analysis of variance table for rrGDM models
#'
#' @description
#' Joint Wald chi-square tests for each \emph{variable} in a fitted
#' \code{rrGDM} model, however many terms (design-matrix columns) the variable
#' contributes. Where \code{\link{wald_test}} gives one test per column of the
#' design matrix, \code{anova()} groups every coefficient that belongs to the
#' same variable into a single test, so that:
#' \itemize{
#'   \item a factor with several levels gets one test rather than one per level
#'   \item terms such as \code{x + I(x^2)}, \code{poly(x, 2)} or \code{splines::ns(x, 3)} are
#'     tested together as "the effect of \code{x}"
#'   \item a variable that enters both the location model (\code{formula}) and
#'     the dispersion model (\code{formula_sd}) can be tested on both at once, or
#'     on each separately.
#' }
#'
#' @param object A fitted model of class \code{rrGDM}, fitted with
#'   \code{se = TRUE}.
#' @param ... Not used. Only a single model can be supplied (comparison of
#'   several models is not supported).
#' @param interactions Logical. Controls how terms involving more than one
#'   variable (e.g., \code{x:f}) are handled. If \code{TRUE} (default), a
#'   variable's test includes its main-effect terms \emph{and} every
#'   interaction term it takes part in, i.e., it tests whether the variable has
#'   any effect on composition (or dispersion) at all. If \code{FALSE},
#'   interaction terms are omitted. Variables that only appear in interactions are then
#'   omitted.
#' @param separate Logical. A variable can appear in both the location model
#'   (\code{formula}) and the dispersion model (\code{formula_sd}). If
#'   \code{FALSE} (default), these are tested jointly in a single row
#'   (\code{effect = "location + dispersion"}), accounting for the covariance
#'   between location and dispersion coefficients. If \code{TRUE}, the
#'   location and dispersion effects of each variable get a row each.
#'
#' @details
#' Intercepts and random effects are excluded.
#'
#' For each row, the null hypothesis is that all of the coefficients involved
#' are jointly zero. As in \code{\link{wald_test}},
#' structural zeros of the loadings matrix under \code{constraint == 1} are
#' excluded from both the statistic and the degrees of freedom. If the joint
#' covariance matrix of the coefficients being tested is singular or
#' ill-conditioned (judged on the correlation scale; reciprocal condition number
#' below \code{1e-8}), the Wald statistic and p-value for that row are set to
#' \code{NA} with a warning, since the test is not meaningful when coefficients
#' are unidentified or nearly collinear. If the Hessian at the optimum was not
#' positive definite, the function stops with an error.
#'
#' With \code{interactions = TRUE} the test compares the model against one where
#' the variable is dropped entirely. With \code{interactions = FALSE} the test of a
#' variable that interacts with others refers to its effect when the interacting
#' variables are at their reference level (factors) or zero (continuous
#' predictors).
#'
#' As with \code{\link{wald_test}}, tests on location effects can be
#' anti-conservative with few samples.
#'
#' @return
#' A \code{data.frame} of class \code{anova_rrGDM}, with one row per variable
#' (or per variable and effect type if \code{separate = TRUE}) and columns:
#' \itemize{
#'   \item \code{variable}: name of the variable.
#'   \item \code{effect}: \code{"location"}, \code{"dispersion"} or
#'     \code{"location + dispersion"}.
#'   \item \code{terms}: the model terms grouped in the test.
#'   \item \code{Df}: degrees of freedom (number of free coefficients tested).
#'   \item \code{Wald}: the Wald chi-squared statistic.
#'   \item \code{Pr_Chisq}: p-value.
#'   \item \code{Sig}: significance code.
#' }
#'
#' @seealso \code{\link{wald_test}} for per-term tests, effect sizes and
#'   dispersion estimates.
#'
#' @examples
#' \dontrun{
#' m <- rrGDM(Y = com,
#'            X = env,
#'            k = 2,
#'            formula = ~ pH + I(pH^2) + Region,
#'            formula_sd = ~ Region,
#'            family = 'binomial')
#'
#' # One row per variable; Region tested on location and dispersion jointly
#' anova(m)
#' }
#'
#' @export
#' @method anova rrGDM
#' @importFrom stats anova pchisq
anova.rrGDM <- function(object,
                        ...,
                        interactions = TRUE,
                        separate = FALSE) {

  # Checks
  if (!inherits(object, "rrGDM")) {
    cli::cli_abort("{.arg object} must be a fitted {.cls rrGDM} model.")
  }

  if (...length() > 0) {
    cli::cli_abort(c(
      "x" = "{.fn anova} for {.cls rrGDM} objects accepts a single model.",
      "i" = "Comparison of several models is not supported."
    ))
  }

  if (!is.logical(interactions) || length(interactions) != 1 || is.na(interactions)) {
    cli::cli_abort("{.arg interactions} must be a single {.code TRUE}/{.code FALSE} value.")
  }
  if (!is.logical(separate) || length(separate) != 1 || is.na(separate)) {
    cli::cli_abort("{.arg separate} must be a single {.code TRUE}/{.code FALSE} value.")
  }

  if (is.null(object$sd_rep)) {
    cli::cli_abort(c(
      "x" = "{.arg object} has no {.field sd_rep}.",
      "i" = "Refit the model with {.code se = TRUE}, or run {.fn include_se} on the fitted object."
    ))
  }

  check_pd_hessian(object)

  if (is.null(object$constraint)) {
    cli::cli_abort("{.arg object} has no {.field constraint} recorded; try refitting the model with {.fn rrGDM}. Please report if the issue persists.")
  }

  constraint <- object$constraint
  k <- object$k
  sdr_vals <- object$sd_rep$value
  sdr_cov <- object$sd_rep$cov

  X_design <- object$X_design
  W_design <- object$W_design
  n_env <- if (is.null(X_design)) 0L else ncol(X_design)

  # ------------------- ESTIMATES AND POSITIONS IN sd_rep -------------------
  lam_idx <- which(names(sdr_vals) == "lambda")
  beta_idx <- which(names(sdr_vals) == "beta")

  if (length(lam_idx) != n_env * k) {
    cli::cli_abort(c(
      "x" = "Could not locate the expected number of {.field lambda} entries in {.field sd_rep}."
    ))
  }
  if (length(beta_idx) != ncol(W_design)) {
    cli::cli_abort(c(
      "x" = "Could not locate the expected number of {.field beta} entries in {.field sd_rep}.",
      "i" = "Expected {ncol(W_design)}, found {length(beta_idx)}."
    ))
  }

  lambda_hat <- matrix(sdr_vals[lam_idx], nrow = n_env, ncol = k)
  beta_hat <- unname(sdr_vals[beta_idx])

  # Check that the values match the model report (guards the indexing above)
  report <- tryCatch(get_report(object), error = function(e) NULL)
  if (!is.null(report)) {
    if (!is.null(report$lambda) && all(dim(report$lambda) == dim(lambda_hat)) &&
        !isTRUE(all.equal(unname(lambda_hat), unname(report$lambda), tolerance = 1e-6))) {
      cli::cli_abort(c(
        "x" = "ADREPORT-derived {.field lambda} does not match {.code get_report(object)$lambda}.",
        "i" = "If the problem persists, please report the issue."
      ))
    }
    if (!is.null(report$beta) && length(report$beta) == length(beta_hat) &&
        !isTRUE(all.equal(beta_hat, unname(report$beta), tolerance = 1e-6))) {
      cli::cli_abort(c(
        "x" = "ADREPORT-derived {.field beta} does not match {.code get_report(object)$beta}.",
        "i" = "If the problem persists, please report the issue."
      ))
    }
  }

  if (constraint == 0 && n_env > 0) {
    cli::cli_warn(c(
      "!" = "Model was fit with {.code constraint = 0} (no identifiability constraint).",
      "i" = "Individual loadings in {.field lambda} are not identified under rotation of the latent space; interpret location tests with caution."
    ))
  }

  # Structural zeros of lambda (constraint == 1) must really be zero
  if (constraint == 1 && n_env > 0 && k > 1) {
    zeros <- lambda_hat[upper.tri(lambda_hat)]
    if (!isTRUE(all.equal(unname(zeros), rep(0, length(zeros)), tolerance = 1e-8))) {
      cli::cli_abort(c(
        "x" = "Unexpected non-zero loading found where a structural zero was expected."
      ))
    }
  }

  # Free (non-structural-zero) latent dimensions for row i of lambda
  free_dims <- function(i) {
    if (constraint == 1) seq_len(min(i, k)) else seq_len(k)
  }

  # VARIABLES AND TERMS
  loc <- anova_term_map(X_design, object$terms_fix)
  disp <- anova_term_map(W_design, object$terms_sd)

  # Terms of a map that belong to variable v (honouring `interactions`)
  terms_of <- function(map, v) {
    if (length(map$vars) == 0) return(character(0))
    keep <- vapply(map$vars, function(vs) {
      v %in% vs && (interactions || length(vs) == 1L)
    }, logical(1))
    names(map$vars)[keep]
  }

  all_vars <- as.character(unique(c(unlist(loc$vars, use.names = FALSE),
                                    unlist(disp$vars, use.names = FALSE))))

  has_interactions <- any(lengths(c(loc$vars, disp$vars)) > 1L)

  # TESTS
  singular_rows <- character(0)

  make_row <- function(v, effect, loc_terms, disp_terms) {
    est <- numeric(0)
    pos <- integer(0)

    if (length(loc_terms) > 0) {
      for (i in loc$cols$col[loc$cols$term %in% loc_terms]) {
        fd <- free_dims(i)
        est <- c(est, lambda_hat[i, fd])
        pos <- c(pos, lam_idx[i + (fd - 1L) * n_env])
      }
    }
    if (length(disp_terms) > 0) {
      j <- disp$cols$col[disp$cols$term %in% disp_terms]
      est <- c(est, beta_hat[j])
      pos <- c(pos, beta_idx[j])
    }

    res <- wald_joint(est, sdr_cov[pos, pos, drop = FALSE])
    if (!res$ok) singular_rows <<- c(singular_rows, paste0(v, " (", effect, ")"))

    term_txt <- if (length(loc_terms) > 0 && length(disp_terms) > 0) {
      paste0("location: ", paste(loc_terms, collapse = ", "),
             "; dispersion: ", paste(disp_terms, collapse = ", "))
    } else {
      paste(c(loc_terms, disp_terms), collapse = ", ")
    }

    data.frame(variable = v,
               effect = effect,
               terms = term_txt,
               Df = res$df,
               Wald = res$W,
               Pr_Chisq = if (is.na(res$W) || res$df == 0) NA_real_ else
                 stats::pchisq(res$W, df = res$df, lower.tail = FALSE),
               stringsAsFactors = FALSE)
  }

  rows <- list()
  for (v in all_vars) {
    lt <- terms_of(loc, v)
    dt <- terms_of(disp, v)

    if (length(lt) == 0 && length(dt) == 0) next   # e.g., only in interactions

    if (length(lt) > 0 && length(dt) > 0 && !separate) {
      rows[[length(rows) + 1]] <- make_row(v, "location + dispersion", lt, dt)
    } else {
      if (length(lt) > 0) rows[[length(rows) + 1]] <- make_row(v, "location", lt, character(0))
      if (length(dt) > 0) rows[[length(rows) + 1]] <- make_row(v, "dispersion", character(0), dt)
    }
  }

  out <- if (length(rows) > 0) {
    do.call(rbind, rows)
  } else {
    data.frame(variable = character(0), effect = character(0), terms = character(0),
               Df = integer(0), Wald = numeric(0), Pr_Chisq = numeric(0),
               stringsAsFactors = FALSE)
  }

  if (length(singular_rows) > 0) {
    cli::cli_warn(c(
      "!" = "Covariance matrix was singular or ill-conditioned for: {.val {singular_rows}}.",
      "i" = "{.field Wald} and {.field Pr_Chisq} were set to {.code NA} for {cli::qty(length(singular_rows))}{?this/these} test{?s}.",
      "i" = "This usually indicates unidentified or nearly collinear coefficients."
    ))
  }

  out$Sig <- cut(
    out$Pr_Chisq,
    breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
    labels = c("***", "**", "*", ".", " "))

  rownames(out) <- NULL

  structure(out,
            class = c("anova_rrGDM", "data.frame"),
            interactions = interactions,
            separate = separate,
            has_interactions = has_interactions,
            constraint = constraint)
}


#' Map design-matrix columns to terms and variables (internal)
#'
#' @param design Design matrix with an \code{"assign"} attribute (as returned
#'   by \code{model.matrix}); intercept column, if any, has \code{assign == 0}.
#' @param tt A \code{terms} object for the same formula.
#' @return A list with \code{cols}, a data frame with one row per non-intercept
#'   design column (\code{col}: column index; \code{term}: term label), and
#'   \code{vars}, a named list (one element per term label) with the raw
#'   variable names each term involves.
#' @noRd
anova_term_map <- function(design, tt) {
  empty <- list(cols = data.frame(col = integer(0), term = character(0),
                                  stringsAsFactors = FALSE),
                vars = list())

  if (is.null(design) || is.null(tt) || ncol(design) == 0) return(empty)

  labels <- attr(tt, "term.labels")
  if (length(labels) == 0) return(empty)

  assign <- attr(design, "assign")
  if (is.null(assign) || length(assign) != ncol(design)) {
    cli::cli_abort(c(
      "x" = "Could not map design-matrix columns to model terms.",
      "i" = "The design matrix has no {.field assign} attribute; try refitting the model with {.fn rrGDM}."
    ))
  }

  vars <- lapply(labels, function(l) unique(all.vars(parse(text = l)[[1]])))
  names(vars) <- labels

  cols <- which(assign > 0)
  list(cols = data.frame(col = cols, term = labels[assign[cols]],
                         stringsAsFactors = FALSE),
       vars = vars)
}


#' Joint Wald statistic (internal)
#'
#' Computes \eqn{W = \hat\theta^\top \Sigma^{-1} \hat\theta}. Conditioning is
#' judged on the correlation scale (the statistic does not depend on the units
#' of the coefficients, so neither should the check). If \code{Sigma} is not
#' finite, not positive definite, or has a reciprocal condition number (of its
#' correlation matrix) below \code{rcond_tol}, no statistic is returned: the
#' coefficients are then unidentified or (nearly) collinear and a Wald test is
#' not meaningful. No pseudoinverse is used.
#'
#' @param est Numeric vector of estimates.
#' @param Sigma Covariance matrix of \code{est}.
#' @param rcond_tol Smallest acceptable reciprocal condition number of the
#'   correlation matrix.
#' @return A list with the statistic \code{W} (\code{NA} if it could not be
#'   computed), degrees of freedom \code{df} (the length of \code{est}) and a
#'   flag \code{ok} (\code{FALSE} if the covariance matrix was unusable).
#' @noRd
wald_joint <- function(est, Sigma, rcond_tol = 1e-8) {
  df <- length(est)
  if (df == 0L) return(list(W = NA_real_, df = 0L, ok = TRUE))

  fail <- list(W = NA_real_, df = as.integer(df), ok = FALSE)

  if (anyNA(est) || any(!is.finite(Sigma))) return(fail)

  Sigma <- (Sigma + t(Sigma)) / 2
  s <- sqrt(diag(Sigma))
  if (any(!is.finite(s)) || any(s <= 0)) return(fail)

  R <- Sigma / outer(s, s)
  ch <- tryCatch(chol(R), error = function(e) NULL)   # fails if not positive definite
  if (is.null(ch) || rcond(R) < rcond_tol) return(fail)

  # W = u' R^{-1} u with u = est / s and R = ch'ch
  z <- backsolve(ch, est / s, transpose = TRUE)
  list(W = sum(z^2), df = as.integer(df), ok = TRUE)
}


#' Stop if the standard errors cannot be used (internal)
#'
#' Stops if the Hessian at the optimum was not positive definite (\code{sdreport}
#' then returns NaN covariances: the model has not converged to a proper maximum
#' or has unidentified parameters), or if the covariance matrix of the reported
#' quantities is missing (\code{sdreport} run with
#' \code{getReportCovariance = FALSE} stores \code{NA} there).
#' @param object A fitted \code{rrGDM} object with an \code{sd_rep}.
#' @noRd
check_pd_hessian <- function(object) {
  if (isFALSE(hess_ok(object))) {
    cli::cli_abort(c(
      "x" = "The Hessian at the optimum is not positive definite (or is nearly singular), so standard errors and Wald tests are not available.",
      "i" = "This usually means the model has not converged or has unidentified parameters (e.g., empty or collinear predictor levels, or too many latent dimensions for the data).",
      "*" = "Try {.code n_starts > 1}, fewer latent dimensions or a simpler model."
    ))
  }
  cov_rep <- object$sd_rep$cov
  if (!is.matrix(cov_rep) || nrow(cov_rep) != length(object$sd_rep$value)) {
    cli::cli_abort(c(
      "x" = "The standard-error report of {.arg object} has no covariance matrix.",
      "i" = "Refit with {.code se = TRUE}, or run {.fn include_se} on the fitted object without {.code getReportCovariance = FALSE}."
    ))
  }
  invisible(TRUE)
}


#' Print method for anova_rrGDM objects
#'
#' @param x An object of class \code{anova_rrGDM}, as returned by
#'   \code{\link{anova.rrGDM}}.
#' @param digits Number of decimal places for the Wald statistics. Default is 3.
#' @param ... Currently unused.
#' @return Invisibly returns \code{x}.
#' @export
#' @method print anova_rrGDM
print.anova_rrGDM <- function(x, digits = 3, ...) {
  inter <- isTRUE(attr(x, "interactions"))
  sep <- isTRUE(attr(x, "separate"))

  cli::cli_h2("Wald tests by variable")
  cli::cat_bullet(c(
    if (inter) "Interaction terms included in each variable's test" else
      "Interaction terms excluded: main-effect terms only",
    if (sep) "Location and dispersion effects tested separately" else
      "Location and dispersion effects of the same variable tested jointly"))
  cat("\n")

  if (nrow(x) == 0) {
    cli::cli_alert_info("No variables to test.")
    return(invisible(x))
  }

  d <- x
  class(d) <- "data.frame"

  # Text columns left-aligned, numeric columns right-aligned (headers included)
  col <- function(header, values, left = FALSE) {
    values <- as.character(values)
    w <- max(nchar(c(header, values), type = "width"))
    format(c(header, values), width = w, justify = if (left) "left" else "right")
  }

  sig <- as.character(d$Sig)
  sig[is.na(sig)] <- ""

  tab <- cbind(
    col("variable", d$variable, left = TRUE),
    col("effect", d$effect, left = TRUE),
    col("Df", d$Df),
    col("Wald", formatC(d$Wald, format = "f", digits = digits)),
    col("Pr_Chisq", format_pval(d$Pr_Chisq)),
    col("", sig, left = TRUE))

  lines <- apply(tab, 1, paste, collapse = "  ")
  cat(paste0(sub("\\s+$", "", lines), "\n"), sep = "")

  # Terms are only listed when they add information beyond the variable name
  trivial <- d$terms == d$variable |
    d$terms == paste0("location: ", d$variable, "; dispersion: ", d$variable)
  if (any(!trivial)) {
    cat("\nTerms grouped in each test:\n")
    cat(paste0("  ", d$variable[!trivial], " (", d$effect[!trivial], "): ",
               d$terms[!trivial], "\n"), sep = "")
  }

  cat("\n---\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  if (anyNA(d$Wald)) {
    cat("NA: covariance matrix of the coefficients was singular or ill-conditioned.\n")
  }
  if (!inter && isTRUE(attr(x, "has_interactions"))) {
    cat("Interactions excluded: tests refer to each variable's effect when the\n")
    cat("variables it interacts with are at their reference level (or zero).\n")
  }
  cat("Tests on location effects can be anti-conservative with few samples.\n")
  invisible(x)
}
