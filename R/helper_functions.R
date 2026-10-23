#' Map dissimilarities in [0, 1] onto a new interval (internal)
#'
#' A linear map of [0, 1] onto \code{new_range}. \code{x} is assumed to lie in
#' [0, 1] (as bounded dissimilarities do); values outside are mapped linearly and
#' fall outside \code{new_range}. \code{new_range} is validated in
#' \code{rrGDM()}.
#'
#' @param x A numeric vector with values in [0, 1].
#' @param new_range Numeric vector of length 2 (the bounds of the target
#'   interval, in any order) or \code{NULL} to leave \code{x} unchanged.
#'   Default is c(0.01, 0.99).
#' @return A numeric vector of the same length as \code{x}.
#' @noRd
scale_dist <- function(x, new_range = c(0.01, 0.99)) {
  if (!is.null(new_range)) {
    min(new_range) + x * (max(new_range) - min(new_range))
  } else {
    x
  }
}

#' Calculate Numerator and Denominator for Bray-Curtis Dissimilarity
#'
#' @description
#' Numerator and denominator of the Bray-Curtis dissimilarity between two sites,
#' with pairwise deletion of NAs. Used by \code{make_Y_diss()} (\code{num_den = TRUE})
#' only when the community matrix contains NAs; otherwise it uses vectorised code.
#'
#' @param x A numeric vector representing species abundances at site 1.
#' @param y A numeric vector representing species abundances at site 2.
#'
#' @return A named numeric vector of length 2 containing the elements \code{num_bra} and \code{den_bra}.
#' @noRd
bray_num_den <- function(x, y) {
  num <- sum(abs(x - y), na.rm = TRUE)
  den <- sum(x + y, na.rm = TRUE)
  return(c(num_bra = num, den_bra = den))
}

#' Build inverse-link function for an rrGDM model
#'
#' @param link Character string: "identity", "negexp", or "neggaus".
#' @return A function mapping latent Euclidean distance to expected dissimilarity.
#' @noRd
make_inv_link <- function(link) {
  switch(link,
         'identity' = function(x) x,
         'negexp' = function(x) -expm1(-x),
         'neggaus'  = function(x) -expm1(-x^2),
         cli::cli_abort("Link {.val {link}} is not supported."))
}

#' Model report at the best-fit parameters (internal)
#'
#' `obj$report()` with no argument uses `obj$env$last.par`, which later calls
#' (e.g. `sdreport`) can overwrite. `best_par` is the snapshot taken after fitting.
#' @noRd
get_report <- function(object) object$obj$report(object$best_par)



#' Collapse a formula into a single string (internal)
#' @param f A formula.
#' @noRd
format_formula <- function(f) paste(trimws(deparse(f)), collapse = " ")




#' Subset a long data frame of dissimilarities to the pairs in a sparse_diss object (internal)
#'
#' Rows follow \code{combn(n, 2)} (\code{s1 < s2}), the same order as
#' \code{as.vector(dist())}. \code{sparse$edges} has \code{s1 < s2} too (as
#' created by \code{sparse_diss()}); the result is in the order of the edges.
#' @noRd
subset_pairs <- function(Y_df, sparse, n, id = NULL, dist_method = NULL) {
  if (!inherits(sparse, "sparse_diss")) {
    cli::cli_abort("{.arg sparse} must be a {.cls sparse_diss} object (see {.fn sparse_diss}), not a {.cls {class(sparse)}}.")
  }
  if (sparse$n_samples != n) {
    cli::cli_abort(c(
      "x" = "{.arg sparse} and {.arg Y} must have the same samples.",
      "i" = "{.arg sparse} was created from {sparse$n_samples} sample{?s}; {.arg Y} has {n}."
    ))
  }
  if (!is.null(id) && !identical(as.character(id), as.character(sparse$id))) {
    cli::cli_abort(c(
      "x" = "Row names of {.arg Y} differ from those of the data used to create {.arg sparse}.",
      "i" = "Samples must be in the same order."
    ))
  }
  if (!is.null(dist_method) && !(sparse$dist_method %in% c(dist_method, "user-supplied dist"))) {
    cli::cli_warn(c(
      "!" = "{.arg sparse} was created with {.code dist_method = {.val {sparse$dist_method}}}, but the model uses {.val {dist_method}}.",
      "i" = "Pairs were selected using a different dissimilarity."
    ))
  }
  if (sparse$n_components > 1) {
    cli::cli_warn(c(
      "!" = "The pairs in {.arg sparse} form {sparse$n_components} disconnected groups of samples.",
      "i" = "Their relative positions in the latent space are not constrained. Consider a larger {.arg n_random} in {.fn sparse_diss}."
    ))
  }

  key_all <- (Y_df$s1 - 1) * n + Y_df$s2
  key_sparse <- (sparse$edges[, 1] - 1) * n + sparse$edges[, 2]
  keep <- match(key_sparse, key_all)

  if (anyNA(keep)) {
    cli::cli_abort("Some pairs in {.arg sparse} could not be found in the dissimilarities. Was it created from a different {.arg Y}?")
  }

  Y_df[keep, , drop = FALSE]
}

#' Observed dissimilarities and latent distances for all pairs (internal)
#'
#' Both vectors follow \code{combn(n, 2)}, the order of \code{as.vector(dist())}.
#' For sparse fits the observed values are \code{object$obs_full} (all pairs);
#' otherwise they are \code{Y_diss / Y_den} (all pairs were fitted).
#'
#' @param object A fitted \code{rrGDM} object.
#' @param L_hat Matrix of latent coordinates (samples x dimensions).
#' @return A list with \code{obs} (observed dissimilarities) and \code{dist}
#'   (Euclidean distances between samples in the latent space).
#' @noRd
fit_pairs <- function(object, L_hat) {
  obs <- if (is.null(object$sparse) || is.null(object$obs_full)) {
    object$Y_diss / object$Y_den
  } else {
    object$obs_full
  }
  list(obs = obs, dist = as.vector(stats::dist(L_hat)))
}


#' Is a Hessian positive definite, with a tolerance? (internal)
#'
#' Checks the Hessian of the negative log-likelihood at an optimum (all
#' eigenvalues should be positive). The Hessian is first rescaled to a
#' correlation-type matrix, so the check does not depend on the units of the
#' parameters.
#'
#' @param H Numeric Hessian matrix (of the negative log-likelihood).
#' @param tol Smallest acceptable ratio between the smallest and the largest
#'   eigenvalue of the rescaled Hessian.
#' @return \code{TRUE} or \code{FALSE}.
#' @noRd
hess_is_pd <- function(H, tol = 1e-9) {
  if (is.null(H) || any(!is.finite(H))) return(FALSE)
  H <- (H + t(H)) / 2
  d <- diag(H)
  if (any(d <= 0)) return(FALSE)
  ev <- eigen(H / sqrt(outer(d, d)), symmetric = TRUE, only.values = TRUE)$values
  min(ev) > tol * max(ev)
}

#' Positive-definiteness flag from an sdreport object (internal)
#'
#' Combines TMB's own \code{pdHess} with \code{\link{hess_is_pd}} applied to the
#' Hessian implied by \code{sd_rep$cov.fixed}.
#'
#' @param sd_rep An \code{sdreport} object (or \code{NULL}).
#' @return \code{NA} if \code{sd_rep} is \code{NULL}, otherwise \code{TRUE}/\code{FALSE}.
#' @noRd
compute_pd_hess <- function(sd_rep) {
  if (is.null(sd_rep)) return(NA)
  isTRUE(sd_rep$pdHess) &&
    isTRUE(tryCatch(hess_is_pd(solve(sd_rep$cov.fixed)), error = function(e) FALSE))
}

#' Is the Hessian of a fitted rrGDM usable? (internal)
#'
#' Uses the stored \code{pd_hess} flag; for models fitted before it existed,
#' falls back to TMB's \code{pdHess}.
#'
#' @param object A fitted \code{rrGDM} object.
#' @return \code{NA} if there is no \code{sd_rep}, otherwise \code{TRUE}/\code{FALSE}.
#' @noRd
hess_ok <- function(object) {
  if (is.null(object$sd_rep)) return(NA)
  if (!is.null(object$pd_hess)) object$pd_hess else isTRUE(object$sd_rep$pdHess)
}
