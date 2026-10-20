#' Given a numeric vector \code{x}, this function re-scales \code{x} into the interval \[\code{min}, \code{max}\].
#'
#' @param x A numeric vector
#' @param new_range Numeric scalar: the upper and lower bound of the target interval. Default is \[\code{0.01}, \code{0.99}\].
#' @return A numeric vector of the same length as \code{x}, with values in the interval \[\code{min(new_range)}, \code{max(new_range)}\].
#' @noRd

scale_dist <- function(x, new_range = c(0.01, 0.99)) {
  if (!is.null(new_range)) {
    stopifnot(
      "x must be numeric" = is.numeric(x),
      "new_range must be a numeric vector of length 2" = is.numeric(new_range) && length(new_range) == 2,
      "new_range[1] and new_range[2] cannot be equal" = new_range[1] != new_range[2])
    min(new_range) + x * (max(new_range) - min(new_range))
  } else {
    x
  }
}

#' Calculate Numerator and Denominator for Bray-Curtis Dissimilarity
#'
#' @description
#' An internal helper function that computes the individual components of the
#' Bray-Curtis distance between two vectors. This is used when \code{num_den = TRUE} in the model fitting function
#'
#' @param x A numeric vector representing species abundances at site 1.
#' @param y A numeric vector representing species abundances at site 2.
#'
#' @return A named numeric vector of length 2 containing the elements \code{num} and \code{den}.
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
         'negexp'   = function(x) 1 - exp(-x),
         'neggaus'  = function(x) 1 - exp(-(x^2)),
         cli::cli_abort("Link {.val {link}} is not supported."))
}



#' Model report at the best-fit parameters (internal)
#' `obj$report()` with no argument uses `obj$env$last.par`, which later calls
#' (e.g. `sdreport`) can overwrite. `best_par` is the snapshot taken after fitting.
#' @noRd
get_report <- function(object) object$obj$report(object$best_par)



#' Collapse a formula into a single string (internal)
#' @param f A formula.
#' @noRd
format_formula <- function(f) paste(trimws(deparse(f)), collapse = " ")




#' Subset a long data frame of dissimilarities to the pairs in a sparse_diss object (internal)
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

#' Observed and modelled distances
#' @noRd
fit_pairs <- function(object, L_hat) {
  obs <- if (is.null(object$sparse) || is.null(object$obs_full)) {
    object$Y_diss / object$Y_den
  } else {
    object$obs_full
  }
  list(obs = obs, dist = as.vector(stats::dist(L_hat)))
}
