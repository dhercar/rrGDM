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
