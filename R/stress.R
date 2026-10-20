#' Calculate Stress for a fitted rrGDM model
#'
#' @description
#' Computes a goodness-of-fit "stress" statistic for a fitted \code{rrGDM}
#' object, comparing the observed dissimilarities used to fit the model
#' against a model-implied quantity computed from the conditional latent
#' coordinates (i.e., including sample-level errors): either the raw Euclidean
#' distance between fitted latent positions (\code{scale = "link"}), or the
#' inverse-link-transformed expected dissimilarity (\code{scale = "response"}).
#'
#' For \code{type = "nonmetric"}, this is Kruskal's stress-1,
#' \deqn{Stress = \sqrt{\frac{\sum (\hat{d} - D)^2}{\sum D^2}}}
#' where \eqn{D} is the chosen model-implied quantity and \eqn{\hat{d}}
#' ("disparities") is the isotonic (monotone) regression of \eqn{D} on the
#' observed dissimilarities.
#'
#' For \code{type = "metric"}, observed and model-implied dissimilarities are
#' compared directly, with no regression step:
#' \deqn{Stress = \sqrt{\frac{\sum (y_{ij} - \hat{y}_{ij})^2}{\sum y_{ij}^2}}}
#' where \eqn{y_{ij}} are the observed dissimilarities and \eqn{\hat{y}_{ij}}
#' the model-implied ones. Because it does not allow for a linear
#' recalibration, it also penalises a poorly calibrated link. This is not
#' directly comparable to the stress reported by NMDS software.
#'
#' Stress is returned as a proportion between 0 (best) and 1 (worst) for
#' bounded dissimilarities.
#'
#' @param object A fitted object of class \code{rrGDM}.
#' @param type Character string: \code{"nonmetric"} (default) computes
#'   Kruskal's stress-1 using isotonic regression, as in standard NMDS.
#'   \code{"metric"} compares observed and model-implied dissimilarities
#'   directly (see Description).
#' @param scale Character string: \code{"link"} (default) uses the raw
#'   Euclidean distance between fitted latent positions in the ordination.
#'   This is generally the only scale on which stress values from RRGDM and
#'   other methods (e.g., \code{vegan::metaMDS}, \code{MASS::isoMDS}) should be
#'   compared directly. \code{"response"} instead uses the inverse-link
#'   transformed expected dissimilarity (\code{object$inv_link(dist)}), which is
#'   on the same scale as the observed dissimilarities used to fit the model.
#'   This is the appropriate scale for \code{type = "metric"}.
#'
#' @return A single numeric value: the stress.
#'
#' @importFrom stats dist isoreg
#' @export
stress <- function(object,
                   type  = c("nonmetric", "metric"),
                   scale = c("link", "response")) {

  type <- match.arg(type)
  scale <- match.arg(scale)

  if (!inherits(object, "rrGDM")) {
    cli::cli_abort("{.arg object} must be a fitted {.cls rrGDM} model.")
  }

  if (scale == "link" && type == "metric") {
    cli::cli_warn(c(
      "!" = "{.code scale = \"link\"} with {.code type = \"metric\"} compares raw latent Euclidean distances directly with observed dissimilarities.",
      "i" = "These are generally on different scales unless {.code link = \"identity\"}, so the resulting stress is typically not meaningful.",
      "*" = "Consider using {.code type = \"nonmetric\"} with {.code scale = \"link\"}, or {.code scale = \"response\"} with {.code type = \"metric\"} to assess fit directly on the dissimilarity scale."
    ))
  }

  L_hat <- get_report(object)$L

  if (is.null(L_hat)) {
    cli::cli_abort(c(
      "x" = "Could not retrieve fitted latent coordinates {.field L} from {.arg object}."
    ))
  }

  fp <- fit_pairs(object, L_hat)
  raw_dist <- fp$dist

  fit_vec <- if (scale == "link") {
    as.vector(raw_dist)
  } else {
    as.vector(object$inv_link(raw_dist))
  }

  obs_vec <- fp$obs

  if (length(obs_vec) != length(fit_vec)) {
    cli::cli_abort(c(
      "x" = "Observed and modelled dissimilarity vectors differ in length."
    ))
  }

  ok <- is.finite(obs_vec) & is.finite(fit_vec)
  if (!all(ok)) {
    n_drop <- sum(!ok)
    cli::cli_warn(c("!" = "{n_drop} non-finite pair{?s} removed before computing stress."))
    obs_vec <- obs_vec[ok]
    fit_vec <- fit_vec[ok]
  }

  if (length(obs_vec) < 2) {
    cli::cli_abort("At least 2 finite pairs are required to compute stress.")
  }

  if (type == "nonmetric") {
    ord <- order(obs_vec, fit_vec)
    obs_sort <- obs_vec[ord]
    fit_sort <- fit_vec[ord]
    dhat  <- stats::isoreg(obs_sort, fit_sort)$yf
    num <- sum((dhat - fit_sort)^2)
    denom <- sum(fit_sort^2)

    if (denom == 0) {
      cli::cli_abort("Sum of squared modelled distances is 0; stress is undefined.")
    }
  } else {
    num <- sum((obs_vec - fit_vec)^2)
    denom <- sum(obs_vec^2)
    if (denom == 0) {
      cli::cli_abort("Sum of squared observed dissimilarities is 0; stress is undefined.")
    }
  }

  sqrt(num / denom)
}
