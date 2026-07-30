#' Calculate Stress for a fitted rrGDM model
#'
#' @description
#' Computes a goodness-of-fit "stress" statistic for a fitted \code{rrGDM}
#' object, comparing the observed dissimilarities used to fit the model
#' against a model-implied quantity: either the raw Euclidean distance
#' between fitted latent positions (\code{scale = "link"}), or the
#' link-transformed expected dissimilarity (\code{scale = "response"}),
#' \deqn{Stress = \sqrt{\frac{\sum (\hat{d} - D)^2}{\sum D^2}}}
#' where \eqn{D} is the chosen model-implied quantity and \eqn{\hat{d}}
#' ("disparities") is a transformation of the observed dissimilarities
#' regressed onto \eqn{D}. Returned as a percentage.
#'
#' @param object A fitted object of class \code{rrGDM}.
#' @param type Character string: \code{"nonmetric"} (default) computes
#'   disparities via isotonic (monotone) regression, as in standard
#'   Kruskal / NMDS stress. \code{"metric"} uses linear regression instead.
#' @param scale Character string: \code{"link"} (default) uses the raw
#'   Euclidean distance between fitted latent positions in ordination,
#'   Generally the only scale on which stress values from RRGDM and other
#'   methods (e.g., (\code{vegan::metaMDS, MASS::isoMDS}))
#'   should be compared directly.
#'   \code{"response"} instead uses the link-transformed expected
#'   dissimilarity (\code{object$inv_link(dist)}), which is on the same
#'   scale as the observed dissimilarities used to fit the model (e.g.
#'   [0,1]-bounded for indices like Bray-Curtis or Jaccard, unbounded for
#'   metrics like Euclidean). This is a useful descriptive goodness-of-fit
#'   on the response/dissimilarity scale.
#'
#' @return A single numeric value: the stress, bounded between 0 (best) and 1 (worst).
#'
#' @importFrom stats dist isoreg lm fitted
#' @export
stress <- function(object,
                        type  = c("nonmetric", "metric"),
                        scale = c("link", "response")) {

  type  <- match.arg(type)
  scale <- match.arg(scale)

  if (!inherits(object, "rrGDM")) {
    cli::cli_abort("{.arg object} must be a fitted {.cls rrGDM} model.")
  }

  if (scale == "link" && type == "metric") {
    cli::cli_warn(c(
      "!" = "{.code scale = \"link\"} with {.code type = \"metric\"} regresses raw latent Euclidean distances linearly against observed dissimilarities.",
      "i" = "These are generally on different scales unless {.code link = \"identity\"}, so the resulting stress is typically not meaningful.",
      "*" = "Consider using {.code type = \"nonmetric\"} with {.code scale = \"link\"}, or {.code scale = \"response\"} with {.code type = \"metric\"} to assess fit directly on the dissimilarity scale."
    ))
  }

  L_hat <- object$obj$report()$L

  if (is.null(L_hat)) {
    cli::cli_abort(c(
      "x" = "Could not retrieve fitted latent coordinates {.field L} from {.arg object}."
    ))
  }

  raw_dist <- stats::dist(L_hat)

  fit_vec <- if (scale == "link") {
    as.vector(raw_dist)
  } else {
    as.vector(object$inv_link(raw_dist))
  }

  obs_vec <- object$Y_diss / object$Y_den

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

  ord <- order(obs_vec)
  obs_sort <- obs_vec[ord]
  fit_sort <- fit_vec[ord]

  dhat <- if (type == "nonmetric") {
    stats::isoreg(obs_sort, fit_sort)$yf
  } else {
    stats::fitted(stats::lm(fit_sort ~ obs_sort))
  }

  denom <- sum(fit_sort^2)
  if (denom == 0) {
    cli::cli_abort("Sum of squared modelled distances is 0; stress is undefined.")
  }

  sqrt(sum((dhat - fit_sort)^2) / denom)
}
