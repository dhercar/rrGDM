#' Conditional correlation-based R2 (corR2) for a fitted rrGDM model
#'
#' @description
#' Computes the conditional \code{corR2} of a fitted \code{rrGDM} object: the
#' squared Pearson correlation between the observed dissimilarities used to
#' fit the model and the dissimilarities implied by the conditional latent
#' coordinates (i.e., including sample-level errors and any other random
#' effects), after applying the inverse link function.
#'
#' @details
#' The corR2 is useful to assess how well the latent space captures the
#' observed dissimilarities (e.g., to choose the number of latent dimensions
#' \code{k}), but it is \strong{not} a measure of the explanatory power of the
#' predictors. Because the conditional coordinates include the sample-level
#' errors, the corR2 can be high even if the predictors explain little. Like
#' \code{\link{stress}}, it increases as \code{k} increases, so it is
#' best used descriptively.
#'
#' @param object A fitted object of class \code{rrGDM}.
#' @param message Logical. Should a message clarifying the interpretation of
#'   the corR2 be printed? Default is \code{TRUE}.
#'
#' @return A single numeric value: the conditional corR2, between 0 and 1.
#'
#' @importFrom stats cor dist sd
#' @export
cor_r2 <- function(object, message = TRUE) {

  if (!inherits(object, "rrGDM")) {
    cli::cli_abort("{.arg object} must be a fitted {.cls rrGDM} model.")
  }
  if (!is.logical(message) || length(message) != 1 || is.na(message)) {
    cli::cli_abort("{.arg message} must be a single {.val TRUE}/{.val FALSE} value.")
  }

  L_hat <- get_report(object)$L

  if (is.null(L_hat)) {
    cli::cli_abort(c(
      "x" = "Could not retrieve fitted latent coordinates {.field L} from {.arg object}."
    ))
  }

  # Conditional fitted dissimilarities (same pair order as Y_diss)
  fp <- fit_pairs(object, L_hat)
  fit_vec <- as.vector(object$inv_link(fp$dist))
  obs_vec <- fp$obs

  if (length(obs_vec) != length(fit_vec)) {
    cli::cli_abort(c(
      "x" = "Observed and modelled dissimilarity vectors differ in length."
    ))
  }

  ok <- is.finite(obs_vec) & is.finite(fit_vec)
  if (!all(ok)) {
    n_drop <- sum(!ok)
    cli::cli_warn(c("!" = "{n_drop} non-finite pair{?s} removed before computing corR2."))
    obs_vec <- obs_vec[ok]
    fit_vec <- fit_vec[ok]
  }

  if (length(obs_vec) < 3) {
    cli::cli_abort("At least 3 finite pairs are required to compute corR2.")
  }

  if (stats::sd(obs_vec) == 0 || stats::sd(fit_vec) == 0) {
    cli::cli_abort("Observed or fitted dissimilarities are constant.")
  }

  r2 <- stats::cor(obs_vec, fit_vec)^2

  if (message) {
    cli::cli_inform(c(
      "i" = "{.fn cor_r2} returns the {.emph conditional} corR2: the squared correlation between observed and fitted dissimilarities, with fitted values computed from the latent coordinates including sample-level errors.",
      "*" = "It is useful to assess how well the latent space captures dissimilarities, but not the explanatory power of the predictors."
    ))
  }
  r2
}
