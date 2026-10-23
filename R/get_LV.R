#' Extract latent coordinates and loadings from an rrGDM model
#'
#' @description
#' Extracts the main components of the fitted latent space of an \code{rrGDM}
#' object as a list of matrices (see also \code{\link{biplot.rrGDM}}).
#'
#' @param object A fitted model of class \code{rrGDM}.
#' @param rotate Logical indicating whether to apply varimax rotation
#'   (via \code{\link[stats]{varimax}}) before returning the results.
#'   Default is \code{TRUE}. The rotation is estimated from the conditional
#'   coordinates (\code{L}).
#'
#' @details
#' Latent positions of the samples are built as
#' \code{L_m = X \%*\% lambda} (marginal: determined by the predictors only) and
#' \code{L = L_m + u_sample (+ Z \%*\% u_re)} (conditional: also including
#' sample-level errors and any other random effects). Rotating all matrices
#' with the same rotation matrix preserves these relationships and all
#' pairwise distances between samples.
#'
#' @return
#' A \code{list} of matrices, with latent dimensions (\code{LV1}, \code{LV2},
#' ...) as columns:
#' \itemize{
#'   \item \code{lambda}: environmental loadings (predictors x latent dimensions).
#'   \item \code{u_sample}: sample-level errors (samples x latent dimensions).
#'   \item \code{u_re}: additional random effects, e.g. \code{(1 | group)}
#'     (groups x latent dimensions). Omitted if the model has none.
#'   \item \code{L}: conditional coordinates of the samples (samples x latent dimensions).
#'   \item \code{L_m}: marginal coordinates of the samples (samples x latent dimensions).
#'   \item \code{rotation}: the rotation matrix applied (the identity matrix if
#'     \code{rotate = FALSE}), useful to rotate other quantities consistently.
#' }
#'
#' @examples
#' set.seed(1)
#' # Simulate data
#' n <- 50 # sites
#' p <- 25 # species
#' env <- data.frame(temp = rnorm(n), moist = rnorm(n))
#' eta <- as.matrix(env) %*% matrix(rnorm(2 * p), 2, p)
#' Y <- matrix(rpois(n * p, exp(1 + 0.8 * eta)), n, p,
#'             dimnames = list(paste0("site", 1:n), paste0("sp", 1:p)))
#'
#' fit <- rrGDM(Y, env, k = 2, formula = ~ temp + moist,
#'              family = "binomial")
#'
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   # basic biplot
#'   biplot(fit)
#'   # add species correlations with latent variables
#'   biplot(fit, type = "marginal", species_cor = TRUE)
#'   # colour by predictor
#'   biplot(fit, col_by = "moist") + ggplot2::scale_colour_viridis_c()
#' }
#' @export
get_LV <- function(object, rotate = TRUE) {

  if (!inherits(object, "rrGDM")) {
    cli::cli_abort("{.arg object} must be a fitted {.cls rrGDM} model.")
  }
  if (!is.logical(rotate) || length(rotate) != 1 || is.na(rotate)) {
    cli::cli_abort("{.arg rotate} must be a single {.val TRUE}/{.val FALSE} value.")
  }

  #  EXTRACT REPORT
  report <- tryCatch(get_report(object), error = function(e) {
    cli::cli_abort(c(
      "x" = "Could not extract model report from {.arg object$obj}.",
      "i" = "Original error: {conditionMessage(e)}"
    ))
  })

  k <- object$k # num LV
  L <- as.matrix(report$L) # conditional coordinates

  # ROTATION
  rot <- diag(k)
  if (rotate && k > 1) {
    rot <- tryCatch(
      stats::varimax(L, normalize = FALSE)$rotmat,
      error = function(e) {
        cli::cli_warn(c(
          "!" = "Varimax rotation failed; returning unrotated scores instead.",
          "i" = "Original error: {conditionMessage(e)}"
        ))
        diag(k)
      })}

  # OUT
  lv_names <- paste0("LV", seq_len(k))
  sample_names <- rownames(object$Y)

  # Rotate and name
  rotate_mat <- function(mat, row_names) {
    out <- as.matrix(mat) %*% rot
    if (length(row_names) != nrow(out)) row_names <- NULL
    dimnames(out) <- list(row_names, lv_names)
    out
  }

  out <- list(
    lambda = rotate_mat(report$lambda, colnames(object$X_design)),
    u_sample = if (!is.null(report$u_sample)) rotate_mat(report$u_sample, sample_names) else NULL,
    u_re = if (!is.null(report$u) && length(report$u) > 0) rotate_mat(report$u, colnames(object$Z_design)) else NULL,
    L = rotate_mat(L, sample_names),
    L_m = rotate_mat(report$L_m, sample_names),
    rotation = rot
  )

  out[!vapply(out, is.null, logical(1))]
}
