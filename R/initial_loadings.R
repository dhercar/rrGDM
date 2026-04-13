#' Initialize Latent Variables via PCoA + RDA
#'
#' @description
#' Provides starting values for site-level latent variables and environmental
#' loadings from a PCoA or dbRDA. If predictors are present, the top k constrained axes (ranked by fitted score variance)
#' are used. If no predictors are present, the top k PCoA axes are used directly.
#' If \code{QR = TRUE}, rotates the site loading matrix to be lower triangular
#' with positive diagonal.
#'
#' @param Y A \code{matrix} of community data (sites x species).
#' @param X_design If provided, a model \code{matrix} of predictors.
#' @param dist_method Character string passed to \code{vegan::vegdist}.
#' @param binary Logical; if \code{TRUE}, presence/absence is used.
#' @param k Integer; number of latent dimensions.
#' @param link Character string for the model's link function.
#' @param QR Logical; if \code{TRUE}, constrains the site loading matrix to
#'   be lower triangular with positive diagonal.
#'
#' @return A \code{list} containing:
#' \itemize{
#'   \item \code{u_sample}: Matrix of site-level latent variables (n x k).
#'   \item \code{lambda_init}: Matrix of environmental loadings (p x k), or
#'     \code{NULL} if no predictors.
#' }
#' @noRd
initial_loadings <- function(Y,
                             X_design  = NULL,
                             dist_method = "bray",
                             binary  = FALSE,
                             k  = 2,
                             link   = "negexp",
                             QR  = TRUE) {

  # Safe sign to avoid 0
  sign2 <- function(x) ifelse(x >= 0, 1, -1)

  # Check if predictors are provided
  has_X <- !is.null(X_design) && (ncol(X_design) > 0)


  # Dissimilarity matrix and link transformation
  if (dist_method == "sorensen") {
    binary <- TRUE
    dist_method <- "bray"
  }

  D <- as.matrix(vegan::vegdist(Y, method = dist_method, binary = binary))

  D_linear <- switch(link,
                     negexp  = { D_safe <- ifelse(D >= 1, 0.99, D); -log(1 - D_safe) },
                     neggaus = { D_safe <- ifelse(D >= 1, 0.99, D); sqrt(-log(1 - D_safe)) },
                     identity = D,
                     sqrt     = sqrt(D),
                     cli::cli_abort("Link {.val {link}} is not supported."))


  # PCoA via double-centred inner product matrix
  n <- nrow(D_linear)
  J <- diag(n) - matrix(1 / n, n, n)
  B <- -0.5 * (J %*% (D_linear^2) %*% J)
  e <- eigen(B, symmetric = TRUE)

  # Retain only axes with strictly positive eigenvalues
  pos_idx <- which(e$values > 0)
  n_pos   <- length(pos_idx)

  # Check if there are sufficient non-negative dimensions
  if (n_pos < k) {
    cli::cli_abort(c(
      "Insufficient positive eigenvalues to initialise {k} latent dimension{?s}.",
      "i" = "{n_pos} positive eigenvalue{?s} {?was/were} found, but {k} {?is/are} required.",
      ">" = "Possible causes:",
      " " = "- {.var k} is large relative to the number of sites.",
      " " = "- The dissimilarity matrix is very sparse (e.g., with many 0s or 1s)).",
      "!" = "Try reducing {.var k}, using more sites if available, or choosing a different {.var dist_method}."
    ))
  }

  eig_pos  <- e$values[pos_idx]
  evec_pos <- e$vectors[, pos_idx, drop = FALSE]

  scale_correction = sqrt(sum(eig_pos) / sum(eig_pos[1:k]))

  # PCOA scores
  U_pcoa <- evec_pos %*% diag(sqrt(eig_pos), nrow = n_pos)

  # RDA ON PCOA scores if predictors are present
  if (has_X) {
    X_mat <- as.matrix(X_design)
    X_centered <- X_mat # Not centring to respect LV scale (?)
    qx <- qr(X_centered)

    # Fitted (constrained) scores for all PCoA axes
    U_fitted <- qr.fitted(qx, U_pcoa)

    # Eigen-decomposition of the fitted scores cross-product
    S   <- crossprod(U_fitted)
    rda <- eigen(S, symmetric = TRUE)
    V_k <- rda$vectors[, 1:k, drop = FALSE] # Top k axes
    u_raw      <- U_fitted %*% V_k  # site scores

    lambda_raw <- qr.coef(qx, u_raw)
    lambda_raw[is.na(lambda_raw)] <- 0


  } else {
    # No predictors: use top k PCoA axes directly
    u_raw      <- U_pcoa[, 1:k, drop = FALSE]
    lambda_raw <- NULL
  }


  # QR rotation: lower-triangular u with positive diagonal contraint
  if (QR && k > 1) {
    qr_decomp <- qr(t(u_raw))
    R_mat     <- qr.R(qr_decomp)
    sign_vec  <- sign2(diag(R_mat))

    u_final <- sweep(t(R_mat), 2, sign_vec, "*")

    if (has_X) {
      Q_rot        <- qr.Q(qr_decomp)
      lambda_final <- sweep(lambda_raw %*% Q_rot, 2, sign_vec, "*")
    } else {
      lambda_final <- NULL
    }

  } else if (QR && k == 1) {
    s            <- sign2(u_raw[1, 1])
    u_final      <- u_raw * s
    lambda_final <- if (has_X) lambda_raw * s else NULL # Positive diagonal

  } else {
    u_final      <- u_raw
    lambda_final <- lambda_raw
  }

  colnames(u_final) <- paste0("LV", 1:k)

  if (!is.null(lambda_final)) lambda_final <- lambda_final * scale_correction

  list(
    u_sample    = u_final*scale_correction,
    lambda_init = lambda_final
  )
}
