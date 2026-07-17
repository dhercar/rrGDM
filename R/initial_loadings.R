#' Initialize Latent Variables via PCoA + RDA
#'
#' @description
#' Provides starting values for site-level latent variables and environmental
#' loadings from a PCoA or dbRDA. The top k constrained axes (ranked by fitted score variance)
#' are used, and \code{u_sample} is initialised from the corresponding *residual*
#' (unexplained) scores -- i.e. what is not captured by the environmental
#' predictors -- since that is what the sample-level error term represents.
#' If no predictors are present, the top k PCoA axes are used directly as
#' \code{u_sample} (there is nothing to partition out).
#' If \code{QR = TRUE}, rotates the site loading matrix to be lower triangular
#' with positive diagonal, and applies the same rotation to the environmental
#' loadings so that \code{X \%*\% lambda + u_sample} is unchanged.
#'
#' @param Y A \code{matrix} of community data (sites x species).
#' @param X_design If provided, a model \code{matrix} of predictors.
#' @param dist_method Character string passed to \code{vegan::vegdist}.
#' @param binary Logical; if \code{TRUE}, presence/absence is used.
#' @param k Integer; number of latent dimensions.
#' @param link Character string for the model's link function.
#' @param QR Logical; if \code{TRUE}, constrains the site loading matrix to
#'   be lower triangular with positive diagonal.
#' @param constraint Integer. Which matrix receives the triangular identifiability constraint
#'
#' @return A \code{list} containing:
#' \itemize{
#'   \item \code{u_sample}: Matrix of site-level latent variables (n x k),
#'     representing the residual (environmentally-unexplained) component when
#'     predictors are supplied.
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
                             constraint = 2,
                             QR  = TRUE) {

  if (constraint == 0) QR <- FALSE

  # Safe sign to avoid exact 0s
  sign2 <- function(x) ifelse(x >= 0, 1, -1)

  # Check if predictors are provided
  has_X <- !is.null(X_design) && (ncol(X_design) > 0)


  if (constraint == 1 && !has_X) {
    cli::cli_abort(c(
      "x" = "{.arg constraint == 1} requires environmental predictors."
    ))
  }

  if (!(constraint %in% c(0, 1, 2))) {
    cli::cli_abort("{.arg constraint} must be either {.val 0} (no constraint), {.val 1} or {.val 2}.")
  }

  # Dissimilarity matrix and link transformation
  if (dist_method == "sorensen") {
    binary <- TRUE
    dist_method <- "bray"
  }

  # Compute fill diss. matrix, ignoring num_den specifications (e.g., for binomial family)
  D <- as.matrix(vegan::vegdist(Y, method = dist_method, binary = binary))

  # Cap 1s to avoid crash
  D_linear <- switch(link,
                     negexp  = {
                       d_max  <- max(D[D < 1], 0)                          # largest non-1 dissimilarity
                       d_cap  <- max(0.99, d_max + 0.5 * (1 - d_max))      # whichever is larger
                       D_safe <- ifelse(D >= 1, d_cap, D)
                       - log(1 - D_safe)},
                     neggaus = {
                       d_max  <- max(D[D < 1], 0)
                       d_cap  <- max(0.99, d_max + 0.5 * (1 - d_max))
                       D_safe <- ifelse(D >= 1, d_cap, D)
                       sqrt(-log(1 - D_safe))},
                     identity = D,
                     cli::cli_abort("Link {.val {link}} is not supported."))

  # ------------------------------------ PCOA ------------------
  n <- nrow(D_linear)
  J <- diag(n) - matrix(1 / n, n, n)
  B <- -0.5 * (J %*% (D_linear^2) %*% J)
  e <- eigen(B, symmetric = TRUE)

  # Retain only positive VAP
  pos_idx <- which(e$values > 0)
  n_pos   <- length(pos_idx)

  if (n_pos < k) {
    cli::cli_abort(c(
      "Insufficient positive eigenvalues to initialise {k} latent dimension{?s}.",
      "i" = "{n_pos} positive eigenvalue{?s} {?was/were} found, but {k} {?is/are} required.",
      ">" = "Possible causes:",
      " " = "- {.var k} is large relative to the number of sites.",
      " " = "- The dissimilarity matrix is very sparse (e.g., with many 0s or 1s)).",
      "!" = "Try reducing {.var k}, using more sites if available, or a different {.var dist_method}."
    ))
  }

  eig_pos  <- e$values[pos_idx]
  evec_pos <- e$vectors[, pos_idx, drop = FALSE]
  scale_correction = sqrt(sum(eig_pos) / sum(eig_pos[1:k])) # correction for reduced var. when using k dims.

  # PCOA scores
  U_pcoa <- evec_pos %*% diag(sqrt(eig_pos), nrow = n_pos) # site-level scores


  # ------------------------------------ RDA ------------------
  if (has_X) {
    X_mat <- as.matrix(X_design)
    qx <- qr(X_mat)

    U_fitted <- qr.fitted(qx, U_pcoa)
    U_resid  <- qr.resid(qx, U_pcoa)

    # Eigendecomposition of the fitted scores
    rda <- eigen(crossprod(U_fitted), symmetric = TRUE)
    V_k <- rda$vectors[, 1:k, drop = FALSE] # take top k axes

    u_env_k   <- U_fitted %*% V_k   # environmentally-explained site scores
    u_raw <- U_resid  %*% V_k   # residual

    lambda_raw <- qr.coef(qx, u_env_k)
    lambda_raw[is.na(lambda_raw)] <- 0

  } else {
    # No predictors, use top k PCoA axes
    u_raw      <- U_pcoa[, 1:k, drop = FALSE]
    lambda_raw <- NULL
  }

  # -------------ROTATE TO MEET CONTRAINTS 1 OR 2 ------------------
  if (QR && k > 1) {
    if (constraint == 1) {
     # Rotate lambda to lower tri, and apply rot. to u_sample
     p_lambda <- nrow(lambda_raw)
     rank_r <- min(k, p_lambda) # deals with the case of p < k

     qr_decomp <- qr(t(lambda_raw))
     Q_full    <- qr.Q(qr_decomp, complete = TRUE)
     R_econ    <- qr.R(qr_decomp)

     # If p < lambda, pad R with 0s
     if (k > p_lambda) {
       R_full <- rbind(R_econ, matrix(0, k - p_lambda, p_lambda))
     } else {
       R_full <- R_econ
     }

     sign_vec <- rep(1, k)
     sign_vec[seq_len(rank_r)] <- sign2(diag(R_econ))

     lambda_final <- sweep(t(R_full), 2, sign_vec, "*")
     u_final      <- sweep(u_raw %*% Q_full, 2, sign_vec, "*")

    } else if (constraint == 2) {
      # Rotate u_sample to lower tri, apply rot. to lambda
      # This part does not check dims, assumes n sample > k
      qr_decomp <- qr(t(u_raw))
      R_mat     <- qr.R(qr_decomp)
      sign_vec  <- sign2(diag(R_mat))

      u_final <- sweep(t(R_mat), 2, sign_vec, "*")

      if (has_X) {
        Q_rot <- qr.Q(qr_decomp)
        lambda_final <- sweep(lambda_raw %*% Q_rot, 2, sign_vec, "*")
      } else {
        lambda_final <- NULL
      }
    }

  } else if (QR && k == 1) {

    if (constraint == 1) {
      s            <- sign2(lambda_raw[1, 1])
      lambda_final <- lambda_raw * s
      u_final      <- u_raw * s
    } else {
      s            <- sign2(u_raw[1, 1])
      u_final      <- u_raw * s
      lambda_final <- if (has_X) lambda_raw * s else NULL
    }

  } else {
    u_final      <- u_raw
    lambda_final <- lambda_raw
  }


  colnames(u_final) <- paste0("LV", 1:k)
  if (!is.null(lambda_final)) lambda_final <- lambda_final * scale_correction

  list(                                                              # FIX: whole return was missing
    u_sample    = u_final * scale_correction,
    lambda_init = lambda_final)
}
