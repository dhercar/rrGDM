#' Initialize Latent Variables (u_sample) via PCoA
#'
#' @description
#' Provides starting values for site-level latent variables and environmental loadings.
#' If \code{QR = TRUE}, rotates site loading matrix (u_sample) to be lower triangular (upper triangle = 0) with positive diagonal
#'
#' @param Y A \code{matrix} of community data (sites x species).
#' @param X_design If provided, a model \code{matrix} passed as \code{data} to vegan:: dbRDA.
#' @param dist_method Character string for \code{vegan::vegdist}.
#' @param binary Logical; if \code{TRUE}, presence/absence is used.
#' @param k Integer; number of latent dimensions.
#' @param link Character string for the model's link function.
#' @param QR Logical; if \code{TRUE}, constrains the Site loading matrix to be lower triangular.
#'
#' @return A \code{list} containing:
#' \itemize{
#'   \item \code{u_sample}: Matrix of site-level latent variables (n x k).
#'   \item \code{lambda_init}: Matrix of environmental loadings (p x k).
#' }

initial_loadings <- function(Y,
                             X_design = NULL,
                             init_method = NULL,
                             dist_method = "bray",
                             binary = FALSE,
                             k = 2,
                             link = "negexp",
                             QR = TRUE,
                             alpha = TRUE) {

  # Safe sign (no 0)
  sign2 <- function(x) ifelse(x >= 0, 1, -1)

  # Check if X_design is provided
  has_X <- !is.null(X_design) && (ncol(X_design) > 0)

  # 1. Distance & link
  if (dist_method == 'sorensen') { binary <- TRUE; dist_method <- 'bray' }
  D <- as.matrix(vegan::vegdist(Y, method = dist_method, binary = binary))

  D_linear <- switch(link,
                     'negexp'   = { D_safe <- ifelse(D >= 1, 0.99, D); -log(1 - D_safe) },
                     'neggaus'  = { D_safe <- ifelse(D >= 1, 0.99, D); sqrt(-log(1 - D_safe)) },
                     'identity' = { D },
                     'sqrt'     = { sqrt(D) },
                     cli::cli_abort("Link not supported"))



  if (init_method == 'PCOA') {

    ######## PCOA ########

  D2 <- D_linear^2
  n <- nrow(D2)
  J <- diag(n) - (1/n) * matrix(1, n, n)
  B <- -0.5 * (J %*% D2 %*% J)
  e <- eigen(B, symmetric = TRUE)
  pos_eig <- pmax(e$values, 0)
  var_explained <- pos_eig[1:k] / sum(pos_eig)

  # This is to scale initial loadings to compensat for missing dimensions
  #load_scale <- sqrt(sum(pos_eig) / sum(pos_eig[1:k]))

  # Check if sufficient positive eigen
  if (k >= sum(pos_eig > 0)) cli::cli_abort('failed to obtain initial values. k is larger than the dimensions positive eigenvalues in PCOA')

  # scale eigenvectors
  scaling_matrix <- if (k == 1) matrix(sqrt(pos_eig[1])) else diag(sqrt(pos_eig[1:k]))

  # Raw PCoA Scores (Sites x LVs)
  u_raw <- e$vectors[, 1:k, drop = FALSE] %*% scaling_matrix #* load_scale

  } else if (init_method == 'NMDS') {

    ######## NMDS ########

    nmds <- vegan::metaMDS(D_linear,
                           k = k,
                           autotransform = FALSE,
                           trace = FALSE)


    mean_D  <- mean(as.dist(D_linear))
    mean_nmds <- mean(dist(nmds$points))
    s       <- mean_D / mean_nmds
    u_raw <- nmds$points * s
  }
  #——————————————————————————————————————————————— RDA (DBRDA) —————————————————————————————————————————————

  if (has_X) {
    X_mat <- as.matrix(X_design)
    X_centered <- scale(X_mat, scale = FALSE)
    qx <- qr(X_centered)
    lambda_raw <- qr.coef(qx, u_raw)[, 1:k, drop = FALSE]
    # Residual Site Scores (u)
    u_res <- qr.resid(qx, u_raw)[, 1:k, drop = FALSE]
  } else {
    u_res <- u_raw
    lambda_raw <- NULL
  }

  # Constrain Matrix to be Lower Triangular
  if (QR && k > 1) {
    # Rotate so upper triangle = 0
    qr_decomp <- qr(t(u_res))
    u_final <- t(qr.R(qr_decomp))
    Q_rot <- qr.Q(qr_decomp) # Extract rotation

    # Change direction to get positive diagonal
    sign_vec <- sign2(diag(qr.R(qr_decomp)))
    u_final <- sweep(u_final, 2, sign_vec, "*")

    # Apply same rotation to lambda
    if (has_X) {
      lambda_final <- lambda_raw %*% Q_rot
      lambda_final <- sweep(lambda_final, 2, sign_vec, "*")
    } else {
      lambda_final <- NULL
    }

  } else if (QR && k == 1) {

    u_final <- u_res * sign2(u_res[1,1])
    if (has_X) lambda_final <- lambda_raw * sign2(u_res[1,1]) else lambda_final <- NULL

  } else { # NO QR
    u_final <- u_res
    lambda_final <- lambda_raw
  }

  #——————————————————————————————————————————————— alpha —————————————————————————————————————————————
  # CURRENTLY NOT USED
  # median residual distance as a robust initial for alpha
  # D_rr <- as.matrix(dist(u_raw))
  # alpha_init <- pmax(median(D_linear[lower.tri(D_linear)] - D_rr[lower.tri(D_rr)]), 1e-10)

  colnames(u_final) <- paste0("LV", 1:ncol(u_final))
  return(list(
    u_sample = u_final,
    lambda_init = lambda_final
 #   var_explained = sum(var_explained)
  ))
}


