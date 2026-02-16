#' Initialize Latent Variables (u_sample) via PCoA
#'
#' @description
#' Provides starting values for site-level latent variables.
#' If \code{QR = TRUE}, rotates site loading matrix (u_sample) to be lower triangular (upper triangle = 0) with positive diagonal
#'
#' @param Y A \code{matrix} of community data (sites x species).
#' @param dist_method Character string for \code{vegan::vegdist}.
#' @param binary Logical; if \code{TRUE}, presence/absence is used.
#' @param k Integer; number of latent dimensions.
#' @param link Character string for the model's link function.
#' @param QR Logical; if \code{TRUE}, constrains the Site matrix to be lower triangular.
#'
#' @return A \code{list} containing \code{u_sample} and \code{var_explained}.
initial_loadings <- function(Y,
                             dist_method = "bray",
                             binary = FALSE,
                             k = 2,
                             link = "negexp",
                             QR = TRUE) {

  # 1. Distance & Linearization
  if (dist_method == 'sorensen') { binary <- TRUE; dist_method <- 'bray' }
  D <- as.matrix(vegan::vegdist(Y, method = dist_method, binary = binary))

  D_linear <- switch(link,
                     'negexp'   = { D_safe <- ifelse(D >= 1, 0.99, D); -log(1 - D_safe) },
                     'neggaus'  = { D_safe <- ifelse(D >= 1, 0.99, D); sqrt(-log(1 - D_safe)) },
                     'identity' = { D },
                     'sqrt'     = { sqrt(D) },
                     cli::cli_abort("Link not supported")
  )

  # 2. PCoA (Double Centering)
  D2 <- D_linear^2
  n <- nrow(D2)
  J <- diag(n) - (1/n) * matrix(1, n, n)
  B <- -0.5 * (J %*% D2 %*% J)
  e <- eigen(B, symmetric = TRUE)

  pos_eig <- pmax(e$values, 0)

  if (k >= length(pos_eig)) cli::cli_abort('failed to obtain initial values. k is larger than the dimensions positive eigenvalues in PCOA')

  scaling_matrix <- if (k == 1) matrix(sqrt(pos_eig[1])) else diag(sqrt(pos_eig[1:k]))

  # Raw PCoA Scores (Sites x LVs)
  u_raw <- e$vectors[, 1:k, drop = FALSE] %*% scaling_matrix
  var_explained <- pos_eig[1:k] / sum(pos_eig)

  # Constrain Site Matrix to be Lower Triangular
  if (QR && k > 1) {
    # Rotate so upper triangle = 0
    qr_res <- qr(t(u_raw))
    u_final <- t(qr.R(qr_res))

    # Positive diagonal
    R_diag <- diag(qr.R(qr_res)) #get signs
    sign_vec <- sign(R_diag)
    u_final <- sweep(u_final, 2, sign_vec, "*") # Flip columns

  } else if (QR && k == 1) {
    u_final <- u_raw * sign(u_raw[1,1])
  } else {
    u_final <- u_raw
  }

  colnames(u_final) <- paste0("LV", 1:ncol(u_final))

  return(list(
    u_sample = u_final,
    var_explained = sum(var_explained)
  ))
}


