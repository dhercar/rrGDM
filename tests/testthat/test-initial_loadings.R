test_that("QR rotation of initial loadings preserves latent distances and triangular structure", {
  set.seed(1)
  Y <- matrix(rpois(15 * 10, 4), 15, 10)
  X <- matrix(rnorm(15 * 3), 15, 3)
  k <- 2

  for (link in c("negexp", "neggaus", "identity")) {
    for (n_env in 1:3) {
      Xi <- X[, seq_len(n_env), drop = FALSE]

      unrot <- initial_loadings(Y, X_design = Xi, k = k, link = link,
                                constraint = 0, QR = FALSE)
      rot   <- initial_loadings(Y, X_design = Xi, k = k, link = link,
                                constraint = 1, QR = TRUE)

      coord_unrot <- Xi %*% unrot$lambda_init + unrot$u_sample
      coord_rot   <- Xi %*% rot$lambda_init + rot$u_sample

      expect_equal(as.matrix(dist(coord_rot)), as.matrix(dist(coord_unrot)),
                   tolerance = 1e-8, info = paste("link =", link, "n_env =", n_env))
      expect_true(all(abs(rot$lambda_init[upper.tri(rot$lambda_init)]) < 1e-8),
                  info = paste("link =", link, "n_env =", n_env))
    }
  }
})
