#' Predict method for rrGDM models
#'
#' @description
#' Obtains predictions from a fitted \code{rrGDM} object. This method can either project
#' new environmental data into the latent ordination space or calculate pairwise
#' dissimilarities between sites. It supports parametric bootstrapping to assess
#' parameter uncertainty and site-level dispersion.
#'
#' @param object An object of class \code{rrGDM}.
#' @param X_new A data frame containing new environmental predictor values.
#' If \code{NULL}, predictions are made for the original data used in model fitting.
#' @param boot Logical. Should predictions account for parameter and site-level uncertainty via simulation? Defaults to \code{TRUE}.
#' @param sample_param Logical. If \code{TRUE}, parameters (\code{lambda} and
#' \code{beta}) are sampled from the variance-covariance matrix.
#' @param n_sim Integer specifying the number of simulation draws. Defaults to 1000.
#' @param site_error Logical. Should site-level random variation
#' be included in the predictions? Defaults to \code{TRUE}.
#' @param ordination Logical. If \code{TRUE}, returns coordinates in the latent
#' space. If \code{FALSE}, returns pairwise dissimilarities. Defaults to \code{FALSE}.
#' @param type Character string; \code{"response"} returns dissimilarities on the
#' observed scale (e.g., 0-1), while \code{"link"} returns Euclidean distances
#' in the latent space. Defaults to \code{"response"}.
#' @param probs Numeric vector of probabilities for calculating quantiles from
#' the simulated distribution. Defaults to \code{c(0.025, 0.5, 0.975)}.
#' @param ... not used
#'
#' @return
#' If \code{ordination = TRUE}, a \code{data.frame} of latent coordinates.
#' If \code{ordination = FALSE}, a \code{data.frame} of pairwise dissimilarities, with columns for the requested quantiles.
#'
#' @method predict rrGDM
#' @importFrom MASS mvrnorm
#' @import Matrix cli
#' @export
predict.rrGDM <- function(object,
                          X_new = NULL,
                          boot = TRUE,
                          sample_param = TRUE,
                          n_sim = 1000,
                          site_error = TRUE,
                          ordination = FALSE,
                          type = 'response',
                          probs = c(0.025, 0.5, 0.975),
                          ...) {

  if ((site_error == FALSE) & (ordination == FALSE)) {
    cli::cli_inform(c(
      "i" = "Returning dissimilarities as distances between {.strong centroids}.",
      "!" = "Note: Distances between centroids are not equivalent to the expected (average) distances or dissimilarities between samples, which can be much larger.",
      "*" = "To calculate expected dissimilarities between samples, set {.code boot = TRUE} and {.code site_error = TRUE}."
    ))
  }

  # Inverse link
  inv_link <- switch(object$link,
                     'negexp' = function(x) 1 - exp(-x),
                     'neggaus' = function(x) 1 - exp(-(x^2)),
                     'identity' = function(x) x)

  if (type == 'link') inv_link <- function(x) x

  # Prep data
  if (is.null(X_new)) {
    # Use original data
    X_new <- object$X
    X_design <- object$X_design
    W_design <- object$W_design

  } else {
    # Use new data
    X_new <- data.frame(X_new)
    W_new <- data.frame(X_new)
    X_design <- model.matrix(object$formula_fix, data = X_new)
    X_design <- X_design[, !colnames(X_design) %in% "(Intercept)", drop = FALSE] # remove intercept
    W_design <- model.matrix(object$formula_sd, data = X_new)
  }

  if (boot) { # Return samples
    if (sample_param) { # Parametric bootstrapping?
      sim <- MASS::mvrnorm(n = n_sim, mu = object$sd_rep$value, Sigma = object$sd_rep$cov) # sample from vcov
      samples <- lapply(1:nrow(sim), function(x){ # Expected coordinates
        sim_i <- sim[x,]
        lambda <- matrix(sim_i[names(sim_i) == 'lambda'], ncol = object$k)
        beta <- sim_i[names(sim_i) == 'beta']
        tau <- exp(W_design %*% beta)
        L <- X_design %*% lambda
        if (site_error)  { # Include site level error (dispersion)? otherwise param uncertainty only
          L = L  + matrix(rnorm(nrow(X_design) * object$k, 0, rep(tau, object$k)), ncol = object$k)
        }
        data.frame(draw = x, L = L)
      })

    } else {
      samples <- lapply(1:n_sim, function(x){
        L <- X_design %*% object$obj$report()$lambda
        tau <- exp(W_design %*% object$obj$report()$beta)
        if (site_error)  { # Include site level error (dispersion)? otherwise expected value
          L = L  + matrix(rnorm(nrow(X_design) * object$k, 0, rep(tau, object$k)), ncol = object$k)
        }
        data.frame(draw = x, L = L)
      })
    }

    if (ordination == TRUE) {
      # Return coordinates directly (all samples)
      samples <- lapply(samples, function(x) cbind(x, sample = 1:nrow(x), X_new))
      samples <- do.call(rbind, samples)
      return(samples)
    } else {
      # Compute expected dissimilarities from coordinates
      comb <- t(combn(nrow(X_new), 2))
      colnames(comb) <- c('s1', 's2')
      samples_dist <- lapply(samples, function(x){
        dist <- as.matrix(dist(x[,-1]))
        inv_link(dist[comb])
      })

      # Return dissimilarities (summarised into quantiles)
      samples_dist <- do.call(cbind, samples_dist)
      samples_data <- t(apply(samples_dist, 1, function(x) quantile(x, probs = probs)))
      return(cbind(comb, samples_data))

    }

  } else { # Return expected values

    # Expected coordinates
    L <- X_design %*% object$obj$report()$lambda

    if (ordination) {
      # Return expected coordinates
      predict_ord <- data.frame(L = L, sample = 1:nrow(X_new), X_new)
      return(predict_ord)
    } else {
      comb <- t(combn(nrow(X_new), 2))
      colnames(comb) <- c('s1', 's2')
      predict_dis <- data.frame(dist = inv_link(as.matrix(dist(L))[comb]), comb)
      return(predict_dis)
    }
  }
}


