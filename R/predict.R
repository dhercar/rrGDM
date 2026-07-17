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
#' @param sim_site Logical. Should site-level random variation
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
#
predict.rrGDM <- function(object,
                          X_new = NULL,
                          boot = TRUE,
                          sample_param = TRUE,
                          n_sim = 1000,
                          sim_site = TRUE,
                          ordination = FALSE,
                          type = 'response',
                          probs = c(0.025, 0.5, 0.975),
                          conditional = FALSE,
                          ...) {

  if (conditional == TRUE) boot <- FALSE

  # -------------------- BASIC CHECKS ---------------------------------

  # Conditional predictions only available for fitted data
  if ((conditional == TRUE) && !is.null(X_new)) {
    cli::cli_abort(c(
      "x" = "{.arg conditional = TRUE} is currently only available for the samples used to fit the model.",
      "*" = "Set {.arg X_new = NULL} to get conditional predictions for the fitted samples.",
      "*" = "Or set {.arg conditional = FALSE} to predict for new data: use {.code sim_site = TRUE} to include population-level error, and {.code sample_param = TRUE} to also include parameter uncertainty."
    ))
  }

  # Warning about biased average dissimilarities
  if ((sim_site == FALSE) & (ordination == FALSE) & (conditional == FALSE)) {
    cli::cli_inform(c(
      "i" = "Returning dissimilarities as distances between {.strong centroids}.",
      "!" = "Note: Distances between centroids are not equivalent to the expected (average) distances or dissimilarities between samples, which can be much larger.",
      "*" = "To calculate expected dissimilarities between samples, set {.code boot = TRUE} and {.code sim_site = TRUE}."
    ))
  }

  # Inverse link
  inv_link <- if (type == 'link') function(x) x else object$inv_link

  # --------------------------- PREP DATA --------------------------------------------

  # If no new data data is provided, use data in model object
  if (is.null(X_new)) {
    X_new <- object$X
    X_design <- object$X_design
    W_design <- object$W_design

  } else {
    X_new <- data.frame(X_new)
    tt_fix <- delete.response(object$terms_fix)

    mf_fix <- model.frame(
      tt_fix,
      data = X_new,
      xlev = object$xlevels_fix,
      na.action = na.pass)

    X_tmp <- model.matrix(
      tt_fix,
      data = mf_fix,
      contrasts.arg = object$contrasts_fix)

    # Remove intercept if needed
    if ("(Intercept)" %in% colnames(X_tmp)) {
      intercept_idx <- which(colnames(X_tmp) == "(Intercept)")
      X_design <- X_tmp[, -intercept_idx, drop = FALSE]
      attr(X_design, "assign") <- attr(X_tmp, "assign")[-intercept_idx]
      attr(X_design, "contrasts") <- object$contrasts_fix
    } else {
      X_design <- X_tmp
    }

    # --- DISPERSION MODEL ---
    tt_sd <- delete.response(object$terms_sd)

    mf_sd <- model.frame(
      tt_sd,
      data = X_new,
      xlev = object$xlevels_sd,
      na.action = na.pass)

    W_design <- model.matrix(
      tt_sd,
      data = mf_sd,
      contrasts.arg = object$contrasts_sd)
  }

  # --------------------------- GET COORDINATES --------------------------------------

  if (!ordination && nrow(X_new) < 2) {
    cli::cli_abort(c(
      "x" = "Pairwise dissimilarity prediction requires at least 2 rows in {.arg X_new}.",
      "i" = "Set {.arg ordination = TRUE} to predict latent coordinates for a single row."
    ))
  }

  if (boot) { # Return samples

    if (sample_param == TRUE) { # Parametric bootstrapping?
      if (is.null(object$sd_rep)) {
        cli::cli_abort(c(
          "x" = "Parameter bootstrapping requires standard errors, but this model has no {.field sd_rep}.",
          "i" = "Refit with {.code se = TRUE}, or set {.code sample_param = FALSE}."
        ))
      }

      idx <- names(object$sd_rep$value) %in% c('lambda', 'beta')
      sim <- MASS::mvrnorm(n = n_sim, mu = object$sd_rep$value[idx], Sigma = object$sd_rep$cov[idx, idx]) # sample from vcov
      samples <- lapply(1:nrow(sim), function(x){ # Expected coordinates
        sim_i <- sim[x,]
        lambda <- matrix(sim_i[names(sim_i) == 'lambda'], ncol = object$k)
        beta <- sim_i[names(sim_i) == 'beta']
        tau <- exp(W_design %*% beta)
        L <- X_design %*% lambda
        if (sim_site)  { # Include site level error (dispersion)? otherwise param uncertainty only
          L = L  + matrix(rnorm(nrow(X_design) * object$k, 0, rep(tau, object$k)), ncol = object$k)
        }
        data.frame(draw = x, L = L)
      })

    } else if (sample_param == FALSE)  {
      samples <- lapply(1:n_sim, function(x){
        L <- X_design %*% object$obj$report()$lambda # Point estimate of lambda
        if (sim_site)  {  # Include site level error (dispersion)? otherwise expected value
          tau <- exp(W_design %*% object$obj$report()$beta)
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

    } else if (ordination == FALSE) {
      # Compute expected dissimilarities from coordinates
      comb <- t(combn(nrow(X_new), 2))
      colnames(comb) <- c('s1', 's2')

      samples_dist <- lapply(samples, function(x){
        dist <- as.matrix(dist(x[,-1]))
        inv_link(dist[comb])})

      # Return dissimilarities (summarised into quantiles)
      samples_dist <- do.call(cbind, samples_dist)
      samples_data <- t(apply(samples_dist, 1, function(x) quantile(x, probs = probs)))
      return(cbind(comb, samples_data))
    }

  } else {
    if (conditional == TRUE) {
     L <- object$obj$report()$L
    } else if (conditional == FALSE) {
     L <- X_design %*% object$obj$report()$lambda
    }
      if (ordination == TRUE) {
        predict_ord <- data.frame(L = L, sample = 1:nrow(X_new), X_new)
        return(predict_ord)

      } else if (ordination == FALSE) {
        comb <- t(combn(nrow(X_new), 2))
        colnames(comb) <- c('s1', 's2')
        predict_dis <- data.frame(dist = inv_link(as.matrix(dist(L))[comb]), comb)
        return(predict_dis)
    }
  }
}

