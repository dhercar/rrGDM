#' Predict method for rrGDM models
#'
#' @description
#' Obtains predictions from a fitted \code{rrGDM} object. This method can either project
#' new environmental data into the latent ordination space or calculate pairwise
#' dissimilarities between sites. It supports parametric bootstrapping to assess
#' parameter uncertainty and site-level dispersion.
#'
#' @details
#' Simulations (\code{boot = TRUE}) use R's random number generator: call
#' \code{set.seed()} beforehand for reproducible results.
#'
#' @param object An object of class \code{rrGDM}.
#' @param newdata A data frame containing new environmental predictor values.
#' If \code{NULL}, predictions are made for the original data used in model fitting.
#' @param boot Logical. Should predictions account for parameter and site-level uncertainty via simulation? Defaults to \code{TRUE}.
#' If \code{FALSE}, a single prediction is returned and \code{sample_param}, \code{sim_site},
#' \code{n_sim} and \code{probs} are ignored.
#' @param sample_param Logical. If \code{TRUE}, parameters (\code{lambda} and
#' \code{beta}) are sampled from the variance-covariance matrix. This requires standard
#' errors (fit the model with \code{se = TRUE} or run \code{\link{include_se}}). Defaults to \code{TRUE}.
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
#' @param conditional Logical. If \code{TRUE}, returns predictions from the fitted (conditional) latent positions of the samples used to fit the model, including sample-level and random-effect deviations. Only available when \code{newdata = NULL}. Sets \code{boot = FALSE}. Defaults to \code{FALSE}.
#' @param ... Ignored.
#'
#' @return
#' A \code{data.frame}.
#' If \code{ordination = TRUE}: the latent coordinates (columns \code{L.1}, ..., \code{L.k}),
#' the sample number (\code{sample}) and the columns of \code{newdata}. With
#' \code{boot = TRUE} there is one row per sample and draw, and a column \code{draw}.
#' If \code{ordination = FALSE}: one row per pair of rows of \code{newdata}, with the row
#' numbers of the pair (\code{s1}, \code{s2}) and either the predicted dissimilarity
#' (\code{dist}, if \code{boot = FALSE}) or one column per value of \code{probs}
#' (named \code{"2.5\%"}, \code{"50\%"}, ..., if \code{boot = TRUE}).
#'
#' @examples
#' set.seed(2)
#' # simulate data
#' n <- 25 # sites
#' p <- 25 # species
#'
#' env <- data.frame(temp = rnorm(n), moist = rnorm(n))
#' eta <- as.matrix(env) %*% matrix(rnorm(2 * p), 2, p)
#' Y <- matrix(rpois(n * p, exp(1 + 0.8 * eta)), n, p,
#'             dimnames = list(paste0("site", 1:n), paste0("sp", 1:p)))
#'
#' fit <- rrGDM(Y, env, k = 2,
#'              formula = ~ temp + moist,
#'              family = "binomial")
#' #  predict
#' pred <- predict(fit, boot = FALSE) # point predictions, all pairs
#' head(pred)
#'
#' # Predict ordination coordinates
#' pred_ord <- predict(fit,
#'                     ordination = TRUE, # latent positions
#'                     conditional = TRUE)
#'
#' head(pred_ord)
#' plot(pred_ord[,c(1,2)])
#'
#' # with newdata
#' new <- data.frame(temp = c(-1, 0, 1), moist = c(0, 0, 0))
#' pred_new <- predict(fit, newdata = new)
#' head(pred_new)
#'
#' # Predict position in ordination
#' pred_new_ord <- predict(fit,
#'                         newdata = new,
#'                         ordination = TRUE,
#'                         n_sim = 500)
#'
#' # Visualise sample-level variability
#' plot(pred_new_ord$L.1,
#'      pred_new_ord$L.2,
#'      xlab = 'L1',
#'      ylab = 'L2',
#'      asp = 1,
#'      col = c('blue', 'grey', 'red')[as.numeric(as.factor(pred_new_ord$temp))])
#'
#' @method predict rrGDM
#' @export
predict.rrGDM <- function(object,
                          newdata = NULL,
                          boot = TRUE,
                          sample_param = TRUE,
                          n_sim = 1000,
                          sim_site = TRUE,
                          ordination = FALSE,
                          type = c('response', 'link'),
                          probs = c(0.025, 0.5, 0.975),
                          conditional = FALSE,
                          ...) {

  # -------------------- BASIC CHECKS ---------------------------------

  type <- match.arg(type)

  flags <- list(boot = boot, sample_param = sample_param, sim_site = sim_site,
                ordination = ordination, conditional = conditional)
  for (nm in names(flags)) {
    if (!(isTRUE(flags[[nm]]) || isFALSE(flags[[nm]]))) {
      cli::cli_abort("{.arg {nm}} must be either {.code TRUE} or {.code FALSE}.")
    }
  }
  if (!is.numeric(n_sim) || length(n_sim) != 1 || !is.finite(n_sim) || n_sim < 1 || n_sim != round(n_sim)) {
    cli::cli_abort("{.arg n_sim} must be a single positive integer.")
  }
  if (!is.numeric(probs) || anyNA(probs) || any(probs < 0 | probs > 1)) {
    cli::cli_abort("{.arg probs} must be numeric, with values between 0 and 1.")
  }

  if (conditional) boot <- FALSE

  # Conditional predictions only available for fitted data
  if (conditional && !is.null(newdata)) {
    cli::cli_abort(c(
      "x" = "{.arg conditional = TRUE} is currently only available for the samples used to fit the model.",
      "*" = "Set {.arg newdata = NULL} to get conditional predictions for the fitted samples.",
      "*" = "Or set {.arg conditional = FALSE} to predict for new data: use {.code sim_site = TRUE} to include population-level error, and {.code sample_param = TRUE} to also include parameter uncertainty."
    ))
  }

  # Warning about biased average dissimilarities
  if ((!boot || !sim_site) && !ordination && !conditional) {
    cli::cli_inform(c(
      "i" = "Returning dissimilarities as distances between {.strong centroids}.",
      "!" = "Note: Distances between centroids are not equivalent to the expected (average) distances or dissimilarities between samples, which can be much larger.",
      "*" = "To calculate expected dissimilarities between samples, set {.code boot = TRUE} and {.code sim_site = TRUE}."
    ))
  }

  if (!conditional && length(object$map_re) > 0) {
    cli::cli_inform(c(
      "i" = "Random-effect terms such as {.code (1 | group)} are not included: only fixed effects and sample-level dispersion are used.",
      "*" = "Expected dissimilarities between samples will be understated if group-level variance is large."
    ))
  }

  inv_link <- if (type == 'link') function(x) x else object$inv_link

  # --------------------------- PREP DATA --------------------------------------------

  # If no new data is provided, use data in model object
  if (is.null(newdata)) {
    newdata <- object$X
    X_design <- object$X_design
    W_design <- object$W_design
    if (is.null(newdata)) newdata <- data.frame(row.names = seq_len(nrow(W_design)))

  } else {
    newdata <- as.data.frame(newdata)
    tt_fix <- stats::delete.response(object$terms_fix)

    mf_fix <- stats::model.frame(
      tt_fix,
      data = newdata,
      xlev = object$xlevels_fix,
      na.action = stats::na.pass)

    X_tmp <- stats::model.matrix(
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
    tt_sd <- stats::delete.response(object$terms_sd)

    mf_sd <- stats::model.frame(
      tt_sd,
      data = newdata,
      xlev = object$xlevels_sd,
      na.action = stats::na.pass)

    W_design <- stats::model.matrix(
      tt_sd,
      data = mf_sd,
      contrasts.arg = object$contrasts_sd)

    if (anyNA(X_design) || anyNA(W_design)) {
      cli::cli_abort("{.arg newdata} has missing values in the variables used by the model.")
    }
  }

  # --------------------------- GET COORDINATES --------------------------------------

  if (!ordination) {
    if (nrow(newdata) < 2) {
      cli::cli_abort(c(
        "x" = "Pairwise dissimilarity prediction requires at least 2 rows in {.arg newdata}.",
        "i" = "Set {.arg ordination = TRUE} to predict latent coordinates for a single row."
      ))
    }
    # Pairs of rows, in the order of as.vector(dist())
    comb <- t(utils::combn(nrow(newdata), 2))
    colnames(comb) <- c('s1', 's2')
  }

  if (boot) {

    if (sample_param) { # Parametric bootstrapping?
      if (is.null(object$sd_rep)) {
        cli::cli_abort(c(
          "x" = "Parameter bootstrapping requires standard errors, but this model has no {.field sd_rep}.",
          "i" = "Refit with {.code se = TRUE}, run {.fn include_se} on the fitted object, or set {.code sample_param = FALSE}."
        ))
      }
      check_pd_hessian(object)

      idx <- names(object$sd_rep$value) %in% c('lambda', 'beta')
      sim <- MASS::mvrnorm(n = n_sim, mu = object$sd_rep$value[idx], Sigma = object$sd_rep$cov[idx, idx]) # sample from vcov
      if (is.null(dim(sim))) sim <- t(sim) # n_sim = 1 returns a vector

      samples <- lapply(seq_len(nrow(sim)), function(x){
        sim_i <- sim[x,]
        lambda <- matrix(sim_i[names(sim_i) == 'lambda'], ncol = object$k)
        beta <- sim_i[names(sim_i) == 'beta']
        tau <- exp(W_design %*% beta)
        L <- X_design %*% lambda
        if (sim_site)  { # Include site level error (dispersion)? otherwise param uncertainty only
          L <- L + matrix(stats::rnorm(nrow(X_design) * object$k, 0, rep(tau, object$k)), ncol = object$k)
        }
        data.frame(draw = x, L = L)
      })

    } else {
      rep_par <- get_report(object)
      L_base <- X_design %*% rep_par$lambda
      tau <- exp(W_design %*% rep_par$beta)
      samples <- lapply(seq_len(n_sim), function(x){
        L <- L_base
        if (sim_site)  {
          L <- L + matrix(stats::rnorm(nrow(X_design) * object$k, 0, rep(tau, object$k)), ncol = object$k)
        }
        data.frame(draw = x, L = L)
      })
    }

    if (ordination) {
      # Return coordinates directly (all samples)
      samples <- lapply(samples, function(x) cbind(x, sample = seq_len(nrow(x)), newdata))
      samples <- do.call(rbind, samples)
      return(samples)

    } else {
      # Compute expected dissimilarities from coordinates
      samples_dist <- lapply(samples, function(x){
        inv_link(as.vector(stats::dist(x[,-1])))})

      # Return dissimilarities (summarised into quantiles)
      samples_dist <- do.call(cbind, samples_dist)
      samples_data <- t(matrix(apply(samples_dist, 1, stats::quantile, probs = probs), nrow = length(probs)))
      colnames(samples_data) <- names(stats::quantile(0, probs = probs))
      return(data.frame(comb, samples_data, check.names = FALSE))
    }

  } else {
    if (conditional) {
      L <- get_report(object)$L
    } else {
      L <- X_design %*% get_report(object)$lambda
    }

    if (ordination) {
      predict_ord <- data.frame(L = L, sample = seq_len(nrow(newdata)), newdata)
      return(predict_ord)

    } else {
      predict_dis <- data.frame(comb, dist = inv_link(as.vector(stats::dist(L))))
      return(predict_dis)
    }
  }
}
