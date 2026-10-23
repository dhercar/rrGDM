#' Fit a reduced-rank generalised dissimilarity model
#' @description Fit a reduced-rank Generalised Dissimilarity Model using TMB
#' @param Y A matrix or data frame containing the response variable (usually species abundance matrix).
#' @param X A data frame containing the predictor variables.
#' @param k An integer specifying the number of latent variables to be used in the model. Default is 2.
#' @param formula Formula specifying the model to be fitted.
#' @param formula_sd Formula specifying the model for the dispersion parameter. Default is ~1 (intercept only).
#' @param constraint Integer. Specifies identifiability constraints. 0: no constraint; 1: Constraint applied to environmental loadings matrix (\eqn{\Lambda}); 2: Constraint applied to site random effects matrix (U). Default is 1.
#' @param dist_method A character string specifying the distance method to be used (e.g., "euclidean", "bray"). Default is "bray".
#' @param trans Function used to transform the response matrix before computing dissimilarities (e.g., log, sqrt). Default is no transformation.
#' @param sparse Optional object of class \code{sparse_diss} (see \code{\link{sparse_diss}}). If supplied, the model is fitted using only the pairs of samples listed in \code{sparse$edges}, instead of all pairs. It must be created from the same samples as \code{Y}. Default is NULL (all pairs).
#' @param binary Logical indicating whether the response variable should be treated as binary (e.g., presence-absence). Default is FALSE.
#' @param scale_diss A numeric vector of length two indicating a new scale for dissimilarities. Useful for distributions bounded between 0 and 1 (e.g., beta family). Ignored, with a warning, for the \code{"binomial"} and \code{"beta.binomial"} families.
#' @param family A character string specifying the family of the response variable: \code{"gaussian"}, \code{"binomial"}, \code{"beta"} or \code{"beta.binomial"}. Default is \code{"gaussian"}.
#' @param link Character string specifying the link function to be used: \code{"identity"}, \code{"negexp"} or \code{"neggaus"}. Default is 'identity' for Gaussian distribution, and 'negexp' for beta and binomial distributions.
#' @param se Logical indicating whether to calculate standard errors (TRUE) or not (FALSE). Default is TRUE.
#' @param sample_error Logical. If TRUE, includes site-level random effects to account for variation not explained by predictors. Default is TRUE.
#' @param seed An integer passed to \code{set.seed()} at the start of the fit, so it changes the state of R's random number generator. Default is NULL (no seed is set).
#' @param n_starts Integer. Number of optimisation starts. The first start uses the values from \code{init_method} (or \code{start_fit}); each further start adds Gaussian noise (SD \code{init_jitter}) to the best solution so far (\code{previous_best = TRUE}) or to the initial values (\code{FALSE}). With more than one start, only fits that converged and have a positive-definite Hessian are eligible as the final fit. Default is 1.
#' @param init_jitter Numeric. Standard deviation of the Gaussian noise added to the starting values of every start after the first. Only used when \code{n_starts > 1}. Default is 0.1.
#' @param previous_best Logical. Controls where the jitter for multiple starts after the first is applied.
#' If \code{TRUE} (default), each start jitters around the best solution found so far (perturb-and-reoptimise). If \code{FALSE}, every start jitters around the initial values (from \code{init_method}), so starts do not depend on earlier fits. Only relevant when \code{n_starts > 1}.
#' @param init_method A character string specifying a method used to obtain initial values for optimisation. Currently, only \code{'PCOA'} and \code{'dbRDA'} are available; \code{NULL} starts from small random values. Default is \code{'dbRDA'}. Ignored if \code{start_fit} is supplied.
#' @param start_fit Optional fitted \code{rrGDM} model used to initialise the optimisation (e.g., a fast sparse fit before refitting with more pairs). It must have the same samples, \code{k}, \code{family}, \code{link}, \code{constraint}, predictors and random-effect structure as the new model; the set of pairs (\code{sparse}) may differ. If supplied, \code{init_method} is ignored. Default is NULL.
#' @param control A list of control parameters passed to \code{\link[stats]{nlminb}}. \code{iter.max} defaults to 1000.
#' @param trace Logical. If \code{TRUE}, TMB prints the optimisation progress and the objective of each start is reported. Default is \code{FALSE}.
#' @return An object of class \code{"rrGDM"}: a list with the components below.
#'
#' \strong{Estimation}
#' \describe{
#'   \item{\code{opt}}{Result of \code{\link[stats]{nlminb}} for the fixed-effect
#'     parameters}
#'   \item{\code{obj}}{The \pkg{TMB} object created by \code{\link[TMB]{MakeADFun}}}
#'   \item{\code{best_par}}{Numeric vector with estimated parameters at the
#'     optimum, including the modes of random effects}
#'   \item{\code{sd_rep}}{An \code{\link[TMB]{sdreport}} object with standard
#'     errors and the covariance matrix of the estimates, or \code{NULL} if
#'     \code{se = FALSE} or the standard errors could not be computed. See
#'     \code{\link{include_se}}}
#'   \item{\code{pd_hess}}{\code{TRUE} if the Hessian at the optimum is positive
#'     definite and not nearly singular, \code{FALSE} if not, and \code{NA} if
#'     there is no \code{sd_rep}}
#' }
#'
#' \strong{Data}
#' \describe{
#'   \item{\code{Y}, \code{X}}{The community matrix and the predictor data as
#'     supplied}
#'   \item{\code{Y_diss}, \code{Y_den}}{Observed values for the fitted pairs.
#'     For the \code{"gaussian"} and \code{"beta"} families, \code{Y_diss} is
#'     the dissimilarity (rescaled to \code{scale_diss} if given) and
#'     \code{Y_den} is 1. For \code{"binomial"} and \code{"beta.binomial"},
#'     \code{Y_diss} is the numerator of the dissimilarity and
#'     \code{Y_den} the denominator, so that
#'     \code{Y_diss / Y_den} is the dissimilarity.}
#'   \item{\code{pairs}}{Two-column matrix with the row numbers of \code{Y} of the
#'     two samples in each fitted pair, in the same order as \code{Y_diss}.}
#'   \item{\code{sparse}}{The \code{\link{sparse_diss}} object supplied in
#'     \code{sparse}, or \code{NULL} if all pairs were used.}
#'   \item{\code{obs_full}}{Only for fits with \code{sparse}: observed
#'     dissimilarities of \emph{all} pairs of samples, in the order of
#'     \code{\link[stats]{dist}}}
#' }
#'
#' \strong{Design matrices}
#' \describe{
#'   \item{\code{X_design}}{Design matrix of the predictors in \code{formula}
#'     (samples x predictors, without intercept; zero columns if there are no
#'     predictors).}
#'   \item{\code{W_design}}{Design matrix of \code{formula_sd} (samples x terms,
#'     with intercept)}
#'   \item{\code{Z_design}}{Sparse indicator matrix for the random intercepts
#'   in \code{formula} (samples x total number of groups).}
#'   \item{\code{map_re}}{Integer vector giving, for each column of
#'     \code{Z_design}, the random-effect term it belongs to (index into
#'     \code{re_vars}). \code{integer(0)} if there are no random effects.}
#'   \item{\code{re_vars}}{Names of the grouping variables with random
#'     intercepts. \code{character(0)} if there are none.}
#' }
#'
#' \strong{Model specification}
#' \describe{
#'   \item{\code{k}, \code{family}, \code{link}}{Number of latent dimensions,
#'     family and link used}
#'   \item{\code{inv_link}}{Function that maps distances in the latent space to
#'     expected dissimilarities (identity, \eqn{1 - e^{-d}} or
#'     \eqn{1 - e^{-d^2}}).}
#'   \item{\code{constraint}}{Identifiability constraint used}
#'   \item{\code{formula}, \code{formula_sd}}{The formulas supplied.}
#'   \item{\code{formula_fix}}{\code{formula} without the random-effect terms
#'     (\code{~0} if there are no predictors).}
#'   \item{\code{terms_fix}, \code{terms_sd}}{\code{\link[stats]{terms}} objects
#'     of the fixed-effect and \code{formula_sd} models.}
#'   \item{\code{xlevels_fix}, \code{xlevels_sd}, \code{contrasts_fix},
#'     \code{contrasts_sd}}{Factor levels and contrasts used to build the design
#'     matrices. Used to build the design matrices for new data.}
#'   \item{\code{call}}{The matched call.}
#' }
#'
#' @details
#' With \code{constraint = 1} (or \code{2}) the upper triangle of the loadings
#' \eqn{\Lambda} (or of the sample-level errors \eqn{U}) is fixed at 0 so that
#' the model is identifiable. A positive diagonal is currently not imposed, so the sign of
#' a latent axis is arbitrary.
#'
#' @seealso \code{\link{summary.rrGDM}}, \code{\link{anova.rrGDM}},
#'   \code{\link{wald_test}}, \code{\link{get_LV}}, \code{\link{predict.rrGDM}},
#'   \code{\link{biplot.rrGDM}}, \code{\link{stress}}, \code{\link{cor_r2}},
#'   \code{\link{include_se}}, \code{\link{sparse_diss}}
#'
#' @examples
#'   # Toy community with 25 sites, 25 species, 2 environmental gradients
#'   set.seed(2)
#'   n <- 25 # sites
#'   p <- 25 # species
#'
#'   env <- data.frame(temp = rnorm(n), moist = rnorm(n))
#'   eta <- as.matrix(env) %*% matrix(rnorm(2 * p), 2, p)
#'   Y <- matrix(rpois(n * p, exp(1 + 0.8 * eta)), n, p)
#'
#'   fit <- rrGDM(Y, env, k = 2,
#'                formula = ~ temp + moist,
#'                family = "binomial")
#'   summary(fit)
#'
#'   # Model the sample-level dispersion too
#'   fit2 <- rrGDM(Y, env, k = 2,
#'                 formula = ~ temp + moist,
#'                 formula_sd = ~ temp,
#'                 family = "binomial")
#'   summary(fit2)
#'
#' @export
rrGDM <- function(Y,
                  X = NULL,
                  k = 2,
                  binary = FALSE,
                  formula = ~.,
                  formula_sd = ~ 1,
                  dist_method = 'bray',
                  scale_diss = NULL,
                  family = 'gaussian',
                  link = NULL,
                  se = TRUE,
                  sample_error = TRUE,
                  constraint = 1,
                  trans = function(x) x,
                  sparse = NULL,
                  seed = NULL,
                  init_jitter = 0.1,
                  n_starts = 1,
                  previous_best = TRUE,
                  init_method = 'dbRDA',
                  start_fit = NULL,
                  trace = FALSE,
                  control = list()) {

  if (anyNA(Y)) cli::cli_abort("{.arg Y} contains missing values.")
  family <- match.arg(family, c('gaussian', 'binomial', 'beta', 'beta.binomial'))
  if (!is.null(link)) link <- match.arg(link, c('identity', 'negexp', 'neggaus'))
  if (!is.null(init_method)) init_method <- match.arg(init_method, c('dbRDA', 'PCOA'))
  if (!is.numeric(k) || length(k) != 1 || k < 1 || k != round(k)) cli::cli_abort("{.arg k} must be a single positive integer.")
  if (!is.numeric(n_starts) || length(n_starts) != 1 || n_starts < 1 || n_starts != round(n_starts)) cli::cli_abort("{.arg n_starts} must be a single positive integer.")
  if (!is.null(seed)) set.seed(seed)
  if (!(constraint %in% c(0,1,2))) cli::cli_abort("Identifiable constraint should be set to either 0, 1, or 2")
  if (!is.null(scale_diss) && (family %in% c('binomial','beta.binomial'))) {
    cli::cli_warn("Ignoring scale_diss for binomial and beta-binomial families ")
    scale_diss <- NULL
  }
  if (!is.null(scale_diss) &&
      !(is.numeric(scale_diss) && length(scale_diss) == 2 && !anyNA(scale_diss) && scale_diss[1] != scale_diss[2])) {
    cli::cli_abort("{.arg scale_diss} must be {.code NULL} or two different numbers, e.g. {.code c(0.01, 0.99)}.")
  }
  if (!is.logical(previous_best) || length(previous_best) != 1 || is.na(previous_best)) {
    cli::cli_abort("{.arg previous_best} must be either TRUE or FALSE.")
  }

  # ---------- DATA PREP ----------
  cli::cli_process_start('Preparing data')

  # Map Family and Link to index
  fam_idx <- switch(family,
                    "gaussian" = 0,
                    "binomial" = 1,
                    "beta" = 2,
                    "beta.binomial" = 3,
                    cli::cli_abort("Family not supported"))

  if (is.null(link)) {
    link <- switch(family,
                   "gaussian" = 'identity',
                   "binomial" = 'negexp',
                   "beta" = 'negexp',
                   "beta.binomial" = 'negexp')
  }

  link_idx <- switch(link,
                     "identity" = 0,
                     "negexp" = 1,
                     "neggaus" = 2,
                     cli::cli_abort("Link not supported"))

  if (link == 'identity' && family != 'gaussian') cli::cli_warn("{.code link = 'identity'} does not keep fitted dissimilarities within (0, 1) for the {.val {family}} family; consider {.code 'negexp'} or {.code 'neggaus'}.")

  # Prep. diss data
  Y_df <- make_Y_diss(Y,
                      trans = trans,
                      binary = binary,
                      method = dist_method,
                      num_den = (family %in% c('binomial', 'beta.binomial')))
  if (anyNA(Y_df)) cli::cli_abort("Dissimilarities computed with {.val {dist_method}} contain missing values.")

  obs_full <- NULL
  if (!is.null(sparse)) {
    obs_full <- scale_dist(Y_df[, 3], scale_diss)
    if (family %in% c('binomial', 'beta.binomial')) obs_full <- obs_full / as.numeric(Y_df[, 4])
    Y_df <- subset_pairs(Y_df, sparse, n = nrow(Y), id = rownames(Y), dist_method = dist_method)
  }


  D_mat <- Y_df[, c(1,2)]  # sample indices of each pair
  Y_val <- Y_df[,3]  # Dissimilarity value (or n unshared species if num_den type of dissimilarity)
  Y_val <- scale_dist(Y_val, scale_diss)  # re-scale between set interval, if defined by user (bypassed if NULL)
  Y_den <- if (family %in% c('binomial', 'beta.binomial')) as.numeric(Y_df[, 4]) else rep(1, length(Y_val)) # Denominator (used in binomial and beta-binomial)

  if (family == 'beta' && any(Y_val <= 0 | Y_val >= 1)) {
    cli::cli_abort(c(
      "x" = "The beta family needs dissimilarities strictly between 0 and 1, but some are exactly 0 or 1.",
      "i" = "Use {.arg scale_diss}, e.g. {.code scale_diss = c(0.01, 0.99)}, or the {.val beta.binomial} family."))
  }

  data_design <- prepare_data(formula = formula,
                              formula_sd = formula_sd,
                              data = X,
                              n = nrow(Y))

  # ---------- IDENTIFIABILITY CHECKS ----------
  n_env <- ncol(data_design$X_design)
  n_samples <- nrow(Y)

  if (n_env == 0 && constraint == 1) {
    constraint <- 2
    cli::cli_warn("No predictors in {.arg formula}: identifiability constraint set to 2.")
  }

  if (constraint == 0) cli::cli_warn("No identifiability constraint. Consider using constraint == 1, or constraint == 2")
  if (constraint == 1 && k > n_env) cli::cli_warn("k is larger than the number of predictors. This may cause problems. Consider 'constraint == 2'")
  if ((constraint == 2) && !sample_error) cli::cli_abort('sample_error needs to be set to TRUE if constraint = 2')
  if ((constraint == 2) && k > n_samples) cli::cli_abort("{.arg k} ({k}) cannot be larger than the number of samples ({n_samples}).")

  # ---------- TMB DATA ----------
  tmb_data <- list(
    Y = as.vector(Y_val),
    Y_den = as.vector(Y_den),
    X = as.matrix(data_design$X_design),
    W = as.matrix(data_design$W_design),
    D = as.matrix(D_mat - 1), # 0-indexing for C++
    Z = data_design$Z_design,
    has_sample_re = as.integer(sample_error),
    has_random = as.integer(data_design$has_random),
    family = fam_idx,
    link = link_idx,
    constraint = constraint,
    map_re = as.integer(data_design$map_re - 1),
    n_lv = k)
  # ---------- BASE TMB PARAMS ----------
  parameters <- list(
    lambda = matrix(stats::rnorm(n_env * k, 0, 1e-5), n_env, k),
    u_sample = matrix(stats::rnorm(n_samples*k,0,  1e-5), n_samples, k),
    log_scale = switch(family,
                       "beta" = 2,
                       "beta.binomial" = log(stats::median(Y_den)),
                       -1),
    beta = rep(0, ncol(data_design$W_design)),
    u = matrix(0, nrow = ncol(data_design$Z_design), ncol = k),
    log_sigma_re = rep(0, data_design$n_re_vars))

  if (!sample_error) parameters$u_sample <- matrix(0, nrow = 1, ncol = 1)

  cli::cli_process_done()

  # ---------- GOOD INITS -------------
  if (!is.null(start_fit)) {

    if (!inherits(start_fit, 'rrGDM')) cli::cli_abort("{.arg start_fit} must be a fitted {.cls rrGDM} model.")
    if (start_fit$k != k || start_fit$family != family || start_fit$link != link || start_fit$constraint != constraint) {
      cli::cli_abort("{.arg start_fit} must have the same {.arg k}, {.arg family}, {.arg link} and {.arg constraint} as the new model.")
    }
    if (nrow(start_fit$Y) != n_samples ||
        (!is.null(rownames(Y)) && !identical(rownames(start_fit$Y), rownames(Y)))) {
      cli::cli_abort("{.arg start_fit} must be fitted to the same samples (same order) as {.arg Y}.")
    }

    # All parameters of the previous fit, including random-effect modes
    start_par <- start_fit$obj$env$parList(start_fit$opt$par, start_fit$best_par)

    for (nm in names(parameters)) {
      if (length(start_par[[nm]]) != length(parameters[[nm]])) {
        cli::cli_abort("{.arg start_fit} has a different number of {.field {nm}} (predictors, {.arg sample_error} or random effects differ).")
      }
      parameters[[nm]][] <- as.vector(start_par[[nm]])  # keep the dimensions of the new object
    }

  } else if (!is.null(init_method) && init_method %in% c('PCOA', 'dbRDA')) {

    cli::cli_process_start(paste0('Finding initial parameters (', init_method, ')'))

    use_X <- (init_method == 'dbRDA') && (n_env > 0)

    loadings <- initial_loadings(
      trans(Y),
      X_design  = if (use_X) data_design$X_design else NULL,
      dist_method = dist_method,
      k  = k,
      constraint  = if (use_X) constraint else 2,
      QR  = use_X || constraint == 2,
      link = link,
      binary = binary)

    if (sample_error) {
      parameters$u_sample <- parameters$u_sample + as.matrix(loadings$u_sample)
      parameters$beta[1]  <- log(stats::sd(loadings$u_sample))
    }

    if (!is.null(loadings$lambda_init)) {
      parameters$lambda <- parameters$lambda + as.matrix(loadings$lambda_init)
    }

    if (family %in% c('beta', 'beta.binomial')) {
      L <- matrix(0, n_samples, k)
      if (sample_error) L <- L + parameters$u_sample
      if (n_env > 0)    L <- L + as.matrix(data_design$X_design) %*% parameters$lambda
      pr <- as.matrix(D_mat)
      d  <- sqrt(rowSums((L[pr[, 1], , drop = FALSE] - L[pr[, 2], , drop = FALSE])^2))
      mu <- pmin(pmax(make_inv_link(link)(d), 1e-3), 1 - 1e-3)

      if (family == 'beta.binomial') {
        ok  <- Y_den > 1
        z   <- (Y_val[ok] - Y_den[ok] * mu[ok])^2 / (Y_den[ok] * mu[ok] * (1 - mu[ok]))
        rho <- sum(z - 1) / sum(Y_den[ok] - 1)
      } else {
        z   <- (Y_val - mu)^2 / (mu * (1 - mu))
        rho <- mean(z)
      }
      if (is.finite(rho)) { # ignore if there are non-finite values
        rho <- min(max(rho, 1e-4), 0.99)
        parameters$log_scale <- log(1 / rho - 1)
      } else {
        cli::cli_warn("Moment matching for the dispersion start value gave a non-finite result; using the default start value instead.")
      }
    }

    cli::cli_process_done()
  }

  # ---------- STRUCTURAL ZEROS AND MAP ----------
  mapping <- list()

  if (family == 'binomial') { # no scale param
    mapping$log_scale <- factor(rep(NA, length(parameters$log_scale)))
  }

  if (constraint == 1) {
    # upper triangle of Lambda starts at exactly 0 and ignored in optimisation.
    parameters$lambda[upper.tri(parameters$lambda)] <- 0
    map_lambda <- matrix(seq_len(n_env * k), n_env, k)
    map_lambda[upper.tri(map_lambda)] <- NA
    mapping$lambda <- factor(as.vector(map_lambda))
  }

  if (!sample_error) { # u_sample is fixed at 0 (set above)
    mapping$u_sample <- factor(rep(NA, length(parameters$u_sample)))

  } else if (constraint == 2) {
    # Corner constraint on U: upper triangle fixed at 0.
    parameters$u_sample[upper.tri(parameters$u_sample)] <- 0
    map_u <- matrix(seq_len(n_samples * k), n_samples, k)
    map_u[upper.tri(map_u)] <- NA
    mapping$u_sample <- factor(as.vector(map_u))
  }

  if (data_design$has_random == 0) {   # ignore other random effects if not provided
    mapping$u <- factor(rep(NA, length(parameters$u)))
    mapping$log_sigma_re <- factor(rep(NA, length(parameters$log_sigma_re)))
  }

  # ------- FIT MODEL ---------
  obj <- TMB::MakeADFun(
    data = tmb_data,
    parameters = parameters,
    DLL = 'rrGDM',
    random = c('u', 'u_sample'),
    random.start = expression(last.par.best[random]),
    silent = !trace,
    map = mapping)

  # ------------- OPTIMISE --------------
  control_def <- list(iter.max = 1000)
  control_def[names(control)] <- control

  cli::cli_process_start('optimisation')

  best_obj <- Inf
  best_fit <- NULL
  ok_obj <- Inf
  ok_fit <- NULL
  ok_re <- NULL
  ok_H <- NULL
  last_error <- "no error message (the objective was not finite)"
  re_idx <- obj$env$random

  # Initial values
  init_fixed <- obj$par
  init_re <- if (length(re_idx) > 0) obj$env$last.par.best[re_idx] else numeric(0)
  base_fixed <- init_fixed
  base_re <- init_re

  for (i in seq_len(n_starts)) {
    # apply jitter: best fit so far, or the initial values
    centre_fixed <- if (previous_best) base_fixed else init_fixed
    centre_re <- if (previous_best) base_re else init_re

    if (length(re_idx) > 0) {
      re_start <- if (i == 1) centre_re else centre_re + stats::rnorm(length(re_idx), 0, init_jitter)
      obj$env$last.par.best[re_idx] <- re_start
      obj$env$last.par[re_idx] <- re_start
    }

    obj$env$value.best <- Inf
    start <- if (i == 1) centre_fixed else centre_fixed + stats::rnorm(length(centre_fixed), 0, init_jitter)
    fit <- tryCatch(
      stats::nlminb(start = start,
                    objective = obj$fn,
                    gradient = obj$gr,
                    control = control_def),
      error = function(e) {
        last_error <<- conditionMessage(e)
        list(objective = Inf, par = start)
      })

    is_best <- is.finite(fit$objective) && fit$objective < best_obj
    if (is_best) {
      best_obj <- fit$objective
      best_fit <- fit
      base_fixed <- fit$par
      if (length(re_idx) > 0) base_re <- obj$env$last.par.best[re_idx]
    }

    if (n_starts > 1 && is.finite(fit$objective) && fit$convergence == 0 && fit$objective < ok_obj) {
      obj$env$value.best <- Inf
      obj$fn(fit$par)  #inner modes at this optimum
      H <- tryCatch(stats::optimHess(fit$par, obj$fn, obj$gr), error = function(e) NULL)
      if (!is.null(H) && all(is.finite(H)) && hess_is_pd(H)) {
        ok_obj <- fit$objective
        ok_fit <- fit
        ok_H <- H
        ok_re <- if (length(re_idx) > 0) obj$env$last.par.best[re_idx] else numeric(0)
      }
    }

    if (trace) cli::cli_alert_info("Start {i}/{n_starts}: objective = {signif(fit$objective, 8)}{if (is_best && i > 1) ' (new best)' else ''}")
  }

  cli::cli_process_done()

  if (is.null(best_fit)) {
    cli::cli_abort(c(
      "x" = "Model optimisation failed for all {.arg n_starts} starting values.",
      "i" = "Last error: {last_error}",
      "i" = "Try increasing {.arg n_starts}, changing {.arg init_method}, or simplifying the model."
    ))
  }

  if (n_starts > 1) {
    if (!is.null(ok_fit)) {
      if (ok_obj > best_obj + 1e-2) {
        cli::cli_alert_info("Reporting the best converged fit (objective {round(ok_obj, 2)}); a lower objective ({round(best_obj, 2)}) came from a fit that did not converge to a positive-definite optimum.")
      }

      best_fit <- ok_fit
      best_obj <- ok_obj
      base_re <- ok_re

    } else {
      cli::cli_warn("No start converged to a positive-definite optimum: estimates and standard errors are not reliable.")
    }
  } else if (best_fit$convergence != 0) {
    cli::cli_warn("The optimiser did not converge (code {best_fit$convergence}: {best_fit$message}).")
  }

  # Update objective environment with fixed and random for best fit
  if (length(re_idx) > 0) {
    obj$env$last.par.best[re_idx] <- base_re
    obj$env$last.par[re_idx]  <- base_re
  }

  obj$env$value.best <- Inf
  obj$fn(best_fit$par) # evaluate once at the optimum so that obj$env holds the matching random-effect modes
  best_par <- obj$env$last.par # snapshot before sdrep

  sd_rep <- NULL

  if (isTRUE(se)) {
    sd_rep <- tryCatch(
      TMB::sdreport(obj,
                    getReportCovariance = TRUE,
                    hessian.fixed = if (n_starts > 1) ok_H),
      error = function(e) {
        cli::cli_warn("Standard errors could not be computed (the Hessian at the optimum might not be positive definite). Original error: {conditionMessage(e)}")
        NULL
      })
  }

  res <- list(
    Y = Y,
    X = X,
    Y_diss = Y_val,
    Y_den = Y_den,
    pairs = as.matrix(D_mat),
    sparse = sparse,
    obs_full = obs_full,
    X_design = data_design$X_design,
    W_design = data_design$W_design,
    Z_design = data_design$Z_design,
    map_re = data_design$map_re,
    re_vars = data_design$re_vars,
    opt = best_fit,
    obj = obj,
    sd_rep = sd_rep,
    pd_hess = compute_pd_hess(sd_rep),
    best_par = best_par,
    k = k,
    family = family,
    link = link,
    inv_link = make_inv_link(link),
    formula = formula,
    formula_fix = data_design$formula_fix,
    formula_sd = formula_sd,
    terms_fix = data_design$terms_fix,
    terms_sd  = data_design$terms_sd,
    xlevels_fix = data_design$xlevels_fix,
    xlevels_sd = data_design$xlevels_sd,
    contrasts_fix = data_design$contrasts_fix,
    contrasts_sd = data_design$contrasts_sd,
    call = match.call(),
    constraint = constraint
  )

  class(res) <- "rrGDM"
  return(res)
}


#' Calculate standard errors
#'
#' @param object A fitted \code{rrGDM} object.
#' @param ... Further arguments passed to \code{\link[TMB]{sdreport}}.
#'
#' @return A model object, including the sd report
#'
#' @examples
#' set.seed(2)
#' # simulate data
#' n <- 25 # sites
#' p <- 25 # species
#' env <- data.frame(temp = rnorm(n), moist = rnorm(n))
#' eta <- as.matrix(env) %*% matrix(rnorm(2 * p), 2, p)
#' Y <- matrix(rpois(n * p, exp(1 + 0.8 * eta)), n, p,
#'             dimnames = list(paste0("site", 1:n), paste0("sp", 1:p)))
#'
#' # fit model
#' fit <- rrGDM(Y, env, k = 2,
#'              formula = ~ temp + moist,
#'              family = "binomial")
#'
#' fit0 <- rrGDM(Y, env, k = 2, formula = ~ temp + moist, family = "binomial", se = FALSE)
#' fit0 <- include_se(fit0)
#' wald_test(fit0)
#' @export
#'
include_se <- function(object, ...) {
  if (!inherits(object, "rrGDM")) cli::cli_abort("{.arg object} must be a fitted {.cls rrGDM} model.")
  object$obj$env$last.par.best <- object$best_par
  sd_rep <- TMB::sdreport(object$obj, par.fixed = object$opt$par, ...)
  object$sd_rep <- sd_rep
  object$pd_hess <- compute_pd_hess(sd_rep)
  object
}

#' Print method for rrGDM models
#'
#' @description
#' Prints a compact overview of a fitted \code{rrGDM} object.
#'
#' @param x A fitted model of class \code{rrGDM}.
#' @param digits Number of significant digits used to print loadings. Default is 3.
#' @param n_u Number of rows to print for the sample-level error matrix
#'   \code{u} and, for each grouping variable, for the random-effect matrix
#'   (these have one row per sample or per group, so they are truncated).
#'   Default is 5.
#' @param ... Currently unused.
#'
#' @details
#' Convergence information includes the convergence code returned by
#' \code{\link[stats]{nlminb}} (0 indicates successful convergence), whether
#' the Hessian at the optimum is positive definite (only available if standard
#' errors were requested) and the largest absolute gradient of the fixed
#' parameters. Loadings are printed as estimated (i.e., unrotated). Use
#' \code{\link{get_LV}} to extract the full matrices.
#'
#' @method print rrGDM
#' @export
print.rrGDM <- function(x, digits = 3, n_u = 5, ...) {

  if (!inherits(x, "rrGDM")) {
    cli::cli_abort("{.arg x} must be a fitted {.cls rrGDM} model.")}
  if (!is.numeric(n_u) || length(n_u) != 1 || is.na(n_u) || n_u < 0) {
    cli::cli_abort("{.arg n_u} must be a single non-negative number.")}

  lv <- get_LV(x, rotate = FALSE)
  has_sample_re <- !is.null(lv$u_sample)
  has_random <- length(x$re_vars) > 0
  n_env <- ncol(x$X_design)
  if (is.null(n_env)) n_env <- 0L

  show <- function(label, value) cat(sprintf("  %-24s %s\n", label, value))

  # Print the first n_u rows
  print_head <- function(mat, what) {
    n_show <- min(nrow(mat), n_u)
    n_more <- nrow(mat) - n_show
    if (n_show > 0) print(mat[seq_len(n_show), , drop = FALSE], digits = digits)
    if (n_more > 0) {
      cat(sprintf("... %d more %s%s; use get_LV() to extract the full matrix.\n",
                  n_more, what, if (n_more == 1) "" else "s"))
    }
  }

  cat("Reduced-rank Generalised Dissimilarity Model (rrGDM)\n\n")
  cat("Call:\n", paste(deparse(x$call), collapse = "\n"), "\n", sep = "")

  # General model info
  cli::cli_h2("Model")

  mean_lab <- if (n_env > 0) format_formula(stats::formula(x$terms_fix)) else "none (no predictors)"
  re_lab <- if (has_random) paste0("(1 | ", x$re_vars, ")", collapse = " + ") else "none"
  sd_lab <- if (has_sample_re) format_formula(stats::formula(x$terms_sd)) else "none (no sample-level errors)"
  con_lab <- switch(as.character(x$constraint),
                    "0" = "none",
                    "1" = "lambda (lower triangular)",
                    "2" = "u (lower triangular)",
                    as.character(x$constraint))

  show("Method:", "maximum likelihood (Laplace approximation, TMB)")
  show("Family:", x$family)
  show("Link:", x$link)
  show("Latent dimensions:", x$k)
  n_all <- choose(nrow(x$X_design), 2)
  show("Samples (pairs):", if (is.null(x$sparse)) {
    sprintf("%d (%d)", nrow(x$X_design), length(x$Y_diss))
  } else {
    sprintf("%d (%d; sparse, %.1f%% of all pairs)", nrow(x$X_design), length(x$Y_diss), 100 * length(x$Y_diss) / n_all)
  })
  show("Fixed effects:", mean_lab)
  show("Random effects:", re_lab)
  show("Dispersion model:", sd_lab)
  show("Identifiability:", con_lab)

  # Convergence info
  cli::cli_h2("Convergence")

  pd_hess <- if (is.null(x$sd_rep)) {
    "not computed"
  } else if (isTRUE(hess_ok(x))) {
    "yes"
  } else {
    "NO"
  }

  max_grad <- x$sd_rep$gradient.fixed
  grad_lab <- if (is.null(max_grad)) "not computed" else formatC(max(abs(max_grad)), digits = 2, format = "e")

  show("nlminb code:", sprintf("%s (%s)", x$opt$convergence, x$opt$message))
  show("Hessian pos. definite:", pd_hess)
  show("Max |gradient|:", grad_lab)

  # Loadings
  cli::cli_h2("Environmental loadings (lambda)")
  if (nrow(lv$lambda) > 0) {
    print(lv$lambda, digits = digits)
  } else {
    cli::cli_alert_info("No predictors in the mean model.")
  }

  # Sample error
  cli::cli_h2("Sample-level errors (u)")
  if (!has_sample_re) {
    cli::cli_alert_info("Sample-level errors not included in the model.")
  } else {
    print_head(lv$u_sample, "sample")
  }

  # Random effects
  if (has_random) {
    cli::cli_h2("Random effects (u)")
    if (is.null(lv$u_re)) {
      cli::cli_alert_info("Random effects not found in the model report.")
    } else {
      u_re <- lv$u_re
      # One block per grouping variable (rows of u_re follow x$map_re)
      grp <- if (length(x$map_re) == nrow(u_re)) x$map_re else rep(1L, nrow(u_re))

      for (v in seq_along(x$re_vars)) {
        u_v <- u_re[grp == v, , drop = FALSE]
        if (nrow(u_v) == 0) next
        cat(sprintf("(1 | %s): %d group%s\n", x$re_vars[v], nrow(u_v), if (nrow(u_v) == 1) "" else "s"))
        print_head(u_v, "group")
      }
    }
  }

  cat("\nUse summary() for Wald tests, dispersion effects and fit statistics.\n")
  invisible(x)
}
