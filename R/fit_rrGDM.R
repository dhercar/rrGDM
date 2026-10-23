#'rrGDM' function
#' @description Fit a reduced-rank Generalised Dissimilarity Model using TMB
#' @param Y A matrix or data frame containing the response variable (usually species abundance matrix).
#' @param X A data frame containing the predictor variables.
#' @param k An integer specifying the number of latent variables to be used in the model. Default is 2.
#' @param formula Formula specifying the model to be fitted.
#' @param formula_sd Formula specifying the model for the dispersion parameter. Default is ~1 (intercept only).
#' @param constraint Integer. Specifies identifiability constraints. 1: Constraint applied to environmental loadings matrix (\eqn{\Lambda}); 2: Constraint applied to site random effects matrix (U). Default is 1.
#' @param dist_method A character string specifying the distance method to be used (e.g., "euclidean", "bray"). Default is "bray".
#' @param trans Function used to transform the response matrix before computing dissimilarities (e.g., log, sqrt). Default is no transformation.
#' @param sparse Optional object of class \code{sparse_diss} (see \code{\link{sparse_diss}}). If supplied, the model is fitted using only the pairs of samples listed in \code{sparse$edges}, instead of all pairs. It must be created from the same samples as \code{Y}. Default is NULL (all pairs).
#' @param binary Logical indicating whether the response variable should be treated as binary (e.g., presence-absence). Default is FALSE.
#' @param scale_diss A numeric vector of length two indicating a new scale for dissimilarities. Useful for distributions bounded between 0 and 1 (e.g., beta family)
#' @param family A character string specifying the family of the response variable (e.g., "beta", "binomial"). Default is "gaussian".
#' @param link Character string specifying the link function to be used (e.g., "negexp"). Default is 'identity' for Gaussian distribution, and 'negexp' for beta and binomial distributions.
#' @param se Logical indicating whether to calculate standard errors (TRUE) or not (FALSE). Default is TRUE.
#' @param sample_error Logical. If TRUE, includes site-level random effects to account for variation not explained by predictors.
#' @param seed An integer specifying the seed for random number generation. Default is NULL.
#' @param init_jitter Numeric. A small standard deviation used to jitter initial parameter values.
#' @param previous_best Logical. Controls where the jitter for multiple starts after the first is applied
#' If \code{TRUE} (default), each start jitters around the best solution found so far (perturb-and-reoptimise). If \code{FALSE}, every start jitters around the initial values (from \code{init_method}), so starts do not depend on earlier fits. Only relevant when \code{n_starts > 1}.
#' @param init_method A character string specifying a method used to obtain initial values for optimisation. Currently, only \code{'PCOA'} and \code{'dbRDA'} available. Ignored if \code{start_fit} is supplied.
#' @param start_fit Optional fitted \code{rrGDM} model used to initialise the optimisation (e.g., a fast sparse fit before refitting with more pairs). It must have the same samples, \code{k}, \code{family}, \code{link}, \code{constraint}, predictors and random-effect structure as the new model; the set of pairs (\code{sparse}) may differ. If supplied, \code{init_method} is ignored. Default is NULL.
#' @param control A list of control parameters passed to stats::nlminb
#' @return An object of class 'rrGDM' containing the fitted model object, the optimized parameters, and other relevant information.
#' @importFrom TMB MakeADFun sdreport
#' @importFrom stats nlminb rnorm median sd
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
  if (!is.null(seed)) set.seed(seed)
  if (!(constraint %in% c(0,1,2))) cli::cli_abort("Identifiable constraint should be set to either 0, 1, or 2")

  if (!is.null(scale_diss) && (family %in% c('binomial','beta.binomial'))) {
    cli::cli_warn("Ignoring scale_diss for binomial and beta-binomial families ")
    scale_diss <- NULL
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


  # Prep. diss data
  Y_df <- make_Y_diss(Y,
                      trans = trans,
                      binary = binary,
                      method = dist_method,
                      num_den = (family %in% c('binomial', 'beta.binomial')))
  if (anyNA(Y_df)) cli::cli_abort("Dissimilarities computed with {.val {dist_method}} contain missing values.")

  obs_full <- NULL
  if (!is.null(sparse)) {
    obs_full <- scale_dist(Y_df[, 3], sort(scale_diss))
    if (family %in% c('binomial', 'beta.binomial')) obs_full <- obs_full / as.numeric(Y_df[, 4])
    Y_df <- subset_pairs(Y_df, sparse, n = nrow(Y), id = rownames(Y), dist_method = dist_method)
  }


  D_mat <- Y_df[, c(1,2)]  # Design matrix (row col pairs)
  Y_val <- Y_df[,3]  # Dissimilarity value (or n unshared species if num_den type of dissimilarity)
  Y_val <- scale_dist(Y_val, sort(scale_diss))  # re-scale between set interval, if defined by user (bypassed if NULL)
  Y_den <- if (family %in% c('binomial', 'beta.binomial')) as.numeric(Y_df[, 4]) else rep(1, length(Y_val)) # Denominator (used in binomial and beta-binomial)

  # Prep predictors
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

  # Check constraints
  if (constraint == 0) cli::cli_warn("No identifiability constraint. Consider using constraint == 1, or constraint == 2")
  if (constraint == 1 && k > n_env) cli::cli_warn("k is larger than the number of predictors. This may cause problems. Consider 'constraint == 2'")
  if ((constraint == 2) && !sample_error) cli::cli_abort('sample_error needs to be set to TRUE if constraint = 2')
  if ((constraint == 2) && k > n_samples) cli::cli_abort(paste0("k (", k, ") cannot be larger than the number of samples (", n_samples, ")."))

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
    lambda = matrix(rnorm(n_env * k, 0, 1e-5), n_env, k),
    u_sample = matrix(rnorm(n_samples*k,0,  1e-5), n_samples, k),
    log_scale = switch(family,
                       "beta" = 2,
                       "beta.binomial" = log(median(Y_den)),
                       -1),
    beta = rep(0, ncol(data_design$W_design)),
    u = matrix(0, nrow = ncol(data_design$Z_design), ncol = k),
    log_sigma_re = rep(0, data_design$n_re_vars))

  if (sample_error == FALSE) parameters$u_sample <- matrix(0, nrow = 1, ncol = 1)

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
      parameters$beta[1]  <- log(sd(loadings$u_sample))
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
      mu <- pmin(pmax(switch(link, identity = d, negexp = 1 - exp(-d), neggaus = 1 - exp(-d^2)), 1e-3), 1 - 1e-3)

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
        cli::cli_text('Moment matching produced non-infinite values')
      }
    }

    cli::cli_process_done()
  }


  # ---------- 6. STRUCTURAL ZEROS AND MAP ----------
  mapping <- list()

  if (family == 'binomial') { # no scale param
    mapping$log_scale <- factor(rep(NA, length(parameters$log_scale)))
  }

  if (constraint == 1) {
    # upper triangle of Lambda starts at exactly 0 and ignored in optimisation.
    # IMPORTANT !!!: Positive diagonal is NOT hard-imposed at the moment, so sign flips are possible
    parameters$lambda[upper.tri(parameters$lambda)] <- 0
    map_lambda <- matrix(seq_len(n_env * k), n_env, k)
    map_lambda[upper.tri(map_lambda)] <- NA
    mapping$lambda <- factor(as.vector(map_lambda))
  }

  if (!sample_error) { # ignore in optimisation
    parameters$u_sample <- matrix(0, nrow = 1, ncol = 1)
    mapping$u_sample <- factor(rep(NA, length(parameters$u_sample)))

  } else if (constraint == 2) {
    # Corner constraint on U: upper triangle fixed at 0.
    # IMPORTANT!!!! : Positive diagonal is NOT hard-imposed at the moment
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
  control_def <- list(iter.max = 1000)   # nlminb higher iter.max
  control_def[names(control)] <- control

  # optimise
  cli::cli_process_start('optimisation')

  best_obj <- Inf
  best_fit <- NULL
  ok_obj <- Inf
  ok_fit <- NULL
  ok_re <- NULL
  ok_H <- NULL
  re_idx <- obj$env$random

  # Initial values
  init_fixed <- obj$par
  init_re <- if (length(re_idx) > 0) obj$env$last.par.best[re_idx] else numeric(0)
  base_fixed <- init_fixed
  base_re <- init_re

  # First fit with higher tolerance for faster convergence (not when warm-starting from a previous fit):
  control_first <- control_def
  if (n_starts > 1 && is.null(control$rel.tol) && is.null(start_fit)) control_first$rel.tol <- 1e-7

  for (i in seq_len(n_starts)) {
    # apply jitter: best fit so far, or the initial values
    centre_fixed <- if (previous_best) base_fixed else init_fixed
    centre_re <- if (previous_best) base_re else init_re

    if (length(re_idx) > 0) {
      re_start <- if (i == 1) centre_re else centre_re + rnorm(length(re_idx), 0, init_jitter)
      obj$env$last.par.best[re_idx] <- re_start
      obj$env$last.par[re_idx]      <- re_start
    }

    obj$env$value.best <- Inf
    start <- if (i == 1) centre_fixed else centre_fixed + rnorm(length(centre_fixed), 0, init_jitter)
    fit <- tryCatch(
      stats::nlminb(start = start,
                    objective = obj$fn,
                    gradient = obj$gr,
                    control = if (i == 1) control_first else control_def),
      error = function(e) list(objective = Inf, par = start))

    if (is.finite(fit$objective) && fit$objective < best_obj) {
      if (i > 1) cat('new best: ')
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

    cat('model ', i, '; ', fit$objective, '\n')
  }

  cli::cli_process_done()

  if (is.null(best_fit)) {
    cli::cli_abort(c(
      "x" = "Model optimisation failed for all {.arg n_starts} starting values.",
      "i" = "Try increasing {.arg n_starts}, changing {.arg init_method}, or simplifying the model."
    ))
  }

  if (n_starts > 1) {
    if (!is.null(ok_fit)) {
      if (ok_obj > best_obj + 1e-3) {
        cli::cli_alert_info("Reporting the best converged fit (objective {round(ok_obj, 2)}); a lower objective ({round(best_obj, 2)}) came from a fit that did not converge to a positive-definite optimum.")
      }

      best_fit <- ok_fit
      best_obj <- ok_obj
      base_re <- ok_re

    } else {
      cli::cli_warn("No start converged to a positive-definite optimum: estimates and standard errors are not reliable. Try more starts, a larger {.arg init_jitter}, or a different {.arg init_method}.")
    }
  }

  # Update objective environment with fixed and random for best fit
  if (length(re_idx) > 0) {
    obj$env$last.par.best[re_idx] <- base_re
    obj$env$last.par[re_idx]  <- base_re
  }

  obj$env$value.best <- Inf
  obj$fn(best_fit$par) # update objective (see * below)
  best_par <- obj$env$last.par # snapshot before sdrep

  # SD report
  sd_rep <- NULL

  if (se == TRUE) {
    sd_rep <- tryCatch(
      TMB::sdreport(obj,
                    getReportCovariance = TRUE,
                    hessian.fixed = if (n_starts > 1) ok_H),
      error = function(e) {
        cli::cli_warn('standard errors could not be computed. Hessian at optimum might not positive definite.')
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
#' @param obj A fitted rrgdm object
#' @param ... Parameters passed to `TMB::sdreport()`
#'
#' @return A model object, including the sd report
#' @export
#' @importFrom TMB sdreport
#'
include_se <- function(obj, ...) {
  sd_rep <- TMB::sdreport(obj$obj, getReportCovariance = TRUE, ...)
  obj$sd_rep <- sd_rep
  obj$pd_hess <- compute_pd_hess(sd_rep)
  return(obj)
}


#' Collapse a formula into a single string (internal)
#' @param f A formula.
#' @noRd
format_formula <- function(f) paste(trimws(deparse(f)), collapse = " ")

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
#' @importFrom stats formula
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
