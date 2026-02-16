#'rrGDM' function
#' @description Fit a reduced-rank Generalised Dissimilarity Model using TMB
#' @param Y A matrix or data frame containing the response variable (usually species abundance matrix).
#' @param X A data frame containing the predictor variables.
#' @param k An integer specifying the number of latent variables to be used in the model. Default is 2.
#' @param formula Formula specifying the model to be fitted.
#' @param formula_sd Formula specifying the model for the dispersion parameter. Default is ~1 (intercept only).
#' @param constraint Integer. Specifies identifiability constraints. 1: Constraint applied to environmental loadings matrix ($\lambda$); 2: Constraint applied to site random effects matrix (U). Default is 1.
#' @param dist_method A character string specifying the distance method to be used (e.g., "euclidean", "bray"). Default is "bray".
#' @param trans Function used to transform the response matrix before computing dissimilarities (e.g., log, sqrt). Default is no transformation.
#' @param binary Logical indicating whether the response variable should be treated as binary (e.g., presence-absence). Default is FALSE.
#' @param scale_diss A numeric vector of length two indicating a new scale for dissimilarities. Useful for distributions bounded between 0 and 1 (e.g., beta family)
#' @param family A character string specifying the family of the response variable (e.g., "beta", "binomial"). Default is "gaussian".
#' @param link Character string specifying the link function to be used (e.g., "negexp"). Default is 'identity' for Gaussian distribution, and 'negexp' for beta and binomial distributions.
#' @param se Logical indicating whether to calculate standard errors (TRUE) or not (FALSE). Default is TRUE.
#' @param sample_error Logical. If TRUE, includes site-level random effects to account for variation not explained by predictors.
#' @param seed An integer specifying the seed for random number generation. Default is NULL.
#' @param init_jitter Numeric. A small standard deviation used to jitter initial parameter values.
#' @param find_initials A character string specifying a method used to obtain initial values for optimisation. Currently, only \code{find_initials = 'PCOA'} available.
#' @param control A list of control parameters passed to stats::nlminb
#' @return An object of class 'rrGDM' containing the fitted model object, the optimized parameters, and other relevant information.
#' @importFrom TMB MakeADFun sdreport
#' @importFrom stats nlminb
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
                  seed = NULL,
                  init_jitter = 0.01,
                  find_initials = 'PCOA',
                  trace = FALSE,
                  control = list()) {

  if (!is.null(seed)) set.seed(seed)

  # ToDo: CHECKS
  # - Y has no NA
  # - X has no NA
  # - formula_sd has *at least* an intercept

  # Already taken care of by other functions:
  #  - binary is logical is checked by make_y_df()
  # - `trans` is null or a function is checked by make_y_df()
  # - `formula` and `formula_sd` are already checked within prep_data()
  # - dist method is already checked within make_Y_diss()

  # ---------- 1. DATA PREP ----------
  cli::cli_process_start('Preparing data')

  # Map Family and Link to index
  fam_idx <- switch(family,
                    "gaussian" = 0,
                    "binomial" = 1,
                    "beta" = 2,
                    cli::cli_abort("Family not supported"))

  if (is.null(link)) {
    link <- switch(family,
                       "gaussian" = 'identity',
                       "binomial" = 'negexp',
                       "beta" = 'negexp')
  }

  link_idx <- switch(link,
                    "identity" = 0,
                    "negexp" = 1,
                    "neggaus" = 2,
                     cli::cli_abort("Link not supported"))



  # Prep. diss data
  Y_df <- make_Y_diss(Y, trans = trans, binary = binary, num_den = (family == 'binomial'))
  D_mat <- Y_df[, c(1,2)]                                                              # Design (row col pairs)
  Y_val <- Y_df[,3]
  Y_val <- scale_dist(Y_val, sort(scale_diss))                                         # scale between set interval, if defined
  Y_den  <- if (family == 'binomial') as.numeric(Y_df[, 4]) else rep(1, length(Y_val)) # Denominator (used in binomial)

  # Prep predictors
  data_design <- prepare_data(formula = formula, formula_sd = formula_sd, data = X, n = nrow(Y))

  # ---------- 2. IDENTIFIABILITY CHECKS ----------

  n_env <- if (!is.null(data_design$X_design)) ncol(data_design$X_design) else 0
  n_samples <- nrow(Y)

  # Check constraints
  if ((constraint == 1) && !is.matrix(data_design$X_design)) cli::cli_abort('Predictors should be provided if constraint = 1')
  if (constraint == 1 && k > n_env) cli::cli_warn("k is larger than the number of predictors. This may cause problems. Consider 'constraint == 2'")
  if ((constraint == 2) && !sample_error) cli::cli_abort('sample_error needs to be set to TRUE if constraint = 2')
  if ((constraint == 2) && k > n_samples) cli::cli_abort(paste0("k (", k, ") cannot be larger than the number of samples (", n_samples, ")."))

  # ---------- 3. SET PARAM SIZE ----------
  # I define diagonal, upper and lower matrix independently to allow constraints

  # BETA: Loading matrix: n_env x k
  n_diag_beta  <- min(n_env, k)
  n_lower_beta <- sum(row(matrix(0, n_env, k)) > col(matrix(0, n_env, k)))
  n_upper_beta <- sum(row(matrix(0, n_env, k)) < col(matrix(0, n_env, k)))

  # ---------- 4. TMB DATA AND PARAM ----------
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
    n_lv = k
  )

  n_diag_lambda  <- min(n_env, k)
  n_lower_lambda <- sum(row(matrix(0, n_env, k)) > col(matrix(0, n_env, k)))
  n_upper_lambda <- sum(row(matrix(0, n_env, k)) < col(matrix(0, n_env, k)))

  parameters <- list(
    lambda_diag = rnorm(n_diag_lambda, 0, init_jitter),
    lambda_lower = rnorm(n_lower_lambda, 0, init_jitter),
    lambda_upper = rnorm(n_upper_lambda, 0, init_jitter),
    u_sample = matrix(rnorm(n_samples*k,0, init_jitter), n_samples, k),
    log_scale = ifelse(family == "beta", 2, -1), # phi starts at ~7.4 for Beta
    beta = rep(0, ncol(data_design$W_design)),
    u = matrix(0, nrow = ncol(data_design$Z_design), ncol = k),
    log_sigma_re = rep(0, data_design$n_re_vars)
  )

  if (sample_error == F) parameters$u_sample <- matrix(1, nrow = 1, ncol = 1)

  # ---------- 5. IGNORE PARAMETERS ----------

  mapping <- list()

  if (family == 'binomial') {
    # Ignore scale parameter if binomial
    mapping$log_scale <- factor(rep(NA, length(parameters$log_scale)))
  }

  if (constraint == 1) {
    # Upper triangle of lambda is zero and ignored in optimisation
    mapping$lambda_upper <- factor(rep(NA, length(parameters$lambda_upper)))
    parameters$lambda_upper <- 0
  }

  if (!sample_error) {
    mapping$u_sample <- factor(rep(NA, length(parameters$u_sample)))
  } else if (constraint == 2) {
    # Lower-triangular constraint (set to 0)
    parameters$u_sample[upper.tri(parameters$u_sample)] <- 0
    # Set to NA in mapping
    map_u <- matrix(seq_len(n_samples * k), n_samples, k)
    map_u[upper.tri(map_u)] <- NA
    mapping$u_sample <- factor(as.vector(map_u))
  }

  #  other random effects
  if (data_design$has_random == 0) {
    mapping$u <- factor(rep(NA, length(parameters$u)))
    mapping$log_sigma_re <- factor(rep(NA, length(parameters$log_sigma_re)))
  }

  cli::cli_process_done()


  if ( find_initials == 'PCOA' ) {
    cli::cli_process_start('Finding initial parameters')
    loadings <- rrGDM:::initial_loadings(Y,
                                         dist_method = dist_method,
                                         k = k,
                                         QR = TRUE,
                                         link = link,
                                         binary = binary)

    parameters$u_sample <- as.matrix(loadings$u_sample)
    parameters$beta[1] <- log(sd(loadings$u_sample))
    cli::cli_progress_done()
  }

  # --- 6. FIT MODEL ---

  obj <- TMB::MakeADFun(
    data = tmb_data,
    parameters = parameters,
    DLL = 'rrGDM',
    random = c('u', 'u_sample'),
    silent = !trace,
    map = mapping,
  )

  # nlminb control default
  control_def <- list(rel.tol = 1e-10, eval.max = 1000, iter.max = 1000, trace = 0)
  control_def[names(control)] <- control # change user preferences


  # optimise
  cli::cli_process_start('optimisation')
  opt <- stats::nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr,
    control = control_def
  )
  cat(opt$message)
  cli::cli_process_done()

  # SD report
  sd_rep <- NULL
  if (se == TRUE) sd_rep <- TMB::sdreport(obj)

  res <- list(
    Y = Y,
    X = X,
    Y_diss = Y_val,
    Y_den = Y_den,
    X_design = data_design$X_design,
    W_design = data_design$W_design,
    Z_design = data_design$Z_design,
    map_re = data_design$map_re,
    opt = opt,
    obj = obj,
    sd_rep = sd_rep,
    k = k,
    family = family,
    link = link,
    formula = formula,
    formula_fix = data_design$formula_fix,
    formula_sd = formula_sd,
    call = match.call())

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
  sd_rep <- sdreport(obj$obj, ...)
  obj$sd_rep <- sd_rep
  return(obj)
}




