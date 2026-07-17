#' Prepare design matrices and data (Internal)
#' @noRd
prepare_data <- function(formula = NULL,
                         formula_sd = NULL,
                         data = NULL,
                         n = NULL) {

  # ---- Input data validation start ----
  if (is.null(n)) {
    cli::cli_abort(c("{.var n} must be supplied",
                     "x" = "The number of observations is required for design matrix dimensions."))
  }

  if (!is.null(formula) && !inherits(formula, "formula")) {
    cli::cli_abort("{.var formula} must be a {.cls formula}, not a {.cls {class(formula)}}.")
  }

  if (!is.null(formula_sd) && !inherits(formula_sd, "formula")) {
    cli::cli_abort("{.var formula_sd} must be a {.cls formula}, not a {.cls {class(formula_sd)}}.")
  }

  if (!is.null(formula_sd)) {
    if (attr(terms(formula_sd), "intercept") == 0) {
      cli::cli_abort(c(
        "{.var formula_sd} must include an intercept.",
        "i" = "Dispersion models require a baseline (intercept)",
        "!" = "Remove {.code 0} or {.code -1} from your {.var formula_sd}."
      ))
    }
  }

  if (!is.null(data) && !(is.matrix(data) || is.data.frame(data))) {
    cli::cli_abort("{.var data} must be a {.cls matrix} or {.cls data.frame}.")
  }
  # ---- Input data validation ends ----

  if (is.null(formula_sd) || is.null(data)) {
    formula_sd  <- ~1
  }

  # Initialise terms to NULL
  fixed_terms   <- NULL
  sd_terms      <- NULL
  xlevels_fix   <- NULL
  xlevels_sd    <- NULL
  contrasts_fix <- NULL
  contrasts_sd  <- NULL

  if (is.null(formula) || is.null(data)) {
    X_design      <- integer(0)
    fixed_formula <- ~0

  } else {
    # --- FIXED EFFECTS ---
    fixed_formula <- lme4::nobars(formula)

    # Standard way to capture data-dependent transformations:
    # model.frame -> terms(mf) -> model.matrix
    mf_fix      <- model.frame(fixed_formula, data = data)
    fixed_terms <- attr(mf_fix, "terms")
    X_tmp       <- model.matrix(fixed_terms, data = mf_fix)

    xlevels_fix <- if (is.data.frame(data)) {
      lapply(data[sapply(data, is.factor)], levels)
    } else {
      NULL
    }

    if ("(Intercept)" %in% colnames(X_tmp)) {
      intercept_idx <- which(colnames(X_tmp) == "(Intercept)")
      X_design <- X_tmp[, -intercept_idx, drop = FALSE]
      attr(X_design, "assign") <- attr(X_tmp, "assign")[-intercept_idx]
      attr(X_design, "contrasts") <- attr(X_tmp, "contrasts")
    } else {
      X_design <- X_tmp
    }

    contrasts_fix <- attr(X_tmp, "contrasts")

    # --- DISPERSION MODEL ---
    if (is.null(data)) {
      W_design <- integer(0)
    } else if (!is.null(lme4::findbars(formula_sd))) {
      cli::cli_abort(c("Random effects not allowed in dispersion model",
                       "x" = "Check {.var formula_sd}: random effects are only supported in the main formula."))
    }

    sd_formula <- update(formula_sd, ~ 1 + .)
    mf_sd      <- model.frame(sd_formula, data = data)
    sd_terms   <- attr(mf_sd, "terms")
    W_design   <- model.matrix(sd_terms, data = mf_sd)

    xlevels_sd <- if (is.data.frame(data)) {
      lapply(data[sapply(data, is.factor)], levels)
    } else {
      NULL
    }

    contrasts_sd <- attr(W_design, "contrasts")
  }

  # Random parts
  random_parts <- if (!is.null(formula)) lme4::findbars(formula) else NULL

  # Return with no re
  if ((length(random_parts) == 0) || is.null(data)) {
    Z_design <- Matrix::sparseMatrix(
      i = integer(0), j = integer(0), x = numeric(0), dims = c(n, 1)
    )

    return(list(
      X_design      = X_design,
      W_design      = W_design,
      Z_design      = Z_design,
      map_re        = integer(0),
      n_re          = 0,
      n_re_vars     = 0,
      has_random    = FALSE,
      formula_fix   = fixed_formula,
      terms_fix     = fixed_terms,
      terms_sd      = sd_terms,
      xlevels_fix   = xlevels_fix,
      xlevels_sd    = xlevels_sd,
      contrasts_fix = contrasts_fix,
      contrasts_sd  = contrasts_sd
    ))
  }

  # Random effects processing (remains same as your original)
  Z_list <- list()
  re_vars <- character()
  n_groups_vars <- integer()

  for (re in random_parts) {
    if (re[[2]] == 1) {
      var <- deparse(re[[3]])
      if (!(var %in% colnames(data))) {
        cli::cli_abort("Variable {.val {var}} not found in the provided data.")
      }

      factor_vec <- factor(data[[var]])
      re_levels <- levels(factor_vec)
      n_groups <- length(re_levels)
      Z_re <- Matrix::sparse.model.matrix(~0 + factor_vec)
      colnames(Z_re) <- paste0(var, "_", re_levels)

      Z_list[[var]] <- Z_re
      re_vars <- c(re_vars, var)
      n_groups_vars <- c(n_groups_vars, n_groups)
    } else {
      cli::cli_abort(c("Random slopes are not supported",
                       "i" = "Only random intercepts in the form {.code (1 | group)} are allowed."))
    }
  }

  Z_design <- do.call(cbind, Z_list)
  map_re <- rep(seq_along(re_vars), n_groups_vars)

  return(list(
    X_design      = X_design,
    W_design      = W_design,
    Z_design      = Z_design,
    map_re        = map_re,
    n_re          = length(map_re),
    n_re_vars     = length(re_vars),
    has_random    = TRUE,
    formula_fix   = fixed_formula,
    terms_fix     = fixed_terms,
    terms_sd      = sd_terms,
    xlevels_fix   = xlevels_fix,
    xlevels_sd    = xlevels_sd,
    contrasts_fix = contrasts_fix,
    contrasts_sd  = contrasts_sd
  ))
}
