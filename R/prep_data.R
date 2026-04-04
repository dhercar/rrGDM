#' Prepare design matrices for TMB (Internal)
#' @noRd
prepare_data <- function(formula = NULL,
                         formula_sd = NULL,
                         data = NULL,
                         n = NULL) {

  # ---- Input Validation using cli ----
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

  # ---- Logic start ----
  if (is.null(formula_sd) || is.null(data)) {
    formula_sd  <- ~1
  }

  if (is.null(formula) || is.null(data)) {
    X_design = integer(0)
    W_design = integer(0)
    fixed_formula = ~0

  } else {
    # Fixed part
    fixed_formula <- lme4::nobars(formula)
    X_tmp <- model.matrix(fixed_formula, data)

    if ('(Intercept)' %in% colnames(X_tmp)) {
      assign_attr <- attr(X_tmp, "assign")[-1]
      contr_attr <- attr(X_tmp, "contrasts")
      X_design <- X_tmp[, !colnames(X_tmp) %in% "(Intercept)", drop = FALSE] # remove intercept
      attr(X_design, "assign") <- assign_attr
      attr(X_design, "contrasts") <- contr_attr
    } else {
      X_design <- X_tmp
    }

    X_design

    # SD eq
    if (!is.null(lme4::findbars(formula_sd))) {
      cli::cli_abort(c("Random effects not allowed in dispersion model",
                       "x" = "Check {.var formula_sd}: random effects are only supported in the main formula."))
    }
    W_design <- model.matrix(update(formula_sd, ~ 1+.), data)
  }

  # Random parts
  random_parts <- lme4::findbars(formula)

  # Return with no re
  if ((length(random_parts) == 0) || is.null(data)) {

    # empty Z matrix (expected by TMB model)
    Z_design <- Matrix::sparseMatrix(
      i = integer(0),
      j = integer(0),
      x = numeric(0),
      dims = c(n, 1)
    )

    return(list(X_design = X_design,
                W_design = W_design,
                Z_design = Z_design,
                map_re = integer(0),
                n_re = 0,
                n_re_vars = 0,
                has_random = FALSE,
                formula_fix = fixed_formula))
  }

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

      # Create sparse model matrix for re
      Z_re <- Matrix::sparse.model.matrix(~0 + factor(data[[var]]))
      colnames(Z_re) <- paste0(var, "_", re_levels)

      # Add to lists
      Z_list[[var]] <- Z_re
      re_vars <- c(re_vars, var)
      n_groups_vars <- c(n_groups_vars, n_groups)
    } else {
      cli::cli_abort(c("Random slopes are not supported",
                       "i" = "Only random intercepts in the form {.code (1 | group)} are allowed."))
    }
  }

  Z_design <- do.call(cbind, Z_list)

  # Create mapping from columns of Z to random effect groups
  map_re <- rep(seq_along(re_vars), n_groups_vars)

  return(list(
    X_design = X_design,
    W_design = W_design,
    Z_design = Z_design,
    map_re = map_re,
    n_re = length(map_re),
    n_re_vars = length(re_vars),
    has_random = TRUE,
    formula_fix = fixed_formula
  ))
}

