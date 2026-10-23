#' Prepare design matrices and data (Internal)
#' @importFrom stats model.frame model.matrix terms update
#' @importFrom reformulas nobars findbars
#' @importFrom Matrix sparseMatrix sparse.model.matrix
#' @noRd
prepare_data <- function(formula = NULL,
                         formula_sd = NULL,
                         data = NULL,
                         n = NULL) {

  if (is.null(n)) {
    cli::cli_abort(c("{.var n} must be supplied",
                     "x" = "The number of observations is required for design matrix dimensions."))
  }

  # check formulas
  if (!is.null(formula) && !inherits(formula, "formula")) {
    cli::cli_abort("{.var formula} must be a {.cls formula}, not a {.cls {class(formula)}}.")
  }

  if (!is.null(formula_sd) && !inherits(formula_sd, "formula")) {
    cli::cli_abort("{.var formula_sd} must be a {.cls formula}, not a {.cls {class(formula_sd)}}.")
  }

  if (!is.null(formula_sd)) {
    if (attr(stats::terms(formula_sd), "intercept") == 0) {
      cli::cli_abort(c(
        "{.var formula_sd} must include an intercept.",
        "!" = "Check your {.var formula_sd} and avoid {.code 0} or {.code -1} from your {.var formula_sd}."
      ))
    }
  }

  if (!is.null(data) && !(is.matrix(data) || is.data.frame(data))) {
    cli::cli_abort("{.var data} must be a {.cls matrix} or {.cls data.frame}.")
  }
  if (is.matrix(data)) data <- as.data.frame(data)

  if (!is.null(data) && nrow(data) != n) {
    cli::cli_abort(c(
      "x" = "{.arg X} and {.arg Y} must have the same number of rows (sites).",
      "i" = "{.arg X} has {nrow(data)} row{?s}; {.arg Y} has {n}."
    ))
  }

  if (is.null(formula_sd) || is.null(data)) {
    formula_sd  <- ~1
  }

  # Initialise terms to NULL / EMPTY
  fixed_terms <- stats::terms(~0)
  sd_terms <- stats::terms(~1)
  xlevels_fix <- NULL
  xlevels_sd <- NULL
  contrasts_fix <- NULL
  contrasts_sd <- NULL

  if (is.null(data)) {
    X_design <- matrix(numeric(0), nrow = n, ncol = 0)
    W_design <- matrix(1, nrow = n, ncol = 1, dimnames = list(NULL, "(Intercept)"))
    fixed_formula <- ~0

  } else {

    if (is.null(formula)) {
      X_design <- matrix(numeric(0), nrow = n, ncol = 0)
      fixed_formula <- ~0

    } else {
      # --- FIXED EFFECTS ---
      fixed_formula <- reformulas::nobars(formula)
      mf_fix <- stats::model.frame(fixed_formula, data = data, na.action = stats::na.pass)
      bad <- names(mf_fix)[vapply(mf_fix, anyNA, logical(1))]
      if (length(bad) > 0) cli::cli_abort("Missing values in {.arg formula} variable{?s}: {.var {bad}}.")
      fixed_terms <- attr(mf_fix, "terms")
      X_tmp <- stats::model.matrix(fixed_terms, data = mf_fix)

      xlevels_fix <- stats::.getXlevels(fixed_terms, mf_fix)


      if ("(Intercept)" %in% colnames(X_tmp)) {
        intercept_idx <- which(colnames(X_tmp) == "(Intercept)")
        X_design <- X_tmp[, -intercept_idx, drop = FALSE]
        attr(X_design, "assign") <- attr(X_tmp, "assign")[-intercept_idx]
        attr(X_design, "contrasts") <- attr(X_tmp, "contrasts")
      } else {
        X_design <- X_tmp
      }

      contrasts_fix <- attr(X_tmp, "contrasts")
    }

    # --- DISPERSION MODEL ---
    if (!is.null(reformulas::findbars(formula_sd))) {
      cli::cli_abort(c("Random effects not allowed in dispersion model",
                       "x" = "Check {.var formula_sd}: random effects are only supported in the main formula."))
    }

    sd_formula <- stats::update(formula_sd, ~ 1 + .)
    mf_sd <- stats::model.frame(sd_formula, data = data, na.action = stats::na.pass)
    bad <- names(mf_sd)[vapply(mf_sd, anyNA, logical(1))]
    if (length(bad) > 0) cli::cli_abort("Missing values in {.arg formula_sd} variable{?s}: {.var {bad}}.")
    sd_terms <- attr(mf_sd, "terms")
    W_design <- stats::model.matrix(sd_terms, data = mf_sd)

    xlevels_sd <- stats::.getXlevels(sd_terms, mf_sd)
    contrasts_sd <- attr(W_design, "contrasts")
  }

  # Random parts
  random_parts <- if (!is.null(formula)) reformulas::findbars(formula) else NULL

  if ((length(random_parts) == 0) || is.null(data)) {
    Z_design <- Matrix::sparseMatrix(
      i = integer(0), j = integer(0), x = numeric(0), dims = c(n, 1)
    )

    return(list(
      X_design = X_design,
      W_design = W_design,
      Z_design = Z_design,
      map_re = integer(0),
      n_re = 0,
      n_re_vars = 0,
      re_vars = character(0),
      has_random = FALSE,
      formula_fix = fixed_formula,
      terms_fix = fixed_terms,
      terms_sd = sd_terms,
      xlevels_fix = xlevels_fix,
      xlevels_sd = xlevels_sd,
      contrasts_fix = contrasts_fix,
      contrasts_sd  = contrasts_sd
    ))
  }

  # Random effects
  Z_list <- list()
  re_vars <- character()
  n_groups_vars <- integer()

  for (re in random_parts) {
    if (re[[2]] == 1) {
      var <- deparse(re[[3]])
      if (!(var %in% colnames(data))) {
        cli::cli_abort(c(
          "Variable {.val {var}} not found in the provided data.",
          if (grepl(":", var, fixed = TRUE))
            c("i" = "Nested or crossed groupings such as {.code (1 | a/b)} are not supported: create a single grouping variable in your data.")
        ))
      }

      if (var %in% re_vars) {
        cli::cli_abort("Grouping variable {.val {var}} appears in more than one random-effect term.")
      }

      if (anyNA(data[[var]])) {
        cli::cli_abort("Grouping variable {.val {var}} has missing values.")
      }

      factor_vec <- factor(data[[var]])
      re_levels <- levels(factor_vec)
      n_groups <- length(re_levels)
      if (n_groups < 2) {
        cli::cli_abort("Grouping variable {.val {var}} must have at least 2 levels.")
      }
      Z_re <- Matrix::sparse.model.matrix(~0 + factor_vec)
      colnames(Z_re) <- paste0(var, "_", re_levels)

      Z_list[[var]] <- Z_re
      re_vars <- c(re_vars, var)
      n_groups_vars <- c(n_groups_vars, n_groups)
    } else {
      cli::cli_abort(c("Random slopes are not supported",
                       "i" = "Only random intercepts in the form {.code (1 | group)} are currently allowed."))
    }
  }

  Z_design <- do.call(cbind, Z_list)
  map_re <- rep(seq_along(re_vars), n_groups_vars)

  return(list(
    X_design = X_design,
    W_design = W_design,
    Z_design = Z_design,
    map_re = map_re,
    n_re = length(map_re),
    n_re_vars = length(re_vars),
    re_vars = re_vars,
    has_random = TRUE,
    formula_fix = fixed_formula,
    terms_fix = fixed_terms,
    terms_sd = sd_terms,
    xlevels_fix = xlevels_fix,
    xlevels_sd = xlevels_sd,
    contrasts_fix = contrasts_fix,
    contrasts_sd  = contrasts_sd
  ))
}
