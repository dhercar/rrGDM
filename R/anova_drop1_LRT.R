#' Likelihood Ratio Test for rrGDM Models
#'
#' @description
#' Compares two fitted `rrGDM` models using a likelihood ratio test (LRT).
#' The LRT is useful to compare nested models to assess the significance of one or more predictors.
#'
#' @param object  A fitted model of class `rrGDM`. **Should be the fuller
#'   (more parameters) model**.
#' @param object2 A fitted model of class `rrGDM`. Should be the reduced
#'   (fewer parameters).
#' @param ... Currently unused.
#'
#' @details
#' The function does **not** verify that the models are truly nested. A
#' warning is issued when the models have equal numbers of parameters, but
#' no structural check is performed.
#'
#' @return
#' A `data.frame` with one row per model.
#'
#' @seealso [drop1.rrGDM()] for sequential single-term deletion tests.
#'
#' @export
#' @method drop1 rrGDM
anova.rrGDM <- function(object, object2, ...) {
  if (missing(object2)) {
    cli::cli_abort(
      "Requires two models to compare."
    )
  }

  name_1 <- deparse(substitute(object))
  name_2 <- deparse(substitute(object2))

  n_par_1 <- length(object$opt$par)
  n_par_2 <- length(object2$opt$par)


  if (n_par_1 < n_par_2) {
    cli::cli_warn(c(
      "!" = "First model has fewer parameters than the second.",
      "i" = "Convention is to pass the fuller model first: {.code anova(full_model, reduced_model)}."
    ))
  }

  nll_1   <- object$opt$objective
  nll_2   <- object2$opt$objective
  lrt     <- 2 * (nll_2 - nll_1)
  ll_dif  <- nll_2 - nll_1
  df_diff <- n_par_1 - n_par_2             # negative when order reversed

  if (lrt < -1e-6) {
    cli::cli_warn(c(
      "!" = "LRT statistic is negative ({round(lrt, 4)})",
      "i" = "The models supplied may ve in the wrong order, or may not have converged."
    ))
  }

  p_val <- pchisq(lrt, df = abs(df_diff), lower.tail = FALSE)

  res <- data.frame(
    Model    = c(name_1, name_2),
    n_par    = c(n_par_1, n_par_2),
    nll      = round(c(nll_1, nll_2), 2),
    Df       = c(NA, df_diff),
    ll_dif   = round(c(NA, ll_dif), 3),
    LRT      = round(c(NA, max(0, lrt)), 3),
    Pr_Chisq = c(NA, p_val)
  )

  return(res)
}


#' Drop-One Likelihood Ratio Tests for Fixed Effects in an rrGDM Model
#'
#' @description
#' Performs sequential single-term deletion tests for all fixed-effect terms
#' in a fitted `rrGDM` model. Each term is removed in turn, the model is
#' refitted, and the change in fit is assessed via [anova.rrGDM()].
#'
#' @param object A fitted model of class `rrGDM`. The full model from which
#'   terms will be dropped one at a time.
#' @param ... Currently unused.
#'
#' @details
#'
#' If a reduced model fails to converge or generates errors during fitting, a warning
#' is issued for that term and a row of `NA` values is returned in its place.
#'
#' @return
#' A `data.frame` with one row per fixed-effect term.
#'
#' @seealso [anova.rrGDM()] for a direct two-model comparison.
#'
#' @export
#' @method drop1 rrGDM
drop1.rrGDM <- function(object, ...) {
  fix_terms     <- attr(terms(object$formula_fix, data = object$data_frame), "term.labels")
  original_call <- object$call

  caller_env <- parent.frame()

  cli::cli_h2("Drop-1 LRT for Fixed Effects")
  cli::cli_inform("Testing {.val {length(fix_terms)}} term{?s}...")

  results <- lapply(fix_terms, function(term) {
    cli::cli_inform("Testing deletion of term: {.val {term}}")

    new_formula      <- update.formula(object$formula, as.formula(paste("~ . -", term)))
    new_call         <- original_call
    new_call$formula <- new_formula

    fit_red <- try(eval(new_call, envir = caller_env), silent = FALSE)

    if (inherits(fit_red, "try-error")) {
      cli::cli_warn("Model failed to converge without {.val {term}}.")
      return(data.frame(
        Term = term, nll = NA, ll_dif = NA,
        LRT = NA, Df = NA, Pr_Chisq = NA
      ))
    }

    # object is always the full model; fit_red is always reduced.
    # Suppress the expected "fewer parameters" warning from anova.rrGDM.
    # Index by position rather than model name since names now vary.
    lrt_res <- suppressWarnings(anova(object, fit_red))
    red_row <- lrt_res[2, ]

    data.frame(
      Term     = term,
      nll      = red_row$nll,
      ll_dif   = red_row$ll_dif,
      LRT      = red_row$LRT,
      Df       = red_row$Df,
      Pr_Chisq = red_row$Pr_Chisq
    )
  })

  res_df <- do.call(rbind, results)
  res_df$Sig <- cut(
    res_df$Pr_Chisq,
    breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
    labels = c("***", "**", "*", ".", " ")
  )

  cli::cli_h2("Single Term Deletions (LRT)")

  return(res_df)
}
