#' Biplot for rrGDM objects
#'
#' @description
#' \code{\link[stats]{biplot}} method for objects of class \code{rrGDM},
#' showing sample positions in two of the fitted latent dimensions,
#' optionally overlaid with environmental loadings and/or species
#' correlations with latent axes.
#'
#' @param x A fitted model of class \code{rrGDM}.
#' @param k Integer vector of length 2 indicating which latent dimensions to
#'   plot. Default is \code{c(1, 2)}.
#' @param rotate Logical indicating whether to apply varimax rotation
#'   (via \code{\link[stats]{varimax}}) before plotting. Default is
#'   \code{TRUE}. The result of varimax rotation is that ordination axes are ordered relative to their
#'   associated importance.
#' @param type Character string: either \code{'conditional'} (site scores
#'   including sample-level error, \code{L}) or \code{'marginal'} (site
#'   scores as determined by environmental loadings only, \code{L_m}). Default is
#'   \code{'conditional'}.
#' @param col_by Optional character string naming a variable
#'   used to colour points. Ignored if \code{add_points = FALSE}.
#' @param shape_by Optional character string naming a variable in
#'   \code{x$X} used to shape points. Ignored if \code{add_points = FALSE}.
#' @param add_points Logical indicating whether to add the default sample
#'   points layer via \code{ggplot2::geom_point()}. Default is \code{TRUE}.
#'   Set to \code{FALSE} to omit this layer entirely and add a fully custom
#'   one instead. The returned plot's default data already contains the
#'   rotated latent coordinates together with all
#'   columns used to fit the model, so
#'   \code{biplot(model, add_points = FALSE) + ggplot2::geom_point(ggplot2::aes(fill = pH), shape = 21)}
#'   works directly, with more control over aesthetics, shapes, and legends.
#' @param env_load Logical indicating whether to add environmental loadings
#'   as arrows. Default is \code{TRUE}.
#' @param species_cor Logical indicating whether to add species labels
#'   positioned by their (Spearman) correlation with each plotted latent
#'   axis. Default is \code{FALSE}.
#' @param species_cor_n Optional integer. If \code{species_cor = TRUE}, only
#'   the \code{species_cor_n} species with the strongest correlation are
#'   labelled. Useful when the community matrix has hundreds or thousands
#'   of species. Default is \code{20}. Set to \code{NULL} to label all species.
#' @param env_load_scale Numeric scale factor applied to environmental
#'   loading arrows. Default is 1.
#' @param sp_load_scale Numeric scale factor applied to species-correlation
#'   labels. Default is 1.
#' @param text_size Numeric text size used for species labels. Default is 3.
#' @param compare Logical indicating whether to overlay marginal scores on
#'   top of conditional scores, connected by a segment. Only valid when
#'   \code{type = 'conditional'}. Default is \code{FALSE}.
#' @param ... Currently unused.
#'
#' @return
#' A \code{ggplot} object. Unlike
#' \code{\link[stats]{biplot.default}}, the object is
#' returned so it can be further customised (e.g.
#' \code{biplot(model) + ggplot2::theme_bw()}).
#'
#' @details
#' The function extracts fitted latent coordinates, optionally rotates them via varimax rotation
#' to improve interpretability, and plots the two latent dimensions
#' requested. Environmental loadings and species correlations,
#' if requested, are rotated using the same rotation matrix so that they
#' remain aligned with the sample scores.
#'
#' @examples
#' \dontrun{
#' com <- gllvm::microbialdata$Y
#' env <- gllvm::microbialdata$Xenv
#' env$pH_s <- scale(env$pH)
#'
#' m1 <- rrGDM(Y = com,
#'             X = env,
#'             k = 2,
#'             seed = 1,
#'             formula = ~ pH_s + as.factor(Region),
#'             family = 'binomial')
#'
#' biplot(m1,
#'        shape_by = 'Region',
#'        col_by = 'pH',
#'        type = 'conditional')
#'
#' # Custom:
#' biplot(m1, add_points = FALSE) +
#'   ggplot2::geom_point(ggplot2::aes(fill = pH), shape = 21, size = 3) +
#'   ggplot2::scale_fill_viridis_c()
#'
#' # Add species correlation with latent axes
#' biplot(m1, species_cor = TRUE, species_cor_n = 20)
#' }
#'
#' @importFrom grid arrow unit
#' @importFrom stats cor varimax biplot
#' @method biplot rrGDM
#' @export
biplot.rrGDM <- function(x,
                         k = c(1, 2),
                         rotate = TRUE,
                         type = 'conditional',
                         col_by = NULL,
                         shape_by = NULL,
                         add_points = TRUE,
                         env_load = TRUE,
                         species_cor = FALSE,
                         species_cor_n = 20,
                         env_load_scale = 1,
                         sp_load_scale = 1,
                         text_size = 3,
                         compare = FALSE,
                         ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort(c(
      "x" = "biplot method for rrGDM objects requires the {.pkg ggplot2} package, which is not installed.",
      "i" = "You can install it via: {.code install.packages(\"ggplot2\")}."
    ))
  }

  # ---------------------- Initial checks ---------------------

  if (!inherits(x, "rrGDM")) {
    cli::cli_abort("{.arg x} must be a fitted {.cls rrGDM} model.")
  }
  if (!is.numeric(k) || length(k) != 2 || anyNA(k) || any(k != as.integer(k))) {
    cli::cli_abort("{.arg k} must be an integer vector of length 2 (e.g. {.code c(1, 2)}).")
  }

  k <- as.integer(k)

  if (k[1] == k[2]) {
    cli::cli_abort("{.arg k} must reference two distinct latent dimensions.")
  }
  if (is.null(x$k) || x$k < 2) {
    cli::cli_abort(c(
      "x" = "A biplot requires a model fitted with at least 2 latent dimensions.",
      "i" = "This model was fitted with {.code k = {x$k}}."
    ))
  }
  if (any(k < 1) || any(k > x$k)) {
    cli::cli_abort("{.arg k} must reference dimensions between 1 and {x$k} (the model's number of latent dimensions).")
  }
  if (!is.logical(rotate) || length(rotate) != 1 || is.na(rotate)) {
    cli::cli_abort("{.arg rotate} must be a single {.val TRUE}/{.val FALSE} value.")
  }

  type <- match.arg(type, c('conditional', 'marginal'))

  if (!is.logical(env_load) || !is.logical(species_cor) || !is.logical(compare) || !is.logical(add_points)) {
    cli::cli_abort("{.arg env_load}, {.arg species_cor}, {.arg compare} and {.arg add_points} must be logical.")
  }
  if (!is.numeric(env_load_scale) || !is.numeric(sp_load_scale) || !is.numeric(text_size)) {
    cli::cli_abort("{.arg env_load_scale}, {.arg sp_load_scale} and {.arg text_size} must be numeric.")
  }

  if (!missing(species_cor_n) && !species_cor) {
    cli::cli_warn("{.arg species_cor_n} is ignored when {.code species_cor = FALSE}.")
  }

  if (species_cor && !is.null(species_cor_n)) {
    if (!is.numeric(species_cor_n) || length(species_cor_n) != 1 ||
        is.na(species_cor_n) || species_cor_n != as.integer(species_cor_n) ||
        species_cor_n < 1) {
      cli::cli_abort("{.arg species_cor_n} must be a single positive integer, or {.val NULL}.")
    }
    species_cor_n <- as.integer(species_cor_n)
  }

  if (compare && type != 'conditional') {
    cli::cli_abort(c(
      "x" = "{.arg compare = TRUE} is only supported when {.code type = 'conditional'}.",
      "i" = "Set {.code type = 'conditional'} to overlay marginal scores on conditional scores."
    ))
  }

  if (!add_points && (!is.null(col_by) || !is.null(shape_by))) {
    cli::cli_warn("{.arg col_by}/{.arg shape_by} are ignored when {.code add_points = FALSE}.")
  }

  if (add_points) {
    # X columns available for col_by / shape_by, robust to X being a matrix or data.frame
    x_names <- if (is.null(x$X)) character(0) else colnames(as.data.frame(x$X, check.names = FALSE))

    if (!is.null(col_by) && !(col_by %in% x_names)) {
      cli::cli_warn("{.arg col_by} ({.val {col_by}}) was not found in {.field x$X}; ignoring.")
      col_by <- NULL
    }
    if (!is.null(shape_by) && !(shape_by %in% x_names)) {
      cli::cli_warn("{.arg shape_by} ({.val {shape_by}}) was not found in {.field x$X}; ignoring.")
      shape_by <- NULL
    }
  }

  # ---------------------------- EXTRACT REPORT ----------------------------

  report <- tryCatch(x$obj$report(), error = function(e) {
    cli::cli_abort(c(
      "x" = "Could not extract model report from {.arg x$obj}.",
      "i" = "Original error: {conditionMessage(e)}"
    ))
  })

  points <- switch(type,
                   conditional = report$L,
                   marginal    = report$L_m)

  points <- as.matrix(points)

  # ---------------------------- ROTATION ----------------------------

  rot <- diag(x$k)

  #  Try to rotate
  if (rotate) {
    rot <- tryCatch(
      stats::varimax(points, normalize = FALSE)$rotmat,
      error = function(e) {
        cli::cli_warn(c(
          "!" = "Varimax rotation failed; plotting unrotated scores instead.",
          "i" = "Original error: {conditionMessage(e)}"
        ))
        diag(x$k)
      }
    )
  }

  pts_mat <- (points %*% rot)[, k, drop = FALSE]
  colnames(pts_mat) <- c('LV.1', 'LV.2')
  points_rot <- data.frame(pts_mat, x$X, check.names = FALSE)

  # ----------------- BASE PLOT --------------------------
  # Empty plot
  plot_obj <- ggplot2::ggplot(points_rot, ggplot2::aes(x = .data$LV.1, y = .data$LV.2)) +
    ggplot2::coord_equal() +
    ggplot2::xlab(paste0('LV', k[1])) +
    ggplot2::ylab(paste0('LV', k[2]))

  # Specified col or shape?
  has_col <- add_points && !is.null(col_by)
  has_shape <- add_points && !is.null(shape_by)

  if (add_points) {
    if (has_col && has_shape) {
      plot_obj <- plot_obj +
        ggplot2::geom_point(ggplot2::aes(color = .data[[col_by]], shape = as.factor(.data[[shape_by]]))) +
        ggplot2::labs(color = col_by, shape = shape_by)
    } else if (has_col) {
      plot_obj <- plot_obj +
        ggplot2::geom_point(ggplot2::aes(color = .data[[col_by]])) +
        ggplot2::labs(color = col_by)
    } else if (has_shape) {
      plot_obj <- plot_obj +
        ggplot2::geom_point(ggplot2::aes(shape = .data[[shape_by]])) +
        ggplot2::labs(shape = shape_by)
    } else {
      plot_obj <- plot_obj + ggplot2::geom_point()
    }
  }

  # ---------------------------- ADD ENV LOADINGS ----------------------------

  if (env_load) {
    if (is.null(x$X_design) || ncol(x$X_design) == 0) {
      cli::cli_warn("{.arg env_load = TRUE} but the model has no predictors in {.field X_design}; skipping loadings.")
    } else if (is.null(report$lambda)) {
      cli::cli_warn("Model report does not contain loadings {.field lambda}; skipping environmental loadings.")
    } else {
      load_mat <- as.matrix(as.matrix(report$lambda) %*% rot * env_load_scale)[, k, drop = FALSE]
      colnames(load_mat) <- c('LV.1', 'LV.2')
      loadings_rot <- data.frame(var = colnames(x$X_design), load_mat, check.names = FALSE)

      plot_obj <- plot_obj +
        ggplot2::geom_segment(data = loadings_rot,
                              ggplot2::aes(x = 0, y = 0, xend = .data$LV.1, yend = .data$LV.2),
                              arrow = grid::arrow(angle = 15, length = grid::unit(0.2, 'cm'))) +
        ggplot2::geom_text(data = loadings_rot,
                           size = text_size,
                           ggplot2::aes(label = .data$var, x = .data$LV.1 * 1.2, y = .data$LV.2 * 1.2))
    }
  }

  # ---------------------------- SP CORRELATIONS ----------------------------

  if (species_cor) {
    if (is.null(x$Y)) {
      cli::cli_abort(c(
        "x" = "{.arg species_cor = TRUE} requires the original community matrix.",
        "i" = "No {.field Y} component was found on {.arg x}."
      ))
    }

    Y_com <- as.matrix(x$Y)

    # Check in case rows have been removed in the data prep
    if (nrow(Y_com) != nrow(points_rot)) {
      cli::cli_abort(c(
        "x" = "The community matrix ({.field x$Y}) does not match the number of fitted samples.",
        "i" = "{.field x$Y} has {nrow(Y_com)} row{?s}, but the model has {nrow(points_rot)} sample{?s}."
      ))
    }

    sp_names <- colnames(Y_com)
    if (is.null(sp_names)) sp_names <- paste0('sp', seq_len(ncol(Y_com)))

    sp_load <- cbind(
      LV.1 = as.numeric(stats::cor(Y_com, points_rot$LV.1, method = 'spearman')),
      LV.2 = as.numeric(stats::cor(Y_com, points_rot$LV.2, method = 'spearman'))
    )

    sp_load_df <- data.frame(sp_load, sp = sp_names, check.names = FALSE)

    if (!is.null(species_cor_n) && nrow(sp_load_df) > species_cor_n) {
      cli::cli_warn(c(
        "!" = "Labelling only the {species_cor_n} species with the strongest correlation ({nrow(sp_load_df)} available).",
        "i" = "Increase {.arg species_cor_n}, or set it to {.val NULL}, to label more."
      ))
      strength <- sqrt(rowSums(sp_load^2))
      keep <- order(strength, decreasing = TRUE)[seq_len(species_cor_n)]
      sp_load_df <- sp_load_df[keep, , drop = FALSE]
    }

    plot_obj <- plot_obj +
      ggplot2::geom_text(data = sp_load_df,
                         colour = 'red',
                         size = text_size,
                         ggplot2::aes(label = .data$sp,
                                      x = .data$LV.1 * sp_load_scale,
                                      y = .data$LV.2 * sp_load_scale))
  }

  # ---------------------------- COND. VS MARG. ----------------------------

  if (compare) {
    if (is.null(report$L_m)) {
      cli::cli_warn("{.arg compare = TRUE} The fitted model does not contain marginal scores; skipping comparison.")
    } else {
      lm_mat <- (as.matrix(report$L_m) %*% rot)[, k, drop = FALSE]
      colnames(lm_mat) <- c('LV_m.1', 'LV_m.2')
      points_m <- data.frame(points_rot, lm_mat, check.names = FALSE)

      plot_obj <- plot_obj +
        ggplot2::geom_point(data = points_m,
                            shape = 21,
                            ggplot2::aes(x = .data$LV_m.1, y = .data$LV_m.2, fill = 'marginal')) +
        ggplot2::geom_segment(data = points_m,
                              ggplot2::aes(x = .data$LV.1, y = .data$LV.2,
                                           xend = .data$LV_m.1, yend = .data$LV_m.2),
                              colour = 'grey')
    }
  }

  return(plot_obj)
}
