#' Select a sparse set of sample pairs for model fitting
#'
#' @description
#' Selects a subset of sample pairs to be used for fitting \code{\link{rrGDM}},
#' such that every sample is in exactly the same number of pairs
#' (\code{knn + n_random}). Each sample is paired with its \code{knn} nearest
#' neighbours, which inform its local position, and with \code{n_random}
#' random samples, which ensure that the global structure (including distant samples)
#' is represented. Selection is based on the observed
#' dissimilarities only.
#'
#' @param Y A community matrix or data frame (samples x species), or an object
#'   of class \code{\link[stats]{dist}}. If a \code{dist} object is supplied,
#'   \code{dist_method}, \code{trans} and \code{binary} are ignored.
#' @param knn Integer. Number of nearest neighbours per sample. Default is 10.
#' @param n_random Integer. Number of random pairs per sample. Default is 3.
#' @param max_zeros Integer or \code{Inf}. Maximum number of pairs with a
#'   dissimilarity of exactly 0 per sample. Default is \code{Inf}.
#' @param dist_method Dissimilarity index: a method of
#'   \code{\link[vegan]{vegdist}} (e.g. \code{'bray'}, \code{'jaccard'},
#'   \code{'euclidean'}); \code{'sorensen'} is Bray-Curtis on presence-absence
#'   data. Ignored if \code{Y} is a \code{dist} object. Default is \code{'bray'}.
#' @param trans Transformation applied to \code{Y} before calculating
#'   dissimilarities. Default is no transformation.
#' @param binary Logical indicating whether \code{Y} should be converted to
#'   presence-absence before calculating dissimilarities. Default is FALSE.
#' @param k_pcoa Integer (at least 2) or \code{NULL}. Number of PCoA axes of the
#'   full dissimilarity matrix stored for 2D representation. Use
#'   \code{NULL} to skip it. Default is 2.
#' @param link Character string: \code{'identity'} (default), \code{'negexp'} or
#'   \code{'neggaus'}. Link function of the model to be fitted (see
#'   \code{\link{rrGDM}}). The PCoA is computed on dissimilarities transformed to
#'   the scale of latent distances (e.g. \code{-log(1 - diss)} for
#'   \code{'negexp'}), after setting dissimilarities above 0.99 to 0.99. This
#'   does not affect which pairs are selected.
#' @param seed An integer passed to \code{set.seed()} at the start, so it changes
#'   the state of R's random number generator. Default is NULL (no seed is set).
#'
#' @return
#' An object of class \code{sparse_diss}: a list including \code{edges}, a
#' two-column integer matrix (\code{s1}, \code{s2}) of the pairs kept, which
#' \code{\link{rrGDM}} uses to subset the dissimilarities (argument
#' \code{sparse}), and \code{pairs}, a data frame with the pairs, their
#' dissimilarity (\code{diss}) and their type (\code{'knn'}, \code{'random'}
#' or \code{'zero'}). Each pair is stored once, with \code{s1 < s2}.
#' \code{degree} gives the number of pairs per sample. The other components are
#' the sample names (\code{id}), the number of samples and pairs (\code{n_samples},
#' \code{n_pairs_full}, \code{n_pairs}, \code{n_zero_pairs}), the settings used
#' (\code{knn}, \code{n_random}, \code{n_per_sample}, \code{max_zeros},
#' \code{link}, \code{dist_method}, \code{binary}), the number of connected groups
#' of samples (\code{n_components}), the PCoA used by \code{\link{plot.sparse_diss}}
#' (\code{pcoa}, \code{NULL} if \code{k_pcoa = NULL}) and the matched call
#' (\code{call}).
#'
#' @details
#' Pairs are sorted by dissimilarity and kept if both samples
#' have fewer than \code{knn} pairs (and fewer than \code{max_zeros} for zero
#' dissimilarities). The remaining pairs up to \code{knn + n_random}, are
#' filled with random pairs among samples with free slots. Because the number of
#' pairs per sample is fixed, a sample may be paired with a sample that is not
#' among its own \code{knn} nearest neighbours. Pairs with a dissimilarity of
#' exactly 0 (identical samples) are only added in the nearest-neighbour step,
#' never in the random step.
#'
#' One sample has one pair fewer if both the number of samples and
#' \code{knn + n_random} are odd. A warning is given if some samples fall short
#' or if the pairs do not connect all samples (increase \code{n_random}).
#'
#' @seealso \code{\link{plot.sparse_diss}}, \code{\link{print.sparse_diss}}
#'
#' @importFrom stats cmdscale as.dist runif
#' @export
sparse_diss <- function(Y,
                        knn = 10,
                        n_random = 3,
                        max_zeros = Inf,
                        dist_method = 'bray',
                        trans = function(x) x,
                        binary = FALSE,
                        k_pcoa = 2,
                        link = 'identity',
                        seed = NULL) {

  #  Initial checks
  is_count <- function(v, allow_inf = FALSE) {
    if (!is.numeric(v) || length(v) != 1 || is.na(v) || v < 0) return(FALSE)
    (is.infinite(v) && allow_inf) || (is.finite(v) && v == round(v))
  }

  if (!(inherits(Y, 'dist') || is.matrix(Y) || is.data.frame(Y))) {
    cli::cli_abort(c("{.arg Y} must be a sample x species matrix, data frame or {.cls dist} object.",
                     "x" = "You have supplied a {.cls {class(Y)}}."))
  }

  if (!is_count(knn)) {
    cli::cli_abort("{.arg knn} must be a single non-negative integer.")
  }
  if (!is_count(n_random)) {
    cli::cli_abort("{.arg n_random} must be a single non-negative integer.")
  }
  if (!is_count(max_zeros, allow_inf = TRUE)) {
    cli::cli_abort("{.arg max_zeros} must be a single non-negative integer (or {.val Inf}).")
  }
  if (!is.null(k_pcoa) && (!is_count(k_pcoa) || k_pcoa < 2)) {
    cli::cli_abort("{.arg k_pcoa} must be {.code NULL} or a single integer of at least 2.")
  }
  n_per_sample <- as.integer(knn + n_random)
  if (n_per_sample < 1) {
    cli::cli_abort("{.arg knn} + {.arg n_random} must be at least 1.")
  }
  link <- match.arg(link, c('identity', 'negexp', 'neggaus'))
  if (!is.null(seed)) set.seed(seed)

  # Dissimilarity matrix
  if (inherits(Y, 'dist')) {
    D <- as.matrix(Y)
    n <- nrow(D)
    id <- rownames(D)
    dist_method <- 'user-supplied dist'
    binary <- NA
    if (is.null(id)) id <- as.character(seq_len(n))
  } else {
    n <- nrow(Y)
    id <- rownames(Y)
    if (is.null(id)) id <- as.character(seq_len(n))

    Y_df <- make_Y_diss(Y, method = dist_method, trans = trans, binary = binary)
    if (max(Y_df$s2) != n) {
      cli::cli_abort("Samples were removed while calculating dissimilarities. {.fn sparse_diss} requires all samples to be retained.")
    }
    D <- matrix(0, n, n)
    D[as.matrix(Y_df[, c('s1', 's2')])] <- Y_df$diss
    D <- D + t(D)
  }

  if (n < 3) cli::cli_abort("At least 3 samples are needed to select pairs.")
  if (anyNA(D) || any(!is.finite(D))) {
    cli::cli_abort("Dissimilarities contain missing or non-finite values.")
  }
  if (n_per_sample > n - 1) {
    cli::cli_abort(c("{.arg knn} + {.arg n_random} ({n_per_sample}) must be smaller than the number of samples ({n}).",
                     "i" = "Each sample can be paired with at most {n - 1} other samples."))
  }

  n_pairs_full <- n * (n - 1) / 2
  is_zero <- (D == 0)
  diag(is_zero) <- FALSE

  #Similar samples
  ij <- which(upper.tri(D), arr.ind = TRUE)
  ord <- order(D[ij], stats::runif(nrow(ij)))
  a <- ij[ord, 1]
  b <- ij[ord, 2]
  zero <- is_zero[ij][ord]

  keep <- logical(length(a))
  free <- rep(as.integer(knn), n)
  zero_free <- rep(max_zeros, n)
  n_open <- if (knn > 0) n else 0L

  if (n_open > 1L) {
    for (e in seq_along(a)) {
      i <- a[e]
      j <- b[e]
      if (free[i] == 0L || free[j] == 0L) next
      if (zero[e]) {
        if (zero_free[i] < 1 || zero_free[j] < 1) next
        zero_free[i] <- zero_free[i] - 1
        zero_free[j] <- zero_free[j] - 1
      }
      keep[e] <- TRUE
      free[i] <- free[i] - 1L
      free[j] <- free[j] - 1L
      n_open <- n_open - (free[i] == 0L) - (free[j] == 0L)
      if (n_open < 2L) break
    }
  }

  # Random pairs
  paired <- matrix(FALSE, n, n)
  paired[cbind(a[keep], b[keep])] <- TRUE
  paired <- paired | t(paired)
  is_near <- paired
  slots <- n_per_sample - rowSums(paired)
  ok_pair <- function(x, y) x != y & !paired[cbind(x, y)] & !is_zero[cbind(x, y)]

  repeat {
    open <- which(slots > 0L)
    if (length(open) == 0L) break
    i <- open[which.max(slots[open] + stats::runif(length(open)))]

    cand <- open[open != i & !paired[i, open] & !is_zero[i, open]]
    m <- min(slots[i], length(cand))
    if (m > 0L) {
      j <- cand[sample.int(length(cand), m, prob = slots[cand])]
      paired[i, j] <- TRUE
      paired[j, i] <- TRUE
      slots[j] <- slots[j] - 1L
      slots[i] <- slots[i] - m
    }
    if (slots[i] == 0L) next

    others <- setdiff(which(slots > 0L), i)
    j <- if (length(others) > 0L) others[sample.int(length(others), 1L)] else i
    E <- which(paired & upper.tri(paired), arr.ind = TRUE)
    E <- E[!is_zero[E], , drop = FALSE]
    fits <- list(ok_pair(rep(i, nrow(E)), E[, 1]) & ok_pair(rep(j, nrow(E)), E[, 2]),
                 ok_pair(rep(i, nrow(E)), E[, 2]) & ok_pair(rep(j, nrow(E)), E[, 1]))
    if (j == i && slots[i] < 2L || !(any(fits[[1]]) || any(fits[[2]]))) {
      slots[i] <- 0L  # cannot be completed
      next
    }
    flip <- if (!any(fits[[1]])) 2L else if (!any(fits[[2]])) 1L else sample.int(2L, 1L)
    r <- which(fits[[flip]])[sample.int(sum(fits[[flip]]), 1L)]
    cc <- E[r, flip]
    dd <- E[r, 3L - flip]
    paired[E[r, 1], E[r, 2]] <- FALSE
    paired[E[r, 2], E[r, 1]] <- FALSE
    is_near[E[r, 1], E[r, 2]] <- FALSE
    is_near[E[r, 2], E[r, 1]] <- FALSE
    paired[i, cc] <- TRUE
    paired[cc, i] <- TRUE
    paired[j, dd] <- TRUE
    paired[dd, j] <- TRUE
    slots[i] <- slots[i] - 1L
    slots[j] <- slots[j] - 1L
  }

  # Assemble pairs
  idx <- which(paired & upper.tri(paired), arr.ind = TRUE)
  idx <- idx[order(idx[, 1], idx[, 2]), , drop = FALSE]
  s1 <- as.integer(idx[, 1])
  s2 <- as.integer(idx[, 2])
  type <- ifelse(is_zero[idx], 'zero', ifelse(is_near[idx], 'knn', 'random'))

  pairs <- data.frame(s1 = s1,
                      s2 = s2,
                      diss = D[cbind(s1, s2)],
                      type = type,
                      stringsAsFactors = FALSE)

  #  Graph summaries
  degree <- tabulate(c(s1, s2), nbins = n)
  names(degree) <- id
  n_components <- max(graph_components(n, s1, s2))

  if (any(degree < n_per_sample)) {
    n_short <- sum(degree < n_per_sample)
    cli::cli_warn(c(
      "!" = "{n_short} sample{?s} ha{?s/ve} fewer than {n_per_sample} pairs.",
      "i" = "This happens if {.code n * (knn + n_random)} is odd, or when few candidate pairs are left (small data sets)."
    ))
  }
  if (n_components > 1) {
    cli::cli_warn(c(
      "!" = "The selected pairs form {n_components} disconnected groups of samples.",
      "i" = "Their relative positions are not constrained by the data. Increase {.arg n_random}."
    ))
  }

  #  PCoA
  pcoa_res <- NULL
  if (!is.null(k_pcoa)) {
    k_pcoa <- min(as.integer(k_pcoa), n - 1L)
    D_link <- switch(link,
                     'identity' = D,
                     'negexp' = -log(1 - pmin(D, 0.99)),
                     'neggaus' = sqrt(-log(1 - pmin(D, 0.99))))
    pc <- tryCatch(suppressWarnings(stats::cmdscale(stats::as.dist(D_link), k = k_pcoa, eig = TRUE)),
                   error = function(e) {
                     cli::cli_warn(c("!" = "PCoA could not be computed. {.fn plot} will not be available.",
                                     "i" = "Original error: {conditionMessage(e)}"))
                     NULL
                   })
    if (!is.null(pc) && ncol(pc$points) >= 2) {
      k_got <- ncol(pc$points)
      if (k_got < k_pcoa) {
        cli::cli_warn("Only {k_got} PCoA axes with positive eigenvalues could be computed ({.arg k_pcoa} = {k_pcoa}).")
      }
      pts <- pc$points
      colnames(pts) <- paste0('PCoA', seq_len(k_got))
      pos_eig <- pc$eig[pc$eig > 0]
      pcoa_res <- list(coords = data.frame(id = id, pts, stringsAsFactors = FALSE),
                       var_exp = pc$eig[seq_len(k_got)] / sum(pos_eig),
                       k = k_got)
    } else if (!is.null(pc)) {
      cli::cli_warn("Fewer than 2 PCoA axes with positive eigenvalues. {.fn plot} will not be available.")
    }
  }

  res <- list(
    edges = cbind(s1 = s1, s2 = s2),
    pairs = pairs,
    id = id,
    n_samples = n,
    n_pairs_full = n_pairs_full,
    n_pairs = nrow(pairs),
    n_zero_pairs = sum(is_zero[upper.tri(is_zero)]),
    knn = knn,
    n_random = n_random,
    n_per_sample = n_per_sample,
    max_zeros = max_zeros,
    link = link,
    dist_method = dist_method,
    binary = binary,
    degree = degree,
    n_components = n_components,
    pcoa = pcoa_res,
    call = match.call()
  )

  class(res) <- 'sparse_diss'
  return(res)
}


#' Connected components of a graph (internal)
#'
#' @param n Number of nodes.
#' @param s1,s2 Integer vectors with the two ends of each edge.
#' @return Integer vector with the component of each node.
#' @noRd
graph_components <- function(n, s1, s2) {
  adj <- split(c(s2, s1), factor(c(s1, s2), levels = seq_len(n)))
  comp <- integer(n)
  n_comp <- 0L

  for (start in seq_len(n)) {
    if (comp[start] > 0L) next
    n_comp <- n_comp + 1L
    queue <- start
    comp[start] <- n_comp
    while (length(queue) > 0) {
      nb <- unlist(adj[queue], use.names = FALSE)
      nb <- unique(nb[comp[nb] == 0L])
      comp[nb] <- n_comp
      queue <- nb
    }
  }
  comp
}


#' Print method for sparse_diss objects
#'
#' @description
#' Prints a compact overview of a \code{sparse_diss} object: the number of
#' samples and pairs (full and included), the selection settings, the pairs
#' kept by type and the number of pairs per sample.
#'
#' @param x An object of class \code{sparse_diss}, as returned by
#'   \code{\link{sparse_diss}}.
#' @param digits Number of significant digits used to print dissimilarities.
#'   Default is 3.
#' @param ... Currently unused.
#'
#' @method print sparse_diss
#' @export
print.sparse_diss <- function(x, digits = 3, ...) {

  if (!inherits(x, 'sparse_diss')) {
    cli::cli_abort("{.arg x} must be a {.cls sparse_diss} object.")
  }

  show <- function(label, value) cat(sprintf("  %-26s %s\n", label, value))
  fmt_int <- function(v) format(v, big.mark = ',', scientific = FALSE, trim = TRUE)

  p <- x$pairs

  cat("Sparse pair selection (sparse_diss)\n\n")
  cat("Call:\n", paste(deparse(x$call), collapse = "\n"), "\n", sep = "")

  # Pairs
  cli::cli_h2("Pairs")
  show("Dissimilarity:", x$dist_method)
  show("Samples:", fmt_int(x$n_samples))
  show("Pairs (full):", fmt_int(x$n_pairs_full))
  show("Pairs included (%):", sprintf("%s (%.1f%%)", fmt_int(x$n_pairs), 100 * x$n_pairs / x$n_pairs_full))

  # Settings
  cli::cli_h2("Selection")
  show("Pairs per sample:", sprintf("%d (%s knn + %s random)", x$n_per_sample, x$knn, x$n_random))
  show("Max zeros per sample:", x$max_zeros)
  show("PCoA axes stored:", if (is.null(x$pcoa)) 'none' else sprintf("%d (link: %s)", x$pcoa$k, x$link))

  # Pairs by type
  cli::cli_h2("Pairs by type")
  types <- c('knn', 'random', 'zero')
  d <- data.frame(type = types,
                  n = vapply(types, function(t) sum(p$type == t), numeric(1)),
                  stringsAsFactors = FALSE)
  d$pct <- sprintf("%.1f%%", 100 * d$n / x$n_pairs)
  d$mean_diss <- vapply(types, function(t) if (any(p$type == t)) mean(p$diss[p$type == t]) else NA_real_, numeric(1))
  d$mean_diss <- ifelse(is.na(d$mean_diss), '-', formatC(d$mean_diss, digits = digits, format = 'fg'))
  d$n <- fmt_int(d$n)
  print(d[d$n != '0', , drop = FALSE], row.names = FALSE, right = TRUE)

  if (x$n_zero_pairs > 0) {
    cat(sprintf("\n  %s pair%s with zero dissimilarity in the full data. %s kept.\n",
                fmt_int(x$n_zero_pairs), if (x$n_zero_pairs == 1) '' else 's',
                fmt_int(sum(p$diss == 0))))
  }

  # Graph
  cli::cli_h2("Graph")
  if (min(x$degree) == max(x$degree)) {
    show("Pairs per sample (degree):", sprintf("%d for all samples", min(x$degree)))
  } else {
    show("Pairs per sample (degree):", sprintf("min %d, median %s, max %d",
                                               min(x$degree), format(stats::median(x$degree)), max(x$degree)))
  }
  show("Connected components:", x$n_components)
  if (x$n_components > 1) {
    cli::cli_alert_warning("The graph is not connected: the relative position of components is not constrained. Consider a larger {.arg n_random}.")
  }

  if (!is.null(x$pcoa)) cat("\nUse plot() to see the selected pairs on a PCoA.\n")
  invisible(x)
}


#' Plot method for sparse_diss objects
#'
#' @description
#' \code{\link[graphics]{plot}} method for objects of class \code{sparse_diss}.
#' Shows the samples on two axes of a PCoA of the full dissimilarity
#' matrix, with the selected pairs drawn as segments between samples.
#'
#' @param x An object of class \code{sparse_diss}, as returned by
#'   \code{\link{sparse_diss}}.
#' @param y Currently unused.
#' @param k Integer vector of length 2 indicating which PCoA axes to plot.
#'   Must be within the \code{k_pcoa} axes stored in \code{x}. Default is
#'   \code{c(1, 2)}.
#' @param edge_types Optional character vector with the type of pairs to draw
#'   (any of \code{'knn'}, \code{'random'}, \code{'zero'}).
#'   Default is \code{NULL} (all).
#' @param edge_col_by Character string: colour segments by pair \code{'type'}
#'   (default) or by dissimilarity (\code{'diss'}).
#' @param edge_alpha Numeric transparency of segments. Default is 0.4.
#' @param edge_linewidth Numeric width of segments. Default is 0.3.
#' @param add_points Logical indicating whether to add the default sample
#'   points layer via \code{ggplot2::geom_point()}. Default is \code{TRUE}.
#'   Set to \code{FALSE} to omit it and add a custom one instead.
#' @param point_size Numeric size of points. Default is 1.5.
#' @param ... Currently unused.
#'
#' @return
#' A \code{ggplot} object, so that it can be further customised (e.g.
#' \code{plot(sp) + ggplot2::theme_bw()}). The default data of the plot
#' contains the PCoA coordinates of each sample, in columns \code{PCoA.1} and
#' \code{PCoA.2} (the plotted axes, in order).
#'
#' @details
#' The PCoA is computed on the full dissimilarity matrix, on the link scale (see
#' \code{k_pcoa} and \code{link} in \code{\link{sparse_diss}}).
#'
#' @method plot sparse_diss
#' @export
plot.sparse_diss <- function(x,
                             y = NULL,
                             k = c(1, 2),
                             edge_types = NULL,
                             edge_col_by = 'type',
                             edge_alpha = 0.4,
                             edge_linewidth = 0.3,
                             add_points = TRUE,
                             point_size = 1.5,
                             ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort(c(
      "x" = "plot method for sparse_diss objects requires the {.pkg ggplot2} package, which is not installed.",
      "i" = "You can install it via: {.code install.packages(\"ggplot2\")}."
    ))
  }

  # checks
  if (!inherits(x, 'sparse_diss')) {
    cli::cli_abort("{.arg x} must be a {.cls sparse_diss} object.")
  }
  if (is.null(x$pcoa)) {
    cli::cli_abort(c("x" = "This object does not contain a PCoA.",
                     "i" = "Re-run {.fn sparse_diss} with {.arg k_pcoa} set to an integer of at least 2."))
  }
  if (!is.numeric(k) || length(k) != 2 || anyNA(k) || any(k != as.integer(k))) {
    cli::cli_abort("{.arg k} must be an integer vector of length 2 (e.g. {.code c(1, 2)}).")
  }

  k <- as.integer(k)

  if (k[1] == k[2]) {
    cli::cli_abort("{.arg k} must reference two distinct PCoA axes.")
  }
  if (any(k < 1) || any(k > x$pcoa$k)) {
    cli::cli_abort(c("{.arg k} must reference axes between 1 and {x$pcoa$k}.",
                     "i" = "Only {x$pcoa$k} PCoA axes were stored: increase {.arg k_pcoa} in {.fn sparse_diss}."))
  }

  all_types <- c('random', 'zero', 'knn')  # drawing order (last on top)
  if (!is.null(edge_types)) {
    if (!is.character(edge_types) || !all(edge_types %in% all_types)) {
      cli::cli_abort("{.arg edge_types} must be {.code NULL} or a subset of {.val {all_types}}.")
    }
  } else {
    edge_types <- all_types
  }
  edge_col_by <- match.arg(edge_col_by, c('type', 'diss'))
  if (!is.logical(add_points) || length(add_points) != 1 || is.na(add_points)) {
    cli::cli_abort("{.arg add_points} must be a single {.val TRUE}/{.val FALSE} value.")
  }
  if (!is.numeric(edge_alpha) || !is.numeric(edge_linewidth) || !is.numeric(point_size)) {
    cli::cli_abort("{.arg edge_alpha}, {.arg edge_linewidth} and {.arg point_size} must be numeric.")
  }

  var_exp <- x$pcoa$var_exp

  coords <- x$pcoa$coords
  coords$PCoA.1 <- coords[[paste0('PCoA', k[1])]]
  coords$PCoA.2 <- coords[[paste0('PCoA', k[2])]]

  #  Segments (edges)
  seg <- x$pairs[x$pairs$type %in% edge_types, , drop = FALSE]
  seg$type <- factor(seg$type, levels = all_types)
  seg <- seg[order(seg$type), , drop = FALSE]
  seg$x1 <- coords$PCoA.1[seg$s1]
  seg$y1 <- coords$PCoA.2[seg$s1]
  seg$x2 <- coords$PCoA.1[seg$s2]
  seg$y2 <- coords$PCoA.2[seg$s2]

  #  Base plot
  plot_obj <- ggplot2::ggplot(coords, ggplot2::aes(x = .data$PCoA.1, y = .data$PCoA.2)) +
    ggplot2::coord_equal() +
    ggplot2::xlab(sprintf('PCoA %d (%.1f%%)', k[1], 100 * var_exp[k[1]])) +
    ggplot2::ylab(sprintf('PCoA %d (%.1f%%)', k[2], 100 * var_exp[k[2]]))

  if (nrow(seg) > 0) {
    if (edge_col_by == 'type') {
      plot_obj <- plot_obj +
        ggplot2::geom_segment(data = seg,
                              ggplot2::aes(x = .data$x1, y = .data$y1,
                                           xend = .data$x2, yend = .data$y2,
                                           colour = .data$type),
                              alpha = edge_alpha,
                              linewidth = edge_linewidth)
    } else {
      plot_obj <- plot_obj +
        ggplot2::geom_segment(data = seg,
                              ggplot2::aes(x = .data$x1, y = .data$y1,
                                           xend = .data$x2, yend = .data$y2,
                                           colour = .data$diss),
                              alpha = edge_alpha,
                              linewidth = edge_linewidth)
    }
  } else {
    cli::cli_warn("No pairs of the requested {.arg edge_types} were found in the object.")
  }

  if (add_points) {
    plot_obj <- plot_obj + ggplot2::geom_point(size = point_size)
  }

  return(plot_obj)
}
