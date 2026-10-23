#' Format community data (internal)
#' @description
#' Formats a site x species table into a long data frame of community
#' dissimilarities between pairs of sites, in the order of \code{combn(n, 2)}
#' (equal to \code{as.vector(dist())}).
#' @param Y_mat community matrix (sites x species)
#' @param id Optional: Unique name for each location
#' @param method Method used to calculate community dissimilarity between pairs of sites. A dissimilarity index supported by \code{vegan::vegdist}
#' @param num_den For some dissimilarity indices, whether function should return the numerator and denominator in two different columns
#' @param trans Transformation to be applied to raw data before calculating distance matrix. A function applied to the community matrix, e.g. \code{sqrt}. Use \code{binary = TRUE} for presence-absence data.
#' @param drop_empty_rows Should empty rows (i.e., sites containing 0s only) be removed?
#' @param drop_empty_cols Should empty columns (i.e., species without observations) be removed?
#' @param binary  Should abundance data be transform into presence absence?
#' @param check_trans Check for non-finite values in transformation
#' @param na.rm Remove rows with NA?
#' @param na.replace Replace NAs with specific value
#' @return A data frame with the following columns
#' \itemize{
#'   \item s1, s2  -  A combination of sites (rows).
#'   \item Additional columns indicating the dissimilarity between sites. If `num_den = TRUE`, the data frame will contain these individual elements as different columns.
#' }
#' @importFrom vegan vegdist
#' @noRd
make_Y_diss <- function(Y_mat = NULL,
                        id = NULL,
                        method = 'bray',
                        num_den = FALSE,
                        trans = function(x) x,
                        drop_empty_rows = FALSE,
                        drop_empty_cols = TRUE,
                        binary = FALSE,
                        check_trans = TRUE,
                        na.rm = FALSE,
                        na.replace = NULL) {
  # BASIC CHECKS
  if (!(isTRUE(binary) || isFALSE(binary))) {
    cli::cli_abort(c(
      "x" = "{.arg binary} must be a single logical value ({.code TRUE} or {.code FALSE}).",
      "i" = "You provided an object of class {.cls {class(binary)}}."
    ))
  }
  if (is.null(trans)) {
    trans <- function(x) x
  }
  if (!is.function(trans)) {
    cli::cli_abort(c(
      "x" = "{.arg trans} must be a {.cls function} (e.g., {.fn sqrt} or {.code function(x) x}).",
      "i" = "You provided an object of class {.cls {class(trans)}}."
    ))
  }
  if (is.null(Y_mat)) {
    cli::cli_abort(c("{.var Y_mat} must be supplied",
                     "x" = "You have not provided a community matrix"))
  }
  if (!(is.matrix(Y_mat) || is.data.frame(Y_mat))) {
    cli::cli_abort(c("{.var Y_mat} must be a site x species matrix or data frame",
                     "x" = "You have supplied a {.cls {class(Y_mat)}}"))
  }

  Y_mat <- as.matrix(Y_mat)
  storage.mode(Y_mat) <- "numeric"

  if (!is.null(id) && length(id) != nrow(Y_mat)) {
    cli::cli_abort("{.arg id} must have one value per row of {.arg Y_mat} ({nrow(Y_mat)}), not {length(id)}.")
  }

  # ---- NA handling ----
  if (!is.null(na.replace)) {
    n_na <- sum(is.na(Y_mat))
    if (n_na > 0) {
      Y_mat[is.na(Y_mat)] <- na.replace
      cli::cli_warn(c("Values replaced!", "i" = "{n_na} NA{?s} replaced with {na.replace}"))
    }
  } else if (na.rm) {
    NA_row <- rowSums(is.na(Y_mat)) > 0
    if (any(NA_row)) {
      if (!is.null(id)) id <- id[!NA_row]
      Y_mat <- Y_mat[!NA_row, , drop = FALSE]
      cli::cli_warn(c("Rows removed!", "i" = "{sum(NA_row)} row{?s} containing NAs removed"))
    }
  }

  # drop rows/cols
  if (drop_empty_rows) {
    row_drop <- rowSums(Y_mat, na.rm = TRUE) == 0
    if (any(row_drop)) {
      site_char <- if (is.null(id)) paste(which(row_drop), collapse = ",") else paste(id[row_drop], collapse = ",")
      Y_mat <- Y_mat[!row_drop, , drop = FALSE]
      if (!is.null(id)) id <- id[!row_drop]
      cli::cli_warn(c("Sites removed!", "i" = "Row{?s}: {site_char}"))
    }
  }

  if (drop_empty_cols) {
    col_drop <- colSums(Y_mat, na.rm = TRUE) == 0
    if (any(col_drop)) {
      nms <- colnames(Y_mat)
      sp_char <- if (is.null(nms)) paste(which(col_drop), collapse = ",") else paste(nms[col_drop], collapse = ",")
      Y_mat <- Y_mat[, !col_drop, drop = FALSE]
      cli::cli_warn(c("Species removed!", "i" = "Column{?s}: {sp_char}"))
    }
  }

  # transformation
  Y_mat <- tryCatch({ trans(Y_mat) }, error = function(e) cli::cli_abort("Transformation failed: {e$message}"))

  if (check_trans && any(!is.finite(Y_mat))) {
    cli::cli_warn(c('Non-finite values created!', 'i' = 'Transformation produced {sum(!is.finite(Y_mat))} non-finite value{?s}'))
  }

  if (method == 'sorensen') binary <- TRUE
  if (binary) Y_mat <- (Y_mat > 0) * 1

  # Compute dissimilarity
  n <- nrow(Y_mat)
  if (n < 2) cli::cli_abort("Less than 2 sites remaining after cleaning. Cannot calculate distance.")

  # All pairs of sites (s1 < s2), in the order of combn(n, 2) and as.vector(dist())
  D <- cbind(s1 = rep.int(seq_len(n - 1L), (n - 1L):1L),
             s2 = unlist(lapply(2:n, function(i) i:n)))

  if (num_den) {
    if (method %in% c('sorensen', 'jaccard')) {
      # Shared species (a) and species found only in the first (b) or second (c) site
      Y_pa <- (Y_mat > 0) * 1
      n_sp <- rowSums(Y_pa)
      a <- tcrossprod(Y_pa)[D]
      b <- n_sp[D[, 1]] - a
      c <- n_sp[D[, 2]] - a

      if (method == 'sorensen') {
        dist_data <- data.frame(D, num_sor = (b + c), den_sor = (2*a + b + c))
      } else {
        dist_data <- data.frame(D, num_jac = (b + c), den_jac = (a + b + c))
      }
    } else if (method == 'bray') {
      if (anyNA(Y_mat)) {
        # pairwise deletion of NAs
        vals <- t(apply(D, 1, function(x) bray_num_den(Y_mat[x[1],], Y_mat[x[2],])))
        num <- vals[, 1]; den <- vals[, 2]
      } else {
        # num = sum|x - y|, den = sum(x + y)
        num <- as.vector(stats::dist(Y_mat, method = "manhattan"))
        n_ab <- rowSums(Y_mat)
        den <- n_ab[D[, 1]] + n_ab[D[, 2]]
      }
      dist_data <- data.frame(D, num_bra = num, den_bra = den)
    } else {
      cli::cli_abort("{.var num_den} only supported for 'bray', 'sorensen' and 'jaccard'.")
    }
  } else {
    meth_val <- if (method == 'sorensen') 'bray' else method
    dist_data <- data.frame(D, diss = as.vector(vegan::vegdist(Y_mat, method = meth_val)))
  }

  if (!is.null(id)) {
    dist_data$s1 <- id[dist_data$s1]
    dist_data$s2 <- id[dist_data$s2]
  }

  dist_data
}
