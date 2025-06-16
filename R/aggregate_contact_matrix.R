#' Aggregate a Square Contact Matrix to Custom Age Groups
#'
#' Aggregates a contact matrix (e.g., 101x101 for ages 0–100) into a smaller square matrix
#' according to custom `age_breaks`. Optionally applies population weighting.
#'
#' @param mat A square numeric matrix where both rows and columns correspond to single-year ages (e.g., 0:100).
#' @param age_breaks A numeric vector of age cut points (e.g. `c(0, 1, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, Inf)`).
#' @param weights Optional numeric vector of length equal to `nrow(mat)` (e.g., population by age).
#'   Required if `method = "weighted.mean"`.
#' @param method Aggregation method. One of `"weighted.mean"` (default) or `"sum"`.
#'
#' @return A square matrix of size `(length(age_breaks) - 1)` aggregated by age group.
#' @export
aggregate_contact_matrix <- function(mat, age_breaks, weights = NULL, method = c("weighted.mean", "sum")) {

  method <- match.arg(method, c("weighted.mean", "sum"))

  if (nrow(mat) != ncol(mat)) stop("Matrix must be square")

  n <- nrow(mat)
  full_ages <- 0:(n - 1)
  age_groups <- cut(full_ages, breaks = age_breaks, right = FALSE, labels = FALSE)
  n_age_new <- length(age_breaks) - 1

  groupings <- split(seq_len(n), age_groups)
  agg_mat <- matrix(NA_real_, nrow = n_age_new, ncol = n_age_new)

  for (i in seq_len(n_age_new)) {
    for (j in seq_len(n_age_new)) {
      rows <- groupings[[i]]
      cols <- groupings[[j]]
      submat <- mat[rows, cols, drop = FALSE]

      if (method == "weighted.mean") {
        if (is.null(weights)) stop("Weights must be provided for method = 'weighted.mean'")
        w_i <- weights[rows] / sum(weights[rows])
        w_j <- weights[cols] / sum(weights[cols])
        agg_mat[i, j] <- sum(outer(w_i, w_j) * submat, na.rm = TRUE)

      } else if (method == "sum") {
        agg_mat[i, j] <- sum(submat, na.rm = TRUE)
      }
    }
  }

  rownames(agg_mat) <- colnames(agg_mat) <- paste0(
    head(age_breaks, -1), "-", head(age_breaks[-1] - 1, -1)
  )
  rownames(agg_mat)[n_age_new] <- colnames(agg_mat)[n_age_new] <- paste0(age_breaks[n_age_new], "+")
  return(agg_mat)
}
