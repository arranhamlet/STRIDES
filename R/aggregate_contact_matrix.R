#' Aggregate a Square Contact Matrix Using Population Weights
#'
#' Aggregates a contact matrix to custom age groups using the mean population over time
#' as weights to preserve transmission dynamics.
#'
#' @param mat Square matrix of contact rates (e.g., 101 x 101 for ages 0–100).
#' @param age_breaks Vector of age group cutpoints (e.g. c(0, 5, 10, ..., Inf)).
#' @param population Long-format data.frame with `time`, `age`, and `value` columns.
#' @param symmetric Logical; if TRUE, enforces reciprocity between age groups.
#'
#' @return Aggregated square contact matrix (n_group x n_group).
aggregate_contact_matrix <- function(mat, age_breaks, population, symmetric = TRUE) {

  # Step 1: Compute average population per age over time
  pop_by_age <- population %>%
    dplyr::group_by(age) %>%
    dplyr::summarise(weight = mean(value, na.rm = TRUE), .groups = "drop")

  # Step 2: Map ages to age groups
  if (!is.infinite(tail(age_breaks, 1))) {
    age_breaks[length(age_breaks)] <- Inf
  }
  pop_by_age$group <- cut(pop_by_age$age - 1, breaks = age_breaks, labels = FALSE, right = FALSE)

  # Step 3: Aggregate contact matrix using weighted means
  age_index <- cut(seq_len(nrow(mat)) - 1, breaks = age_breaks, labels = FALSE, right = FALSE)
  n_groups <- max(age_index)
  agg_mat <- matrix(0, nrow = n_groups, ncol = n_groups)

  for (i in seq_len(n_groups)) {
    for (j in seq_len(n_groups)) {
      rows <- which(age_index == i)
      cols <- which(age_index == j)

      pop_weights <- pop_by_age$weight[match(rows, pop_by_age$age)]
      pop_weights <- pop_weights / sum(pop_weights, na.rm = TRUE)
      sub_mat <- mat[rows, cols, drop = FALSE]

      # Weighted average over rows
      agg_mat[i, j] <- sum(pop_weights %*% sub_mat, na.rm = TRUE)
    }
  }

  # Step 4: Optional reciprocity enforcement
  if (symmetric) {
    group_weights <- pop_by_age %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(pop = sum(weight, na.rm = TRUE), .groups = "drop") %>%
      dplyr::pull(pop)

    for (i in seq_len(n_groups)) {
      for (j in seq_len(n_groups)) {
        agg_mat[i, j] <- (agg_mat[i, j] * group_weights[j] + agg_mat[j, i] * group_weights[i]) /
          (group_weights[i] + group_weights[j])
      }
    }
  }

  return(agg_mat)
}
