#' Aggregate Age Structure of a Model Output Object
#'
#' Aggregates a long-format data frame with an age dimension (`dim1`) into custom age groups
#' based on specified `age_breaks`. Supports sum, mean, or population-weighted.mean aggregation.
#'
#' @param obj A data frame with at least the columns `dim1` and `value`. Typically model output like vaccination coverage or seeding.
#' @param age_breaks A numeric vector of breakpoints defining age bands (e.g. c(0, 1, 5, 10, 20, 30, 40, 50, 60, 70, 80, Inf)).
#' @param method Character. One of `"sum"`, `"mean"`, or `"weighted.mean"`. Controls how values are aggregated within each age band.
#' @param weights Optional numeric vector of weights (e.g., population) for `method = "weighted.mean"`. Must match the full set of model ages (e.g., 0:100 or 1:101).
#'
#' @return A data frame with aggregated `dim1` values matching the number of age groups defined by `age_breaks`.
#' @export
aggregate_age_structure <- function(obj, age_breaks, method = c("sum", "mean", "weighted.mean"), weights = NULL) {

  method <- match.arg(method, c("sum", "mean", "weighted.mean"))

  # Identify all age indices present in the model
  all_model_ages <- sort(unique(obj$dim1))
  n_age_model <- max(all_model_ages)

  # Create full age vector assuming dim1 starts at 1 and represents ages 0:(n-1)
  full_ages <- 0:(max(n_age_model, length(weights)) - 1)

  # Determine age group index for each possible model age
  age_group_index <- cut(full_ages, breaks = age_breaks, right = FALSE, labels = FALSE)
  n_age_new <- length(age_breaks) - 1

  # Map only those dim1 values present in the object
  obj_ages <- obj$dim1
  obj$age_group <- age_group_index[obj_ages]

  if (method == "sum") {
    out <- obj %>%
      dplyr::group_by(across(-c(value, dim1))) %>%
      dplyr::summarise(
        dim1 = unique(age_group),
        value = sum(value, na.rm = TRUE),
        .groups = "drop"
      )

  } else if (method == "mean") {
    out <- obj %>%
      dplyr::group_by(across(-c(value, dim1))) %>%
      dplyr::summarise(
        dim1 = unique(age_group),
        value = mean(value, na.rm = TRUE),
        .groups = "drop"
      )

  } else if (method == "weighted.mean") {
    if (is.null(weights)) stop("Weights must be provided for method = 'weighted.mean'")
    if (length(weights) != length(full_ages)) stop("Length of weights must match full age vector (e.g., 0:100)")

    # Create named weight lookup
    weight_map <- setNames(weights, as.character(seq_along(full_ages)))

    # Match weights to only the dim1 values present
    obj$weight <- weight_map[as.character(obj$dim1)]

    out <- obj %>%
      dplyr::group_by(across(-c(value, weight, dim1))) %>%
      dplyr::summarise(
        dim1 = unique(age_group),
        value = weighted.mean(value, weight, na.rm = TRUE),
        .groups = "drop"
      )
  }

  # Final formatting
  out$dim1 <- as.integer(out$dim1)
  out <- dplyr::filter(out, dim1 <= n_age_new)

  return(out)
}
