#' Aggregate Age Structure of a Model Output Object
#'
#' Aggregates a long-format data frame with an age dimension (`dim1`) into custom age groups
#' based on specified `age_breaks`. Supports sum, mean, or population-weighted.mean aggregation.
#'
#' @param obj A data frame with at least the columns `dim1` and `value`. If using time-varying weights, must also include a time column (default: `dim4`).
#' @param age_breaks A numeric vector of breakpoints defining age bands (e.g. c(0, 1, 5, 10, 20, 30, 40, 50, 60, 70, 80, Inf)).
#' @param method One of `"sum"`, `"mean"`, or `"weighted.mean"`. Controls aggregation method.
#' @param weights Optional. Either a numeric vector of length 101 (ages 0–100) or a matrix/data.frame of dimension (time × age).
#' @param time_var Name of the column in `obj` representing time (used only if `weights` is time-varying). Default is `"dim4"`.
#'
#' @return A data frame with aggregated `dim1` values and other dimensions preserved.
#'
#' @importFrom dplyr group_by summarise mutate filter across
#' @importFrom data.table first
#' @export
aggregate_age_structure <- function(obj,
                                    age_breaks,
                                    method = c("sum", "mean", "weighted.mean"),
                                    weights = NULL,
                                    time_var = "dim4") {

  method <- match.arg(method, c("sum", "mean", "weighted.mean"))
  full_ages <- 1:101
  age_group_index <- cut(full_ages, breaks = age_breaks + 1, right = FALSE, labels = FALSE)
  n_age_new <- length(age_breaks) - 1

  # Assign age group to each row based on dim1
  obj$age_group <- age_group_index[obj$dim1]

  if (method == "sum") {
    out <- obj %>%
      dplyr::group_by(dplyr::across(-c(value, dim1))) %>%
      dplyr::summarise(
        dim1 = data.table::first(age_group),
        value = sum(value, na.rm = TRUE),
        .groups = "drop"
      )

  } else if (method == "mean") {
    out <- obj %>%
      dplyr::group_by(dplyr::across(-c(value, dim1))) %>%
      dplyr::summarise(
        dim1 = data.table::first(age_group),
        value = mean(value, na.rm = TRUE),
        .groups = "drop"
      )

  } else if (method == "weighted.mean") {
    if (is.null(weights)) stop("Weights must be provided for method = 'weighted.mean'")

    if (is.matrix(weights) || is.data.frame(weights)) {
      if (!(time_var %in% names(obj))) stop("Time variable not found in data and required for time-varying weights")

      obj$weight <- mapply(function(age, time) {
        weights[time, age]  # age 0 maps to column 1
      }, obj$dim1, obj[[time_var]])

    } else {
      # Static vector weights
      weight_map <- setNames(weights, as.character(0:100))
      obj$weight <- weight_map[as.character(obj$dim1)]
    }

    # Group, normalize weights, and aggregate
    out <- obj %>%
      dplyr::group_by(dplyr::across(-c(value, weight, dim1))) %>%
      dplyr::mutate(norm_weight = weight / sum(weight, na.rm = TRUE)) %>%
      dplyr::summarise(
        dim1 = data.table::first(age_group),
        value = sum(value * norm_weight, na.rm = TRUE),
        .groups = "drop"
      )
  }

  out$dim1 <- as.integer(out$dim1)
  out <- dplyr::filter(out, dim1 <= n_age_new)

  return(out)
}
