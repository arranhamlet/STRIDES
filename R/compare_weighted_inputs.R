#' Compare Two Inputs Using Population Weights or Sum
#'
#' Computes population-weighted averages or total sums for two input objects and compares them.
#'
#' @param default_input A data frame or matrix containing the baseline input values.
#' @param default_weight A data frame or matrix of weights (e.g., N0) corresponding to `default_input`.
#' @param adapted_input A data frame or matrix containing the transformed input values (e.g., aggregated).
#' @param adapted_weight A data frame or matrix of weights corresponding to `adapted_input`.
#' @param method Either "weighted" (default) or "sum".
#'
#' @return A list containing summary statistics for each input comparison.
#' @export
compare_weighted_inputs <- function(default_input,
                                    default_weight,
                                    adapted_input,
                                    adapted_weight,
                                    method = "weighted") {
  # Check type consistency
  if (!inherits(default_input, "data.frame") || !inherits(adapted_input, "data.frame")) {
    stop("Inputs must be data frames.")
  }
  if (!inherits(default_weight, "data.frame") || !inherits(adapted_weight, "data.frame")) {
    stop("Weights must be data frames.")
  }

  if (method == "weighted") {

    # Merge value and weight
    default <- dplyr::left_join(default_input, default_weight %>%
                                  dplyr::group_by(dim1) %>%
                                  dplyr::mutate(weight = value / sum(value)) %>%
                                  dplyr::select(dim1, weight),
                                by = "dim1")
    adapted <- dplyr::left_join(adapted_input, adapted_weight %>%
                                  dplyr::group_by(dim1) %>%
                                  dplyr::mutate(weight = value / sum(value)) %>%
                                  dplyr::select(dim1, weight),
                                by = "dim1")

    wm_default <- sum(default$value * default$weight, na.rm = TRUE) / sum(default$weight, na.rm = TRUE)
    wm_adapted <- sum(adapted$value * adapted$weight, na.rm = TRUE) / sum(adapted$weight, na.rm = TRUE)

  } else if (method == "sum") {
    wm_default <- sum(default_input$value, na.rm = TRUE)
    wm_adapted <- sum(adapted_input$value, na.rm = TRUE)

  } else {
    stop("Invalid method: must be 'weighted' or 'sum'.")
  }

  abs_diff <- abs(wm_default - wm_adapted)
  rel_diff <- abs_diff / (abs(wm_default) + 1e-9)

  return(list(
    method = method,
    default = wm_default,
    adapted = wm_adapted,
    absolute_difference = abs_diff,
    relative_difference = rel_diff
  ))
}
