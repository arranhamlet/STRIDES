#' Aggregate Age Structure of a Model Output Object
#'
#' Aggregates a long-format data frame with an age dimension (dim1) into custom age groups.
#' Supports sum, mean, weighted.mean, and rate-based aggregation with optional plotting check.
#'
#' @param obj A data frame with at least the columns dim1 and value. If using time-varying weights, must also include a time column (default: dim4).
#' @param age_breaks A numeric vector of breakpoints defining age bands (e.g. c(0, 1, 5, ..., Inf)).
#' @param method One of "sum", "mean", "weighted.mean", or "rate". "rate" returns population-weighted rate while preserving total events.
#' @param weights A numeric vector or data frame of population weights, required for "rate" and "weighted.mean".
#' @param time_var The name of the time column. Default is "dim4".
#' @param tol Numerical tolerance for preservation check. Default is 1e-6.
#'
#' @return A data frame with aggregated dim1 values. If method = "rate", a plot is attached as an attribute.
#'
#' @importFrom dplyr group_by summarise mutate left_join rename select across
#' @importFrom data.table first
#' @importFrom ggplot2 ggplot aes geom_line theme_minimal labs
#' @keywords internal
aggregate_age_structure <- function(obj,
                                    age_breaks,
                                    method = c("sum", "mean", "weighted.mean", "rate"),
                                    weights = NULL,
                                    time_var = "dim4",
                                    tol = 1e-6) {

  original <- obj

  #Set the final one as inf to catch anything above it
  age_breaks[length(age_breaks)] <- Inf

  method <- match.arg(method, c("sum", "mean", "weighted.mean", "rate"))
  full_ages <- 1:101
  age_group_index <- cut(full_ages, breaks = age_breaks + 1, right = FALSE, labels = FALSE)
  obj$age_group <- age_group_index[obj$dim1]

  if (method == "sum") {
    out <- obj %>%
      dplyr::group_by(dplyr::across(-c(value, dim1))) %>%
      dplyr::summarise(
        dim1 = data.table::first(age_group),
        value = sum(value, na.rm = TRUE),
        .groups = "keep"
      )

  } else if (method == "mean") {
    out <- obj %>%
      dplyr::group_by(dplyr::across(-c(value, dim1))) %>%
      dplyr::summarise(
        dim1 = data.table::first(age_group),
        value = mean(value, na.rm = TRUE),
        .groups = "drop"
      )

  } else if (method %in% c("weighted.mean", "rate")) {
    if (is.null(weights)) stop("Weights must be provided for method = 'weighted.mean' or 'rate'")

    # Attach population
    if (is.data.frame(weights) || is.matrix(weights)) {
      weights <- as.data.frame(weights)
      if (is.null(time_var) || !(time_var %in% names(obj))) {
        static_weights <- weights %>%
          dplyr::group_by(age) %>%
          dplyr::summarise(population = mean(value, na.rm = TRUE), .groups = "drop")
        obj <- obj %>%
          dplyr::left_join(static_weights, by = c("dim1" = "age"))
      } else {
        obj <- obj %>%
          dplyr::left_join(
            weights %>%
              dplyr::rename(population = value) %>%
              dplyr::select(age, time = !!sym("time"), population),
            by = setNames(c("age", "time"), c("dim1", time_var))
          )
      }
    } else {
      weight_map <- setNames(weights, as.character(0:100))
      obj$population <- weight_map[as.character(obj$dim1 - 1)]
    }

    group_vars <- intersect(names(obj), c("dim2", "dim3", time_var, "age_group"))

    if (method == "rate") {
      obj <- obj %>%
        dplyr::mutate(events = value * population)

      out <- obj %>%
        dplyr::group_by(dplyr::across(all_of(group_vars))) %>%
        dplyr::summarise(
          dim1 = first(age_group),
          population = sum(population, na.rm = TRUE),
          events = sum(events, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::mutate(value = events / population)

      # Check preservation
      total_events_input <- sum(obj$events, na.rm = TRUE)
      total_events_output <- sum(out$value * out$population, na.rm = TRUE)
      if (abs(total_events_input - total_events_output) > tol) {
        warning(sprintf(
          "Total events not preserved: input = %.6f, output = %.6f (diff = %.6f)",
          total_events_input, total_events_output,
          abs(total_events_input - total_events_output)
        ))
      }

      # Create comparison plot
      time_col <- time_var
      original_check <- obj %>%
        dplyr::group_by(!!sym(time_col)) %>%
        dplyr::summarise(events = sum(events, na.rm = TRUE), .groups = "drop") %>%
        dplyr::mutate(source = "original")

      aggregated_check <- out %>%
        dplyr::mutate(events = value * population) %>%
        dplyr::group_by(!!sym(time_col)) %>%
        dplyr::summarise(events = sum(events, na.rm = TRUE), .groups = "drop") %>%
        dplyr::mutate(source = "aggregated")

      plot_obj <- dplyr::bind_rows(original_check, aggregated_check) %>%
        ggplot2::ggplot(ggplot2::aes(x = !!sym(time_col), y = events, color = source,
                                     linetype = source), alpha = 0.5) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = "Event Preservation Check (Rate Method)",
          x = time_col, y = "Total Events", color = "Source", linetype = "Source"
        )

      print(plot_obj)

    } else {
      # Weighted mean (simple proportional)
      out <- obj %>%
        dplyr::group_by(dplyr::across(all_of(group_vars))) %>%
        dplyr::mutate(prop = population / sum(population, na.rm = TRUE)) %>%
        dplyr::mutate(value = value * prop) %>%
        dplyr::summarise(
          dim1 = first(age_group),
          value = sum(value, na.rm = TRUE),
          .groups = "drop"
        )
    }
  }

  out$dim1 <- as.integer(out$dim1)
  return(out)
}
