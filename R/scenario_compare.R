#' Compare Two Transmission Scenarios
#'
#' This function compares two transmission scenarios by:
#' 1. Combining case output data and generating comparative line and ribbon plots.
#' 2. Producing summary statistics including peak and total cases.
#' 3. Generating cumulative case trajectories.
#' 4. Creating boxplots of peak and total case distributions (requires simulated `run` column).
#' 5. Comparing susceptible proportions over time between the scenarios.
#'
#' @param scenario_1 A named list of outputs from scenario 1, including:
#'   - `select_state_plot$data`: Time-series summary data (`data.table`)
#'   - `box_data`: Peak/total case data by `run`
#'   - `susceptibility_data`: Long-format susceptibility proportions by time
#' @param scenario_2 Same structure as `scenario_1`, for the second scenario.
#'
#' @return Invisibly returns a named list of plots and summary statistics:
#' \item{case_plot}{Line and ribbon plot comparing median cases over time.}
#' \item{cumulative_plot}{Cumulative case trajectories.}
#' \item{boxplot_peak}{Boxplot comparing peak cases.}
#' \item{boxplot_total}{Boxplot comparing total cases.}
#' \item{susceptibility_plot}{Line plot comparing susceptible proportions over time.}
#' \item{summary_table}{A `gt` object showing formatted peak/total case summaries.}
#'
#' @importFrom data.table rbindlist data.table :=
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_boxplot labs
#'   scale_y_continuous scale_x_continuous theme_bw theme facet_wrap
#' @importFrom gt gt fmt_number cols_label tab_options tab_style
#'   cell_text cells_column_labels opt_table_font google_font opt_row_striping
#' @importFrom dplyr mutate select
#' @importFrom tidyr pivot_longer
#' @importFrom scales comma percent_format pretty_breaks
#' @export
scenario_compare <- function(scenario_1, scenario_2) {

  # --- Combine time series case output (ribbons + lines) ---
  case_combo <- data.table::rbindlist(list(
    dplyr::mutate(
      subset(scenario_1$select_state_plot$data, state_plot == "Cases"),
      type = 1
    ),
    dplyr::mutate(
      subset(scenario_2$select_state_plot$data, state_plot == "Cases"),
      type = 2
    )
  ))

  # Filter for time range where cases are non-zero
  case_time_range <- case_combo[median != 0, range(time)]

  subset_case_plot <- ggplot2::ggplot(
    case_combo[time %in% case_time_range[1]:case_time_range[2]],
    ggplot2::aes(x = time, y = median, ymin = lower, ymax = upper,
                 color = factor(type), fill = factor(type))
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(alpha = 0.25, show.legend = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Time", y = "Cases", color = "Scenario") +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks())

  # --- Summary statistics ---
  case_summary <- case_combo[
    , .(
      peak_median_cases = max(median, na.rm = TRUE),
      time_of_peak       = time[which.max(median)],
      total_median_cases = sum(median, na.rm = TRUE),
      total_lower_cases  = sum(lower, na.rm = TRUE),
      total_upper_cases  = sum(upper, na.rm = TRUE)
    ),
    by = type
  ]

  diff_row <- data.table::data.table(
    type               = "Difference",
    peak_median_cases = diff(case_summary$peak_median_cases),
    time_of_peak       = diff(case_summary$time_of_peak),
    total_median_cases = diff(case_summary$total_median_cases),
    total_lower_cases  = diff(case_summary$total_lower_cases),
    total_upper_cases  = diff(case_summary$total_upper_cases)
  )

  case_summary_full <- data.table::rbindlist(list(case_summary, diff_row), fill = TRUE)

  case_summary_format <- gt::gt(case_summary_full) %>%
    gt::fmt_number(
      columns = where(is.numeric),
      decimals = 0,
      use_seps = TRUE
    ) %>%
    gt::cols_label(
      type = "Scenario",
      peak_median_cases   = "Peak Median",
      time_of_peak        = "Time of Peak",
      total_median_cases  = "Total Median",
      total_lower_cases   = "Total Lower",
      total_upper_cases   = "Total Upper"
    ) %>%
    gt::tab_options(
      table.font.size = "small",
      data_row.padding = gt::px(3)
    ) %>%
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_column_labels(gt::everything())
    ) %>%
    gt::opt_table_font(font = gt::google_font("Roboto")) %>%
    gt::opt_row_striping()

  # --- Cumulative cases ---
  cumulative_cases <- case_combo[
    order(time),
    .(
      time       = unique(time),
      cum_median = cumsum(median),
      cum_lower  = cumsum(lower),
      cum_upper  = cumsum(upper)
    ),
    by = type
  ]

  cumulative_case_plot <- ggplot2::ggplot(
    cumulative_cases[time %in% case_time_range[1]:case_time_range[2]],
    ggplot2::aes(x = time, y = cum_median, ymin = cum_lower, ymax = cum_upper,
                 color = factor(type), fill = factor(type))
  ) +
    ggplot2::geom_ribbon(alpha = 0.3) +
    ggplot2::geom_line(size = 0.9) +
    ggplot2::labs(x = "Time", y = "Cumulative Cases", color = "Scenario", fill = "Scenario") +
    ggplot2::theme_bw() +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks())

  # --- Boxplots for peak and total cases ---
  case_box_data <- data.table::rbindlist(list(
    dplyr::mutate(scenario_1$box_data, type = "1"),
    dplyr::mutate(scenario_2$box_data, type = "2")
  ))

  peak_plot <- ggplot2::ggplot(case_box_data, ggplot2::aes(x = type, y = peak_cases, fill = type)) +
    ggplot2::geom_boxplot(outlier.alpha = 0.2) +
    ggplot2::labs(x = NULL, y = "Peak Cases") +
    ggplot2::theme_bw() +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::theme(legend.position = "none")

  total_plot <- ggplot2::ggplot(case_box_data, ggplot2::aes(x = type, y = total_cases, fill = type)) +
    ggplot2::geom_boxplot(outlier.alpha = 0.2) +
    ggplot2::labs(x = "Scenario", y = "Cumulative Cases") +
    ggplot2::theme_bw() +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::theme(legend.position = "none")

  # --- Susceptible proportion comparison ---
  scenario_1_susc <- dplyr::mutate(
    subset(scenario_1$susceptibility_data, status == "Susceptible"),
    type = 1
  )

  scenario_2_susc <- dplyr::mutate(
    subset(scenario_2$susceptibility_data, status == "Susceptible"),
    type = 2
  )

  baseline_susc <- scenario_1_susc |>
    subset(time == min(time)) |>
    dplyr::select(prop, prop_lower, prop_upper) |>
    tidyr::pivot_longer(cols = dplyr::everything())

  susceptibility_plot <- ggplot2::ggplot(
    data = data.table::rbindlist(list(scenario_1_susc, scenario_2_susc)),
    ggplot2::aes(
      x = time,
      y = prop,
      ymin = prop_lower,
      ymax = prop_upper,
      color = factor(type))
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(alpha = 0.25) +
    ggplot2::geom_hline(yintercept = baseline_susc$value[1], linetype = "dashed") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Time", y = "Susceptible Population", color = "Scenario") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks())

  # --- Return ---
  invisible(list(
    case_plot           = subset_case_plot,
    cumulative_plot     = cumulative_case_plot,
    boxplot_peak        = peak_plot,
    boxplot_total       = total_plot,
    susceptibility_plot = susceptibility_plot,
    summary_table       = case_summary_format
  ))
}
