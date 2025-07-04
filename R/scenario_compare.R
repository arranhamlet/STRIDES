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
#'   - `cumulative_case_plot$data`: Cumulative case trajectories
#'   - `susceptibility_data`: Long-format susceptibility proportions by time
#' @param scenario_2 Same structure as `scenario_1`, for the second scenario.
#'
#' @return A named list containing:
#' \item{summary_plot}{A combined patchwork plot of all summary graphics.}
#' \item{summary_table}{A formatted `gt` table of case summary statistics.}
#' \item{individual_plots}{A list of individual `ggplot2` plots for further use.}
#'
#' @importFrom data.table rbindlist data.table :=
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon labs
#'   scale_y_continuous scale_x_continuous theme_bw theme
#' @importFrom gt gt fmt_number cols_label tab_options tab_style
#'   cell_text cells_column_labels opt_table_font opt_row_striping
#' @importFrom dplyr mutate select
#' @importFrom tidyr pivot_longer
#' @importFrom scales comma percent_format pretty_breaks
#' @importFrom patchwork wrap_table plot_layout
#' @export
scenario_compare <- function(scenario_1, scenario_2) {

  # --- Combine time-series case data ---
  case_combo <- data.table::rbindlist(list(
    dplyr::mutate(subset(scenario_1$select_state_plot$data, state_plot == "Cases"), type = 1),
    dplyr::mutate(subset(scenario_2$select_state_plot$data, state_plot == "Cases"), type = 2)
  ))

  # --- Identify time window with cases ---
  case_time_range <- case_combo[median != 0, range(time)]

  # --- Plot: Weekly case time-series ---
  subset_case_plot <- ggplot2::ggplot(
    case_combo[time %in% case_time_range[1]:case_time_range[2]],
    ggplot2::aes(x = time, y = median, ymin = lower, ymax = upper,
                 color = factor(type), fill = factor(type))
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(alpha = 0.25, show.legend = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Time", y = "", subtitle = "Cases", color = "Scenario") +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks())

  # --- Summary statistics (per scenario and differences) ---
  case_summary <- case_combo[
    , .(
      peak_median_cases   = max(median, na.rm = TRUE),
      time_of_peak        = time[which.max(median)],
      total_median_cases  = sum(median, na.rm = TRUE),
      total_lower_cases   = sum(lower, na.rm = TRUE),
      total_upper_cases   = sum(upper, na.rm = TRUE)
    ),
    by = type
  ]

  diff_row <- data.table::data.table(
    type                = "Difference",
    peak_median_cases   = diff(case_summary$peak_median_cases),
    time_of_peak        = diff(case_summary$time_of_peak),
    total_median_cases  = diff(case_summary$total_median_cases),
    total_lower_cases   = diff(case_summary$total_lower_cases),
    total_upper_cases   = diff(case_summary$total_upper_cases)
  )

  case_summary_full <- data.table::rbindlist(list(case_summary, diff_row), fill = TRUE)

  # --- Susceptibility comparison ---
  scenario_1_susc <- dplyr::mutate(
    subset(scenario_1$susceptibility_data, status == "Susceptible"),
    type = 1
  )
  scenario_2_susc <- dplyr::mutate(
    subset(scenario_2$susceptibility_data, status == "Susceptible"),
    type = 2
  )

  # Calculate max susceptible proportion for each scenario
  max_susc_1 <- max(scenario_1_susc$prop, na.rm = TRUE)
  max_susc_2 <- max(scenario_2_susc$prop, na.rm = TRUE)
  percent_increase_susc <- (max_susc_2 - max_susc_1) * 100

  # Add new column for Susceptible % Increase (only in difference row)
  case_summary_full[, susceptible_pct_increase := c(sprintf("%.1f%%", max_susc_1 * 100), sprintf("%.1f%%", max_susc_2 * 100), sprintf("%.1f%%", percent_increase_susc))]

  # --- Format summary statistics as a gt table ---
  case_summary_format <- gt::gt(case_summary_full) |>
    gt::fmt_number(columns = where(is.numeric), decimals = 0, use_seps = TRUE) |>
    gt::cols_label(
      type                       = "Scenario",
      peak_median_cases          = "Peak Median",
      time_of_peak               = "Time of Peak",
      total_median_cases         = "Total Median",
      total_lower_cases          = "Total Lower",
      total_upper_cases          = "Total Upper",
      susceptible_pct_increase   = "Susceptible % Increase"
    ) |>
    gt::tab_options(data_row.padding = gt::px(3)) |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_column_labels(gt::everything())
    ) |>
    gt::opt_table_font(font = "sans") |>
    gt::opt_row_striping()

  # --- Combine cumulative trajectories ---
  cumulative_cases <- data.table::rbindlist(list(
    dplyr::mutate(scenario_1$cumulative_case_plot$data, type = "1"),
    dplyr::mutate(scenario_2$cumulative_case_plot$data, type = "2")
  ))

  cumulative_case_plot <- ggplot2::ggplot(
    cumulative_cases[time %in% case_time_range[1]:case_time_range[2]],
    ggplot2::aes(x = time, y = cum_median, ymin = cum_lower, ymax = cum_upper,
                 color = factor(type), fill = factor(type))
  ) +
    ggplot2::geom_ribbon(alpha = 0.3) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::labs(x = "Time", y = "", subtitle = "Cumulative Cases", color = "Scenario", fill = "Scenario") +
    ggplot2::theme_bw() +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks())

  # --- Plot baseline susceptible proportion ---
  baseline_susc <- scenario_1_susc %>%
    subset(time == min(time)) %>%
    dplyr::select(prop, prop_lower, prop_upper) %>%
    tidyr::pivot_longer(cols = dplyr::everything())

  susceptibility_plot <- ggplot2::ggplot(
    data = data.table::rbindlist(list(scenario_1_susc, scenario_2_susc)),
    ggplot2::aes(x = time, y = prop, ymin = prop_lower, ymax = prop_upper,
                 fill = factor(type), color = factor(type))
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(alpha = 0.25) +
    ggplot2::geom_hline(yintercept = baseline_susc$value[1], linetype = "dashed") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Time", y = "", subtitle = "Susceptible Population", fill = "Scenario", color = "Scenario") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks())

  # --- Combine all plots and table ---
  summary_linegraphs <- (
    (subset_case_plot + ggplot2::theme(legend.position = "none")) +
      cumulative_case_plot + ggplot2::theme(legend.position = "none") +
      susceptibility_plot +
      patchwork::plot_layout(guides = "collect")
  ) /
    patchwork::wrap_table(case_summary_format, panel = "full")

  # --- Return named list of results ---
  list(
    summary_plot = summary_linegraphs,
    summary_table = case_summary_format,
    individual_plots = list(
      case_plot = subset_case_plot,
      cumulative_plot = cumulative_case_plot,
      susceptibility_plot = susceptibility_plot
    )
  )
}
