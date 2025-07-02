#' Generate Summary Plots and Susceptibility Profiles (Optimized)
#'
#' Processes model output and generates:
#' 1. A time-series plot of selected states (weekly only).
#' 2. A stacked bar chart of age-specific susceptibility profiles (annual only).
#' 3. A long-format data.table with stratified susceptibility proportions.
#'
#' @param model_run A `data.frame` or `data.table` in long format, including columns `time`, `state`, `value`, `age`, `risk`, and `vaccination`.
#' @param params A named list of model input parameters.
#'
#' @return A named list of ggplot objects and summary data.
#'
#' @importFrom data.table setDT fifelse fcase
#' @importFrom collapse fgroup_by fsummarise fquantile
#' @importFrom ggplot2 ggplot aes geom_line geom_bar facet_wrap scale_y_continuous labs theme_bw
#' @export
summary_plots <- function(model_run, params) {

  year_start <- if (params$input_data$year_start == "") 1950 else params$input_data$year_start
  age_groups <- as.numeric(unlist(strsplit(params$input_data$age_breaks, ";")))

  data.table::setDT(model_run)

  # Fast filtering
  state_filter <- c("total_pop", "S", "E", "I", "R", "Is", "Rc", "new_case")
  model_run_plot <- model_run[
    state %chin% state_filter & age == "All" & time %% 7L == 1L
  ]

  # Fast and safe state labeling
  model_run_plot[, state_plot := data.table::fcase(
    state == "S",         "Susceptible",
    state == "E",         "Exposed",
    state == "I",         "Infectious",
    state == "R",         "Recovered",
    state == "Is",        "Infectious, severe",
    state == "Rc",        "Recovered, complications",
    state == "total_pop", "Total population",
    state == "new_case",  "Cases"
  )]

  model_run_plot[, state_plot := factor(state_plot, levels = c(
    "Susceptible", "Exposed", "Infectious", "Infectious, severe",
    "Recovered", "Recovered, complications", "Total population", "Cases"
  ))]

  # Collapse summarisation
  summary_plot_data <- collapse::fgroup_by(model_run_plot, time, state_plot) %>%
    collapse::fsummarise(
      median = collapse::fquantile(value, 0.5),
      lower  = collapse::fquantile(value, 0.025),
      upper  = collapse::fquantile(value, 0.975)
    )

  # Susceptibility summary: annual only
  susceptibility_data <- model_run[
    state %chin% c("S", "E", "I", "R", "Is", "Rc") & age != "All" & time %% 365 == 1L
  ][
    , year := year_start + floor(time / 365)
  ][
    , .SD[.N], by = .(year, age, risk, vaccination, state)
  ][
    , status := data.table::fcase(
      state == "S" & vaccination == 1, "Susceptible",
      state == "S", "Vaccine protected",
      vaccination == 1, "Exposure protected",
      default = "Vaccine and exposure protected"
    )
  ][
    , status := factor(status, levels = c(
      "Susceptible", "Vaccine protected", "Exposure protected", "Vaccine and exposure protected"
    ))
  ][
    , .(value = sum(value, na.rm = TRUE)), by = .(year, age, risk, vaccination, status)
  ][
    , prop := value / sum(value), by = .(year, age)
  ]

  # Label age groups
  n_labels <- length(age_groups) - 1
  age_group_labels <- paste0(age_groups[-length(age_groups)], "–", age_groups[-1] - 1)
  age_group_labels[n_labels] <- paste0(age_groups[n_labels], "+")
  susceptibility_data[, age_group_label := factor(age, levels = seq_len(n_labels), labels = age_group_labels)]

  # Plot: State dynamics
  plot_data <- summary_plot_data[time > 365]
  select_state_plot <- ggplot2::ggplot(plot_data,
                                       ggplot2::aes(x = time / 365 + year_start - 1,
                                                    y = median, ymin = lower, ymax = upper)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(alpha = 0.5) +
    ggplot2::facet_wrap(~state_plot, scales = "free_y") +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "", y = "")

  # Plot: Susceptibility bar chart
  latest_year <- max(susceptibility_data$year)

  susceptibility_plot <- ggplot2::ggplot(
    susceptibility_data[year == latest_year],
    ggplot2::aes(x = age_group_label, y = prop, fill = status)
  ) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "", y = "Proportion", fill = "") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1))

  list(
    select_state_plot = select_state_plot,
    susceptibility_plot = susceptibility_plot,
    susceptibility_data = susceptibility_data
  )
}
