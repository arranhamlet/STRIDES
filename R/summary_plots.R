#' Generate Summary Plots and Susceptibility Profiles
#'
#' This function processes model output and generates:
#' 1. A time-series plot of selected states (weekly only).
#' 2. A stacked bar chart of age-specific susceptibility profiles (annual only).
#' 3. A line chart of vaccination coverage over time.
#' 4. A long-format `data.table` with stratified susceptibility proportions.
#' 5. A cumulative case time-series plot.
#' 6. A summary `box_data` table with peak and total case statistics per run.
#'
#' @param model_run A `data.frame` or `data.table` in long format, including columns
#'   `time`, `state`, `value`, `age`, `risk`, `vaccination`, and `run`.
#' @param params A named list of model input parameters, including:
#'   - `input_data$year_start`: Start year of the simulation ("" uses 1950).
#'   - `input_data$age_breaks`: Age breakpoints as a semicolon-separated string.
#'
#' @return A named list containing:
#' \item{select_state_plot}{A ggplot2 object: time-series of key model states.}
#' \item{susceptibility_plot}{A ggplot2 object: latest age-specific susceptibility profile.}
#' \item{vaccination_coverage_plot}{A ggplot2 object: vaccination coverage over time.}
#' \item{cumulative_case_plot}{A ggplot2 object: cumulative case trajectory with uncertainty.}
#' \item{susceptibility_data}{A `data.table`: weekly susceptibility proportions.}
#' \item{box_data}{A `data.table`: peak and total cases per run.}
#'
#' @importFrom data.table setDT fcase
#' @importFrom collapse fgroup_by fsummarise fquantile fsum
#' @importFrom scales comma breaks_pretty percent_format
#' @importFrom ggplot2 ggplot aes geom_line geom_bar geom_ribbon facet_wrap
#'   scale_y_continuous scale_x_continuous labs theme_bw coord_cartesian
#' @importFrom dplyr mutate case_when group_by
#' @export
summary_plots <- function(model_run, params) {

  # --- Setup and Input Parsing ---
  year_start <- if (params$input_data$year_start == "") 1950 else params$input_data$year_start
  age_groups <- as.numeric(unlist(strsplit(params$input_data$age_breaks, ";")))
  data.table::setDT(model_run)

  # --- Filter and label relevant model states ---
  state_filter <- c("total_pop", "S", "E", "I", "R", "Is", "Rc", "new_case")
  model_run_plot <- model_run[state %chin% state_filter & age == "All"]

  model_run_plot[, state_plot := data.table::fcase(
    state == "S",         "Susceptible",
    state == "E",         "Exposed",
    state == "I",         "Infectious",
    state == "R",         "Recovered",
    state == "Is",        "Infectious, severe",
    state == "Rc",        "Recovered, complications",
    state == "total_pop", "Total population",
    state == "new_case",  "Cases"
  )][, state_plot := factor(state_plot, levels = c(
    "Susceptible", "Exposed", "Infectious", "Infectious, severe",
    "Recovered", "Recovered, complications", "Total population", "Cases"
  ))]

  # --- Cumulative case trajectory with uncertainty ---
  cumulative_cases <- model_run_plot[state == "new_case" & age == "All"][
    order(time)
  ][
    , .(
      time = unique(time),
      value = cumsum(value)
      ), by = .(run)
  ][
    , .(
      cum_median = collapse::fquantile(value, 0.5),
      cum_lower  = collapse::fquantile(value, 0.025),
      cum_upper  = collapse::fquantile(value, 0.975)
    ), by = time
  ]

  cumulative_case_plot <- ggplot2::ggplot(
    cumulative_cases,
    ggplot2::aes(x = time, y = cum_median, ymin = cum_lower, ymax = cum_upper)
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(alpha = 0.25) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::scale_x_continuous(breaks = scales::breaks_pretty()) +
    ggplot2::labs(x = "", y = "") +
    ggplot2::theme_bw()

  # --- Weekly state summaries for plotting ---
  summary_plot_data <- model_run_plot[time %% 7L == 1L] %>%
    collapse::fgroup_by(time, state_plot) %>%
    collapse::fsummarise(
      median = collapse::fquantile(value, 0.5),
      lower  = collapse::fquantile(value, 0.025),
      upper  = collapse::fquantile(value, 0.975)
    )

  select_state_plot <- ggplot2::ggplot(
    summary_plot_data,
    ggplot2::aes(x = time / 365 + year_start - 1, y = median, ymin = lower, ymax = upper)
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(alpha = 0.5) +
    ggplot2::facet_wrap(~state_plot, scales = "free_y") +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::scale_x_continuous(breaks = scales::breaks_pretty()) +
    ggplot2::labs(x = "", y = "") +
    ggplot2::theme_bw()

  # --- Susceptibility snapshot by age/year ---
  susceptibility_data_all <- model_run[
    state %chin% c("S", "E", "I", "R", "Is", "Rc") & age != "All" & time %% 365 == 1L
  ][, year := year_start + floor(time / 365)][
    , .SD[.N], by = .(year, age, risk, vaccination, run, state)
  ][, status := data.table::fcase(
    state == "S" & vaccination == 1, "Susceptible",
    state == "S", "Vaccine derived protection",
    vaccination == 1, "Infection derived protection",
    default = "Vaccine and infection derived protection"
  )][, status := factor(status, levels = c(
    "Susceptible", "Vaccine derived protection",
    "Infection derived protection", "Vaccine and infection derived protection"
  ))]

  # --- Detailed weekly susceptibility proportions ---
  detailed_susceptibility <- model_run[
    state %chin% c("S", "E", "I", "R", "Is", "Rc") & age != "All" & time %% 7 == 1L
  ][, status := data.table::fcase(
    state == "S" & vaccination == 1, "Susceptible",
    state == "S", "Vaccine derived protection",
    vaccination == 1, "Infection derived protection",
    default = "Vaccine and infection derived protection"
  )][, status := factor(status, levels = c(
    "Susceptible", "Vaccine derived protection",
    "Infection derived protection", "Vaccine and infection derived protection"
  ))][
    , .(value = sum(value, na.rm = TRUE)), by = .(run, time, status)
  ][
    , .(
      median = collapse::fquantile(value, 0.5),
      lower  = collapse::fquantile(value, 0.025),
      upper  = collapse::fquantile(value, 0.975)
    ), by = .(time, status)
  ][
    , total_median := sum(median), by = time
  ][
    , `:=`(
      prop        = median / total_median,
      prop_lower  = lower  / total_median,
      prop_upper  = upper  / total_median
    )
  ][, total_median := NULL]

  # --- Annual susceptibility proportions by age group ---
  susceptibility_data <- susceptibility_data_all[
    , .(value = sum(value, na.rm = TRUE)), by = .(time, year, age, status)
  ][
    , .(value = sum(value, na.rm = TRUE)), by = .(year, age, status)
  ][
    , prop := value / sum(value), by = .(year, age)
  ]

  n_labels <- length(age_groups) - 1
  age_group_labels <- paste0(age_groups[-length(age_groups)], "–", age_groups[-1] - 1)
  age_group_labels[n_labels] <- paste0(age_groups[n_labels], "+")
  susceptibility_data[, age_group_label := factor(age, levels = seq_len(n_labels), labels = age_group_labels)]

  latest_year <- max(susceptibility_data$year)
  susceptibility_plot <- ggplot2::ggplot(
    susceptibility_data[year == latest_year],
    ggplot2::aes(x = age_group_label, y = prop, fill = status)
  ) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(x = "", y = "Proportion", fill = "") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::theme_bw()

  # --- Vaccination coverage plot ---
  vac_cov <- dplyr::mutate(susceptibility_data,
                           status_simple = dplyr::case_when(
                             grepl("Vaccine", status) ~ "Vaccinated",
                             TRUE ~ "Unvaccinated"
                           )
  ) %>%
    collapse::fgroup_by(year, status_simple) %>%
    collapse::fsummarise(value = collapse::fsum(value)) %>%
    dplyr::group_by(year) %>%
    dplyr::mutate(prop = value / sum(value)) %>%
    dplyr::filter(status_simple == "Vaccinated")

  vaccination_coverage_plot <- ggplot2::ggplot(
    vac_cov,
    ggplot2::aes(x = year, y = prop)
  ) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "", y = "Vaccination coverage (total population)") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::scale_x_continuous(breaks = scales::breaks_pretty()) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::theme_bw()

  # --- Return structured output ---
  list(
    select_state_plot = select_state_plot,
    susceptibility_plot = susceptibility_plot,
    vaccination_coverage_plot = vaccination_coverage_plot,
    cumulative_case_plot = cumulative_case_plot,
    susceptibility_data = detailed_susceptibility
  )
}
