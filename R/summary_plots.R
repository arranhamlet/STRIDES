#' Generate Summary Plots and Susceptibility Profiles
#'
#' Processes model output and generates:
#' 1. A time-series plot of selected states.
#' 2. A stacked bar chart of age-specific susceptibility profiles.
#' 3. A long-format data.table with stratified susceptibility proportions.
#'
#' @param model_run A `data.frame` or `data.table` in long format, including columns `time`, `state`, `value`, `age`, `risk`, and `vaccination`.
#' @param params A named list of model input parameters, including `n_age`, `timestep`, and `input_data$year_start` and `input_data$age_breaks`.
#'
#' @return A named list containing:
#' \item{select_state_plot}{A `ggplot` object showing state dynamics over time.}
#' \item{susceptibility_plot}{A `ggplot` object showing age-specific susceptibility in the latest year.}
#' \item{susceptibility_data}{A `data.table` summarizing protection composition by age and year.}
#'
#' @importFrom data.table setDT fifelse CJ
#' @importFrom dplyr case_when
#' @importFrom ggplot2 ggplot aes geom_line geom_bar facet_wrap scale_y_continuous labs theme_bw
#' @export
summary_plots <- function(model_run, params) {

  year_start <- if (params$input_data$year_start == "") 1950 else params$input_data$year_start
  timestep <- params$input_data$timestep
  age_groups <- as.numeric(unlist(strsplit(params$input_data$age_breaks, ";")))

  data.table::setDT(model_run)

  model_run_plot <- model_run[
    state %in% c("total_pop", "S", "E", "I", "R", "Is", "Rc", "new_case") &
      age == "All" & time >= 365 / 1
  ][
    , state_plot := dplyr::case_when(
      state == "S" ~ "Susceptible",
      state == "E" ~ "Exposed",
      state == "I" ~ "Infectious",
      state == "R" ~ "Recovered",
      state == "Is" ~ "Infectious, severe",
      state == "Rc" ~ "Recovered, complications",
      state == "total_pop" ~ "Total population",
      state == "new_case" ~ "Cases",
      TRUE ~ state
    )
  ][
    , state_plot := factor(state_plot, levels = c(
      "Susceptible", "Exposed", "Infectious", "Infectious, severe",
      "Recovered", "Recovered, complications", "Total population", "Cases"
    ))
  ]

  susceptibility_data <- model_run[
    state %in% c("S", "E", "I", "R", "Is", "Rc") & age != "All"
  ][
    , year := year_start + floor(as.numeric(time) / (365 / 1))
  ][
    order(time)
  ][
    , .SD[.N], by = .(year, age, risk, vaccination, state)
  ][
    , status := data.table::fifelse(state == "S" & vaccination == 1, "Susceptible",
                                    data.table::fifelse(state == "S" & vaccination > 1, "Vaccine protected",
                                                        data.table::fifelse(state != "S" & vaccination == 1, "Exposure protected",
                                                                            "Vaccine and exposure protected")))
  ][
    , status := factor(status, levels = c(
      "Susceptible", "Vaccine protected", "Exposure protected", "Vaccine and exposure protected"
    ))
  ][
    , .(value = sum(value, na.rm = TRUE)), by = .(year, age, risk, vaccination, status)
  ][
    , prop := value / sum(value), by = .(year, age)
  ]

  # Create age group labels based on age index (1 = 0-5, etc.)
  n_labels <- length(age_groups) - 1
  age_group_labels <- paste0(
    age_groups[-length(age_groups)],
    "–",
    age_groups[-1] - 1
  )
  age_group_labels[n_labels] <- paste0(age_groups[n_labels], "+")

  susceptibility_data[, age_group_label := factor(age, levels = seq_len(n_labels), labels = age_group_labels)]

  select_state_plot <- ggplot2::ggplot(
    model_run_plot,
    ggplot2::aes(x = time / (365 / 1) + year_start - 1, y = value)
  ) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~state_plot, scales = "free_y") +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "", y = "")

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
