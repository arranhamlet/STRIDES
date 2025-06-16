#' Generate Summary Plots and Susceptibility Profiles
#'
#' This function processes model output and produces a list of plots including
#' selected states over time and age-specific susceptibility composition for the latest year.
#'
#' @param model_run A `data.frame` or `data.table` of model output in long format,
#'   including columns `time`, `state`, `value`, `age`, `risk`, and `vaccination`.
#' @param params A named list of model input parameters, including `n_age`, `timestep`,
#'   and `input_data$year_start`.
#'
#' @return A named list containing:
#' \item{select_state_plot}{A ggplot object showing key compartment values over time}
#' \item{susceptibility_plot}{A ggplot object showing age-specific susceptibility composition}
#' \item{susceptibility_data}{A data.table summarizing protection categories by age and year}
#'
#' @import data.table
#' @import ggplot2
#' @export
summary_plots <- function(model_run, params) {

  # Set up default values
  year_start <- if (params$input_data$year_start == "") 1950 else params$input_data$year_start
  timestep <- params$input_data$timestep

  age_groups <- if (params$n_age == 101) 0:100 else c(
    "<1", "1–4", "5–9", "10–14", "15–19", "20–29",
    "30–39", "40–49", "50–59", "60–69", "70–79", "80+"
  )

  time_correct <- as.numeric(dplyr::case_when(
    timestep == "day"     ~ 1,
    timestep == "week"    ~ 7,
    timestep == "month"   ~ 30,
    timestep == "quarter" ~ 91.25,
    timestep == "year"    ~ 365,
    TRUE                  ~ NA_real_
  ))

  # Convert to data.table
  data.table::setDT(model_run)

  # Plot-friendly compartment names
  model_run_plot <- model_run[
    state %in% c("total_pop", "S", "E", "I", "R", "Is", "Rc", "new_case") &
      age == "All" & time >= 365 / time_correct
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

  # Susceptibility classification
  susceptibility_data <- model_run[
    state %in% c("S", "E", "I", "R", "Is", "Rc") & age != "All"
  ][
    , year := year_start + floor(as.numeric(time) / (365 / time_correct))
  ][
    , .(value = mean(value, na.rm = TRUE)),
    by = .(year, age, risk, vaccination, state)
  ][
    , status := fifelse(state == "S" & vaccination == 1, "Susceptible",
                        fifelse(state == "S" & vaccination > 1, "Vaccine protected",
                                fifelse(state != "S" & vaccination == 1, "Exposure protected",
                                        "Vaccine and exposure protected")))
  ][
    , status := factor(status, levels = c(
      "Susceptible", "Vaccine protected", "Exposure protected", "Vaccine and exposure protected"
    ))
  ][
    , .(value = mean(value, na.rm = TRUE)),  # Aggregate over state
    by = .(year, age, risk, vaccination, status)
  ][
    , prop := value / sum(value), by = .(year, age)
  ]

  # Map age group labels
  age_group_map <- setNames(age_groups, as.character(seq_along(age_groups)))
  susceptibility_data[, age_group_label := age_group_map[as.character(age)]]

  # Plot: Time series of selected states
  select_state_plot <- ggplot2::ggplot(
    model_run_plot,
    ggplot2::aes(x = time / (365 / time_correct) + year_start - 1, y = value)
  ) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~state_plot, scales = "free_y") +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "", y = "")

  # Plot: Susceptibility bar plot for final year
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
