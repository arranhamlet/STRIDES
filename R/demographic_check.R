#' Compare Model Population to Population Reference Estimates
#'
#' This function compares the population from the model output (`model_run`) to
#' UN World Population Prospects (reference population) data. It returns:
#'
#' 1. A line plot comparing total population over time.
#' 2. A bar plot comparing the age distribution in the final year.
#'
#' @param reference_population A matrix or data frame of reference population data,
#'   where each column is a year and each row is an age group. Values should be in thousands.
#'   Provided by `data_load_process_wrapper()`.
#' @param model_run A data frame of stratified model output,
#'   including `time`, `age`, `state`, `value`, `vaccination`, and `risk`. Typically from `summary_plots()`.
#' @param age_breaks A numeric vector of age group breakpoints (e.g., `c(0, 5, 10, ..., 80, Inf)`).
#'
#' @return A named list with:
#' \describe{
#'   \item{population_plot}{A `ggplot` object comparing total population by year.}
#'   \item{age_breakdown_plot}{A `ggplot` object comparing age distribution in the final year.}
#' }
#'
#' @importFrom collapse fgroup_by fsummarise fsum
#' @importFrom dplyr filter mutate case_when group_by slice ungroup bind_rows
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line geom_bar labs theme_bw scale_y_continuous position_dodge theme
#' @importFrom scales comma percent_format
#' @export
demographic_check <- function(reference_population,
                              model_run,
                              age_breaks) {

  # -------------------------
  # Aggregate model population
  # -------------------------
  population_aggregate <- model_run %>%
    subset(state %in% c("S", "E", "I", "R", "Is", "Rc") & age != "All") %>%
    mutate(year = time %/% 365) %>%
    group_by(state, year, age, vaccination, risk) %>%
    slice(which.max(time)) %>%
    ungroup() %>%
    collapse::fgroup_by(year, age) %>%
    collapse::fsummarise(value = collapse::fsum(value)) %>%
    mutate(age = as.numeric(age))

  # Assign age group labels
  age_groups <- age_breaks_to_labels(paste(age_breaks, collapse = ";"))
  population_aggregate$age_group <- factor(age_groups[population_aggregate$age], levels = age_groups)

  # -------------------------
  # Total population by year
  # -------------------------
  population_year <- population_aggregate %>%
    collapse::fgroup_by(year) %>%
    collapse::fsummarise(population = collapse::fsum(value))

  # Build comparison table
  population_comparison <- data.frame(
    year = seq_len(ncol(reference_population)),
    reference = colSums(reference_population),
    model = population_year$population
  ) %>%
    dplyr::filter(reference != 0) %>%
    tidyr::pivot_longer(cols = c("reference", "model"),
                        names_to = "label", values_to = "value") %>%
    mutate(
      year = as.numeric(year),
      label = case_when(
        label == "reference" ~ "Population reference",
        label == "model" ~ "Model estimate"
      )
    )

  # -------------------------
  # Plot: population over time
  # -------------------------
  population_plot <- ggplot2::ggplot(
    data = population_comparison,
    mapping = aes(x = year, y = value, group = label, fill = label)
  ) +
    geom_line(aes(color = label), show.legend = TRUE) +  # add color here for line, but map fill too
    theme_bw() +
    labs(x = "Year", y = "Population", fill = "", color = "") +
    scale_y_continuous(labels = comma)


  # -------------------------
  # Final year age breakdown
  # -------------------------
  final_year <- max(population_aggregate$year)

  final_year_age_model <- population_aggregate %>%
    filter(year == final_year) %>%
    collapse::fgroup_by(age_group) %>%
    collapse::fsummarise(value = collapse::fsum(value)) %>%
    mutate(label = "Model estimate")

  clean_reference <- reference_population[, colSums(reference_population != 0) > 0]

  final_year_age_reference <- data.frame(
    age_group = final_year_age_model$age_group,
    value = clean_reference[, ncol(clean_reference)],
    label = "Population reference"
  )

  final_age <- bind_rows(final_year_age_model, final_year_age_reference) %>%
    group_by(label) %>%
    mutate(prop = value / sum(value))

  # -------------------------
  # Plot: age structure in final year
  # -------------------------
  age_breakdown_plot <- ggplot2::ggplot(
    data = final_age,
    mapping = aes(x = age_group, y = prop * 100, fill = label)
  ) +
    geom_bar(stat = "identity", position = position_dodge(), show.legend = TRUE) +
    theme_bw() +
    labs(x = "Age group", y = "Percent of population", fill = "") +
    scale_y_continuous(labels = comma)

  # -------------------------
  # Return plots
  # -------------------------
  list(
    population_plot = population_plot,
    age_breakdown_plot = age_breakdown_plot
  )
}
