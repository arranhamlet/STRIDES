#' Generate Aggregated Crude Death Output
#'
#' This function computes the aggregated crude death rates from a preprocessed demographic dataset,
#' using weighted average rates across new age bands. It also computes the total number of deaths before
#' and after aggregation to check for consistency.
#'
#' @param preprocessed A list-like object containing processed demographic data, including `population_data` and `crude_death`.
#' @param new_age_breaks A numeric vector defining the new age group breaks (e.g., c(0, 1, 5, 15, 50, Inf)).
#' @param time_factor A numeric value used to convert annual crude death rate to model timestep (e.g., 365 for daily).
#'
#' @return A list with aggregated crude death data and comparison of total deaths.
#'
#' @importFrom dplyr mutate left_join summarise pull group_by
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @keywords internal
generate_crude_death_outputs <- function(preprocessed, weight_reformatted, age_breaks) {

  #Set up raw data and join the population weights
  raw_d <- preprocessed$processed_demographic_data$crude_death
  comb_d <- raw_d %>%
    left_join(
      weight_reformatted %>%
        dplyr::rename(population = value), by = c("dim3" = "time", "dim1" = "age")
    )

  age_group_index <- cut(1:101, breaks = age_breaks + 1, right = FALSE, labels = FALSE)
  n_age_new <- length(age_breaks) - 1

  # Assign age group to each row based on dim1
  comb_d$age_group <- age_group_index[comb_d$dim1]


  #Work out proportion per age group and then
  agg_combo <- comb_d %>%
    mutate(deaths = population * value) %>%
    group_by(age_group, dim3) %>%
    summarise(
      deaths = sum(deaths),
      population = sum(population)
    ) %>%
    mutate(
      crude_death = deaths
    )

  death_upd <- comb_d %>%
    mutate(
      death = value * population
    ) %>%
    group_by(
      dim3
    ) %>%
    summarise(death = sum(death)) %>%
    mutate(type = "default")

  agg_upd <- agg_combo %>%
    select(-deaths) %>%
    dplyr::rename(death = crude_death,
                  dim1 = age_group) %>%
    group_by(
      dim3
    ) %>%
    summarise(death = 1000 * sum(death)) %>%
    mutate(type = "agg")


  #Why no line up
  ggplot(data = rbind(death_upd, agg_upd),
         mapping = aes(x = dim3, y = death, color = type)) +
    geom_line() +
    theme_bw() +
    labs(x = "", y = "")

  sum(agg_combo$crude_death) == sum(comb_d$value * comb_d$population)




}
