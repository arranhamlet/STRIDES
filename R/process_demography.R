#' Prepare Demographic Inputs for a Transmission Model
#'
#' This function processes raw demographic data—covering migration, fertility, mortality,
#' population counts, and age-specific contact matrices—into a structured format suitable for use
#' in an age-, vaccination-, and risk-structured infectious disease transmission model.
#'
#' The function uses tidyverse conventions and assumes consistent
#' naming of age columns (e.g., `x0` to `x100` for all population counts).
#'
#' @param migration A data frame with net migration rates by year and country.
#' @param fertility A data frame with fertility rates by single year of age (15–49) and year.
#' @param mortality A data frame with mortality counts by single year of age (0–100) and year.
#' @param population_all A data frame of total population counts by single year of age (0–100) and year.
#' @param population_female A data frame of female population counts by single year of age (15–49) and year.
#' @param contact_matricies A named list of 16x16 contact matrices by country ISO3 code.
#' @param iso A three-letter ISO country code.
#' @param year_start First year to include.
#' @param year_end Final year to include.
#' @param n_age Number of collapsed age bins to use in the model.
#' @param number_of_vaccines Number of vaccine doses considered.
#' @param n_risk Number of risk strata in the model.
#'
#' @return A named list of processed demographic inputs.
#'
#' @keywords internal
#' @importFrom dplyr filter select mutate bind_rows across case_when left_join all_of
#' @importFrom tidyr expand_grid pivot_longer
#' @importFrom purrr map_dfr
#' @importFrom data.table melt
#' @importFrom tibble tibble rownames_to_column
#' @importFrom stats setNames
#' @importFrom reshape2 melt

process_demography <- function(
    migration,
    fertility,
    mortality,
    population_all,
    population_female,
    contact_matricies,
    iso,
    year_start = "",
    year_end = "",
    n_age = 1,
    number_of_vaccines = 0,
    n_risk = 1
) {

  n_vacc <- if (number_of_vaccines == 0) 1 else number_of_vaccines * 2 + 1

  years_all <- get_years(migration$year, start = year_start, end = year_end)
  time_run_for <- length(years_all)
  time_all <- 0:(time_run_for - 1)

  migration <- dplyr::filter(migration, iso3 == iso, year %in% years_all)
  fertility <- dplyr::filter(fertility, iso3 == iso, year %in% years_all)
  mortality <- dplyr::filter(mortality, iso3 == iso, year %in% years_all)
  population_all <- dplyr::filter(population_all, iso3 == iso, year %in% years_all)
  population_female <- dplyr::filter(population_female, iso3 == iso, year %in% years_all)

  all_ages <- 0:100
  female_ages <- 15:49
  pop_cols <- paste0("x", all_ages)
  fem_cols <- paste0("x", female_ages)

  pop_all_raw <- population_all %>% dplyr::select(dplyr::all_of(pop_cols)) %>% as.matrix()
  pop_all <- collapse_age_bins(pop_all_raw, n_age)

  age_groups <- sapply(split(all_ages, sort(all_ages %% n_age)), min)

  country_contact <- if (!iso %in% names(contact_matricies)) {
    Reduce(`+`, contact_matricies) / length(contact_matricies)
  } else {
    contact_matricies[[iso]]
  }

  reformatted_contact_matrix <- reformat_contact_matrix(country_contact, age_groups)
  reformatted_contact_matrix <- symmetrize_contact_matrix(reformatted_contact_matrix, pop = pop_all[nrow(pop_all), ])
  reformatted_contact_matrix <- project_to_symmetric_doubly_stochastic(reformatted_contact_matrix)

  mort_mat <- mortality %>% dplyr::select(dplyr::all_of(pop_cols)) %>% as.matrix()
  mort_mat <- collapse_age_bins(mort_mat, n_age)
  mortality_rate <- pmin(mort_mat / pop_all, 1)
  mortality_rate[!is.finite(mortality_rate)] <- 1

  mortality_df <- reshape2::melt(t(mortality_rate)) %>%
    stats::setNames(c("dim1", "dim3", "value")) %>%
    dplyr::mutate(dim2 = 1)

  fert_mat <- fertility %>% dplyr::select(dplyr::all_of(fem_cols)) %>% as.data.frame()

  fem_mat <- population_female %>%
    dplyr::select(dplyr::all_of(fem_cols)) %>%
    tibble::rownames_to_column(var = "dim2") %>%
    tidyr::pivot_longer(cols = -dim2, names_to = "dim1", values_to = "population") %>%
    dplyr::mutate(dim1 = as.integer(gsub("x", "", dim1)) + 1)

  births_by_year <- fert_mat %>%
    tibble::rownames_to_column(var = "dim2") %>%
    tidyr::pivot_longer(cols = -dim2, names_to = "dim1", values_to = "value") %>%
    dplyr::mutate(dim1 = as.integer(gsub("x", "", dim1)) + 1) %>%
    dplyr::left_join(fem_mat, by = c("dim1", "dim2")) %>%
    dplyr::mutate(value = value * population)

  mig_rates <- migration[["migration_rate_1000"]]

  migration_in_number <- purrr::map_dfr(seq_len(nrow(pop_all)), function(i) {
    mig_vals <- round(pop_all[i, ] * mig_rates[i])
    chunk <- split_and_sum(mig_vals, n_age)
    tibble::tibble(
      dim1 = seq_len(n_age),
      dim2 = 1,
      dim3 = 1,
      dim4 = i,
      value = chunk
    )
  })

  init_vals <- round(pop_all[1, ] * 1000)
  init_chunk <- split_and_sum(init_vals, n_age)

  N0_df <- tibble::tibble(
    dim1 = seq_len(n_age),
    dim2 = 1,
    dim3 = 1,
    value = init_chunk
  )

  total_population_df <- purrr::map_dfr(seq_len(nrow(pop_all)), function(i) {
    chunk <- split_and_sum(round(pop_all[i, ] * 1000), n_age)
    tibble::tibble(
      dim1 = seq_len(n_age),
      dim2 = 1,
      dim3 = 1,
      dim4 = i,
      value = chunk
    )
  })

  migration_distribution_values <- tidyr::expand_grid(
    dim1 = 1,
    dim2 = seq_along(time_all)
  ) %>% dplyr::mutate(value = 1)

  list(
    N0 = N0_df,
    crude_birth = births_by_year,
    crude_death = mortality_df,
    tt_migration = time_all,
    migration_in_number = migration_in_number,
    migration_distribution_values = migration_distribution_values,
    population_data = pop_all,
    female_population = fem_mat,
    contact_matrix = reformatted_contact_matrix,
    input_data = tibble::tibble(
      iso = iso,
      year_start = min(years_all),
      year_end = max(years_all),
      n_age = n_age,
      n_vacc = n_vacc,
      n_risk = n_risk
    )
  )
}
