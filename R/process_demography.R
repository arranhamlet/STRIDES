#' Prepare Demographic Inputs for a Transmission Model
#'
#' This function processes raw demographic data—covering migration, fertility, mortality,
#' population counts, and age-specific contact matrices—into a structured format suitable for use
#' in an age-, vaccination-, and risk-structured infectious disease transmission model.
#'
#' The function uses tidyverse conventions (`dplyr`, `purrr`, `tidyr`) and assumes consistent
#' naming of age columns (e.g., `x0` to `x100` for all population counts).
#'
#' @param migration A data frame with net migration rates by year and country. Must include columns
#'   `iso3`, `year`, and `migration_rate_1000`.
#' @param fertility A data frame with fertility rates by single year of age (15–49) and year.
#'   Must include columns `iso3`, `year`, and `x15` to `x49`.
#' @param mortality A data frame with mortality counts by single year of age (0–100) and year.
#'   Must include columns `iso3`, `year`, and `x0` to `x100`.
#' @param population_all A data frame of total population counts by single year of age (0–100) and year.
#'   Must include columns `iso3`, `year`, and `x0` to `x100`.
#' @param population_female A data frame of female population counts by single year of age (15–49) and year.
#'   Must include columns `iso3`, `year`, and `x15` to `x49`.
#' @param contact_matricies A named list of 16x16 contact matrices by country ISO3 code. If a matrix for the
#'   specified country is not available, the mean of all provided matrices is used.
#' @param iso A three-letter ISO country code (character) used to subset the demographic inputs.
#' @param year_start First year to include (numeric or empty string ""). If `""`, uses earliest year in data.
#' @param year_end Final year to include (numeric or empty string ""). If `""`, uses latest year in data.
#' @param n_age Integer. Number of collapsed age bins to use in the model (e.g., 16 for 5-year bins).
#' @param number_of_vaccines Integer. Number of vaccine doses considered. Determines number of strata.
#' @param n_risk Integer. Number of risk strata in the model.
#'
#' @return A named list containing:
#' \describe{
#'   \item{N0}{Tibble of initial population counts by age, risk, and vaccination state for model time 0.}
#'   \item{crude_birth}{Tibble of crude birth rates by time and risk.}
#'   \item{crude_death}{Tibble of annual mortality rates by age.}
#'   \item{tt_migration}{Vector of 0-indexed model time points corresponding to migration years.}
#'   \item{migration_in_number}{Tibble of estimated net migration counts by time, age, and strata.}
#'   \item{migration_distribution_values}{Tibble of ones across all time, for uniform migration distribution.}
#'   \item{population_data}{Matrix of collapsed population counts by age and time (scaled to thousands).}
#'   \item{contact_matrix}{Symmetric, doubly stochastic contact matrix matching collapsed age bins.}
#'   \item{input_data}{Tibble summarizing the scenario metadata (country, years, stratification counts).}
#' }
#'
#' @examples
#' # process_demography(migration, fertility, mortality, population_all, population_female,
#' #                    contact_matricies, iso = "KEN", year_start = 2000, year_end = 2020,
#' #                    n_age = 16, number_of_vaccines = 2, n_risk = 2)
#'
#' @keywords internal
#' @import dplyr
#' @import purrr
#' @import tidyr
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

  migration <- migration %>% filter(iso3 == iso, year %in% years_all)
  fertility <- fertility %>% filter(iso3 == iso, year %in% years_all)
  mortality <- mortality %>% filter(iso3 == iso, year %in% years_all)
  population_all <- population_all %>% filter(iso3 == iso, year %in% years_all)
  population_female <- population_female %>% filter(iso3 == iso, year %in% years_all)

  all_ages <- 0:100
  female_ages <- 15:49
  pop_cols <- paste0("x", all_ages)
  fem_cols <- paste0("x", female_ages)

  pop_all_raw <- population_all %>% select(all_of(pop_cols)) %>% as.matrix()
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

  mort_mat <- mortality %>% select(all_of(pop_cols)) %>% as.matrix()
  mort_mat <- collapse_age_bins(mort_mat, n_age)
  mortality_rate <- pmin(mort_mat / pop_all, 1)
  mortality_rate[!is.finite(mortality_rate)] <- 1

  mortality_df <- reshape2::melt(t(mortality_rate)) %>%
    setNames(c("dim1", "dim3", "value")) %>%
    mutate(dim2 = 1)

  fert_mat <- fertility %>% select(all_of(fem_cols)) %>% as.matrix()
  pop_fem <- population_female %>% select(all_of(fem_cols)) %>% as.matrix()
  denom <- rowSums(pop_fem)
  denom[denom == 0] <- NA

  fertility_by_year <- tibble(
    dim1 = n_risk,
    dim2 = time_all + 1,
    value = pmin(rowSums((fert_mat / 1000) * pop_fem) / denom, 1)
  )

  mig_rates <- migration[["migration_rate_1000"]]

  migration_in_number <- map_dfr(seq_len(nrow(pop_all)), function(i) {
    mig_vals <- round(pop_all[i, ] * mig_rates[i])
    chunk <- split_and_sum(mig_vals, n_age)
    tibble(
      dim1 = seq_len(n_age),
      dim2 = 1,
      dim3 = 1,
      dim4 = i,
      value = chunk
    )
  })

  init_vals <- round(pop_all[1, ] * 1000)
  init_chunk <- split_and_sum(init_vals, n_age)

  N0_df <- tibble(
    dim1 = seq_len(n_age),
    dim2 = 1,
    dim3 = 1,
    value = init_chunk
  )

  total_population_df <- map_dfr(seq_len(nrow(pop_all)), function(i) {
    chunk <- split_and_sum(round(pop_all[i, ] * 1000), n_age)
    tibble(
      dim1 = seq_len(n_age),
      dim2 = 1,
      dim3 = 1,
      dim4 = i,
      value = chunk
    )
  })

  migration_distribution_values <- expand_grid(
    dim1 = 1,
    dim2 = seq_along(time_all)
  ) %>% mutate(value = 1)

  list(
    N0 = N0_df,
    crude_birth = fertility_by_year,
    crude_death = mortality_df,
    tt_migration = time_all,
    migration_in_number = migration_in_number,
    migration_distribution_values = migration_distribution_values,
    population_data = pop_all,
    contact_matrix = reformatted_contact_matrix,
    input_data = tibble(
      iso = iso,
      year_start = min(years_all),
      year_end = max(years_all),
      n_age = n_age,
      n_vacc = n_vacc,
      n_risk = n_risk
    )
  )
}
