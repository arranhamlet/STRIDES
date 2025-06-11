#' Prepare Model Inputs for an Infectious Disease Transmission Model
#'
#' This wrapper function prepares demographic, case, and vaccination data for use in an
#' age-, vaccination-, and risk-structured infectious disease model.
#' It processes each input type using internal helper functions and returns a fully structured list.
#'
#' @param migration A data frame of net migration rates by year and country.
#' @param fertility A data frame of fertility rates by age and year.
#' @param mortality A data frame of mortality counts by age and year.
#' @param population_all A data frame of total population counts by age and year.
#' @param population_female A data frame of female population counts by age and year.
#' @param iso A 3-letter ISO country code for filtering the data.
#' @param disease_data A data frame of WHO-reported disease cases across countries and years.
#' @param disease A string indicating the disease of interest (e.g., "measles").
#' @param vaccination_data_routine Routine immunization data including coverage, population, antigen, etc.
#' @param vaccination_data_sia Supplementary immunization activity (SIA) coverage data.
#' @param vaccine A string specifying the vaccine of interest.
#' @param contact_matricies A named list of contact matrices by ISO code.
#' @param year_start First year to include (default: from migration data).
#' @param year_end Final year to include (default: from migration data).
#' @param n_age Number of age groups in the model.
#' @param number_of_vaccines Number of vaccine doses to simulate (used to compute strata).
#' @param n_risk Number of risk groups in the model.
#'
#' @return A named list containing:
#' \describe{
#'   \item{processed_demographic_data}{Output from `process_demography()`.}
#'   \item{processed_case_data}{Filtered disease case data.}
#'   \item{processed_vaccination_data}{Cleaned and summarized routine vaccination data.}
#'   \item{processed_vaccination_sia}{Filtered SIA data for the specified vaccine.}
#'   \item{demographic_plots}{(Placeholder) Plot object to visualize processed inputs.}
#' }
#'
#' @details
#' This function calls internal processing utilities: `process_demography()`,
#' `process_prior_cases()`, `process_vaccination_routine()`, `process_vaccination_sia()`.
#'
#' @keywords internal
#' @importFrom dplyr filter
prepare_model_inputs <- function(
  migration,
  fertility,
  mortality,
  population_all,
  population_female,
  iso,
  disease_data,
  disease,
  vaccination_data_routine,
  vaccination_data_sia,
  vaccine,
  contact_matricies,
  year_start = "",
  year_end = "",
  n_age = 1,
  number_of_vaccines = 0,
  n_risk = 1
) {

  # Step 1: Process demographic data
  demographic <- process_demography(
    migration = migration,
    fertility = fertility,
    mortality = mortality,
    population_all = population_all,
    population_female = population_female,
    contact_matricies = contact_matricies,
    iso = iso,
    year_start = year_start,
    year_end = year_end,
    number_of_vaccines = number_of_vaccines,
    n_risk = n_risk,
    n_age = 101
  )

  # Filter external inputs to final model year
  vaccination_data_routine <- vaccination_data_routine %>%
    filter(year <= demographic$input_data$year_end)

  disease_data <- disease_data %>%
    filter(year <= demographic$input_data$year_end)

  # Step 2: Process disease data
  cases <- process_prior_cases(
    disease_data = disease_data,
    iso = demographic$input_data$iso,
    disease = disease,
    year_start = demographic$input_data$year_start,
    year_end = demographic$input_data$year_end
  )

  # Step 3: Routine vaccination data
  routine <- process_vaccination_routine(
    vaccination_data = vaccination_data_routine,
    iso = demographic$input_data$iso,
    vaccine = vaccine,
    year_start = demographic$input_data$year_start,
    year_end = demographic$input_data$year_end
  )

  # Step 4: SIA vaccination data
  sia <- process_vaccination_sia(
    vaccination_data = vaccination_data_sia,
    iso = demographic$input_data$iso,
    vaccine = vaccine,
    year_start = demographic$input_data$year_start,
    year_end = demographic$input_data$year_end
  )

  # Step 5: Compile output
  list(
    processed_demographic_data = demographic,
    processed_case_data = cases,
    processed_vaccination_data = routine,
    processed_vaccination_sia = sia
  )
}
