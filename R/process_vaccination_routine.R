#' Process Routine Vaccination Coverage Data
#'
#' Filters and processes routine immunization data for a specified country, time range, and optionally vaccine type.
#' Aggregates records by year, antigen, and antigen description, summarising key indicators.
#'
#' @param vaccination_data A data frame or data.table with columns including CODE, YEAR, COVERAGE, ANTIGEN, ANTIGEN_DESCRIPTION, Disease, target_number, and doses.
#' @param iso A 3-letter ISO country code to filter the data by.
#' @param vaccine A specific vaccine name or substring to filter by (default is "All" which returns all vaccines).
#' @param year_start First year to include (default: first year in the data).
#' @param year_end Final year to include (default: last year in the data).
#'
#' @return A data frame summarising median target population, doses delivered, dose order, and coverage by year and antigen.
#' @keywords internal
#'
#' @import data.table
#' @import dplyr
#' @importFrom janitor clean_names
process_vaccination_routine <- function(
    vaccination_data,
    iso,
    vaccine = "All",
    year_start = "",
    year_end = ""
){
  setDT(vaccination_data)
  years <- get_years(vaccination_data$year, year_start, year_end)

  filtered <- vaccination_data[iso3 == iso & year %in% years] %>%
    filter(!is.na(coverage))

  if (vaccine != "All") {
    filtered <- filtered[
      grepl(globalenv()$vaccine, vaccine_description, ignore.case = TRUE) |
        grepl(globalenv()$vaccine, disease, ignore.case = TRUE)
    ]
  }

  filtered %>%
    janitor::clean_names() %>%
    group_by(year, vaccine, vaccine_description, dose_order) %>%
    summarise(
      coverage = median(coverage, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(year) %>%
    mutate(dose_order = case_when(
      grepl("4th", vaccine_description) ~ 4,
      grepl("5th", vaccine_description) ~ 5,
      grepl("6th", vaccine_description) ~ 6,
      TRUE ~ dose_order
    ))
}
