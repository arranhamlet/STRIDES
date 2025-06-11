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
#'
#' @import data.table
#' @keywords internal
process_vaccination_routine <- function(
    vaccination_data,
    iso,
    vaccine = "All",
    year_start = "",
    year_end = ""
) {
  # Ensure data.table format
  setDT(vaccination_data)

  # Standardize column names like janitor::clean_names()
  setnames(vaccination_data, tolower(gsub("[^[:alnum:]]+", "_", names(vaccination_data))))
  # Standardize variable name to match expectations
  setnames(vaccination_data, old = "vaccine", new = "vaccination_name", skip_absent = TRUE)

  years <- get_years(vaccination_data$year, year_start, year_end)

  # Basic filtering
  filtered <- vaccination_data[
    iso3 == iso & year %in% years & !is.na(coverage)
  ]

  # Optional vaccine string filter
  if (vaccine != "All") {
    filtered <- filtered[
      grepl(vaccine, vaccine_description, ignore.case = TRUE) |
        grepl(vaccine, disease, ignore.case = TRUE)
    ]
  }

  # Infer correct dose_order if necessary
  filtered[
    grepl("4th", vaccine_description, ignore.case = TRUE), dose_order := 4
  ][
    grepl("5th", vaccine_description, ignore.case = TRUE), dose_order := 5
  ][
    grepl("6th", vaccine_description, ignore.case = TRUE), dose_order := 6
  ]

  # Aggregate by year, vaccination_name, vaccine_description, dose_order
  result <- filtered[
    ,
    .(coverage = median(coverage, na.rm = TRUE)),
    by = .(year, vaccination_name, vaccine_description, dose_order)
  ][
    order(year)
  ]

  setnames(result, "vaccination_name", "vaccine")

  return(result[])
}
