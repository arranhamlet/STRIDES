#' Process Routine Vaccination Coverage Data
#'
#' Filters and processes routine immunization data for a specified country, time range, and optionally vaccine type.
#' Aggregates records by year, antigen, and antigen description, summarising key indicators such as median coverage.
#'
#' @param vaccination_data A `data.frame` or `data.table` with columns including `ISO3`, `YEAR`, `COVERAGE`, `ANTIGEN`,
#' `ANTIGEN_DESCRIPTION`, `Disease`, `target_number`, and `doses`.
#' @param iso A 3-letter ISO country code to filter the data by.
#' @param vaccine A specific vaccine name or substring to filter by (default is `"All"`, which returns all vaccines).
#' @param year_start First year to include (default: earliest year in data).
#' @param year_end Final year to include (default: latest year in data).
#'
#' @return A `data.table` with columns `year`, `vaccine`, `vaccine_description`, `dose_order`, and median `coverage`.
#'
#' @importFrom data.table setDT setnames
#' @keywords internal
process_vaccination_routine <- function(
    vaccination_data,
    iso,
    vaccine = "All",
    year_start = "",
    year_end = ""
) {
  # Ensure data.table format
  data.table::setDT(vaccination_data)

  # Standardize column names
  data.table::setnames(vaccination_data, tolower(gsub("[^[:alnum:]]+", "_", names(vaccination_data))))

  # Standardize variable name to match expectations
  data.table::setnames(vaccination_data, old = "vaccine", new = "vaccination_name", skip_absent = TRUE)

  years <- get_years(vaccination_data$year, year_start, year_end)

  # Filter by ISO3 and year, ensure coverage is non-missing
  filtered <- vaccination_data[
    iso3 == iso & year %in% years & !is.na(coverage)
  ]

  # Optional vaccine string match (vaccine or disease)
  if (vaccine != "All") {
    filtered <- filtered[
      grepl(vaccine, vaccine_description, ignore.case = TRUE) |
        grepl(vaccine, disease, ignore.case = TRUE)
    ]
  }

  # Infer dose_order if not already present
  filtered[
    grepl("4th", vaccine_description, ignore.case = TRUE), dose_order := 4
  ][
    grepl("5th", vaccine_description, ignore.case = TRUE), dose_order := 5
  ][
    grepl("6th", vaccine_description, ignore.case = TRUE), dose_order := 6
  ]

  # Aggregate: median coverage by year and vaccine descriptors
  result <- filtered[
    ,
    .(coverage = median(coverage, na.rm = TRUE)),
    by = .(year, vaccination_name, vaccine_description, dose_order)
  ][
    order(year)
  ]

  data.table::setnames(result, "vaccination_name", "vaccine")

  return(result[])
}
