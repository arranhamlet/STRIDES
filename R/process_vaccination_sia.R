#' Process Supplementary Immunization Activity (SIA) Data
#'
#' Filters and processes SIA vaccination data for a specified country, time range, and vaccine.
#'
#' @param vaccination_data A `data.frame` or `data.table` containing SIA records. Must include columns `area`, `year`, `coverage`, and `vaccine`.
#' @param iso A 3-letter ISO country code to filter the data by (matched against the `area` column).
#' @param vaccine A specific vaccine name or substring to filter by (default is `"All"` which returns all vaccines).
#' @param year_start First year to include (default: earliest year in data).
#' @param year_end Final year to include (default: latest year in data).
#'
#' @return A filtered `data.table` with records matching the selected country, vaccine, and time period.
#'
#' @importFrom data.table setDT setnames
#' @keywords internal
process_vaccination_sia <- function(
    vaccination_data,
    iso,
    vaccine = "All",
    year_start = "",
    year_end = ""
) {
  # Ensure data.table format
  data.table::setDT(vaccination_data)

  # Standardize column name if needed
  data.table::setnames(vaccination_data, old = "vaccine", new = "vaccination_name", skip_absent = TRUE)

  # Determine time window
  years <- get_years(vaccination_data$year, year_start, year_end)

  # Filter
  result <- vaccination_data[
    !is.na(coverage) &
      area == iso &
      year %in% years &
      (vaccine == "All" | grepl(vaccine, vaccination_name, ignore.case = TRUE))
  ]

  return(result[])
}
