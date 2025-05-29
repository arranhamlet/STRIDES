#' Process Supplementary Immunization Activity (SIA) Data
#'
#' Filters and processes SIA vaccination data for a specified country, time range, and vaccine.
#'
#' @param vaccination_data A data frame or data.table containing SIA records. Must include columns `country`, `year`, `coverage`, and `vaccine`.
#' @param iso A 3-letter ISO country code to filter the data by.
#' @param vaccine A specific vaccine name or substring to filter by (default is "All" which returns all vaccines).
#' @param year_start First year to include (default: first year in the data).
#' @param year_end Final year to include (default: last year in the data).
#'
#' @return A filtered `data.table` with records matching the selected country, vaccine, and time period.
#' @keywords internal
#'
#' @import data.table
#' @import dplyr
process_vaccination_sia <- function(
    vaccination_data,
    iso,
    vaccine = "All",
    year_start = "",
    year_end = ""
){
  setDT(vaccination_data)
  years <- get_years(vaccination_data$year, year_start, year_end)

  vaccination_data <- vaccination_data %>%
    rename(vaccination_name = vaccine)

  vaccination_data[!is.na(coverage) &
                   country == iso &
                   year %in% years &
                   grepl(vaccine, vaccination_name, ignore.case = TRUE)]
}
