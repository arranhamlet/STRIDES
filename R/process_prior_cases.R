#' Process Disease Surveillance Data
#'
#' Filters and subsets disease surveillance data for a specified country, time window, and disease type.
#'
#' @param disease_data A data frame or data.table of disease records with columns for ISO3 code, year, and disease.
#' @param iso A 3-letter ISO country code to filter the data by.
#' @param disease A specific disease to filter by (default is "All" which returns all diseases).
#' @param year_start First year to include (default: first year in data).
#' @param year_end Final year to include (default: last year in data).
#'
#' @return A filtered `data.table` containing only records for the specified country, time range, and disease.
#' @keywords internal
#'
process_prior_cases <- function(
    disease_data,
    iso,
    disease = "All",
    year_start = "",
    year_end = ""
) {

  years <- get_years(disease_data$year, year_start, year_end)

  setDT(disease_data)

  filtered <- disease_data %>%
    dplyr::filter(iso3 == iso, year %in% years, !is.na(cases), cases != 0)

  if (disease != "All") {
    filtered <- filtered %>%
      dplyr::filter(disease_short == disease)
  }

  return(filtered)
}
