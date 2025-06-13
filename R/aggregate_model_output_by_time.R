#' Aggregate Model Output by Time Units (data.table version)
#'
#' Converts model `time` into calendar dates from a reference start date, and aggregates
#' values to a target time unit (e.g., weekly, monthly, yearly), grouped by all other variables.
#'
#' @param df A data.frame or data.table with `time`, `value`, and stratifier columns.
#' @param timestep Character. One of "day", "week", "month", or "year". Time step used in the model.
#' @param target_unit Character. One of "week", "month", or "year". Unit to aggregate into.
#' @param start_date A string date (e.g. "2020-01-01") indicating the model time origin.
#'
#' @return A data.table with time-aggregated results.
#' @import data.table
#' @importFrom lubridate floor_date
#' @export
aggregate_model_output_by_time <- function(df, timestep = "day", target_unit = "week", start_date = "2000-01-01") {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install the 'data.table' package.")
  if (!requireNamespace("lubridate", quietly = TRUE)) stop("Please install the 'lubridate' package.")

  timestep <- tolower(timestep)
  target_unit <- tolower(target_unit)

  step_days <- switch(timestep,
                      day     = 1,
                      week    = 7,
                      month   = 30,
                      quarter = 91.25,
                      year    = 365,
                      stop("Unsupported timestep.")
  )

  date0 <- as.Date(start_date)

  # Ensure data.table
  dt <- data.table::as.data.table(df)

  # Convert model time to actual date
  dt[, date := date0 + time * step_days]

  # Determine time period for aggregation
  dt[, period := switch(target_unit,
                        week  = as.IDate(lubridate::floor_date(date, unit = "week")),
                        month = as.IDate(lubridate::floor_date(date, unit = "month")),
                        year  = as.IDate(lubridate::floor_date(date, unit = "year")),
                        stop("Unsupported target_unit. Use 'week', 'month', or 'year'.")
  )]

  # Identify group columns
  group_cols <- setdiff(names(dt), c("value", "time", "date"))

  # Aggregate
  out <- dt[, .(value = sum(value, na.rm = TRUE)), by = group_cols]

  return(out[])
}
