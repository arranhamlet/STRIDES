#' Aggregate Model Output by Time Units (State-Specific Aggregation)
#'
#' Converts model `time` into calendar dates from a reference start date, and aggregates
#' values to a target time unit (e.g., weekly, monthly, yearly), grouped by all other variables.
#' Aggregation is determined based on the value of the `state` column:
#'   - Stocks (S, E, I, R, Is, Rc, total_pop): take last value
#'   - Flows (new_case): sum
#'   - Rates (Reff, Reff_age): mean
#'
#' @param df A data.frame or data.table with `time`, `value`, `state`, and other stratifiers.
#' @param timestep Character. One of "day", "week", "month", or "year". Time step used in the model.
#' @param target_unit Character. One of "week", "month", or "year". Unit to aggregate into.
#' @param start_date A string date (e.g. "2020-01-01") indicating the model time origin.
#'
#' @return A data.table with time-aggregated results.
#' @importFrom data.table as.data.table
#' @importFrom lubridate floor_date
#' @export
aggregate_model_output_by_time <- function(df, timestep = "day", target_unit = "week", start_date = "2000-01-01") {

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
  dt <- data.table::as.data.table(df)

  # Convert model time to actual date and period
  dt[, date := date0 + time * step_days]
  dt[, period := switch(target_unit,
                        week  = as.IDate(lubridate::floor_date(date, unit = "week")),
                        month = as.IDate(lubridate::floor_date(date, unit = "month")),
                        year  = as.IDate(lubridate::floor_date(date, unit = "year")),
                        stop("Unsupported target_unit. Use 'week', 'month', or 'year'.")
  )]

  # Check required column
  if (!"state" %in% names(dt)) stop("Input data must include a 'state' column.")

  # Define aggregation method per state
  flow_vars <- c("new_case", "t_seeded_al")
  stock_vars <- c("S", "E", "I", "R", "Is", "Rc", "total_pop")
  rate_vars <- c("Reff", "Reff_age")

  # Apply aggregation based on state type
  group_cols <- setdiff(names(dt), c("value", "time", "date"))

  # Split by aggregation type
  out_flow <- dt[state %in% flow_vars, .(value = sum(value, na.rm = TRUE)), by = group_cols]
  out_stock <- dt[state %in% stock_vars, .SD[.N], by = group_cols]  # take last entry
  out_rate <- dt[state %in% rate_vars, .(value = mean(value, na.rm = TRUE)), by = group_cols]

  # Combine all
  out <- data.table::rbindlist(list(out_flow, out_stock, out_rate), use.names = TRUE, fill = TRUE)

  return(out[])
}
