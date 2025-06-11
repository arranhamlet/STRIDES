#' Expand Pre-1980 Routine Vaccination Coverage
#'
#' Uses historical assumptions to estimate vaccination coverage prior to the first observed year.
#' Applies linear interpolation from an assumed introduction year up to the first observed coverage year.
#'
#' @param processed_vaccination Routine coverage data (to determine first year and vaccine).
#' @param vaccination_pre1980 Data frame of historical vaccine assumptions (intro year, starting coverage).
#' @param disease Disease name (lowercase) to match.
#' @param iso 3-letter ISO country code.
#'
#' @return A data frame of additional rows to prepend to coverage time series.
#' @import data.table
#' @keywords internal
expand_pre1980_vaccination <- function(processed_vaccination, vaccination_pre1980, disease, iso) {
  # Ensure input is data.table
  vac_pre1980_sub <- data.table::as.data.table(vaccination_pre1980)
  data.table::setDT(vac_pre1980_sub)

  vac_pre1980_sub[, rowname := .I]  # preserve row identity
  vac_pre1980_sub[, income_group := paste0(paste0(substr(unlist(strsplit(income_group, " |-")), 1, 1), collapse = ""), "C"),
                  by = rowname]
  vac_pre1980_sub[, disease_n := disease]
  vac_pre1980_sub <- vac_pre1980_sub[
    tolower(disease_n) == disease &
      who_region == paste0(get_WHO_region(iso), "O") &
      income_group == as.character(get_income_group(iso))
  ]

  if (nrow(vac_pre1980_sub) == 0 || nrow(processed_vaccination) == 0) return(processed_vaccination)

  # Remove special populations
  processed_vaccination <- processed_vaccination[!grepl("birth|neonatal|pregnant|maternal",
                                                        processed_vaccination$vaccine_description, ignore.case = TRUE), ]

  first_year_vac <- processed_vaccination[processed_vaccination$year == min(processed_vaccination$year), ]

  year_diff <- min(first_year_vac$year) - vac_pre1980_sub$introduction_year
  if (year_diff <= 0) return(processed_vaccination)

  # Precompute interpolation
  vac_prop <- lapply(first_year_vac$coverage, function(e) {
    seq(
      from = vac_pre1980_sub$starting_coverage_percent,
      to = max(vac_pre1980_sub$starting_coverage_percent, e),
      length.out = 6
    )
  })

  # Efficient row expansion
  pre_1980 <- rbindlist(lapply(seq_len(year_diff), function(x) {
    rbindlist(lapply(seq_len(nrow(first_year_vac)), function(k) {
      row <- copy(first_year_vac[k])
      row$year <- vac_pre1980_sub$introduction_year + (x - 1)
      row$coverage <- vac_prop[[k]][min(x, length(vac_prop[[k]]))]
      row
    }))
  }))

  data.table::rbindlist(list(pre_1980, processed_vaccination))
}



#' Format Case and Vaccination Parameters for Model Input (Including Pre-1980)
#'
#' Wrapper function that prepares all time-varying parameters for use in an age-, vaccine-,
#' and risk-structured infectious disease model, including extrapolated vaccination coverage
#' prior to 1980 if available.
#'
#' @param demog_data Output list from `process_demography()`.
#' @param processed_vaccination Data frame of routine vaccine coverage.
#' @param processed_vaccination_sia Data frame of SIA coverage.
#' @param processed_case WHO case reports processed for the country of interest.
#' @param vaccination_schedule Full WHO vaccine schedule data frame.
#' @param vaccination_pre1980 Data frame of historical vaccination introduction assumptions.
#'
#' @return A named list of parameters including coverage arrays and seeding inputs.
#' @import data.table
#' @keywords internal
case_vaccine_to_param <- function(
    demog_data,
    processed_vaccination,
    processed_vaccination_sia,
    processed_case,
    vaccination_schedule,
    vaccination_pre1980
) {

  iso <- demog_data$input_data$iso
  n_age <- demog_data$input_data$n_age
  ages <- 0:(n_age - 1)
  years <- demog_data$input_data$year_start:demog_data$input_data$year_end

  # Use data.table for combined extraction
  vaccination_sources <- c(
    processed_case$disease_description,
    processed_vaccination_sia$vaccination_name,
    processed_vaccination_sia$disease,
    processed_vaccination$vaccine,
    processed_vaccination$vaccine_description
  )
  vaccination_type <- paste(unique(vaccination_sources), collapse = "|")

  if (grepl("Diphtheria|Pertussis", vaccination_type, ignore.case = TRUE)) {
    vaccination_type <- paste0(vaccination_type, "|DTPCV1|DTPCV3|DTaP|DT|DTwP")
  }

  # Expand routine coverage pre-1980
  processed_vaccination <- expand_pre1980_vaccination(
    processed_vaccination,
    vaccination_pre1980,
    disease = tolower(unique(processed_case$disease_description)[1]),
    iso = iso
  )

  # Filter and build vaccination input
  schedule <- filter_vaccine_schedule(vaccination_schedule, vaccination_type, iso)
  routine_df <- build_routine_vaccination_param(processed_vaccination, schedule, ages, years)

  # Build SIA input if needed
  sia_df <- if (nrow(processed_vaccination_sia) > 0) {
    build_sia_vaccination_param(processed_vaccination_sia, ages, years)
  } else {
    NULL
  }

  # Combine and deduplicate
  vacc_df <- combine_vaccination_params(routine_df, sia_df)
  data.table::setDT(vacc_df)

  zero_row <- data.table::data.table(dim1 = 1, dim2 = 1, dim3 = 1, dim4 = 1, value = 0)
  vacc_df <- rbindlist(list(zero_row, vacc_df))
  vacc_df <- vacc_df[, .(value = max(value)), by = .(dim1, dim2, dim3, dim4)]

  # Build seeded case input
  if (nrow(processed_case) > 0) {
    case_df <- build_seeded_case_param(processed_case, demog_data, years, ages)
    tt_seeded <- c(0, match(unique(processed_case$year), years) - 1)
  } else {
    case_df <- data.table::data.table(dim1 = 1, dim2 = 1, dim3 = 1, dim4 = 1, value = 0)
    tt_seeded <- 0
  }

  list(
    tt_vaccination = c(0, match(unique(routine_df$dim4), seq_along(years))),
    vaccination_coverage = vacc_df,
    tt_seeded = tt_seeded,
    seeded = case_df
  )
}
