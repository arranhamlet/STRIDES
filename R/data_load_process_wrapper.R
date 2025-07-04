#' Load and Process All Model Input Data
#'
#' This wrapper function loads internal datasets and processes them into structured model input
#' parameters for a given country, disease, and vaccine. It integrates demographic, vaccination,
#' and disease data, and formats it for use in the dynamic transmission model.
#'
#' @param iso A 3-letter ISO country code.
#' @param disease Name of the disease (e.g., "measles").
#' @param R0 A numeric value or vector specifying the basic reproduction number.
#' @param year_start Optional. Start year for simulation window.
#' @param year_end Optional. End year for simulation window.
#' @param WHO_seed_switch Logical. If `TRUE`, applies custom seeding to replicate data reported to the WHO from 1980.
#' @param aggregate_age Logical. If `TRUE`, aggregates from single-year age groups to custom age groups.
#' @return A named list of structured parameters for the transmission model.
#' @keywords external
#' @export
data_load_process_wrapper <- function(
    iso,
    disease,
    R0,
    year_start = "",
    year_end = "",
    WHO_seed_switch = TRUE,
    aggregate_age = TRUE,
    new_age_breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, Inf)
) {

  # ---- Load Package Data ----
  datasets <- list(
    migration              = PREVAIL::UN_WPP_migration,
    fertility              = PREVAIL::UN_WPP_fertility,
    mortality              = PREVAIL::UN_WPP_mortality,
    population_all         = PREVAIL::UN_WPP_population_total,
    population_female      = PREVAIL::UN_WPP_population_female,
    contact_matricies      = PREVAIL::contact_matricies,
    routine_vaccination    = PREVAIL::WHO_vaccination_routine,
    sia_vaccination        = PREVAIL::VIMC_vaccination_sia,
    disease_data           = PREVAIL::WHO_disease_reports,
    schedule               = PREVAIL::WHO_vaccination_schedule,
    pre1980                = PREVAIL::vaccination_pre1980,
    disease_parameters     = PREVAIL::disease_parameters %>% dplyr::filter(disease == !!disease),
    vaccine_parameters     = PREVAIL::vaccine_parameters %>% dplyr::filter(disease == !!disease)
  )

  number_of_vaccines <- switch(
    disease,
    "measles"    = 2,
    "diphtheria" = 3,
    "pertussis"  = 3,
    0
  )

  # ---- Prepare Inputs ----
  preprocessed <- prepare_model_inputs(
    iso = iso, disease = disease, vaccine = disease,
    n_age = 101, number_of_vaccines = number_of_vaccines,
    migration = datasets$migration, fertility = datasets$fertility, mortality = datasets$mortality,
    population_all = datasets$population_all,
    population_female = datasets$population_female,
    contact_matricies = datasets$contact_matricies,
    disease_data = datasets$disease_data,
    vaccination_data_routine = datasets$routine_vaccination,
    vaccination_data_sia = datasets$sia_vaccination,
    year_start = year_start, year_end = year_end
  )

  # ---- Generate Parameters from Case and Vaccine Data ----
  cv_params <- case_vaccine_to_param(
    demog_data = preprocessed$processed_demographic_data,
    processed_vaccination = preprocessed$processed_vaccination_data,
    processed_vaccination_sia = preprocessed$processed_vaccination_sia,
    processed_case = preprocessed$processed_case_data,
    vaccination_schedule = datasets$schedule %>% dplyr::filter(iso3 == iso),
    vaccination_pre1980 = datasets$pre1980
  )

  # ---- Build Age-Vaccination Modifier Structure ----
  n_vacc <- preprocessed$processed_demographic_data$input_data$n_vacc
  vacc_order <- seq(2, n_vacc, by = 2)

  age_vaccination_beta_modifier <- purrr::map_dfr(vacc_order, function(j) {
    dose_details <- datasets$vaccine_parameters %>%
      dplyr::mutate(order = abs(j - dose)) %>%
      dplyr::filter(order == min(order))

    n_age <- preprocessed$processed_demographic_data$input_data$n_age

    short_term <- base::expand.grid(
      dim1 = seq_len(n_age),
      dim2 = j,
      dim3 = 1,
      value = dose_details$value[dose_details$parameter == "short_term_protection"]
    )

    long_term <- base::expand.grid(
      dim1 = seq_len(n_age),
      dim2 = j + 1,
      dim3 = 1,
      value = dose_details$value[dose_details$parameter == "long_term_protection"]
    )

    dplyr::bind_rows(short_term, long_term)
  })

  # ---- Time Scaling and Change Points ----
  times <- list(
    mig  = base::sort(with(preprocessed$processed_demographic_data, base::floor(c(tt_migration, base::max(tt_migration) + 1) * 365))),
    vac  = base::sort(with(cv_params, base::floor(c(tt_vaccination, base::max(tt_vaccination) + 1) * 365))),
    seed = base::sort(base::floor(cv_params$tt_seeded * 365))
  )

  R0_switch_time <- times$seed[2]

  # ---- Natural Immunity Waning ----
  nat_waning <- datasets$disease_parameters %>%
    dplyr::filter(parameter == "natural immunity waning") %>%
    dplyr::pull(value) %>%
    tidyr::replace_na(0) %>%
    base::gsub("NA", 0, .) %>%
    base::as.numeric() * 365

  # ---- WHO Seeding ----
  if (WHO_seed_switch) {

    # Base seed entries: all except dim4 = 1, and doubled for WHO style
    base_seed <- cv_params$seeded %>%
      dplyr::filter(dim4 != 1) %>%
      dplyr::mutate(dim4 = (dim4 * 2) - 2)

    # Duplicate with value = 0 and dim4 shifted forward
    zero_seed <- base_seed %>%
      dplyr::mutate(value = 0, dim4 = dim4 + 1)

    # Insert "patch" value for initial seeding
    zero_seed <- dplyr::bind_rows(
      zero_seed,
      zero_seed %>%
        dplyr::filter(dim4 == 3) %>%
        dplyr::mutate(dim4 = 1)
    ) %>%
      dplyr::mutate(value = dplyr::case_when(
        dim4 == 1 & dim1 == 1 ~ 10,
        TRUE ~ value
      ))

    seed_data <- dplyr::bind_rows(base_seed, zero_seed) %>%
      dplyr::arrange(dim4)

    # Adjust timepoints to match WHO-style seeding schedule
    original_times <- times$seed
    replicated <- base::unlist(base::lapply(original_times[original_times != 0], function(e) c(e, e + 1)))
    times$seed <- c(0, base::sort(replicated), base::max(replicated) + 364, base::max(replicated) + 365)

  } else {

    # Minimal fallback if WHO seed switch is off
    times$seed <- c(base::min(times$seed), base::max(times$seed) + 1)

    seed_data <- base::data.frame(
      dim1 = 1, dim2 = 1, dim3 = 1,
      dim4 = seq_along(times$seed),
      value = 10
    )
  }

  # ---- Optional Aggregation to New Age Structure ----
  param_inputs <- aggregate_inputs(
    preprocessed = preprocessed,
    cv_params = cv_params,
    seed_data = seed_data,
    new_age_breaks = new_age_breaks,
    aggregate_age = aggregate_age,
    age_vaccination_beta_modifier = age_vaccination_beta_modifier
  )

  inputs <- param_inputs$inputs
  default_inputs <- param_inputs$default_inputs

  # ---- Update Aging Rate and Reproductive Age Bounds ----
  aging_rate <- if (aggregate_age) {
    1 / (365 * base::diff(new_age_breaks))
  } else {
    1 / 365
  }

  if (length(unique(inputs$age_beta_mod$dim1)) == 101) {
    repro_low <- 16
    repro_high <- 50
  } else {
    age_lowers <- new_age_breaks[-base::length(new_age_breaks)]
    age_uppers <- new_age_breaks[-1]
    repro_low  <- base::min(which(age_uppers > 15 & age_lowers < 50))
    repro_high <- base::max(which(age_lowers < 50 & age_uppers > 15))
  }

  # ---- Return Packaged Parameters ----
  packed_params <- param_packager(
    n_age                        = length(unique(inputs$age_beta_mod$dim1)),
    n_vacc                       = n_vacc,
    n_risk                       = preprocessed$processed_demographic_data$input_data$n_risk,
    short_term_waning            = 1 / (max(datasets$vaccine_parameters$value[datasets$vaccine_parameters$parameter == "short_term_waning"]) * 365),
    long_term_waning             = 1 / (max(datasets$vaccine_parameters$value[datasets$vaccine_parameters$parameter == "long_term_waning"]) * 365),
    incubation_rate              = 1 / as.numeric(datasets$disease_parameters$value[datasets$disease_parameters$parameter == "incubation period"]),
    recovery_rate                = 1 / as.numeric(datasets$disease_parameters$value[datasets$disease_parameters$parameter == "infectious period"]),
    severe_recovery_rate         = 1 / as.numeric(datasets$disease_parameters$value[datasets$disease_parameters$parameter == "infectious period"]),
    natural_immunity_waning      = if (nat_waning == 0) 0 else 1 / nat_waning,
    R0                           = if (WHO_seed_switch) c(R0, 0) else R0,
    tt_R0                        = if (WHO_seed_switch) c(0, R0_switch_time) else 0,
    vaccination_coverage         = inputs$vacc_cov,
    contact_matrix               = inputs$contact_matrix,
    age_vaccination_beta_modifier = inputs$age_beta_mod,
    S0                           = inputs$N0,
    Rpop0                        = 0,
    I0                           = 0,
    tt_birth_changes             = times$mig,
    tt_death_changes             = times$mig,
    tt_migration                 = times$mig,
    tt_vaccination_coverage      = times$vac,
    tt_seeded                    = if (WHO_seed_switch) times$seed else c(0, max(times$seed)),
    crude_birth                  = inputs$crude_birth,
    crude_death                  = inputs$crude_death,
    aging_rate                   = aging_rate,
    migration_in_number          = inputs$mig_in,
    migration_distribution_values = inputs$mig_dist,
    seeded                       = inputs$seeded,
    repro_low                    = repro_low,
    repro_high                   = repro_high,
    age_maternal_protection_ends = 1,
    protection_weight_vacc       = 0,
    protection_weight_rec        = 0,
    migration_represent_current_pop = 1,
    population                   = inputs$population,
    female_population            = inputs$female_population,
    new_age_breaks               = new_age_breaks
  )

  # Attach input metadata
  packed_params$input_data <- data.frame(
    iso              = iso,
    disease          = disease,
    R0               = R0,
    year_start       = year_start,
    year_end         = year_end,
    WHO_seed_switch  = WHO_seed_switch,
    aggregate_age    = aggregate_age,
    age_breaks       = if (aggregate_age) {
      paste(new_age_breaks, collapse = ";")
    } else {
      paste(seq(1, 101), collapse = ";")
    }
  )

  packed_params

}
