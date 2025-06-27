#' Load and Process All Model Input Data
#'
#' This wrapper function loads internal datasets and processes them into structured model input
#' parameters for a given country, disease, and vaccine. It integrates demographic, vaccination,
#' and disease data, and formats it for use in the dynamic transmission model.
#'
#' @param iso A 3-letter ISO country code.
#' @param disease Name of the disease (e.g., "measles").
#' @param vaccine Vaccine name or substring to match in schedule and data.
#' @param R0 A numeric value or vector specifying the basic reproduction number.
#' @param timestep One of "day", "week", "month", "quarter", or "year".
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
    vaccine,
    R0,
    timestep = "day",
    year_start = "",
    year_end = "",
    WHO_seed_switch = FALSE,
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
    iso = iso, disease = disease, vaccine = vaccine,
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

    short_term <- expand.grid(
      dim1 = seq_len(preprocessed$processed_demographic_data$input_data$n_age),
      dim2 = j,
      dim3 = 1,
      value = dose_details$value[dose_details$parameter == "short_term_protection"]
    )

    long_term <- expand.grid(
      dim1 = seq_len(preprocessed$processed_demographic_data$input_data$n_age),
      dim2 = j + 1,
      dim3 = 1,
      value = dose_details$value[dose_details$parameter == "long_term_protection"]
    )

    dplyr::bind_rows(short_term, long_term)
  })

  # ---- Time Scaling and Change Points ----
  time_factor <- dplyr::case_when(
    timestep == "day"     ~ 1,
    timestep == "week"    ~ 7,
    timestep == "month"   ~ 30,
    timestep == "quarter" ~ 91.25,
    timestep == "year"    ~ 365
  )

  times <- list(
    mig = sort(with(preprocessed$processed_demographic_data, floor(c(tt_migration, max(tt_migration) + 1) * 365 / time_factor))),
    vac = sort(with(cv_params, floor(c(tt_vaccination, max(tt_vaccination) + 1) * 365 / time_factor))),
    seed = sort(floor(cv_params$tt_seeded * 365 / time_factor))
  )

  R0_switch_time <- times$seed[2]

  # ---- Natural Immunity Waning ----
  nat_waning <- datasets$disease_parameters %>%
    dplyr::filter(parameter == "natural immunity waning") %>%
    dplyr::pull(value) %>%
    tidyr::replace_na(0) %>%
    gsub("NA", 0) %>%
    as.numeric() * 365

  seed_value <- dplyr::case_when(
    timestep == "day"     ~ 1,
    timestep == "week"    ~ 7,
    timestep == "month"   ~ 30,
    timestep == "quarter" ~ 91,
    timestep == "year"    ~ 365
  )

  # ---- WHO Seeding ----
  if (WHO_seed_switch) {

    #Set seeds up
    base_seed <- cv_params$seeded %>%
      subset(dim4 != 1 ) %>%
      mutate(dim4 = (dim4 * 2) - 2)

    zero_seed <- base_seed %>%
      mutate(value = 0) %>%
      mutate(dim4 = dim4 + 1)

    zero_seed <- rbind(zero_seed,
                       zero_seed %>%
                         subset(dim4 == 3) %>%
                         mutate(dim4 = 1)) %>%
      mutate(value = case_when(
        dim4 == 1 & dim1 == 18 ~ 10,
        TRUE ~ value
      ))

    seed_data <- rbind(base_seed, zero_seed) %>%
      arrange(dim4)

    #Set times up
    original_times <- times$seed
    new_times <- c(0, sort(unlist(lapply(original_times[which(original_times != 0)], \(e) c(e, e + 1)))))
    new_times <- c(new_times, max(new_times + 364), max(new_times + 365))
    times$seed <- new_times

  } else {

    # Default fallback if no WHO seeding
    times$seed <- c(min(times$seed), max(times$seed) + 1)

    seed_data <- data.frame(
      dim1 = 18, dim2 = 1, dim3 = 1,
      dim4 = 1:length(times$seed),
      value = seed_value
    )
  }

  # ---- Optional Aggregation to New Age Structure ----
  n_age <- preprocessed$processed_demographic_data$input_data$n_age
  aging_rate <- preprocessed$processed_demographic_data$aging_rate

  weight_reformatted <- preprocessed$processed_demographic_data$population_data %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "time") %>%
    tidyr::pivot_longer(names_to = "age", values_to = "value", cols = -c("time")) %>%
    group_by(time) %>%
    mutate(time = as.integer(time),
           age = as.integer(age))

  default_inputs <- list(
    age_beta_mod   = age_vaccination_beta_modifier,
    vacc_cov       = cv_params$vaccination_coverage,
    N0             = preprocessed$processed_demographic_data$N0,
    crude_death    = preprocessed$processed_demographic_data$crude_death %>% dplyr::mutate(value = value / (365 / time_factor)),
    crude_birth    = preprocessed$processed_demographic_data$crude_birth %>%
      mutate(value = value / (365 / time_factor)),
    mig_in         = preprocessed$processed_demographic_data$migration_in_number %>% dplyr::mutate(value = value / (365 / time_factor)),
    mig_dist       = preprocessed$processed_demographic_data$migration_distribution_values,
    seeded         = seed_data,
    contact_matrix = preprocessed$processed_demographic_data$contact_matrix,
    population = weight_reformatted %>% dplyr::rename("dim1" = "age",
                                                      "dim2" = "time")
  )

  if(aggregate_age) {

    aging_vector <- diff(new_age_breaks) / (365 / time_factor)
    aging_vector[is.infinite(aging_vector)] <- 0

    aggregate_inputs <- function(obj, method = "weighted.mean", weights = NULL, time_var = "dim4") {
      wts <- if (method == "sum") NULL else weights
      aggregate_age_structure(obj, age_breaks = new_age_breaks, method = method, weights = wts, time_var = time_var)
    }

    n_age = length(new_age_breaks) - 1

    death_upd <- preprocessed$processed_demographic_data$crude_death %>% dplyr::mutate(value = value/(365/time_factor))
    death_upd$value <- death_upd$value * weight_reformatted$value

    # ---- Aggregate Crude Birth with Population Weighting ----

    # Step 1: Add fertility = crude_birth for ages 15–49 only
    fertility_df <- preprocessed$processed_demographic_data$crude_birth %>%
      dplyr::mutate(value = value / (365 / time_factor)) %>%
      dplyr::left_join(weight_reformatted %>%
                         mutate(time = as.character(time)), by = c("dim2" = "time", "dim1" = "age")) %>%
      dplyr::mutate(
        fertility = dplyr::if_else(dim1 >= 15 & dim1 <= 49, value.x, 0),
        births = fertility * value.y
      ) %>%
      select(dim1, dim2, value = births)

    inputs <- list(
      age_beta_mod   = aggregate_inputs(age_vaccination_beta_modifier, weights = weight_reformatted, time_var = NULL),
      vacc_cov       = aggregate_inputs(cv_params$vaccination_coverage, weights = weight_reformatted),
      N0             = aggregate_inputs(preprocessed$processed_demographic_data$N0, method = "sum"),
      crude_death    = aggregate_inputs(death_upd, method = "sum"),
      crude_birth    = aggregate_inputs(fertility_df, method = "sum", time_var = "dim2"),
      mig_in         = aggregate_inputs(preprocessed$processed_demographic_data$migration_in_number %>% dplyr::mutate(value = value / (365 / time_factor)), method = "sum"),
      mig_dist       = aggregate_inputs(preprocessed$processed_demographic_data$migration_distribution_values, weights = weight_reformatted, time_var = "dim2"),
      seeded         = aggregate_inputs(seed_data, method = "sum"),
      contact_matrix = aggregate_contact_matrix(preprocessed$processed_demographic_data$contact_matrix, age_breaks = new_age_breaks, method = "sum"),
      population = aggregate_inputs(weight_reformatted %>% dplyr::rename("dim1" = "age",
                                                                  "dim2" = "time"),
                                    time_var = "dim2", method = "sum")
    )

    #Crude birth updated
    crude_birth_updated <- inputs$crude_birth %>%
      mutate(dim2 = as.numeric(dim2),
             dim1 = as.numeric(dim1)) %>%
      arrange(dim1, dim2) %>%
      left_join(
        inputs$population %>%
          select(dim1, dim2, population = value),
        by = c("dim1", "dim2")
      ) %>%
      mutate(value = value/population)

    inputs$crude_birth <- crude_birth_updated
    inputs$crude_death$value <- inputs$crude_death$value / inputs$population$value

  } else {
    inputs <- default_inputs
  }

  if (aggregate_age) {
    aging_vector <- diff(new_age_breaks)
    aging_rate <- time_factor / (365 * aging_vector)
  } else {
    aging_rate <- time_factor / 365
  }

  if (n_age == 101) {
    repro_low <- 16
    repro_high <- 50
  } else {
    age_lowers <- new_age_breaks[-length(new_age_breaks)]
    age_uppers <- new_age_breaks[-1]
    repro_low  <- min(which(age_uppers > 15 & age_lowers < 50))
    repro_high <- max(which(age_lowers < 50 & age_uppers > 15))
  }


  # ---- Return Packaged Parameters ----
  packed_params <- param_packager(
    n_age = n_age,
    n_vacc = n_vacc,
    n_risk = preprocessed$processed_demographic_data$input_data$n_risk,
    short_term_waning = 1 / (max(datasets$vaccine_parameters$value[datasets$vaccine_parameters$parameter == "short_term_waning"]) * 365 / time_factor),
    long_term_waning = 1 / (max(datasets$vaccine_parameters$value[datasets$vaccine_parameters$parameter == "long_term_waning"]) * 365 / time_factor),
    age_vaccination_beta_modifier = inputs$age_beta_mod,
    R0 = if (WHO_seed_switch) c(R0, 0) else R0,
    tt_R0 = if (WHO_seed_switch) c(0, R0_switch_time) else 0,
    cfr_normal = 0,
    cfr_severe = 0,
    incubation_rate = 1 / (as.numeric(datasets$disease_parameters$value[datasets$disease_parameters$parameter == "incubation period"]) * time_factor),
    recovery_rate = 1 / (as.numeric(datasets$disease_parameters$value[datasets$disease_parameters$parameter == "infectious period"]) * time_factor),
    severe_recovery_rate = 1 / (as.numeric(datasets$disease_parameters$value[datasets$disease_parameters$parameter == "infectious period"]) * time_factor),
    natural_immunity_waning = ifelse(nat_waning == 0, 0, 1 / nat_waning * time_factor),
    vaccination_coverage = inputs$vacc_cov,
    contact_matrix = inputs$contact_matrix,
    S0 = inputs$N0,
    Rpop0 = 0,
    I0 = 0,
    tt_birth_changes = times$mig,
    tt_death_changes = times$mig,
    tt_migration = times$mig,
    tt_vaccination_coverage = times$vac,
    crude_birth = inputs$crude_birth,
    crude_death = inputs$crude_death,
    aging_rate = aging_rate,
    migration_in_number = inputs$mig_in,
    migration_distribution_values = inputs$mig_dist,
    tt_seeded = if (WHO_seed_switch) times$seed else c(0, max(times$seed)),
    seeded = inputs$seeded,
    repro_low = repro_low,
    repro_high = repro_high,
    age_maternal_protection_ends = 1,
    protection_weight_vacc = 0,
    protection_weight_rec = 0,
    migration_represent_current_pop = 1,
    population = inputs$population,
    new_age_breaks = new_age_breaks
  )

  packed_params$input_data <- data.frame(
    iso,
    disease,
    vaccine,
    R0,
    timestep ,
    year_start,
    year_end,
    WHO_seed_switch,
    aggregate_age
  )

  packed_params

}
