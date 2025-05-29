
data_load_process_wrapper <- function(
    iso,
    disease,
    vaccine,
    R0,
    timestep = "day",
    year_start = "",
    year_end = "",
    WHO_seed_switch = F
){

  # --------------------------------------------
  # Load Package-Internal Datasets
  # --------------------------------------------
  migration <- PREVAIL::UN_WPP_migration
  fertility <- PREVAIL::UN_WPP_fertility
  mortality <- PREVAIL::UN_WPP_mortality
  population_all <- PREVAIL::UN_WPP_population_total
  population_female <- PREVAIL::UN_WPP_population_female
  contact_matricies <- PREVAIL::contact_matricies
  routine_vaccination_data <- PREVAIL::WHO_vaccination_routine
  sia_vaccination <- PREVAIL::VIMC_vaccination_sia
  full_disease_df <- PREVAIL::WHO_disease_reports
  vaccination_schedule <- PREVAIL::WHO_vaccination_schedule
  vaccination_pre1980 <- PREVAIL::vaccination_pre1980
  vaccine_parameters <- PREVAIL::vaccine_parameters %>% dplyr::filter(disease == !!disease)

  # --------------------------------------------
  # Process Model Input Data
  # --------------------------------------------
  number_of_vaccines <- switch(
    disease,
    "measles" = 2,
    "diphtheria" = 3,
    "pertussis" = 3,
    0
  )

  model_data_preprocessed <- prepare_model_inputs(
    iso = iso,
    disease = disease,
    vaccine = vaccine,
    n_age = 101,
    number_of_vaccines = number_of_vaccines,
    migration = migration,
    fertility = fertility,
    mortality = mortality,
    population_all = population_all,
    population_female = population_female,
    contact_matricies = contact_matricies,
    disease_data = full_disease_df,
    vaccination_data_routine = routine_vaccination_data,
    vaccination_data_sia = sia_vaccination,
    year_start = year_start,
    year_end = year_end
  )

  # --------------------------------------------
  # Integrate Pre-1980 Vaccination Coverage
  # --------------------------------------------
  total_vac <- expand_pre1980_vaccination(
    processed_vaccination = model_data_preprocessed$processed_vaccination_data,
    vaccination_pre1980 = vaccination_pre1980,
    iso = iso,
    disease = disease
  )

  # --------------------------------------------
  # Prepare Case & Vaccination Parameters
  # --------------------------------------------
  case_vaccination_ready <- case_vaccine_to_param(
    demog_data = model_data_preprocessed$processed_demographic_data,
    processed_vaccination = total_vac,
    processed_vaccination_sia = model_data_preprocessed$processed_vaccination_sia,
    processed_case = model_data_preprocessed$processed_case_data,
    vaccination_schedule = vaccination_schedule %>% dplyr::filter(ISO_3_CODE == iso)
  )

  # --------------------------------------------
  # Calculate Vaccine-Derived Protection
  # --------------------------------------------
  n_vacc <- model_data_preprocessed$processed_demographic_data$input_data$n_vacc
  vacc_order <- seq(2, n_vacc, by = 2)

  age_vaccination_beta_modifier <- purrr::map_dfr(vacc_order, function(j) {
    dose_details <- vaccine_parameters %>%
      dplyr::mutate(order = abs(j - dose)) %>%
      dplyr::filter(order == min(order))

    dplyr::bind_rows(
      expand.grid(dim1 = 1:101, dim2 = j,     dim3 = 1, value = dose_details$value[dose_details$parameter == "short_term_protection"]),
      expand.grid(dim1 = 1:101, dim2 = j + 1, dim3 = 1, value = dose_details$value[dose_details$parameter == "long_term_protection"])
    )
  })

  # --------------------------------------------
  # Convert Time Step for Events
  # --------------------------------------------
  time_adjust <- dplyr::case_when(
    timestep == "day"     ~ 1,
    timestep == "week"    ~ 7,
    timestep == "month"   ~ 30,
    timestep == "quarter" ~ 91.25,
    timestep == "year"    ~ 365
  )

  time_changes_mig <- with(model_data_preprocessed$processed_demographic_data, floor(c(tt_migration, max(tt_migration) + 1) * 365 / time_adjust))
  time_changes_vac <- with(case_vaccination_ready, floor(c(tt_vaccination, max(tt_vaccination) + 1) * 365 / time_adjust))
  time_changes_seeded <- floor(case_vaccination_ready$tt_seeded * 365 / time_adjust)
  R0_switch_time <- time_changes_seeded[2]




}
