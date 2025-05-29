
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
  disease_parameters <- PREVAIL::disease_parameters %>% dplyr::filter(disease == !!disease)
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
  # Prepare Case & Vaccination Parameters
  # --------------------------------------------
  case_vaccination_ready <- case_vaccine_to_param(
    demog_data = model_data_preprocessed$processed_demographic_data,
    processed_vaccination = model_data_preprocessed$processed_vaccination_data,
    processed_vaccination_sia = model_data_preprocessed$processed_vaccination_sia,
    processed_case = model_data_preprocessed$processed_case_data,
    vaccination_schedule = vaccination_schedule %>% dplyr::filter(iso3 == iso),
    vaccination_pre1980 = vaccination_pre1980
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

  waning_immunity <- disease_parameters %>%
    mutate(value = replace_na(0)) %>%
    subset(parameter == "natural immunity waning") %>% pull(value) %>% as.numeric() * 365

  if(WHO_seed_switch == T){

    time_changes_seeded <- c(sapply(time_changes_seeded, function(e){
      c(e, e + 1)
    })) %>%
      floor

    seeded_WHO <- case_vaccination_ready$seeded %>%
      subset(dim4 != 1) %>%
      mutate(dim4 = case_when(
        dim4 == 2 ~ dim4 + 1,
        dim4 != 2 ~ dim4 * 2 - 1
      ))

    zero_data <- seeded_WHO %>%
      mutate(dim4 = dim4 - 1,
             value = 0)

    seed <- rbind(seeded_WHO,
                  rbind(subset(zero_data, dim4 == 2) %>%
                          mutate(dim4 = 1),
                        zero_data)) %>%
      arrange(dim4) %>%
      mutate(
        value = case_when(
          dim4 == 2 & value == 0 ~ 10,
          TRUE ~ value
        )
      )

  } else{
    time_changes_seeded <- sort(floor(c(time_changes_seeded, max(time_changes_seeded) + 1)))
  }


  #Finalise parameters
  params <- param_packager(

    # Demographic parameters
    n_age = model_data_preprocessed$processed_demographic_data$input_data$n_age,
    n_vacc = model_data_preprocessed$processed_demographic_data$input_data$n_vacc,
    n_risk = model_data_preprocessed$processed_demographic_data$input_data$n_risk,

    # Vaccine parameters
    short_term_waning = 1/((vaccine_parameters %>% subset(parameter == "short_term_waning") %>% pull(value) %>% as.numeric() %>% max() * 365)/time_adjust),
    long_term_waning = 1/((vaccine_parameters %>% subset(parameter == "long_term_waning") %>% pull(value) %>% as.numeric() %>% max() * 365)/time_adjust),
    age_vaccination_beta_modifier = age_vaccination_beta_modifier,

    # Disease parameters
    R0 = if(WHO_seed_switch == T) c(R0, 0) else R0,
    tt_R0 = if(WHO_seed_switch == T) c(0, R0_switch_time) else 0,

    #Disease parameters
    cfr_normal = 0,
    cfr_severe = 0,
    incubation_rate = 1/subset(disease_parameters, parameter == "incubation period") %>% pull(value) %>% as.numeric() * time_adjust,
    recovery_rate = 1/subset(disease_parameters, parameter == "infectious period") %>% pull(value) %>% as.numeric() * time_adjust,
    severe_recovery_rate = 1/subset(disease_parameters, parameter == "infectious period") %>% pull(value) %>% as.numeric() * time_adjust,

    natural_immunity_waning = if(waning_immunity == 0) 0 else  1/waning_immunity * time_adjust,

    #Setting up vaccination
    vaccination_coverage = case_vaccination_ready$vaccination_coverage,

    #Demographic parameters
    contact_matrix = model_data_preprocessed$processed_demographic_data$contact_matrix,
    S0 = model_data_preprocessed$processed_demographic_data$N0,
    Rpop0 = 0,
    I0 = 0,

    #Time of changes
    tt_birth_changes = time_changes_mig,
    tt_death_changes = time_changes_mig,
    tt_migration = time_changes_mig,
    tt_vaccination_coverage = time_changes_vac,

    #List of when birth_death_changes
    crude_birth = model_data_preprocessed$processed_demographic_data$crude_birth %>%
      mutate(value = value/(365/time_adjust)),
    crude_death = model_data_preprocessed$processed_demographic_data$crude_death %>%
      mutate(value = value/(365/time_adjust)),
    aging_rate = time_adjust/365,
    migration_in_number = model_data_preprocessed$processed_demographic_data$migration_in_number %>%
      mutate(value = value/(365/time_adjust)),
    migration_distribution_values = model_data_preprocessed$processed_demographic_data$migration_distribution_values,

    tt_seeded = if(WHO_seed_switch == T) time_changes_seeded else c(0, max(time_changes_seeded)),
    seeded = if(WHO_seed_switch == T) seed else expand.grid(dim1 = 18, dim2 = 1, dim3 = 1, dim4 = 1, dim5 = 1:2, value = 10),

    #Birth ages
    repro_low = 15,
    repro_high = 49,
    age_maternal_protection_ends = 1,
    protection_weight_vacc = 1,
    protection_weight_rec = 1,
    migration_represent_current_pop = 1

  )


}
