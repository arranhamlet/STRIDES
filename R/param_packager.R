## Refactored param_packager

#' Package Parameters for Transmission Model
#'
#' This function prepares and formats input parameters for an age-, vaccination-, and risk-structured
#' infectious disease transmission model. It handles user-specified arrays or default values and returns
#' a validated, formatted list ready for simulation.
#'
#' @param n_age Number of age groups
#' @param n_vacc Number of vaccination strata
#' @param n_risk Number of risk groups
#' @param S0 Initial susceptible population (array or scalar)
#' @param I0 Initial infections (array or scalar)
#' @param Rpop0 Initial recovered population (array or scalar)
#' @param contact_matrix Contact matrix by age
#' @param incubation_rate Incubation rate (1/duration)
#' @param recovery_rate Recovery rate for mild cases
#' @param severe_recovery_rate Recovery rate for severe cases
#' @param prop_severe Proportion of cases that are severe
#' @param prop_complications Proportion of severe cases with complications
#' @param natural_immunity_waning Waning rate of natural immunity
#' @param vaccination_coverage Array or DF of vaccination coverage
#' @param tt_vaccination_coverage Time points for vaccination coverage
#' @param age_vaccination_beta_modifier Transmission modifier by age/vaccination
#' @param short_term_waning Short-term waning rate by vaccine stratum
#' @param long_term_waning Long-term waning rate by vaccine stratum
#' @param R0 Time-varying basic reproduction number
#' @param tt_R0 Time points at which R0 changes
#' @param seeded Array or DF of seeded infections
#' @param tt_seeded Time points for seeding
#' @param aging_rate Aging rates between age groups
#' @param simp_birth_death Flag to simplify births and deaths (0/1)
#' @param tt_birth_changes Time points for birth rate changes
#' @param tt_death_changes Time points for death rate changes
#' @param crude_birth Birth rates (array or DF)
#' @param crude_death Death rates (array or DF)
#' @param protection_weight_vacc Weight of maternal protection (vaccine-derived)
#' @param protection_weight_rec Weight of maternal protection (natural infection)
#' @param age_maternal_protection_ends Age at which maternal protection ends
#' @param repro_low Lower bound of reproductive age
#' @param repro_high Upper bound of reproductive age
#' @param tt_migration Time points for migration events
#' @param migration_in_number Number of migrants
#' @param migration_distribution_values Migration destination distribution
#' @param migration_represent_current_pop If 1, distribution reflects current population
#' @param cfr_normal Case fatality rate for normal cases
#' @param cfr_severe Case fatality rate for severe cases
#' @param user_specified_foi Flag to override FOI
#' @param initial_FOI Initial force of infection
#' @param foi_turn_off_when_vaccinating Binary flag to zero FOI when vaccinating
#'
#' @return A named list of structured model inputs
#'
#' @keywords internal
param_packager <- function(
    n_age, n_vacc, n_risk,
    S0, I0 = 0, Rpop0 = 0,
    contact_matrix = NULL,
    incubation_rate, recovery_rate, severe_recovery_rate = NULL,
    prop_severe = 0, prop_complications = 0, natural_immunity_waning = 0,
    vaccination_coverage, tt_vaccination_coverage,
    age_vaccination_beta_modifier = 0,
    short_term_waning = 0, long_term_waning = 0,
    R0, tt_R0,
    seeded, tt_seeded,
    aging_rate, simp_birth_death = 0,
    tt_birth_changes, tt_death_changes, crude_birth, crude_death,
    protection_weight_vacc = 0, protection_weight_rec = 0,
    age_maternal_protection_ends = 1, repro_low = 1, repro_high = NULL,
    tt_migration = 0, migration_in_number = 0, migration_distribution_values = 0,
    migration_represent_current_pop = 0,
    cfr_normal = 0, cfr_severe = 0
) {

  format_array <- function(x, dims) {
    if (length(x) == 1) array(x, dim = dims)
    else do.call(array_from_df, c(setNames(as.list(dims), paste0("dim", seq_along(dims))), list(updates = x)))
  }

  if (is.null(contact_matrix)) {
    contact_matrix <- matrix(1 / (n_age^2), nrow = n_age, ncol = n_age)
  }
  if (is.null(severe_recovery_rate)) {
    severe_recovery_rate <- recovery_rate
  }
  if (is.null(repro_high)) repro_high <- n_age

  dims_3d <- c(n_age, n_vacc, n_risk)
  dims_4d_vac <- c(n_age, n_vacc, n_risk, length(tt_vaccination_coverage))
  dims_4d_seed <- c(n_age, n_vacc, n_risk, length(tt_seeded))

  params <- list(
    n_age = n_age, n_vacc = n_vacc, n_risk = n_risk,
    S0 = format_array(S0, dims_3d),
    I0 = format_array(I0, dims_3d),
    Rpop0 = format_array(Rpop0, dims_3d),
    contact_matrix = contact_matrix,
    incubation_rate = incubation_rate,
    recovery_rate = recovery_rate,
    severe_recovery_rate = severe_recovery_rate,
    prop_severe = format_array(prop_severe, dims_3d),
    prop_complications = format_array(prop_complications, n_age),
    natural_immunity_waning = natural_immunity_waning,
    vaccination_coverage = format_array(vaccination_coverage, dims_4d_vac),
    tt_vaccination_coverage = tt_vaccination_coverage,
    age_vaccination_beta_modifier = format_array(age_vaccination_beta_modifier, dims_3d),
    short_term_waning = format_array(short_term_waning, n_vacc),
    long_term_waning = format_array(long_term_waning, n_vacc),
    R0 = R0,
    tt_R0 = tt_R0,
    seeded = format_array(seeded, dims_4d_seed),
    tt_seeded = tt_seeded,
    aging_rate = {
      x <- format_array(aging_rate, n_age)
      x[n_age] <- 0
      x
    },
    simp_birth_death = simp_birth_death,
    tt_birth_changes = tt_birth_changes,
    tt_death_changes = tt_death_changes,
    crude_birth = format_array(crude_birth, c(n_risk, length(tt_birth_changes))),
    crude_death = format_array(crude_death, c(n_age, n_risk, length(tt_death_changes))),
    protection_weight_vacc = protection_weight_vacc,
    protection_weight_rec = protection_weight_rec,
    age_maternal_protection_ends = age_maternal_protection_ends,
    repro_low = repro_low,
    repro_high = repro_high,
    tt_migration = tt_migration,
    migration_in_number = format_array(migration_in_number, c(n_age, n_vacc, n_risk, length(tt_migration))),
    migration_distribution_values = format_array(migration_distribution_values, c(6, length(tt_migration))),
    migration_represent_current_pop = migration_represent_current_pop,
    cfr_normal = format_array(cfr_normal, n_age),
    cfr_severe = format_array(cfr_severe, n_age)
  )

  return(params)
}
