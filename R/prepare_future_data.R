#' Prepare Future Scenario Parameters for Transmission Model
#'
#' This function generates an updated parameter list for projecting future scenarios in a transmission model.
#' It modifies the vaccination coverage, R0 values, and seeding of new cases based on a user-specified set
#' of future events. It also resets migration, birth, and death inputs to their final known values.
#'
#' @param params A named list of model parameters as used in the core transmission model.
#' @param future_events A `data.frame` specifying future event timepoints and modifications. Must contain:
#'   \itemize{
#'     \item \code{year}: Future year (numeric, in calendar years).
#'     \item \code{relative_coverage}: Multiplicative modifier to vaccination coverage.
#'     \item \code{R0}: Basic reproduction number for each future year.
#'     \item \code{introduced_cases}: Number of seed infections to introduce at each timepoint.
#'   }
#' @param current_susceptibility A named list with updated `S0`, `Rpop0` values reflecting the final simulation state.
#'
#' @return A modified version of the `params` list ready to run forward simulations under future scenarios.
#'
#' @export
prepare_future_data <- function(params, future_events, current_susceptibility) {

  ## --- Initialize ---
  new_params <- params
  n_future <- nrow(future_events)
  future_event_tt <- future_events$year * 365
  new_params$input_data$year_start <- if(new_params$input_data$year_start == "") 2024 else new_params$input_data$year_start

  ## --- Extend Vaccination Coverage ---
  prior_vacc <- params$vaccination_coverage
  final_vacc <- prior_vacc[, , , last(dim(prior_vacc))]
  dim(final_vacc) <- c(dim(prior_vacc)[1:3], 1)

  vacc_array <- abind::abind(
    replicate(n = n_future, expr = final_vacc, simplify = FALSE),
    along = 4
  )

  for (i in seq_len(n_future)) {
    vacc_array[, , , i] <- vacc_array[, , , i] * future_events$relative_coverage[i]
  }

  new_params$vaccination_coverage <- vacc_array
  new_params$tt_vaccination_coverage <- future_event_tt
  new_params$no_vacc_changes <- n_future

  ## --- Update R0 ---
  new_params$R0 <- future_events$R0
  new_params$tt_R0 <- future_event_tt
  new_params$no_R0_changes <- n_future

  ## --- Disable Migration ---
  new_params$tt_migration <- 0
  new_params$no_migration_changes <- 1

  final_mig <- params$migration_in_number[, , , last(dim(params$migration_in_number)) - 1]
  final_mig[] <- 0
  dim(final_mig) <- c(dim(params$migration_in_number)[1:3], 1)
  new_params$migration_in_number <- final_mig

  final_mig_dist <- params$migration_distribution_values[, last(dim(params$migration_distribution_values)) - 1]
  final_mig_dist[] <- 0
  dim(final_mig_dist) <- c(dim(params$migration_distribution_values)[1], 1)
  new_params$migration_distribution_values <- final_mig_dist

  ## --- Lock Birth and Death Rates to Final Values ---
  final_birth <- params$crude_birth[, last(dim(params$crude_birth)) - 1]
  dim(final_birth) <- c(dim(params$crude_birth)[1], 1)

  final_death <- params$crude_death[, , last(dim(params$crude_death)) - 1]
  dim(final_death) <- c(dim(params$crude_death)[1:2], 1)

  new_params$crude_birth <- final_birth
  new_params$crude_death <- final_death
  new_params$tt_birth_changes <- 0
  new_params$no_birth_changes <- 1
  new_params$tt_death_changes <- 0
  new_params$no_death_changes <- 1

  ## --- Seed Future Infections by Population Weight ---
  final_pop <- params$population[, ncol(params$population) - 1]
  pop_weights <- final_pop / sum(final_pop)

  seed_df <- Reduce(rbind, lapply(seq_len(n_future), function(i) {
    n_cases <- future_events$introduced_cases[i]
    if (n_cases == 0) return(NULL)
    sampled_ages <- sample(seq_along(pop_weights), n_cases, prob = pop_weights)
    data.frame(dim1 = sampled_ages, dim2 = 1, dim3 = 1, dim4 = i, value = 1)
  }))

  new_params$seeded <- generate_array_df(
    dim1 = params$n_age,
    dim2 = params$n_vacc,
    dim3 = params$n_risk,
    dim4 = n_future,
    updates = seed_df
  ) %>% df_to_array()

  new_params$tt_seeded <- future_event_tt
  new_params$no_seeded_changes <- n_future

  ## --- Set Updated Initial Conditions ---
  new_params$S0 <- current_susceptibility$S0
  new_params$Rpop0 <- current_susceptibility$Rpop0
  new_params$I0[] <- 0

  ## --- Return Future Scenario Parameters ---
  return(new_params)
}
