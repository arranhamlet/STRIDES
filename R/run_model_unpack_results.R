#' Run and Unpack Dust Model Simulation Results
#'
#' Initializes and runs a Dust model using a provided odin model and parameter set, then unpacks the simulation
#' output into a long-format `data.table` using `unpack_dust2()`. Suitable for analysis and plotting.
#'
#' @param params A named list of model parameters including `n_age`, `n_vacc`, `n_risk`, and `tt_migration`.
#' @param no_runs Integer. Number of stochastic particles to simulate (default is 10).
#'
#' @return A `data.table` in long format containing simulation results by state, time, and strata.
#'
#' @details
#' The function uses `dust2::dust_system_create()` to instantiate and run the simulation, and `unpack_dust2()`
#' to reshape the result. The output includes states such as `S`, `E`, `I`, `R`, `Is`, `Rc`, `new_case`, `Reff_age`, etc.
#'
#' @export
run_model_unpack_results <- function(params, no_runs = 10) {

  # Create and initialize the Dust system
  sys <- dust2::dust_system_create(transmission_model(), params, n_particles = no_runs)
  dust2::dust_system_set_state_initial(sys)

  # Define simulation time vector (starts at 0)
  full_time_vector <- 0:(ncol(params$population) * 365 - 1)

  # Run simulation
  y <- dust2::dust_system_simulate(sys, full_time_vector)

  # Unpack model outputs into long format
  clean_df <- unpack_dust2(
    model_system = sys,
    model_object = y,
    dimension_names = list(
      age = list(as.character(1:params$n_age)),
      vaccination = list(as.character(1:params$n_vacc)),
      risk = list(as.character(1:params$n_risk)),
      time = list(full_time_vector),
      no_migration_changes = list(params$tt_migration + 1)
    ),
    which_state_dimensions = list(
      S = c("age", "vaccination", "risk", "time"),
      E = c("age", "vaccination", "risk", "time"),
      I = c("age", "vaccination", "risk", "time"),
      R = c("age", "vaccination", "risk", "time"),
      Is = c("age", "vaccination", "risk", "time"),
      Rc = c("age", "vaccination", "risk", "time"),
      Reff_age = c("age", "time"),
      new_case = c("age", "vaccination", "risk", "time"),
      seropositive = c("age", "time")
    )
  )

  return(clean_df)
}
