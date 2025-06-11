#' Run and Unpack Dust Model Simulation Results
#'
#' This function initializes a dust model using an odin model and provided parameters, runs the simulation
#' for a specified time and number of stochastic particles, and returns a long-format data frame of the
#' unpacked and labeled model state outputs using the `unpack_dust2()` helper.
#'
#' @param odin_model An unevaluated odin model function used to construct the dust system.
#' @param params A named list of model parameters including dimensions (`n_age`, `n_vacc`, `n_risk`),
#'   time points (`tt_migration`), and any other inputs required by the model.
#' @param time Integer. Total number of time steps to simulate (default is 1000). Time starts at 0.
#' @param no_runs Integer. Number of stochastic particles to simulate (default is 10).
#'
#' @return A `data.table` in long format containing the unpacked and labeled simulation results across
#'   dimensions such as time, age, vaccination stratum, risk group, and model state.
#'
#' @details
#' The function uses `dust2::dust_system_create()` to create and simulate the model, and then formats
#' the output using `unpack_dust2()` to return results ready for visualization or analysis.
#'
#' @export

run_model_unpack_results <- function(odin_model, params, time = 1000, no_runs = 10) {

  # Define the dust system and initialize it with given parameters
  sys <- dust2::dust_system_create(odin_model(), params, n_particles = no_runs)

  # Set the initial state for the dust system
  dust2::dust_system_set_state_initial(sys)

  # Define the time vector for simulation (starting from 0)
  full_time_vector <- 0:(time - 1)

  # Run the dust system simulation over the defined time period
  y <- dust2::dust_system_simulate(sys, full_time_vector)

  # Process and clean output data by unpacking the results
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
      Reff_age_prop = c("age", "time"),
      Reff_age = c("age", "time"),
      new_case = c("age", "vaccination", "risk", "time")
    )
  )

  clean_df

}
