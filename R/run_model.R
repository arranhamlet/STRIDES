#' Run Dust Model Simulation
#'
#' Initializes a dust model system using the provided odin model and parameters, runs the simulation
#' over a defined time period with a specified number of stochastic particles, and returns the simulated system.
#'
#' @param odin_model An odin model function (not yet evaluated) used to construct the dust system.
#' @param params A named list of model parameters, including initial conditions, rates, and dimensions.
#' @param time Integer. The number of time steps to simulate (default is 1000). Time starts at 0.
#' @param no_runs Integer. Number of stochastic particles (simulations) to run in parallel (default is 10).
#'
#' @return A `dust_system` object containing the full system state and particle trajectories.
#'   Use `unpack_dust2()` or similar functions to extract outputs.
#'
#' @importFrom dust2 dust_system_create dust_system_set_state_initial dust_system_simulate
#' @export
run_model <- function(odin_model, params, time = 1000, no_runs = 10) {

  # Define and initialize the dust system
  sys <- dust2::dust_system_create(odin_model(), params, n_particles = no_runs)

  # Set the system's initial state
  dust2::dust_system_set_state_initial(sys)

  # Define time vector (starts at 0)
  full_time_vector <- 0:(time - 1)

  # Simulate the model
  dust2::dust_system_simulate(sys, full_time_vector)

  # Return the system object
  return(sys)
}
