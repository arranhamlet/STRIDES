#' Run Dust Model Simulation
#'
#' Initializes a dust model system using the provided odin model and parameters, runs the simulation
#' over a defined time period with a specified number of stochastic particles, and returns the simulated system.
#'
#' @param odin_model An odin model function (not evaluated) used to construct the dust system.
#' @param params A named list of model parameters, including initial conditions, rates, and dimensions.
#' @param time Integer. The number of time steps to simulate (default is 1000). Time starts at 0.
#' @param no_runs Integer. Number of particles (stochastic simulations) to run in parallel (default is 10).
#'
#' @return A `dust_system` object after simulation, which contains particle trajectories and system state.
#'         To extract results, use `unpack_dust2()` or a compatible output processing function.
#'
#' @export

run_model <- function(odin_model, params, time = 1000, no_runs = 10) {

  #Define the dust system and initialize it with given parameters
  sys <- dust2::dust_system_create(odin_model(), params, n_particles = no_runs)

  #Set the initial state for the dust system
  dust2::dust_system_set_state_initial(sys)

  #Define the time vector for simulation (starting from 0)
  full_time_vector <- 0:(time - 1)

  #Run the dust system simulation over the defined time period
  dust2::dust_system_simulate(sys, full_time_vector)

}
