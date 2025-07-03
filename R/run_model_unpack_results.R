#' Run and Unpack Dust Model Simulation Results
#'
#' Initializes a Dust model simulation using an `odin.dust` model (compiled with `transmission_model()`),
#' runs the simulation using specified parameters and particle count, and returns the unpacked results
#' as a long-format `data.table` suitable for downstream analysis.
#'
#' @param params A named list of model parameters, which must include:
#'   \itemize{
#'     \item \code{n_age}: Integer. Number of age groups.
#'     \item \code{n_vacc}: Integer. Number of vaccination strata.
#'     \item \code{n_risk}: Integer. Number of risk groups.
#'     \item \code{tt_migration}: Integer. Number of timepoints with migration changes (used for unpacking migration dimension).
#'     \item Other model-specific parameters required by `transmission_model()`.
#'   }
#' @param simulation_length Integer. Total length of the simulation in model timesteps. Defines the number of steps run.
#' @param no_runs Integer. Number of stochastic particles to simulate (default = 10). Higher values yield more stable outputs.
#'
#' @return A `data.table` in long format with the following columns:
#'   \itemize{
#'     \item \code{state}: The name of the state variable (e.g., "S", "E", "I", etc.).
#'     \item \code{time}: Model time step (starting at 0).
#'     \item \code{age}, \code{vaccination}, \code{risk}: Stratification dimensions.
#'     \item \code{value}: Simulated value for the given state and strata.
#'     \item (Optional) \code{particle}: Particle identifier, if applicable.
#'   }
#'
#' @details
#' Internally uses:
#' \itemize{
#'   \item `dust2::dust_system_create()` to build the Dust model system.
#'   \item `dust2::dust_system_simulate()` to run the stochastic simulation.
#'   \item `unpack_dust2()` to reshape model outputs into long-format `data.table`.
#' }
#'
#' The `dimension_names` list passed to `unpack_dust2()` includes:
#' \itemize{
#'   \item Age groups: 1 to `n_age`
#'   \item Vaccination strata: 1 to `n_vacc`
#'   \item Risk strata: 1 to `n_risk`
#'   \item Time: from 0 to `simulation_length`
#'   \item no_migration_changes: derived as `tt_migration + 1`
#' }
#'
#' State dimensions unpacked include: `"S"`, `"E"`, `"I"`, `"R"`, `"Is"`, `"Rc"`, `"new_case"`, `"Reff_age"`, and `"seropositive"`.
#'
#' @examples
#' \dontrun{
#' params <- list(n_age = 17, n_vacc = 3, n_risk = 2, tt_migration = 5, ...)
#' sim_results <- run_model_unpack_results(params, simulation_length = 365)
#' }
#'
#' @export

run_model_unpack_results <- function(params, simulation_length, no_runs = 10) {

  # Create and initialize the Dust system
  sys <- dust2::dust_system_create(transmission_model(), params, n_particles = no_runs)
  dust2::dust_system_set_state_initial(sys)

  # Define simulation time vector (starts at 0)
  full_time_vector <- 0:simulation_length

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
  ) %>%
    dplyr::filter(time <= max(full_time_vector))

  return(clean_df)
}
