% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_model.R
\name{run_model}
\alias{run_model}
\title{Run Dust Model Simulation}
\usage{
run_model(odin_model, params, time = 1000, no_runs = 10)
}
\arguments{
\item{odin_model}{An odin model function (not yet evaluated) used to construct the dust system.}

\item{params}{A named list of model parameters, including initial conditions, rates, and dimensions.}

\item{time}{Integer. The number of time steps to simulate (default is 1000). Time starts at 0.}

\item{no_runs}{Integer. Number of stochastic particles (simulations) to run in parallel (default is 10).}
}
\value{
A `dust_system` object containing the full system state and particle trajectories.
  Use `unpack_dust2()` or similar functions to extract outputs.
}
\description{
Initializes a dust model system using the provided odin model and parameters, runs the simulation
over a defined time period with a specified number of stochastic particles, and returns the simulated system.
}
