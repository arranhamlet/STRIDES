% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unpack_dust2.R
\name{unpack_dust2}
\alias{unpack_dust2}
\title{Unpack Dust Model State Data}
\usage{
unpack_dust2(
  model_system,
  model_object,
  dimension_names,
  which_state_dimensions
)
}
\arguments{
\item{model_system}{A `dust2` system object with packing info and particle count.}

\item{model_object}{A model output object returned from `dust2::dust_system_simulate()`.}

\item{dimension_names}{A named list of dimension label vectors (e.g., `time`, `age`, `vaccination`, etc.).}

\item{which_state_dimensions}{A named list mapping state names to dimension label sets.}
}
\value{
A `data.table` with unpacked model output, including `value`, `state`, `time`, and relevant stratifiers.
}
\description{
Extracts and unpacks state data from a dust model object, handling both compartmental and non-compartmental states.
Compartmental states are unpacked via `process_obj()` and include time and stratifier dimensions;
non-compartmental states (e.g., derived metrics) are returned as flat records.
}
