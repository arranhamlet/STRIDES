#' Process Dust Model Data
#'
#' Transforms a multi-dimensional dust model output array into a long-format data.table,
#' assigning appropriate names to dimensions and unpacking by particle if needed.
#'
#' @param this_obj The model output array to process.
#' @param present_dimensions Character vector of dimension labels for this state (e.g., `c("age", "vaccination", "risk", "time")`).
#' @param colnames Column names for the output data table.
#' @param time_length Length of the time dimension.
#' @param x Index of the current state in the model.
#' @param model_system The model system object containing metadata like number of particles.
#' @param dimension_names Named list of dimension label values for each dimension.
#' @param dust_state List of all unpacked model state arrays.
#'
#' @return A list of two `data.table`s: the full stratified long-format output and a reduced time-by-state (and optionally run) aggregation.
#'
#' @importFrom data.table as.data.table setnames setcolorder
#' @keywords internal
process_obj <- function(this_obj, present_dimensions, colnames, time_length, x, model_system, dimension_names, dust_state) {
  state_name <- names(dust_state)[x]
  has_particles <- model_system$n_particles > 1
  obj_dims <- length(dim(this_obj))
  expect_dims <- length(present_dimensions)

  # Create dimension names
  dim_names <- vector("list", obj_dims)
  for (i in seq_along(present_dimensions)) {
    dim_label <- present_dimensions[i]
    dim_names[[i]] <- dimension_names[[dim_label]][[1]]
  }

  # Add particle names if needed
  if (has_particles || obj_dims != expect_dims) {
    particle_dim_index <- obj_dims - 1
    dim_names[[particle_dim_index]] <- paste0("run_", seq_len(model_system$n_particles))
  }

  dimnames(this_obj) <- dim_names

  # Melt using base R and convert to data.table
  melted_df <- data.table::as.data.table(as.data.frame.table(this_obj, responseName = "value"))

  if (has_particles || obj_dims != expect_dims) {
    # For multi-particle or higher-dimensional models
    data.table::setnames(melted_df, c(colnames[1:(which(colnames == "time") - 1)], "run", "time", "value"))
    melted_df[, state := state_name]
    data.table::setcolorder(melted_df, c("time", "state", colnames[1:(which(colnames == "time") - 1)], "run", "value"))
    aggregate_df <- melted_df[, .(value = sum(value)), by = .(state, time, run)]
  } else {
    # For single-particle models with expected dimensions
    data.table::setnames(melted_df, c(present_dimensions, "value"))
    melted_df[, state := state_name]
    data.table::setcolorder(melted_df, c("time", "state", setdiff(present_dimensions, "time"), "value"))
    aggregate_df <- melted_df[, .(value = sum(value)), by = .(state, time)]
  }

  list(melted_df, aggregate_df)
}
