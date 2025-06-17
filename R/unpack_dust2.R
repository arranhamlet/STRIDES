#' Unpack Dust Model State Data
#'
#' Extracts and unpacks state data from a dust model object, handling both compartmental and non-compartmental states.
#' Compartmental states are unpacked via `process_obj()` and include time and stratifier dimensions;
#' non-compartmental states (e.g., derived metrics) are returned as flat records.
#'
#' @param model_system A `dust2` system object with packing info and particle count.
#' @param model_object A model output object returned from `dust2::dust_system_simulate()`.
#' @param dimension_names A named list of dimension label vectors (e.g., `time`, `age`, `vaccination`, etc.).
#' @param which_state_dimensions A named list mapping state names to dimension label sets.
#'
#' @return A `data.table` with unpacked model output, including `value`, `state`, `time`, and relevant stratifiers.
#'
#' @importFrom data.table data.table rbindlist fifelse setDT
#' @importFrom dust2 dust_unpack_state
#' @export
unpack_dust2 <- function(model_system, model_object, dimension_names, which_state_dimensions) {

  dust_state <- dust2::dust_unpack_state(model_system, model_object)
  is_compartment <- lengths(model_system$packing_state) != 0
  time_vector <- dimension_names$time[[1]]
  n_particles <- model_system$n_particles

  # --- Compartmental states ---
  unpacked_compartments <- data.table::rbindlist(
    lapply(which(is_compartment), function(x) {
      this_obj <- dust_state[[x]]
      state_name <- names(dust_state)[x]
      dims <- which_state_dimensions[[state_name]]
      processed <- process_obj(this_obj, dims, dims, time_vector, x, model_system, dimension_names, dust_state)
      data.table::rbindlist(processed, fill = TRUE)
    }), fill = TRUE
  )

  # --- Non-compartmental states ---
  unpacked_noncompartments <- data.table::rbindlist(
    lapply(which(!is_compartment), function(x) {
      this_obj <- dust_state[[x]]
      state_name <- names(dust_state)[x]

      if (n_particles > 1) {
        data.table::data.table(
          run = paste0("run_", seq_len(n_particles)),
          value = as.numeric(this_obj),
          time = time_vector,
          state = state_name
        )
      } else {
        data.table::data.table(
          value = as.numeric(this_obj),
          time = time_vector,
          state = state_name
        )
      }
    }), fill = TRUE
  )

  # --- Combine ---
  combined <- data.table::rbindlist(
    list(unpacked_compartments, unpacked_noncompartments), fill = TRUE
  )

  # --- Fill stratifier columns ---
  stratifiers <- setdiff(unique(unlist(which_state_dimensions)), "time")
  for (col in stratifiers) {
    if (!col %in% names(combined)) {
      combined[[col]] <- "All"
    } else {
      combined[[col]] <- as.character(combined[[col]])
      combined[[col]] <- data.table::fifelse(is.na(combined[[col]]), "All", combined[[col]])
    }
  }

  # --- Clean types and state factor ---
  combined[, time := as.numeric(time)]
  combined[, value := as.numeric(value)]

  all_states <- names(dust_state)
  combined[, state := factor(state, levels = all_states)]

  return(combined[])
}
