#' Unpack Dust Model State Data
#'
#' This function extracts and unpacks data from a dust model state object, separating compartments and non-compartments.
#'
#' @param model_system The model system object containing packing information and particle details.
#' @param model_object The dust model object containing state data.
#' @param dimension_names List of dimension names for states, time, and other variables.
#' @param which_state_dimensions List mapping state names to their associated dimensions.
#'
#' @return A data frame containing the unpacked dust model state data.
#' @import data.table
unpack_dust2 <- function(model_system, model_object, dimension_names, which_state_dimensions) {

  dust_state <- dust_unpack_state(model_system, model_object)
  these_are_compartments <- lengths(model_system$packing_state) != 0
  time_length <- dimension_names$time[[1]]
  n_particles <- model_system$n_particles

  # --- Compartmental states ---
  unpacked_compartments <- rbindlist(
    lapply(which(these_are_compartments), function(x) {
      this_obj <- dust_state[[x]]
      state_name <- names(dust_state)[x]
      present_dimensions <- which_state_dimensions[[state_name]]
      colnames <- present_dimensions

      processed <- process_obj(this_obj, present_dimensions, colnames, time_length, x, model_system, dimension_names, dust_state)
      rbindlist(processed, fill = TRUE)
    }), fill = TRUE
  )

  # --- Non-compartmental states ---
  unpacked_noncompartments <- rbindlist(
    lapply(which(!these_are_compartments), function(x) {
      this_obj <- dust_state[[x]]
      state_name <- names(dust_state)[x]

      if (n_particles > 1) {
        # vector of length n_particles
        dt <- data.table(run = paste0("run_", seq_len(n_particles)),
                         value = as.numeric(this_obj),
                         time = time_length,
                         state = state_name)
      } else {
        dt <- data.table(value = as.numeric(this_obj),
                         time = time_length,
                         state = state_name)
      }
      return(dt)
    }), fill = TRUE
  )

  # --- Combine and format ---
  combined <- rbindlist(list(unpacked_compartments, unpacked_noncompartments), fill = TRUE)

  colnames_use <- setdiff(unique(unlist(which_state_dimensions)), "time")

  # Convert to character and fill NAs for stratifiers
  for (col in colnames_use) {
    if (!col %in% names(combined)) {
      combined[[col]] <- "All"
    } else {
      combined[[col]] <- as.character(combined[[col]])
      combined[[col]] <- fifelse(is.na(combined[[col]]), "All", combined[[col]])
    }
  }

  combined[, time := as.numeric(time)]
  combined[, value := as.numeric(value)]
  combined[, state := factor(state, levels = names(these_are_compartments))]

  return(combined[])
}
