#' Calculate Final Susceptibility and Recovered States
#'
#' Extracts the final timepoint values for each compartment and converts them
#' into arrays suitable for model initialisation. This function produces the initial
#' susceptible and recovered populations from the output of a model run.
#'
#' @param model_run A long-format `data.frame` or `data.table` of model output
#'   including columns `time`, `state`, `age`, `vaccination`, `risk`, and `value`.
#'
#' @return A named list with:
#' \describe{
#'   \item{S0}{An array of susceptible individuals by age, vaccination, and risk group.}
#'   \item{Rpop0}{An array of recovered individuals (all non-susceptible states) by age, vaccination, and risk group.}
#' }
#'
#' @importFrom dplyr rename select mutate across
#' @importFrom collapse fgroup_by fsummarise fsum
#' @export
calculate_current_susceptibility <- function(model_run) {

  # -------------------------
  # Final timepoint values
  # -------------------------
  final_values <- model_run %>%
    subset(time == max(time)) %>%
    subset(state %in% c("S", "E", "I", "R", "Is", "Rc") & age != "All") %>%
    dplyr::mutate(dplyr::across(.cols = age:value, .fns = as.numeric))

  # -------------------------
  # Susceptible array (S)
  # -------------------------
  susceptible <- final_values %>%
    subset(state == "S") %>%
    dplyr::rename(dim1 = age, dim2 = vaccination, dim3 = risk) %>%
    dplyr::select(dim1, dim2, dim3, value) %>%
    df_to_array()

  # -------------------------
  # Recovered array (all non-S)
  # -------------------------
  recovered <- final_values %>%
    subset(state != "S") %>%
    collapse::fgroup_by(age, vaccination, risk) %>%
    collapse::fsummarise(value = collapse::fsum(value)) %>%
    dplyr::rename(dim1 = age, dim2 = vaccination, dim3 = risk) %>%
    dplyr::select(dim1, dim2, dim3, value) %>%
    df_to_array()

  # -------------------------
  # Return initial conditions
  # -------------------------
  list(
    S0 = susceptible,
    Rpop0 = recovered
  )
}
