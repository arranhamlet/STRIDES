#' Convert Age Break String to Age Group Labels
#'
#' Converts a semicolon-separated string of numeric age breakpoints into a character
#' vector of age group labels. For example, the input `"0;5;10;...;80;Inf"` will produce
#' the labels `"0–4"`, `"5–9"`, ..., `"80+"`.
#'
#' @param age_string A single character string of age breakpoints, separated by semicolons.
#'   Must include at least two breakpoints. The last breakpoint can be `Inf` to indicate
#'   an open-ended upper age group.
#'
#' @return A character vector of formatted age group labels.
#'
#' @examples
#' age_breaks_to_labels("0;5;10;15;20;25;30;35;40;45;50;55;60;65;70;75;80;Inf")
#' # Returns:
#' # [1] "0–4"   "5–9"   "10–14" "15–19" "20–24" ... "75–79" "80+"
#'
#' @export
age_breaks_to_labels <- function(age_string) {
  # Split string on semicolons and convert to numeric, preserving Inf
  parts <- strsplit(age_string, ";", fixed = TRUE)[[1]]
  breaks <- suppressWarnings(as.numeric(parts))

  if (length(breaks) < 2 || anyNA(breaks)) {
    stop("Input must be a semicolon-separated string of numeric values, e.g., '0;5;10;...;Inf'")
  }

  # Construct labels
  labels <- character(length(breaks) - 1)
  for (i in seq_along(labels)) {
    from <- breaks[i]
    to <- breaks[i + 1]
    labels[i] <- if (is.infinite(to)) {
      paste0(from, "+")
    } else {
      paste0(from, "\u2013", to - 1)  # en dash
    }
  }

  return(labels)
}
