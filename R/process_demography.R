#' Prepare Demographic Inputs for a Transmission Model (with Custom Age Groups)
#'
#' This function processes raw demographic data into a structured format using
#' custom age groups, suitable for transmission modeling.
#'
#' @inheritParams process_demography
#'
#' @return A named list of formatted demographic and structural inputs.
#' @import dplyr
#' @import purrr
#' @import tidyr
#' @importFrom reshape2 melt
process_demography <- function(
    migration,
    fertility,
    mortality,
    population_all,
    population_female,
    contact_matricies,
    iso,
    year_start = "",
    year_end = "",
    number_of_vaccines = 0,
    n_risk = 1
) {
  age_breaks <- c(0, 1, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, Inf)
  age_labels <- paste0(head(age_breaks, -1), "-", head(age_breaks[-1] - 1, -1))
  age_labels[length(age_labels)] <- "80+"
  n_age <- length(age_labels)
  n_vacc <- if (number_of_vaccines == 0) 1 else number_of_vaccines * 2 + 1

  years_all <- get_years(migration$year, start = year_start, end = year_end)
  time_run_for <- length(years_all)
  time_all <- 0:(time_run_for - 1)

  migration <- migration %>% filter(iso3 == iso, year %in% years_all)
  fertility <- fertility %>% filter(iso3 == iso, year %in% years_all)
  mortality <- mortality %>% filter(iso3 == iso, year %in% years_all)
  population_all <- population_all %>% filter(iso3 == iso, year %in% years_all)
  population_female <- population_female %>% filter(iso3 == iso, year %in% years_all)

  all_ages <- 0:100
  female_ages <- 15:49
  pop_cols <- paste0("x", all_ages)
  fem_cols <- paste0("x", female_ages)

  collapse_age_custom <- function(mat, age_breaks) {
    all_ages <- 0:(ncol(mat) - 1)
    age_groups <- cut(all_ages, breaks = age_breaks, right = FALSE, labels = FALSE)
    collapsed <- t(apply(mat, 1, function(row) tapply(row, age_groups, sum, na.rm = TRUE)))
    if (nrow(collapsed) == 1) matrix(collapsed, nrow = 1) else collapsed
  }

  pop_all_raw <- population_all %>% select(all_of(pop_cols)) %>% as.matrix()
  pop_all <- collapse_age_custom(pop_all_raw, age_breaks)

  # Contact matrix handling
  contact_group_bounds <- seq(5, 80, by = 5)  # 0–4, 5–9, ..., 80–84
  contact_group_midpoints <- (head(contact_group_bounds, -1) + tail(contact_group_bounds, -1)) / 2
  contact_to_model_index <- cut(
    contact_group_midpoints,
    breaks = age_breaks,
    right = FALSE,
    labels = FALSE
  )

  country_contact <- if (!iso %in% names(contact_matricies)) {
    Reduce(`+`, contact_matricies) / length(contact_matricies)
  } else {
    contact_matricies[[iso]]
  }

  # Define midpoints of custom age groups
  custom_midpoints <- (head(age_breaks, -1) + tail(age_breaks, -1)) / 2
  custom_midpoints[is.infinite(custom_midpoints)] <- max(age_breaks[is.finite(age_breaks)])

  # Define midpoints of 5-year bins in contact matrix: 0–4, 5–9, ..., 80–84
  contact_group_bounds <- seq(0, 85, by = 5)
  contact_midpoints <- (head(contact_group_bounds, -1) + tail(contact_group_bounds, -1)) / 2

  # Map each custom midpoint to closest original contact matrix row/col
  nearest_idx <- function(x, from) which.min(abs(from - x))
  row_map <- vapply(custom_midpoints, nearest_idx, from = contact_midpoints, FUN.VALUE = integer(1))
  col_map <- vapply(custom_midpoints, nearest_idx, from = contact_midpoints, FUN.VALUE = integer(1))

  # Build reduced contact matrix by nearest value
  reduced_contact_matrix <- matrix(0, nrow = n_age, ncol = n_age)
  for (i in 1:n_age) {
    for (j in 1:n_age) {
      reduced_contact_matrix[i, j] <- country_contact[row_map[i], col_map[j]]
    }
  }

  reformatted_contact_matrix <- symmetrize_contact_matrix(reduced_contact_matrix, pop = pop_all[nrow(pop_all), ])
  reformatted_contact_matrix <- project_to_symmetric_doubly_stochastic(reformatted_contact_matrix)

  mort_mat <- mortality %>% select(all_of(pop_cols)) %>% as.matrix()
  mort_mat <- collapse_age_custom(mort_mat, age_breaks)
  mortality_rate <- pmin(mort_mat / pop_all, 1)
  mortality_rate[!is.finite(mortality_rate)] <- 1

  mortality_df <- reshape2::melt(t(mortality_rate)) %>%
    setNames(c("dim1", "dim3", "value")) %>%
    mutate(dim2 = 1)

  fert_mat <- fertility %>% select(all_of(fem_cols)) %>% as.matrix()
  pop_fem <- population_female %>% select(all_of(fem_cols)) %>% as.matrix()
  denom <- rowSums(pop_fem)
  denom[denom == 0] <- NA

  fertility_by_year <- tibble(
    dim1 = n_risk,
    dim2 = time_all + 1,
    value = pmin(rowSums((fert_mat / 1000) * pop_fem) / denom, 1)
  )

  mig_rates <- migration[["migration_rate_1000"]]

  migration_in_number <- map_dfr(seq_len(nrow(pop_all)), function(i) {
    mig_vals <- round(pop_all[i, ] * mig_rates[i])
    tibble(
      dim1 = seq_len(n_age),
      dim2 = 1,
      dim3 = 1,
      dim4 = i,
      value = mig_vals
    )
  })

  init_vals <- round(pop_all[1, ] * 1000)
  N0_df <- tibble(
    dim1 = seq_len(n_age),
    dim2 = 1,
    dim3 = 1,
    value = init_vals
  )

  total_population_df <- map_dfr(seq_len(nrow(pop_all)), function(i) {
    chunk <- round(pop_all[i, ] * 1000)
    tibble(
      dim1 = seq_len(n_age),
      dim2 = 1,
      dim3 = 1,
      dim4 = i,
      value = chunk
    )
  })

  migration_distribution_values <- expand_grid(
    dim1 = 1,
    dim2 = seq_along(time_all)
  ) %>% mutate(value = 1)

  list(
    N0 = N0_df,
    crude_birth = fertility_by_year,
    crude_death = mortality_df,
    tt_migration = time_all,
    migration_in_number = migration_in_number,
    migration_distribution_values = migration_distribution_values,
    population_data = pop_all,
    contact_matrix = reformatted_contact_matrix,
    aging = 1/diff(age_breaks),
    input_data = tibble(
      iso = iso,
      year_start = min(years_all),
      year_end = max(years_all),
      n_age = n_age,
      n_vacc = n_vacc,
      n_risk = n_risk
    )
  )
}
