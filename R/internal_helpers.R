#' Reformat Contact Matrix to Match Age Vector
#'
#' This function reformats a raw contact matrix (which uses age groups in 5-year intervals)
#' to match a given age vector, mapping each age group in the contact matrix to the nearest
#' value in the age vector.
#'
#' @param contact_matrix_raw A square matrix of contact rates, with rows and columns corresponding to age groups in 5-year intervals.
#' @param age_vector A vector of age categories that the contact matrix will be reformatted to match.
#'
#' @return A matrix of contact rates with dimensions matching the provided age vector.
#' @keywords internal
reformat_contact_matrix <- function(contact_matrix_raw, age_vector) {

  # Define age groups based on 5-year intervals
  age_group <- seq(0, 80, by = 5)
  age_groups <- age_group[2:length(age_group)]

  # Rename rows and columns of the raw contact matrix to match age groups
  rownames(contact_matrix_raw) <- age_groups
  colnames(contact_matrix_raw) <- age_groups

  # Initialize a new age matrix with the desired age_vector dimensions
  new_age_matrix <- matrix(0,
                           ncol = length(age_vector),
                           nrow = length(age_vector),
                           dimnames = list(age_vector, age_vector))

  # Map contact values from raw matrix to new matrix based on closest matching age group
  for (i in age_vector) {
    for (j in age_vector) {
      # Find the nearest matching age group from the raw matrix
      i_idx <- which.min(abs(age_groups - i))
      j_idx <- which.min(abs(age_groups - j))

      # Assign the contact value to the new matrix
      new_age_matrix[colnames(new_age_matrix) == i, rownames(new_age_matrix) == j] <- contact_matrix_raw[i_idx, j_idx]
    }
  }

  new_age_matrix
}


#' Symmetrize a Contact Matrix Using Population Weights
#'
#' This function symmetrizes a square contact matrix using population sizes,
#' ensuring that total contacts between age groups are reciprocal.
#'
#' The symmetrized matrix is computed using the formula:
#' \deqn{C^\text{sym}_{i,j} = \frac{n_i C_{i,j} + n_j C_{j,i}}{n_i + n_j}}
#' where \eqn{C_{i,j}} is the contact rate from age group \eqn{i} to \eqn{j},
#' and \eqn{n_i} and \eqn{n_j} are the population sizes of groups \eqn{i} and \eqn{j}.
#'
#' @param C A square numeric matrix of contact rates (e.g., POLYMOD contact matrix),
#'   where each element \code{C[i, j]} represents the average number of daily contacts
#'   made by an individual in age group \code{i} with individuals in age group \code{j}.
#' @param pop A numeric vector of population sizes for each age group. Must be the same
#'   length as the number of rows and columns in \code{C}.
#'
#' @return A symmetric matrix of the same dimensions as \code{C}, with contact rates
#'   adjusted to satisfy reciprocal contact symmetry based on population sizes.
#'
#' @examples
#' C <- matrix(c(10, 2, 1,
#'               2, 8, 2,
#'               1, 2, 6), nrow = 3, byrow = TRUE)
#' rownames(C) <- colnames(C) <- c("0-4", "5-17", "18+")
#' pop <- c(3000000, 7000000, 50000000)
#' symmetrize_contact_matrix(C, pop)
#'
#' @keywords internal
symmetrize_contact_matrix <- function(C, pop) {
  if (!is.matrix(C)) stop("C must be a matrix")
  if (length(pop) != nrow(C)) stop("Population vector must match matrix dimensions")
  if (nrow(C) != ncol(C)) stop("Matrix must be square")

  n <- length(pop)
  C_sym <- matrix(0, n, n)

  for (i in 1:n) {
    for (j in 1:n) {
      C_sym[i, j] <- (pop[i] * C[i, j] + pop[j] * C[j, i]) / (pop[i] + pop[j])
    }
  }

  rownames(C_sym) <- rownames(C)
  colnames(C_sym) <- colnames(C)

  return(C_sym)
}


#' Project Matrix to Symmetric Doubly Stochastic Form
#'
#' This function transforms a non-negative matrix into a symmetric, doubly stochastic matrix
#' using iterative projections. It first symmetrizes the matrix, then applies Sinkhorn-Knopp
#' normalization to enforce row and column sums of 1.
#'
#' @param mat A non-negative square matrix.
#' @param max_iter Maximum number of iterations (default = 1000).
#' @param tol Convergence tolerance (default = 1e-8).
#'
#' @return A symmetric, doubly stochastic matrix.
#' @keywords internal
project_to_symmetric_doubly_stochastic <- function(mat, max_iter = 1000, tol = 1e-8) {
  if (any(mat < 0)) stop("Matrix must be non-negative.")

  # Initial guess
  X <- mat / sum(mat)  # Normalize to help convergence

  for (i in 1:max_iter) {
    X_old <- X

    # Step 1: Project to symmetric matrix
    X <- (X + t(X)) / 2

    # Step 2: Project to doubly stochastic matrix using Sinkhorn-Knopp
    for (j in 1:50) {  # Inner loop to normalize
      X <- X / rowSums(X)
      X <- t(t(X) / colSums(X))
    }

    # Check convergence
    if (max(abs(X - X_old)) < tol) {
      return(X)
    }
  }

  warning("Did not converge within max_iter iterations.")
  return(X)
}


#' Generate a Sequence of Years from a Vector and Optional Start/End Bounds
#'
#' This helper function returns a sequence of years based on a full vector of years (`yrs`)
#' and optional start and end bounds. If either `start` or `end` is an empty string,
#' the minimum or maximum value in `yrs` is used, respectively.
#'
#' @param yrs A numeric vector of years from which to derive default start or end values.
#'
#' @param start A numeric or character value specifying the start year. If `""`, the minimum year in `yrs` is used.
#'
#' @param end A numeric or character value specifying the end year. If `""`, the maximum year in `yrs` is used.
#'
#' @return A numeric vector of consecutive years from `start` to `end` (inclusive).
#'
#' @examples
#' get_years(2000:2020, start = 2005, end = 2010)
#' get_years(2000:2020, start = "", end = 2010)
#' get_years(2000:2020, start = 2005, end = "")
#' get_years(2000:2020, start = "", end = "")
#'
#' @keywords internal
get_years <- function(yrs, start, end) {
  start <- if (start == "") min(yrs, na.rm = TRUE) else as.numeric(start)
  end <- if (end == "") max(yrs, na.rm = TRUE) else as.numeric(end)
  start:end
}


#' Collapse Age Bins in a Matrix
#'
#' This function reduces the number of age bins in a matrix by aggregating adjacent columns (age groups)
#' into a specified number of broader bins. Each row is processed independently.
#'
#' @param mat A numeric matrix where rows represent observations (e.g., time points or locations)
#'   and columns represent age groups.
#'
#' @param n_bins An integer specifying the number of broader age bins to collapse into.
#'
#' @return A numeric matrix with the same number of rows as `mat` and `n_bins` columns,
#'   where each column represents the sum of values from grouped original age bins.
#'
#' @examples
#' # Collapse a 5x10 matrix into 5x4 matrix with 4 age bins
#' mat <- matrix(1:50, nrow = 5, ncol = 10)
#' collapse_age_bins(mat, 4)
#'
#' @keywords internal
collapse_age_bins <- function(mat, n_bins) {
  go <- t(apply(mat, 1, function(row) {
    group_size <- ceiling(length(row) / n_bins)
    groups <- (seq_along(row) - 1) %/% group_size + 1
    tapply(row, groups, sum)
  }))
  if(n_bins == 1) t(go) else go
}

#' Split a Vector into Bins and Sum Each Bin
#'
#' This function divides a numeric vector into approximately equal-sized groups (bins),
#' then computes the sum of values within each bin.
#'
#' @param vec A numeric vector to be grouped and summed.
#'
#' @param n_bins An integer specifying the number of bins to split the vector into.
#'
#' @return A numeric vector of length `n_bins` (or fewer if `vec` has fewer elements),
#'   where each element is the sum of values in a bin.
#'
#' @examples
#' # Split and sum a vector of length 10 into 3 bins
#' split_and_sum(1:10, 3)
#'
#' @keywords internal
split_and_sum <- function(vec, n_bins) {
  group_size <- ceiling(length(vec) / n_bins)
  groups <- (seq_along(vec) - 1) %/% group_size + 1
  tapply(vec, groups, sum)
}


#' Filter Routine Vaccination Schedule for a Country and Vaccine
#'
#' Filters a WHO vaccination schedule data frame to extract age-specific routine vaccination
#' rounds for a specific vaccine and country.
#'
#' @param vaccination_schedule A data frame of WHO vaccination schedules.
#' @param vaccination_type A character string for matching VACCINE_DESCRIPTION.
#' @param iso A 3-letter ISO country code.
#'
#' @return A filtered and arranged data frame with an added `age_years` numeric column.
#' @keywords internal
filter_vaccine_schedule <- function(vaccination_schedule, vaccination_type, iso) {
  vaccination_schedule %>%
    dplyr::filter(
      !is.na(age_administered),
      grepl(vaccination_type, vaccine_description, ignore.case = TRUE),
      iso3 == iso,
      target_pop_description != "Travellers",
      !grepl("ADULTS|CATCHUP_C", target_pop, ignore.case = TRUE),
      !grepl("contact", age_administered),
      vaccine_code != "TD_S"
    ) %>%
    dplyr::mutate(age_years = dplyr::case_when(
      grepl("Y", age_administered) ~ as.numeric(gsub("[^0-9.-]", "", age_administered)),
      grepl("M", age_administered) ~ as.numeric(gsub("[^0-9.-]", "", age_administered)) / 12,
      grepl("W", age_administered) ~ as.numeric(gsub("[^0-9.-]", "", age_administered)) / 52,
      grepl("D", age_administered) ~ as.numeric(gsub("[^0-9.-]", "", age_administered)) / 365,
      TRUE ~ NA_real_
    )) %>%
    dplyr::arrange(age_years)
}

#' Build Routine Vaccination Parameters from Schedule and Coverage Data
#'
#' Uses filtered coverage data and country-specific schedule to map routine vaccination
#' parameters to model input dimensions.
#'
#' @param vaccination_data A data frame with annual coverage estimates.
#' @param schedule A filtered and age-parsed WHO schedule object from `filter_vaccine_schedule`.
#' @param ages A numeric vector of model age groups.
#' @param years A numeric vector of model years.
#'
#' @return A data frame with dimensions dim1 (age), dim2 (dose group), dim3, dim4 (year), and value.
#' @keywords internal
build_routine_vaccination_param <- function(vaccination_data, schedule, ages, years) {
  df <- fill_missing_years_general(vaccination_data, "year", "coverage") %>%
    dplyr::group_by(vaccine, vaccine_description, dose_order, year) %>%
    dplyr::summarise(coverage = max(coverage), .groups = "drop") %>%
    dplyr::filter(coverage > 0)

  purrr::map_dfr(seq_len(nrow(df)), function(i) {
    row <- df[i, ]
    dose <- min(row$dose_order, max(schedule$schedulerounds, na.rm = TRUE))
    timing <- schedule %>%
      dplyr::filter(schedulerounds == dose) %>%
      dplyr::slice(1) %>%
      dplyr::pull(age_years)
    targets <- if (row$dose_order == 1) 1 else (2 * row$dose_order - 2):(2 * row$dose_order - 1)
    expand.grid(
      dim1 = match(floor(timing), ages),
      dim2 = targets,
      dim3 = 1,
      dim4 = match(row$year, years),
      value = row$coverage / 100
    )
  })
}

#' Build Seeded Case Parameter Array from Case and Demographic Data
#'
#' Constructs a time- and age-structured seeding input array for model initialization using
#' observed or assumed case distributions.
#'
#' @param processed_case A data frame with year and cases.
#' @param demog_data A named list including `population_data`, a matrix by year and age.
#' @param years Vector of model years.
#' @param ages Vector of model age groups.
#'
#' @return A data frame with dimensions dim1 (age), dim2, dim3, dim4 (year offset), and value.
#' @keywords internal
build_seeded_case_param <- function(processed_case, demog_data, years, ages) {
  df <- fill_missing_years_general(processed_case, "year", "cases") %>%
    dplyr::filter(year <= max(years))
  purrr::map_dfr(seq_len(nrow(df)), function(i) {
    row <- df[i, ]
    popdist <- demog_data$population_data[which(row$year == years), ]
    popdist <- popdist / sum(popdist)
    names(popdist) <- ages + 1
    if (row$cases == 0) {
      data.frame(
        dim1 = ages + 1,
        dim2 = 1,
        dim3 = 1,
        dim4 = i + 1,
        value = 0
      )
    } else {
      sampled <- table(sample(names(popdist), row$cases, prob = popdist, replace = TRUE))
      output <- as.data.frame(sampled)
      names(output) <- c("dim1", "value")
      output$dim1 <- as.integer(as.character(output$dim1))
      output$dim2 <- 1
      output$dim3 <- 1
      output$dim4 <- i + 1
      output
    }
  }) %>%
    dplyr::bind_rows(data.frame(dim1 = 1, dim2 = 1, dim3 = 1, dim4 = 1, value = 0))
}

#' Build SIA Vaccination Parameter Array from Age- and Year-Structured Data
#'
#' Maps SIA (supplementary immunization activity) coverage data to model parameter dimensions.
#'
#' @param sia_data A data frame with columns: age, year, coverage.
#' @param ages A numeric vector of model age groups.
#' @param years A numeric vector of model years.
#'
#' @return A data frame with dimensions dim1 (age), dim2, dim3, dim4 (year), and value (coverage).
#' @keywords internal
build_sia_vaccination_param <- function(sia_data, ages, years) {
  df <- fill_missing_years_general(sia_data, "year", "coverage") %>%
    dplyr::group_by(disease, vaccination_name, age, year) %>%
    dplyr::summarise(coverage = max(coverage), .groups = "drop")
  data.frame(
    dim1 = match(df$age, ages),
    dim2 = 1,
    dim3 = 1,
    dim4 = match(df$year, years),
    value = df$coverage
  )
}

#' Fill Missing Years by Group (Tidyverse)
#'
#' For each combination of grouping columns (all columns except year_col and value_col),
#' this function fills in missing consecutive years, setting the value_col to 0,
#' and replicates the grouping information. One additional year is added after the last.
#'
#' @param df A data frame with at least a year and value column.
#' @param year_col The name of the year column (as a string).
#' @param value_col The name of the value column (as a string).
#'
#' @return A tidy data frame with missing years filled for each group.
#'
#' @examples
#' df <- tibble::tibble(
#'   region = c("A", "A", "A", "B", "B"),
#'   source = c("X", "X", "X", "Y", "Y"),
#'   year = c(2000, 2005, 2010, 2000, 2010),
#'   value = c(1, 25, 29, 3, 12)
#' )
#' fill_missing_years_by_group(df, year_col = "year", value_col = "value")
#'
#' @export
fill_missing_years_general <- function(df, year_col, value_col) {

  year_sym <- sym(year_col)
  value_sym <- sym(value_col)

  # Grouping columns are all others
  group_cols <- setdiff(names(df), c(year_col, value_col))

  df_filled <- df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      min_year = min(.data[[year_col]]),
      max_year = max(.data[[year_col]]),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(data = list(tibble(!!year_sym := seq(min_year, max_year + 1)))) %>%
    select(-min_year, -max_year) %>%
    unnest(data) %>%
    left_join(df, by = c(setNames(group_cols, group_cols), year_col)) %>%
    mutate(
      !!value_sym := if_else(is.na(.data[[value_col]]), 0, .data[[value_col]])
    ) %>%
    arrange(across(all_of(c(group_cols, year_col))))

  return(df_filled)
}

#' Combine Routine and SIA Vaccination Coverage Parameters
#'
#' Merges routine and SIA (supplementary immunization activity) vaccination parameter data frames
#' and ensures values are capped at 1 (100% coverage) for each unique parameter combination.
#'
#' @param routine_df A data frame of routine vaccination parameters. Must contain columns:
#'   `dim1`, `dim2`, `dim3`, `dim4`, and `value`.
#' @param sia_df A data frame of SIA vaccination parameters (optional). If `NULL`, only routine data is used.
#'
#' @return A data frame containing combined vaccination parameters with coverage values capped at 1.
#'   The returned data frame contains columns: `dim1`, `dim2`, `dim3`, `dim4`, and `value`.
#'
#' @details
#' This function is typically used in the pipeline preparing model inputs to ensure
#' vaccination coverage parameters are bounded and consolidated before being passed
#' into simulation.
#'
#' @examples
#' # combine_vaccination_params(routine_df, sia_df)
#'
#' @keywords internal
combine_vaccination_params <- function(routine_df, sia_df = NULL) {
  combined <- dplyr::bind_rows(routine_df, sia_df)
  combined %>%
    dplyr::group_by(dim1, dim2, dim3, dim4) %>%
    dplyr::summarise(value = pmin(sum(value), 1), .groups = "drop")
}

#' Convert Long-Format Data Frame to Multi-Dimensional Array
#'
#' This internal function converts a tidy data frame containing dimension columns and a `value` column
#' into a numeric array. Each row specifies the value for a particular position in the output array.
#'
#' @param df A data frame containing dimension columns named `dim1`, `dim2`, ..., and a `value` column.
#'
#' @return A numeric array with shape inferred from the maximum index in each dimension column.
#'
#' @keywords internal
df_to_array <- function(df) {
  dim_cols <- grep("^dim\\d+$", names(df), value = TRUE)
  dims <- vapply(dim_cols, function(col) as.integer(max(df[[col]])), integer(1))
  strides <- c(1L, cumprod(head(dims, -1L)))

  idx_matrix <- as.matrix(df %>% select(all_of(dim_cols)))
  lin_idx <- rowSums(sweep(idx_matrix - 1L, 2, strides, `*`)) + 1L

  arr <- array(0, dim = dims)
  arr[lin_idx] <- df$value
  arr
}



#' Generate a Long-Format Data Frame Representing a Multi-Dimensional Array
#'
#' This internal helper function creates a tidy data frame representing a multi-dimensional array
#' up to six dimensions. Each row corresponds to a unique combination of dimension indices.
#' Optionally, specific values can be overridden using the `updates` argument.
#'
#' @param dim1,dim2,dim3,dim4,dim5,dim6 Integers specifying the size of each dimension.
#'   Only `dim1` is required. Others may be `NULL`.
#' @param default_value Numeric. Default value to assign to each array element (default is 0).
#' @param updates Optional `data.frame` with columns `dim1`, `dim2`, ..., and `value` to override specific positions.
#'
#' @return A `data.frame` with columns `ID`, `dim1`, ..., and `value`, representing all array coordinates.
#'
#' @keywords internal
generate_array_df <- function(dim1, dim2 = NULL, dim3 = NULL, dim4 = NULL, dim5 = NULL, dim6 = NULL,
                              default_value = 0, updates = NULL) {
  dims <- list(dim1, dim2, dim3, dim4, dim5, dim6)
  dims <- dims[!vapply(dims, is.null, logical(1))]
  n_dims <- length(dims)

  index_grid <- do.call(expand.grid, lapply(dims, seq_len))
  names(index_grid) <- paste0("dim", seq_len(n_dims))
  index_grid$value <- default_value

  if (!is.null(updates) && nrow(updates) > 0) {
    key_cols <- paste0("dim", seq_len(n_dims))

    # Coerce both updates and grid to integer
    updates[, (key_cols) := lapply(.SD, as.integer), .SDcols = key_cols]
    index_grid[key_cols] <- lapply(index_grid[key_cols], as.integer)

    # Ensure all combinations are filled with default unless updated
    data.table::setDT(index_grid)
    data.table::setDT(updates)

    full_grid <- merge(index_grid, updates, by = key_cols, all.x = TRUE, suffixes = c("", ".update"))
    full_grid[, value := fifelse(!is.na(value.update), value.update, default_value)]
    full_grid[, value.update := NULL]

    index_grid <- full_grid
  }

cbind(ID = seq_len(nrow(index_grid)), index_grid)
}


#' Convert Long-Format Updates into an Array
#'
#' This internal wrapper function generates a tidy grid using `generate_array_df()` and reshapes it into an array
#' using `df_to_array()`. It simplifies the pipeline for transforming partial long-format input into a structured array.
#'
#' @param ... Arguments passed to `generate_array_df()` (including `updates`, `dim1`, `dim2`, etc.).
#'
#' @return A multi-dimensional array constructed from the provided dimensions and update values.
#'
#' @keywords internal
array_from_df <- function(...) df_to_array(generate_array_df(...))
