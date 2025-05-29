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
