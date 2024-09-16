source(file = "./UsefulFunctions_Homogeneous.R")

bfun_time <- function(Mx, times, b_spline_basis,
                      Inf_bound, Sup_bound, K, M) {
  ################################################################################
  ################################################################################
  #### Creates the B-Spline matrix B = (B^j(X_i))_{i,j} in R^{N*n, K+M}
  ####
  #### Mx : Matrix of observed trajectories
  #### Inf_bound : Lower bound
  #### Sup_bound : Upper bound
  #### K : Number of non zero knots
  #### M : Such that K+M = desired number of B-Spline functions
  ##############################################################################
  ##############################################################################

  N <- ncol(Mx)
  n <- nrow(Mx)

  # Create B-spline basis matrix for the time grid
  matrix_B_time <- bspline(times,b_spline_basis,  0, 1)

  # Initialize the B-spline matrix
  matrix_B <- matrix(0, nrow = N * n, ncol = K + M)

  # Loop over each column (observation) to fill the matrix
  for (i in 1:N) {
    bspline_matrix <- bspline(Mx[, i], b_spline_basis, Inf_bound, Sup_bound)
    matrix_B[((i - 1) * n + 1):(i * n), ] <- bspline_matrix * matrix_B_time
  }

  return(matrix_B)
}

################################################################################
################################################################################

selectdimdrift_time <- function(X, times, iZ, c, SetKspline, Inf_bound,
                                Sup_bound, delta, M, Lconst) {
  ################################################################################
  ################################################################################
  #### Returns the value of K that minimizes the penalized criterion
  ####
  #### X : Matrix of observed trajectories
  #### iZ : Matrix of successive differences
  #### c : Coefficient of the penalization
  #### SetKspline : Possible values for dimension K
  #### M : Such that K+M = desired number of B-Spline functions
  #### Inf_bound : Lower bound of the interval considered
  #### Sup_bound : Upper bound of the interval considered
  #### delta : TimeStep of discretization
  #### Lconst : Constant defining the Upper Bound of Approximation
  ################################################################################
  ################################################################################

  N <- ncol(X)
  n <- nrow(X)
  Z <- iZ

  # List to store B-spline basis for each K
  B_ls <- list()
  a_ls <- list()
  gpen_vec <- numeric(length(SetKspline))

  # Loop through each value of K in SetKspline
  for (i in seq_along(SetKspline)) {
    K <- SetKspline[i]
    cat("K value for selection:", K, "\n")
    # The number of basis functions
    nbasis <- K + M

    # Call the Python function to create the B-spline basis
    bspline_basis <- create_bspline_basis(rangeval = c(Inf_bound, Sup_bound), nbasis = nbasis, degree = M)

    B <- bfun_time(X, times, bspline_basis, Inf_bound, Sup_bound, K, M)
    B_ls[[i]] <- B

    # Optimize coefficients for each B-spline basis
    a_ls[[i]] <- optimfun(B, Z, K, M, Lconst)

    # Compute the penalized criterion for each set of coefficients
    gpen_vec[i] <- gamma_pen(a_ls[[i]], c, K, M, Z, B, n, N)
  }

  # Find the index of the minimum penalized criterion
  i_min <- which.min(gpen_vec)

  # Select the corresponding K
  K_ch <- SetKspline[i_min]

  return(K_ch)
}

driftspline_time <- function(x, times, a_hat, b_spline_basis,
                             Inf_bound, Sup_bound, Lconst) {
  ################################################################################
  ################################################################################
  #### Computes the Drift estimator
  ####
  #### x : Observed trajectory
  #### times : Time grid for the observed trajectory
  #### a_hat : Previously computed coefficients of the B-Spline Basis
  #### Inf_bound : Lower bound of the interval considered
  #### Sup_bound : Upper bound of the interval considered
  #### K : Number of non-zero knots
  #### M : Such that K+M = desired number of B-Spline functions
  #### Lconst : Constant defining the Upper Bound of Approximation
  ################################################################################
  ################################################################################

  # Apply labfun to check if values are within bounds
  lab <- labfun(x, Inf_bound, Sup_bound)
  idmat <- diag(as.numeric(lab))

  # Create B-spline basis for the space and time components
  B_space <- bspline(x, b_spline_basis, Inf_bound, Sup_bound)
  B_time <- bspline(times, b_spline_basis, 0, 1)

  # Element-wise multiplication of the space and time B-spline matrices
  B <- B_space * B_time

  # Calculate the drift estimator
  if (length(idmat) != 0) {
    Bs <- idmat %*% B
    b_hat <- Bs %*% a_hat
  } else {
    b_hat <- rep(0, length(a_hat))
  }

  # Estimate the bound for b_hat
  b_hat <- estimbound(b_hat, Lconst)

  return(b_hat)
}

################################################################################
################################################################################

selectdimdiff_time <- function(X, times, U_, c, SetKspline, Inf_bound,
                               Sup_bound, delta, M, Lconst) {
  ################################################################################
  ################################################################################
  #### Returns the value of K that minimizes the penalized criterion
  ####
  #### X : Matrix of observed trajectories
  #### iZ : Matrix of successive differences
  #### c : Coefficient of the penalization
  #### SetKspline : Possible values for dimension K
  #### M : Such that K+M = desired number of B-Spline functions
  #### Inf_bound : Lower bound of the interval considered
  #### Sup_bound : Upper bound of the interval considered
  #### delta : TimeStep of discretization
  #### Lconst : Constant defining the Upper Bound of Approximation
  ################################################################################
  ################################################################################

  N <- ncol(X)
  n <- nrow(X)
  U <- U_

  # List to store B-spline basis for each K
  B_ls <- list()
  a_ls <- list()
  gpen_vec <- numeric(length(SetKspline))

  # Loop through each value of K in SetKspline
  for (i in seq_along(SetKspline)) {
    K <- SetKspline[i]
    cat("K value for selection:", K, "\n")
    # The number of basis functions
    nbasis <- K + M

    # Call the Python function to create the B-spline basis
    bspline_basis <- create_bspline_basis(rangeval = c(Inf_bound, Sup_bound), nbasis = nbasis, degree = M)

    B <- bfun_time(X, times, bspline_basis, Inf_bound, Sup_bound, K, M)
    B_ls[[i]] <- B

    # Optimize coefficients for each B-spline basis
    a_ls[[i]] <- optimfun(B, U, K, M, Lconst)

    # Compute the penalized criterion for each set of coefficients
    gpen_vec[i] <- gamma_pen(a_ls[[i]], c, K, M, U, B, n, N)
  }

  # Find the index of the minimum penalized criterion
  i_min <- which.min(gpen_vec)

  # Select the corresponding K
  K_ch <- SetKspline[i_min]

  return(K_ch)
}

diffspline_time <- function(x, times, a_hat, b_spline_basis,
                            Inf_bound, Sup_bound, Lconst) {
  ################################################################################
  ################################################################################
  #### Computes the Drift estimator
  ####
  #### x : Observed trajectory
  #### times : Time grid for the observed trajectory
  #### a_hat : Previously computed coefficients of the B-Spline Basis
  #### Inf_bound : Lower bound of the interval considered
  #### Sup_bound : Upper bound of the interval considered
  #### K : Number of non-zero knots
  #### M : Such that K+M = desired number of B-Spline functions
  #### Lconst : Constant defining the Upper Bound of Approximation
  ################################################################################
  ################################################################################

  # Apply labfun to check if values are within bounds
  lab <- labfun(x, Inf_bound, Sup_bound)
  idmat <- diag(as.numeric(lab))

  # Create B-spline basis for the space and time components
  B_space <- bspline(x, b_spline_basis, Inf_bound, Sup_bound)
  B_time <- bspline(times, b_spline_basis, 0, 1)

  # Element-wise multiplication of the space and time B-spline matrices
  B <- B_space * B_time

  # Calculate the drift estimator
  if (length(idmat) != 0) {
    Bs <- idmat %*% B
    b_hat <- Bs %*% a_hat
  } else {
    b_hat <- rep(0, length(a_hat))
  }

  # Estimate the bound for b_hat
  b_hat <- estimbound(b_hat, Lconst)

  return(b_hat)
}









