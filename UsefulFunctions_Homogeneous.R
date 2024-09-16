library(splines)
library(pbapply)
library(reticulate)
library(MASS)
library(rootSolve)
library(fda)
library(data.table)
py_config()

zfun <- function(Mx, delta) {
  ################################################################################
  ################################################################################
  #### Computes the successive differences of order 2
  ####
  #### Mx : Data matrix, containing sampled trajectories along rows.
  #### delta : TimeStep used to normalize the differences
  ################################################################################
  ################################################################################

  # Compute the squared successive differences of order 1 along the rows
  diffs <- diff(Mx, differences = 1)

  # Flatten the result in column-major order (similar to 'F' order in Python)
  result <- as.vector(t(diffs)) / delta

  return(result)
}

zfun_authors <-function(Mx,delta){
  if(ncol(Mx)==1) result <- diff(Mx[,1])[-1]/delta
  if(ncol(Mx)>1){
    result <- diff(Mx[,1])[-1]/delta
    for (i in 2:dim(Mx)[2]){
      result <- c(result, (diff(Mx[,i])[-1])/delta)
    }
  }
  return(result)
}

ufun <- function(Mx, delta) {
  ################################################################################
  ################################################################################
  #### Computes the successive differences of order 2
  ####
  #### Mx : Data matrix, containing sampled trajectories along rows.
  #### delta : TimeStep used to normalize the differences
  ################################################################################
  ################################################################################

  # Compute the squared successive differences of order 1 along the rows
  diffs <- diff(Mx, differences = 1)^2

  # Flatten the result in column-major order (similar to 'F' order in Python)
  result <- as.vector(t(diffs)) / delta

  return(result)
}

ufun_authors <-function(Mx,delta){
  if(ncol(Mx)==1) result <- diff(Mx[,1])[-1]/delta
  if(ncol(Mx)>1){
    result <- diff(Mx[,1])[-1]/delta
    for (i in 2:dim(Mx)[2]){
      result <- c(result, (diff(Mx[,i])[-1])/delta)
    }
  }

  return(delta*(result**2))
}

Gaussian_kernel <- function(u,h,power){
  dnorm(u^power, mean = 0, sd = h)
}

zfun_bis <- function(Mx, delta) {
  ################################################################################
  ################################################################################
  #### Function used to compute the successive differences of order 1
  ####
  #### Mx : Data matrix, containing sampled trajectories along rows.
  #### delta : TimeStep used to normalize the differences
  ################################################################################
  ################################################################################
  # Compute the successive differences of order 1
  N <- ncol(Mx)
  n <- nrow(Mx)
  #M <- as.integer(log(N))
  #h <- 1/log(N)
  power <- 1
  Z <- diff(Mx, differences = 1) / delta
  result <- matrix(0, nrow = n-1, ncol = N)
  for (k in pblapply(1:(n-1), identity)) {
    # Compute pairwise differences for the k-th row of X
    # diffs <- outer(Mx[k, 1:M], Mx[k, ], "-") # (N, N) matrix of pairwise differences
    diffs <- outer(Mx[k, ], Mx[k, ], "-") # (N, N) matrix of pairwise differences

    X_bar_k <- mean(Mx[k, ])
    sigma_hat <- sqrt(sum((Mx[k, ] - X_bar_k)^2) / (N - 1))
    sigma_hat <- ifelse(sigma_hat > 1e-30, sigma_hat, 1e-30)
    h_k <-1.06 * sigma_hat * N^(-1/5)
    # Apply the kernel function to each element in the difference matrix
    K_values <- Gaussian_kernel(diffs, h_k, power) # Apply kernel element-wise

    # Compute the numerator and denominator
    numer <- K_values %*% Z[k, ] # (N, N) %*% (N, 1) -> (N, 1)
    denom <- rowSums(K_values)   # Sum over rows to get (N, 1) vector

    # Element-wise division to get the desired (N, 1) vector
    result[k,] <- numer / denom
  }


  return(as.vector(result))
}


################################################################################
################################################################################

labfun <- function(x, Inf_bound, Sup_bound) {
  ##############################################################################
  ##############################################################################
  #### Checks if every coordinate of a trajectory is within
  #### a specified interval [Inf_bound, Sup_bound]
  ####
  #### x : Observed trajectory
  #### Inf_bound : Lower bound
  #### Sup_bound : Higher bound
  ##############################################################################
  ##############################################################################

  # Check if each element of x is within the specified interval
  result <- (Inf_bound <= x) & (x <= Sup_bound)

  # Convert the logical vector to numeric (1 for TRUE, 0 for FALSE)
  return(as.numeric(result))
}

################################################################################
################################################################################

boundfun <- function(x, Inf_bound, Sup_bound) {
  ################################################################################
  ################################################################################
  #### Replaces every coordinate of x that is outside of [Inf_bound, Sup_bound]
  #### by Inf_bound
  ####
  #### x : Observed trajectory
  #### Inf_bound : Lower bound
  #### Sup_bound : Higher bound
  ################################################################################
  ################################################################################

  # Replace values outside the specified bounds with Inf_bound
  result <- (x - Inf_bound) * ((Inf_bound <= x) & (x <= Sup_bound)) + Inf_bound

  return(result)
}

py_run_string("
import numpy as np

from scipy.interpolate import BSpline

def create_bspline_basis(rangeval, nbasis, degree):
    num_points = int(nbasis - degree + 1)  # Ensure integer for linspace
    knots = np.linspace(rangeval[0], rangeval[1], num_points)
    knots = np.concatenate(([rangeval[0]] * int(degree), knots, [rangeval[1]] * int(degree)))
    return BSpline(knots, np.eye(int(nbasis)), int(degree))
")

# Access the Python function from the R environment
create_bspline_basis <- py$create_bspline_basis

# # Function to apply B-spline and boundary checks
bspline <- function(x, bs_basis, InfInt, SupInt) {
  y <- boundfun(x, InfInt, SupInt)
  id_matrix <- diag(as.numeric(labfun(x, InfInt, SupInt)))

  # Convert y to a Python object using reticulate
  y_py <- r_to_py(y)

  # Evaluate the B-spline basis for y using the Python call syntax
  v_basis <- bs_basis(y_py)

  # Convert the result back to an R matrix
  v_basis <- py_to_r(v_basis)

  # Compute the result using matrix multiplication
  if (length(id_matrix) != 0) {
    result <- id_matrix %*% v_basis
  } else {
    result <- matrix(0, nrow = nrow(v_basis), ncol = ncol(v_basis))
  }

  return(result)
}

bspline_authors <- function(x, InfInt, SupInt, K, M){
  y <- boundfun(x, InfInt, SupInt)
  id <- diag(labfun(x, InfInt, SupInt))
  dm <- K + M
  bs_basis <- create.bspline.basis(rangeval=c(InfInt, SupInt),nbasis = dm)
  v.basis <- eval.basis(y, basisobj = bs_basis)
  if(length(id)!=0) result <- id%*%v.basis
  if(length(id)==0) result <- 0*v.basis

  return(result)
}

bfun <- function(Mx, bs_basis, InfInt, SupInt, K, M){
  N = ncol(Mx);  n = nrow(Mx)
  df.X <- as.data.frame(Mx)
  l.X <- as.list(df.X)
  B. <- lapply(l.X,function(x) bspline(x, bs_basis, InfInt, SupInt, K, M))
  B.. <- lapply(B.,function(x) as.data.frame(x))
  B <- rbindlist(B..)

  return(as.matrix(B))
}

bfun_authors <- function(Mx, InfInt, SupInt, K, M){
  N = ncol(Mx);  n = nrow(Mx)
  df.X <- as.data.frame(Mx[2:(n-1),])
  l.X <- as.list(df.X)
  B. <- lapply(l.X,function(x) bspline_authors(x, InfInt, SupInt, K, M))
  B.. <- lapply(B.,function(x) as.data.frame(x))
  B <- rbindlist(B..)

  return(as.matrix(B))
}
# Xtrain <- matrix(0, nrow = 10, ncol = 10)
# n <- nrow(Xtrain)
# Z <- zfun(Xtrain, 1/n)
# B <- bfun(Xtrain[1:(n-1),], -4, 4, K = 2, M = 3)
# a_hat <- optimfun(B, Z, 2, 3, 1)
# driftspline(Xtrain,a_hat,-4,4,2,3,1)

################################################################################
################################################################################

ridge <- function(x, matrix.B, vector.Z, K, M, L){
  upper_bd = (K+M)*L;
  u = ginv(t(matrix.B)%*%matrix.B + x*diag(K+M))%*%t(matrix.B)%*%vector.Z
  result = sum(u^2) - upper_bd

  return(result)
}


optimfun <- function(matrix.B, vector.Z, K, M, L){
  #cat("K value in optimfun:", K, "\n")
  # cat('Dim matrix B:', dim(matrix.B))
  # cat('Dim vector Z:', length(vector.Z))
  upper_bound <- sqrt(L*(K+M))
  estimator <- ginv(t(matrix.B)%*%matrix.B)%*%t(matrix.B)%*%vector.Z;
  norm_estimator <- sqrt(sum(estimator^2))
  Px = t(matrix.B) %*% matrix.B
  if(det(Px) != 0 & norm_estimator <= upper_bound){
    result <- estimator
  }else{
    ridge.bis <- function(x) ridge(x,matrix.B,vector.Z,K,M,L)
    root = multiroot(ridge.bis, c(0.01));
    lambda <- root$root
    result <- ginv(t(matrix.B)%*%matrix.B+lambda*diag(K+M))%*%t(matrix.B)%*%vector.Z;
  }

  return(result)
}

################################################################################
################################################################################

gamma_pen <- function(a, c, K, M, Z, B, n, N) {
  ################################################################################
  ################################################################################
  #### Computes the penalized criterion in order to select a dimension K.
  ####
  #### a : Coefficients of the B-Spline basis
  #### B : B-Spline basis
  #### Z : Matrix of successive differences
  #### K : Number of non-zero knots
  #### M : Such that K+M = desired number of B-Spline functions
  #### c : Coefficient of the penalization
  #### n : Number of discretization
  #### N : Number of observations considered
  ################################################################################
  ################################################################################

  # Compute the penalized criterion
  penalty_term <- c * (log(N) * (K + M) / N)
  residuals <- Z - B %*% a
  criterion <- (1 / (n * N)) * sum(residuals^2) + penalty_term

  return(criterion)
}

################################################################################
################################################################################

selectdimdrift <- function(X, iZ, c, SetKspline, Inf_bound, Sup_bound,
                           M, Lconst) {
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

    B <- bfun(X, Inf_bound, Sup_bound, K, M)
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

selectdimdrift_authors <- function(X, c, SetKspline, M, lower, upper, delta, Lconst){
  N <- ncol(X); n <- nrow(X) - 2; N. <- ncol(D); Z <- zfun_authors(X, delta)
  B.ls<-lapply(SetKspline, function(x) bfun_authors(X, lower, upper, x, M))
  a.ls<-lapply(1:length(SetKspline), function(x) optimfun(B.ls[[x]], Z, SetKspline[x], M, Lconst))
  gpen.vec<-sapply(1:length(SetKspline), function(x) gamma_pen(a.ls[[x]], c, SetKspline[x], M, Z, B.ls[[x]], n, N))
  i.min<-which.min(gpen.vec)
  K.ch<-SetKspline[i.min]

  result=K.ch

  return(result)
}

# Xtrain <- matrix(0, nrow = 10, ncol = 10)
# n <- nrow(Xtrain)
# Z <- zfun(Xtrain, 1/n)
# B <- bfun(Xtrain, -4, 4, K = 2, M = 3)
# a_hat <- optimfun(B, Z, 2, 3, 1)
# driftspline(Xtrain,a_hat,-4,4,2,3,1)

################################################################################
################################################################################

estimbound <- function(x, Lconst) {
  ################################################################################
  ################################################################################
  #### Truncates x in the interval [-sqrt(Lconst), sqrt(Lconst)]
  ####
  #### x : Observed trajectory
  #### Lconst : Upper bound of the interval
  ################################################################################
  ################################################################################

  # Calculate the bound as the square root of Lconst
  bound <- sqrt(Lconst)

  # Apply truncation using vectorized operations
  result <- ifelse(abs(x) <= bound, x, sign(x) * bound)

  return(result)
}

################################################################################
################################################################################

driftspline <- function(x, a_hat, Inf_bound, Sup_bound, K, M, Lconst) {
  ################################################################################
  ################################################################################
  #### Computes the Drift estimator
  ####
  #### x : Observed trajectory
  #### a_hat : Previously computed coefficients of the B-Spline Basis
  #### Inf_bound : Lower bound of the interval considered
  #### Sup_bound : Upper bound of the interval considered
  #### K : Number of non-zero knots
  #### M : Such that K+M = desired number of B-Spline functions
  #### Lconst : Constant defining the Upper Bound of Approximation
  ################################################################################
  ################################################################################
  # Apply labfun to determine which elements are within the bounds
  lab <- labfun(x, Inf_bound, Sup_bound)
  idmat <- diag(lab)

  # Compute the B-spline basis matrix
  B <- bspline(x, Inf_bound, Sup_bound, K, M)

  # Compute the drift estimator
  if (length(idmat) != 0) {
    Bs <- idmat %*% B
    b_hat <- Bs %*% a_hat
  } else {
    b_hat <- rep(0, length(a_hat))
  }

  # Apply the estimation bounds
  b_hat <- estimbound(b_hat, Lconst)

  return(b_hat)
}

driftspline_authors <- function(x, ach, lower, upper, K, M, Lconst){
  lab <- labfun(x, lower, upper)
  idmat <- diag(lab)
  bMat = bspline_authors(x, lower, upper, K, M)
  if(length(idmat)!=0){
    bsMat <- idmat%*%bMat;
    bvalues <- bsMat%*%ach
  }
  if(length(idmat) == 0) bvalues <- 0*ach;
  bvalues<-estimbound(bvalues,Lconst)

  return(as.vector(bvalues))
}

# Xtrain <- matrix(0, nrow = 10, ncol = 10)
# n <- nrow(Xtrain)
# Z <- zfun(Xtrain, 1/n)
# B <- bfun(Xtrain, -4, 4, K = 2, M = 3)
# a_hat <- optimfun(B, Z, 2, 3, 1)
# driftspline(Xtrain,a_hat,-4,4,2,3,1)

################################################################################
################################################################################

gamma_pen2 <- function(a, c, K, M, U, B, n, N) {
  ################################################################################
  ################################################################################
  #### Computes the penalized criterion in order to select a dimension K.
  ####
  #### a : Coefficients of the B-Spline basis
  #### B : B-Spline basis
  #### U : Matrix of successive differences
  #### K : Number of non-zero knots
  #### M : Such that K+M = desired number of B-Spline functions
  #### c : Coefficient of the penalization
  #### n : Number of discretization
  #### N : Number of observations considered
  ################################################################################
  ################################################################################

  # Compute residuals
  residuals <- U - B %*% a

  # Compute the penalized criterion
  penalty_term <- c * log(N) * (K + M) / (N * n)
  criterion <- (1 / (n * N)) * sum(residuals^2) + penalty_term

  return(criterion)
}

################################################################################
################################################################################

selectdimdiff <- function(X, U_, c, SetKspline, Inf_bound, Sup_bound,
                          M, Lconst) {
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

    B <- bfun(X, Inf_bound, Sup_bound, K, M)
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

selectdimdiff_authors <- function(X,c,set.K,M,lower,upper,delta,L.N){
  N<-ncol(X); n <- nrow(X) - 2; N.<-ncol(D); U<-ufun_authors(X,delta)
  B.ls<-lapply(set.K,function(x) bfun_authors(X,lower,upper,x,M))
  L.N<-log(N)
  a.ls<-lapply(1:length(set.K),function(x) optimfun(B.ls[[x]],U,set.K[x],M,L.N))
  gpen.vec<-sapply(1:length(set.K),function(x) gamma_pen2(a.ls[[x]],c,set.K[x],M,U,B.ls[[x]],n,N))
  i.min<-which.min(gpen.vec)
  K.ch<-set.K[i.min]

  result=K.ch

  return(result)
}

################################################################################
################################################################################

zerofun <- function(x) {
  ################################################################################
  ################################################################################
  #### Ensures that every coordinate of x remains non-negative
  ####
  #### x : Observed Trajectory
  ################################################################################
  ################################################################################

  # Replace negative values with 0.01
  result <- ifelse(x > 0, x, 0.01)

  return(result)
}

################################################################################
################################################################################

diffspline <- function(x, a_hat, Inf_bound, Sup_bound, K, M, Lconst) {
  ################################################################################
  ################################################################################
  #### Computes the Drift estimator
  ####
  #### x : Observed trajectory
  #### a_hat : Previously computed coefficients of the B-Spline Basis
  #### Inf_bound : Lower bound of the interval considered
  #### Sup_bound : Upper bound of the interval considered
  #### K : Number of non-zero knots
  #### M : Such that K+M = desired number of B-Spline functions
  #### Lconst : Constant defining the Upper Bound of Approximation
  ################################################################################
  ################################################################################
  # Apply labfun to determine which elements are within the bounds
  lab <- labfun(x, Inf_bound, Sup_bound)
  idmat <- diag(lab)

  # Compute the B-spline basis matrix
  B <- bspline(x, Inf_bound, Sup_bound, K, M)

  # Compute the drift estimator
  if (length(idmat) != 0) {
    Bs <- idmat %*% B
    b_hat <- Bs %*% a_hat
  } else {
    b_hat <- rep(0, length(a_hat))
  }

  # Apply the estimation bounds
  b_hat <- estimbound(b_hat, Lconst)

  return(b_hat)
}

diffspline_authors <- function(x, alpha.ch, lower, upper, K, M, L.N){
  lab <- labfun(x, lower, upper)
  id.mat <- diag(lab)
  bsMat = bspline_authors(x, lower, upper, K, M)                # Use of the bspline functions
  if(length(id.mat)!=0) bs.Mat <- id.mat%*%bsMat
  if(length(id.mat)==0) bs.Mat <- 0*bsMat
  s.values = bs.Mat%*%alpha.ch
  result <- zerofun(s.values); result <- estimbound(result,L.N)

  return(as.vector(result))
}
