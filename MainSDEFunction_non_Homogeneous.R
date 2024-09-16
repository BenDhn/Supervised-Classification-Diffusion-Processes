source(file = "./UsefulFunctions_Homogeneous.R")
source(file = "./UsefulFunctions_non_Homogeneous.R")

sde_classif_time <- function(Xtrain, ytrain, Xtest, ytest, SetK, labels, M = 3,
                        plot = TRUE, overfit_check = FALSE, norm_fun = FALSE) {
  # Estimation of drift functions
  NbClass <- length(labels)
  Drift_a_hat <- list()  # List of estimated coefficients for the drift per class
  Kdrift <- integer(NbClass)  # List for selected value of dimension K per class
  N <- ncol(Xtrain)  # Number of observed trajectories
  n <- nrow(Xtrain)  # Number of observations per trajectory
  TimeStep <- 1 / n  # Time stamp
  L <- log(N)  # Bound corresponding in the paper to ANi.log(Ni)
  classified_data_list <- list()  # List of iXtrain so that we don't need to compute them again
  b_spline_basis_drift_list <- list()
  times <- seq(0, 1, length.out = n)
  cat('\n-----Beginning Drift estimation-----\n')

################################################################################
###################### Estimation of drift coefficients ########################
################################################################################

  for (i in seq_len(NbClass)) {
    cat('\n-----Class number', i, '-----\n')
    iXtrain <- Xtrain[, which(ytrain == labels[i])]
    classified_data_list[[i]] <- iXtrain
    N_i <- ncol(iXtrain)
    SupInt <- log(N_i)
    InfInt <- -SupInt
    # iZ <- zfun_bis(iXtrain, TimeStep)
    iZ <- zfun(iXtrain, TimeStep)
    cat('\n-----K value selection-----\n')
    Kdrift[i] <- selectdimdrift_time(iXtrain[1:(n-1), ], times[1:(n-1)], iZ,
                                     0.1, SetK, InfInt, SupInt, TimeStep, M, L)
    b_spline_basis_drift_list[[i]] <- create_bspline_basis(rangeval = c(InfInt, SupInt),
                                                           nbasis = Kdrift[i] + M,
                                                           degree = M)
    cat('\n-----Selected value K : ', Kdrift[i], '-----\n')
    iB <- bfun_time(iXtrain[1:(n-1), ], times[1:(n-1)],
                    b_spline_basis_drift_list[[i]], InfInt,
                    SupInt, Kdrift[i], M)
    Drift_a_hat[[i]] <- optimfun(iB, iZ, Kdrift[i], M, L)
  }

################################################################################
################### Estimation of the diffusion coefficient ####################
################################################################################
  cat('\n-----Beginning Diffusion estimation-----\n')
  U <- ufun(Xtrain, TimeStep)
  SupInt_ <- log(N)
  InfInt_ <- -SupInt_
  cat('\n-----K value selection-----\n')
  Kdiff <- selectdimdiff_time(Xtrain[1:(n-1), ], times[1:(n-1)], U, 5,
                              SetK, InfInt_, SupInt_, TimeStep, M, L)
  b_spline_basis_diff <- create_bspline_basis(rangeval = c(InfInt, SupInt),
                                              nbasis = Kdiff + M,
                                              degree = M)
  cat('\n-----Selected value K : ', Kdiff, '-----\n')
  B <- bfun_time(Xtrain[1:(n-1), ], times[1:(n-1)],
                 b_spline_basis_diff, InfInt_, SupInt_, Kdiff, M)
  Diff_a_hat <- optimfun(B, U, Kdiff, M, L)

################################################################################
########################## Prediction of labels ################################
################################################################################
  cat('\n-----Computation of probabilities (test)-----\n')
  n <- nrow(Xtest)
  ExpFspline <- matrix(0, nrow = ncol(Xtest), ncol = NbClass)
  VectorPi <- matrix(0, nrow = ncol(Xtest), ncol = NbClass)
  VectorProba <- sapply(labels, function(k) sum(ytrain == k) / length(ytrain))

  #for (j in pblapply(seq_len(ncol(Xtest)), identity)) {
  for (j in seq_len(ncol(Xtest))) {
    for (k in seq_len(NbClass)) {
      kXtrain <- classified_data_list[[k]]
      SupInt <- log(ncol(kXtrain))
      InfInt <- -SupInt
      drift_estim <- driftspline_time(Xtest[1:(n-1), j],times[1:(n-1)], Drift_a_hat[[k]],
                                      b_spline_basis_drift_list[[k]],
                                      InfInt, SupInt, L)
      diff_estim <- diffspline_time(Xtest[1:(n-1), j], times[1:(n-1)], Diff_a_hat,
                                    b_spline_basis_diff, InfInt_,
                                    SupInt_, L)
      result <- ((drift_estim / diff_estim) * diff(Xtest[, j]) - (TimeStep / 2) * (drift_estim^2 / diff_estim))
      scale_factor <- 100
      ExpFspline[j, k] <- VectorProba[k] * exp(sum(result / scale_factor))
    }
    VectorPi[j, ] <- ExpFspline[j, ] / sum(ExpFspline[j, ])
  }


  cat('\n-----Prediction of labels-----\n')
  argmax_vector <- max.col(VectorPi)
  PredClass <- labels[argmax_vector]

  cat('\n-----Computation of probabilities (train)-----\n')
  # Same thing on the training set if
  if (overfit_check) {
    n <- nrow(Xtrain)
    ExpFspline <- matrix(0, nrow = ncol(Xtrain), ncol = NbClass)
    VectorPi <- matrix(0, nrow = ncol(Xtrain), ncol = NbClass)
    VectorProba <- sapply(labels, function(k) sum(ytrain == k) / length(ytrain))
    norm_list_drift <- array(0, dim = c(ncol(Xtrain), NbClass))
    norm_list_diff <- list()

    for (j in pblapply(seq_len(ncol(Xtrain)), identity)) {
      for (k in seq_len(NbClass)) {
        kXtrain <- classified_data_list[[k]]
        SupInt <- log(ncol(kXtrain))
        InfInt <- -SupInt
        drift_estim <- driftspline_time(Xtrain[1:(n-1), j], times[1:(n-1)], Drift_a_hat[[k]],
                                        b_spline_basis_drift_list[[k]],
                                        InfInt, SupInt, L)
        diff_estim <- diffspline_time(Xtrain[1:(n-1), j], times[1:(n-1)], Diff_a_hat,
                                      b_spline_basis_diff, InfInt_,
                                      SupInt_, L)

        real_drift <- sapply(Xtrain[1:(n - 1), j], function(x) b(x, theta[k]))
        real_diff <- sapply(Xtrain[1:(n - 1), j], function(x) sigma(x))

        drift_norm2 <- sqrt(sum((drift_estim - real_drift)^2))
        diff_norm2 <- sqrt(sum((diff_estim - real_diff)^2))

        norm_list_drift[j,k] <- drift_norm2
        norm_list_diff <- c(norm_list_diff, diff_norm2)

        result <- ((drift_estim / diff_estim) * diff(Xtrain[, j]) - (TimeStep / 2) * (drift_estim^2 / diff_estim))
        scale_factor <- 100
        ExpFspline[j, k] <- VectorProba[k] * exp(sum(result / scale_factor))
      }
      VectorPi[j, ] <- ExpFspline[j, ] / sum(ExpFspline[j, ])
    }


    cat('\n-----Prediction of labels-----\n')
    argmax_vector_train <- max.col(VectorPi)
    PredClass_train <- labels[argmax_vector_train]
  }
################################################################################
################################ Plotting ######################################
################################################################################
  if (plot) {
    cat('\n-----Plotting results-----\n')
    time_grid <- seq(0, 1, length.out = n)

    # Prepare the test data for ggplot
    test_data <- data.frame(Time = rep(time_grid, ncol(Xtest)),
                            Value = as.vector(Xtest),
                            Class = factor(rep(ytest, each = length(time_grid))),
                            Group = rep(1:ncol(Xtest), each = length(time_grid)))  # Group to differentiate lines

    # Plot testing data using ggplot2
    p1 <- ggplot(test_data, aes(x = Time, y = Value, color = Class, group = interaction(Class, Group))) +
      geom_line() +
      scale_color_manual(values = rainbow(length(labels))) +
      labs(title = "Test Data Simulation for Each Class", x = "Time Steps", y = "Values of the Mixed process") +
      theme_minimal() +
      theme(legend.position = "right",
            panel.background = element_rect(fill = "gray90", color = NA),  # Change panel background to gray
            plot.background = element_rect(fill = "gray95", color = NA)) +
      guides(color = guide_legend(title = "Testing Class"))

    # Add horizontal lines for U values
    for (i in seq_len(NbClass)) {
      U <- seq(InfInt, SupInt, length.out = Kdrift[i] + M)
      p1 <- p1 + geom_hline(yintercept = U, color = rainbow(NbClass)[i], linetype = "dashed")
    }

    print(p1)

    # Prepare data for predicted classes on test set
    pred_test_data <- test_data
    pred_test_data$Predicted_Class <- factor(rep(PredClass, each = length(time_grid)), levels = labels) # Ensure factor levels match

    # Plot predictions using ggplot2
    p2 <- ggplot(pred_test_data, aes(x = Time, y = Value, color = Predicted_Class, group = interaction(Predicted_Class, Group))) +
      geom_line() +
      scale_color_manual(values = rainbow(length(labels))) +
      labs(title = "Predictions", x = "Time Steps", y = "Value") +
      theme_minimal() +
      theme(legend.position = "right",
            panel.background = element_rect(fill = "gray90", color = NA),  # Change panel background to gray
            plot.background = element_rect(fill = "gray95", color = NA)) +
      guides(color = guide_legend(title = "Predicted Class"))

    print(p2)

    if (overfit_check) {
      # Prepare the training data for ggplot
      train_data <- data.frame(Time = rep(time_grid, ncol(Xtrain)),
                               Value = as.vector(Xtrain),
                               Class = factor(rep(ytrain, each = length(time_grid))),
                               Group = rep(1:ncol(Xtrain), each = length(time_grid)))  # Group to differentiate lines

      # Plot training data using ggplot2
      p3 <- ggplot(train_data, aes(x = Time, y = Value, color = Class, group = interaction(Class, Group))) +
        geom_line() +
        scale_color_manual(values = rainbow(NbClass)) +
        labs(title = "Train Data Simulation for Each Class", x = "Time Steps", y = "Values of the Mixed process") +
        theme_minimal() +
        theme(legend.position = "right",
              panel.background = element_rect(fill = "gray90", color = NA),  # Change panel background to gray
              plot.background = element_rect(fill = "gray95", color = NA)) +
        guides(color = guide_legend(title = "Training Class"))

      # Add horizontal lines for U values
      for (i in seq_len(NbClass)) {
        U <- seq(InfInt, SupInt, length.out = Kdiff + M)
        p3 <- p3 + geom_hline(yintercept = U, color = rainbow(NbClass)[i], linetype = "dashed")
      }

      print(p3)

      # Prepare data for predicted classes on train set
      pred_train_data <- train_data
      pred_train_data$Predicted_Class <- factor(rep(PredClass_train, each = length(time_grid)), levels = labels) # Ensure factor levels match

      # Plot training predictions using ggplot2
      p4 <- ggplot(pred_train_data, aes(x = Time, y = Value, color = Predicted_Class, group = interaction(Predicted_Class, Group))) +
        geom_line() +
        scale_color_manual(values = rainbow(NbClass)) +
        labs(title = "Train Predictions", x = "Time Steps", y = "Value") +
        theme_minimal() +
        theme(legend.position = "right",
              panel.background = element_rect(fill = "gray90", color = NA),  # Change panel background to gray
              plot.background = element_rect(fill = "gray95", color = NA)) +
        guides(color = guide_legend(title = "Predicted Train Class"))

      print(p4)
    }
  }

################################################################################
############################### Return values ##################################
################################################################################

  if (overfit_check) {
    if (norm_fun) {
      return(list(PredClass = PredClass, PredClass_train = PredClass_train,
                  NormsDrift = norm_list_drift, NormsDiff = norm_list_diff)
      )
    }
    return(list(PredClass = PredClass, PredClass_train = PredClass_train))
  }
  else {
    return(PredClass)
  }

}


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################



sde_classif_time_K <- function(Xtrain, ytrain, Xtest, ytest, K, labels, M = 3,
                             plot = TRUE, overfit_check = FALSE, norm_fun = FALSE) {
  # Estimation of drift functions
  NbClass <- length(labels)
  Drift_a_hat <- list()  # List of estimated coefficients for the drift per class
  Kdrift <- integer(NbClass)  # List for selected value of dimension K per class
  N <- ncol(Xtrain)  # Number of observed trajectories
  n <- nrow(Xtrain)  # Number of observations per trajectory
  TimeStep <- 1 / n  # Time stamp
  L <- log(N)  # Bound corresponding in the paper to ANi.log(Ni)
  classified_data_list <- list()  # List of iXtrain so that we don't need to compute them again
  b_spline_basis_drift_list <- list()
  times <- seq(0, 1, length.out = n)
  cat('\n-----Beginning Drift estimation-----\n')

  ################################################################################
  ###################### Estimation of drift coefficients ########################
  ################################################################################

  for (i in seq_len(NbClass)) {
    cat('\n-----Class number', i, '-----\n')
    cat('\n-----K value :', K, '-----\n')
    iXtrain <- Xtrain[, ytrain == labels[i]]
    classified_data_list[[i]] <- iXtrain
    N_i <- ncol(iXtrain)
    SupInt <- log(N_i)
    InfInt <- -SupInt
    iZ <- zfun(iXtrain, TimeStep)
    b_spline_basis_drift_list[[i]] <- create_bspline_basis(rangeval = c(InfInt, SupInt),
                                                           nbasis = K + M,
                                                           degree = M)
    #iB <- bfun(iXtrain[1:(n-1), ], InfInt, SupInt, Kdrift[i], M)
    iB <- bfun_time(iXtrain[1:(n-1), ], times[1:(n-1)], b_spline_basis_drift_list[[i]], InfInt, SupInt, K, M)
    #cat('NA : ', sum(is.na(iB)))
    Drift_a_hat[[i]] <- optimfun(iB, iZ, K, M, L)
    #cat('NA : ', sum(is.na(Drift_a_hat[[i]])))
  }

  ################################################################################
  ################### Estimation of the diffusion coefficient ####################
  ################################################################################
  cat('\n-----Beginning Diffusion estimation-----\n')
  cat('\n-----K value :', K, '-----\n')
  U <- ufun(Xtrain, TimeStep)
  SupInt_ <- log(N)
  InfInt_ <- -SupInt_
  b_spline_basis_diff <- create_bspline_basis(rangeval = c(InfInt, SupInt),
                                              nbasis = K + M,
                                              degree = M)
  B <- bfun_time(Xtrain[1:(n-1), ], times[1:(n-1)], b_spline_basis_diff, InfInt_, SupInt_, K, M)
  Diff_a_hat <- optimfun(B, U, K, M, L)

  ################################################################################
  ########################## Prediction of labels ################################
  ################################################################################
  cat('\n-----Computation of probabilities (test)-----\n')
  n <- nrow(Xtest)
  ExpFspline <- matrix(0, nrow = ncol(Xtest), ncol = NbClass)
  VectorPi <- matrix(0, nrow = ncol(Xtest), ncol = NbClass)
  VectorProba <- sapply(labels, function(k) sum(ytrain == k) / length(ytrain))

  #for (j in pblapply(seq_len(ncol(Xtest)), identity)) {
  for (j in seq_len(ncol(Xtest))) {
    for (k in seq_len(NbClass)) {
      kXtrain <- classified_data_list[[k]]
      SupInt <- log(ncol(kXtrain))
      InfInt <- -SupInt
      drift_estim <- driftspline_time(Xtest[1:(n-1), j],times[1:(n-1)], Drift_a_hat[[k]],
                                      b_spline_basis_drift_list[[k]], InfInt, SupInt, L)
      diff_estim <- diffspline_time(Xtest[1:(n-1), j], times[1:(n-1)], Diff_a_hat,
                                    b_spline_basis_diff, InfInt_, SupInt_, L)
      result <- ((drift_estim / diff_estim) * diff(Xtest[, j]) - (TimeStep / 2) * (drift_estim^2 / diff_estim))
      #scale_factor <- 100
      #ExpFspline[j, k] <- VectorProba[k] * exp(sum(result / scale_factor))
      ExpFspline[j, k] <- VectorProba[k] * exp(sum(result))
    }
    VectorPi[j, ] <- ExpFspline[j, ] / sum(ExpFspline[j, ])
    #cat('Number of NAs for obs',j,':', sum(is.na(ExpFspline[j, ])),'\n')
  }


  cat('\n-----Prediction of labels-----\n')
  argmax_vector <- max.col(VectorPi)
  PredClass <- labels[argmax_vector]

  cat('\n-----Computation of probabilities (train)-----\n')
  # Same thing on the training set if
  if (overfit_check) {
    n <- nrow(Xtrain)
    ExpFspline <- matrix(0, nrow = ncol(Xtrain), ncol = NbClass)
    VectorPi <- matrix(0, nrow = ncol(Xtrain), ncol = NbClass)
    VectorProba <- sapply(labels, function(k) sum(ytrain == k) / length(ytrain))
    norm_list_drift <- array(0, dim = c(ncol(Xtrain), NbClass))
    norm_list_diff <- list()

    for (j in pblapply(seq_len(ncol(Xtrain)), identity)) {
      for (k in seq_len(NbClass)) {
        kXtrain <- classified_data_list[[k]]
        SupInt <- log(ncol(kXtrain))
        InfInt <- -SupInt
        drift_estim <- driftspline_time(Xtrain[1:(n-1), j], times[1:(n-1)], Drift_a_hat[[k]],
                                        b_spline_basis_drift_list[[k]], InfInt, SupInt, L)
        diff_estim <- diffspline_time(Xtrain[1:(n-1), j], times[1:(n-1)], Diff_a_hat,
                                      b_spline_basis_diff, InfInt_, SupInt_, L)

        real_drift <- sapply(Xtrain[1:(n - 1), j], function(x) b(x, theta[k]))
        real_diff <- sapply(Xtrain[1:(n - 1), j], function(x) sigma(x))

        drift_norm2 <- sqrt(sum((drift_estim - real_drift)^2))
        diff_norm2 <- sqrt(sum((diff_estim - real_diff)^2))

        norm_list_drift[j,k] <- drift_norm2
        norm_list_diff <- c(norm_list_diff, diff_norm2)

        result <- ((drift_estim / diff_estim) * diff(Xtrain[, j]) - (TimeStep / 2) * (drift_estim^2 / diff_estim))
        scale_factor <- 100
        ExpFspline[j, k] <- VectorProba[k] * exp(sum(result / scale_factor))
      }
      VectorPi[j, ] <- ExpFspline[j, ] / sum(ExpFspline[j, ])
    }


    cat('\n-----Prediction of labels-----\n')
    argmax_vector_train <- max.col(VectorPi)
    PredClass_train <- labels[argmax_vector_train]
  }
  ################################################################################
  ################################ Plotting ######################################
  ################################################################################
  if (plot) {
    cat('\n-----Plotting results-----\n')
    time_grid <- seq(0, 1, length.out = n)

    # Prepare the test data for ggplot
    test_data <- data.frame(Time = rep(time_grid, ncol(Xtest)),
                            Value = as.vector(Xtest),
                            Class = factor(rep(ytest, each = length(time_grid))),
                            Group = rep(1:ncol(Xtest), each = length(time_grid)))  # Group to differentiate lines

    # Plot testing data using ggplot2
    p1 <- ggplot(test_data, aes(x = Time, y = Value, color = Class, group = interaction(Class, Group))) +
      geom_line() +
      scale_color_manual(values = rainbow(length(labels))) +
      labs(title = "Test Data Simulation for Each Class", x = "Time Steps", y = "Values of the Mixed process") +
      theme_minimal() +
      theme(legend.position = "right",
            panel.background = element_rect(fill = "gray90", color = NA),  # Change panel background to gray
            plot.background = element_rect(fill = "gray95", color = NA)) +
      guides(color = guide_legend(title = "Testing Class"))

    # Add horizontal lines for U values
    for (i in seq_len(NbClass)) {
      U <- seq(InfInt, SupInt, length.out = K + M)
      p1 <- p1 + geom_hline(yintercept = U, color = rainbow(NbClass)[i], linetype = "dashed")
    }

    print(p1)

    # Prepare data for predicted classes on test set
    pred_test_data <- test_data
    pred_test_data$Predicted_Class <- factor(rep(PredClass, each = length(time_grid)), levels = labels) # Ensure factor levels match

    # Plot predictions using ggplot2
    p2 <- ggplot(pred_test_data, aes(x = Time, y = Value, color = Predicted_Class, group = interaction(Predicted_Class, Group))) +
      geom_line() +
      scale_color_manual(values = rainbow(length(labels))) +
      labs(title = "Predictions", x = "Time Steps", y = "Value") +
      theme_minimal() +
      theme(legend.position = "right",
            panel.background = element_rect(fill = "gray90", color = NA),  # Change panel background to gray
            plot.background = element_rect(fill = "gray95", color = NA)) +
      guides(color = guide_legend(title = "Predicted Class"))

    print(p2)

    if (overfit_check) {
      # Prepare the training data for ggplot
      train_data <- data.frame(Time = rep(time_grid, ncol(Xtrain)),
                               Value = as.vector(Xtrain),
                               Class = factor(rep(ytrain, each = length(time_grid))),
                               Group = rep(1:ncol(Xtrain), each = length(time_grid)))  # Group to differentiate lines

      # Plot training data using ggplot2
      p3 <- ggplot(train_data, aes(x = Time, y = Value, color = Class, group = interaction(Class, Group))) +
        geom_line() +
        scale_color_manual(values = rainbow(NbClass)) +
        labs(title = "Train Data Simulation for Each Class", x = "Time Steps", y = "Values of the Mixed process") +
        theme_minimal() +
        theme(legend.position = "right",
              panel.background = element_rect(fill = "gray90", color = NA),  # Change panel background to gray
              plot.background = element_rect(fill = "gray95", color = NA)) +
        guides(color = guide_legend(title = "Training Class"))

      # Add horizontal lines for U values
      for (i in seq_len(NbClass)) {
        U <- seq(InfInt, SupInt, length.out = K + M)
        p3 <- p3 + geom_hline(yintercept = U, color = rainbow(NbClass)[i], linetype = "dashed")
      }

      print(p3)

      # Prepare data for predicted classes on train set
      pred_train_data <- train_data
      pred_train_data$Predicted_Class <- factor(rep(PredClass_train, each = length(time_grid)), levels = labels) # Ensure factor levels match

      # Plot training predictions using ggplot2
      p4 <- ggplot(pred_train_data, aes(x = Time, y = Value, color = Predicted_Class, group = interaction(Predicted_Class, Group))) +
        geom_line() +
        scale_color_manual(values = rainbow(NbClass)) +
        labs(title = "Train Predictions", x = "Time Steps", y = "Value") +
        theme_minimal() +
        theme(legend.position = "right",
              panel.background = element_rect(fill = "gray90", color = NA),  # Change panel background to gray
              plot.background = element_rect(fill = "gray95", color = NA)) +
        guides(color = guide_legend(title = "Predicted Train Class"))

      print(p4)
    }
  }

  ################################################################################
  ############################### Return values ##################################
  ################################################################################

  if (overfit_check) {
    if (norm_fun) {
      return(list(PredClass = PredClass, PredClass_train = PredClass_train,
                  NormsDrift = norm_list_drift, NormsDiff = norm_list_diff)
      )
    }
    return(list(PredClass = PredClass, PredClass_train = PredClass_train))
  }
  else {
    return(PredClass)
  }

}
