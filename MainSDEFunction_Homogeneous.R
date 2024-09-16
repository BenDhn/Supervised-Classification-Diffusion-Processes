# sde_classif <- function(Xtrain, ytrain, Xtest, ytest, SetK, labels, M = 3, plot = TRUE, overfit_check = FALSE) {
#   # Estimation of drift functions
#   NbClass <- length(labels)
#   Drift_a_hat <- list()  # List of estimated coefficients for the drift per class
#   Kdrift <- integer(NbClass)  # List for selected value of dimension K per class
#   N <- ncol(Xtrain)  # Number of observed trajectories
#   n <- nrow(Xtrain)  # Number of observations per trajectory
#   TimeStep <- 1 / n  # Time stamp
#   L <- log(N)  # Bound corresponding in the paper to ANi.log(Ni)
#   classified_data_list <- list()  # List of iXtrain so that we don't need to compute them again
#   cat('\n-----Beginning Drift estimation-----\n')
#
#   for (i in seq_len(NbClass)) {
#     cat('\n-----Class number', i, '-----\n')
#     iXtrain <- Xtrain[, ytrain == labels[i]]
#     classified_data_list[[i]] <- iXtrain
#     N_i <- ncol(iXtrain)
#     SupInt <- log(N_i)
#     InfInt <- -SupInt
#     iZ <- zfun(iXtrain, TimeStep)
#     Kdrift[i] <- selectdimdrift(iXtrain[1:(n-1), ], iZ, 0.1, SetK, M, InfInt, SupInt, TimeStep, L)
#     cat('\n-----Selected value K : ', Kdrift[i], '-----\n')
#     #iB <- bfun(iXtrain[1:(n-1), ], InfInt, SupInt, Kdrift[i], M)
#     iB <- bfun(iXtrain[1:(n-1), ], InfInt, SupInt, Kdrift[i], M)
#     Drift_a_hat[[i]] <- optimfun(iB, iZ, Kdrift[i], M, L)
#   }
#
#   # Estimation of the diffusion coefficient
#   cat('\n-----Beginning Diffusion estimation-----\n')
#   U <- ufun(Xtrain, TimeStep)
#   SupInt_ <- log(N)
#   InfInt_ <- -SupInt_
#   Kdiff <- selectdimdiff(Xtrain[1:(n-1), ], U, 5, SetK, M, InfInt_, SupInt_, TimeStep, L)
#   cat('\n-----Selected value K : ', Kdiff, '-----\n')
#   B <- bfun(Xtrain[1:(n-1), ], InfInt_, SupInt_, Kdiff, M)
#   Diff_a_hat <- optimfun(B, U, Kdiff, M, L)
#
#   # Computation of probabilities
#   cat('\n-----Computation of probabilities-----\n')
#   n <- nrow(Xtest)
#   ExpFspline <- matrix(0, nrow = ncol(Xtest), ncol = NbClass)
#   VectorPi <- matrix(0, nrow = ncol(Xtest), ncol = NbClass)
#   VectorProba <- sapply(labels, function(k) sum(ytrain == k) / length(ytrain))
#
#   for (j in seq_len(ncol(Xtest))) {
#     for (k in seq_len(NbClass)) {
#       kXtrain <- classified_data_list[[k]]
#       SupInt <- log(ncol(kXtrain))
#       InfInt <- -SupInt
#       drift_estim <- driftspline(Xtest[1:(n-1), j], Drift_a_hat[[k]],
#                                  InfInt, SupInt, Kdrift[k], M, L)
#       diff_estim <- diffspline(Xtest[1:(n-1), j], Diff_a_hat,
#                                InfInt_, SupInt_, Kdiff, M, L)
#       result <- ((drift_estim / diff_estim) * diff(Xtest[, j]) - (TimeStep / 2) * (drift_estim^2 / diff_estim))
#       scale_factor <- 100
#       ExpFspline[j, k] <- VectorProba[k] * exp(sum(result / scale_factor))
#     }
#     VectorPi[j, ] <- ExpFspline[j, ] / sum(ExpFspline[j, ])
#   }
#
#   # Prediction of labels
#   cat('\n-----Prediction of labels-----\n')
#   argmax_vector <- max.col(VectorPi)
#   PredClass <- labels[argmax_vector]
#
#   if (overfit_check) {
#     n <- nrow(Xtrain)
#     ExpFspline <- matrix(0, nrow = ncol(Xtrain), ncol = NbClass)
#     VectorPi <- matrix(0, nrow = ncol(Xtrain), ncol = NbClass)
#     VectorProba <- sapply(labels, function(k) sum(ytrain == k) / length(ytrain))
#
#     for (j in seq_len(ncol(Xtrain))) {
#       for (k in seq_len(NbClass)) {
#         kXtrain <- classified_data_list[[k]]
#         SupInt <- log(ncol(kXtrain))
#         InfInt <- -SupInt
#         drift_estim <- driftspline(Xtrain[1:(n-1), j], Drift_a_hat[[k]], InfInt, SupInt, Kdrift[k], M, L)
#         diff_estim <- diffspline(Xtrain[1:(n-1), j], Diff_a_hat, InfInt_, SupInt_, Kdiff, M, L)
#         result <- ((drift_estim / diff_estim) * diff(Xtrain[, j]) - (TimeStep / 2) * (drift_estim^2 / diff_estim))
#         scale_factor <- 100
#         ExpFspline[j, k] <- VectorProba[k] * exp(sum(result / scale_factor))
#       }
#       VectorPi[j, ] <- ExpFspline[j, ] / sum(ExpFspline[j, ])
#     }
#
#     cat('\n-----Prediction of labels (training)-----\n')
#     argmax_vector_train <- max.col(VectorPi)
#     PredClass_train <- labels[argmax_vector_train]
#   }
#
#   # Plotting
#   if (plot) {
#     cat('\n-----Plotting results-----\n')
#     time_grid <- seq(0, 1, length.out = n)
#
#     par(mfrow = c(1, 2))
#
#     # Plot testing data
#     for (i in seq_len(NbClass)) {
#       class_data <- Xtest[, ytest == labels[i]]
#       color_i <- rainbow(NbClass)[i]
#       plot(time_grid, class_data[, 1], type = 'l', col = color_i, ylim = range(Xtest),
#            xlab = 'Time Steps', ylab = 'Values of the Mixed process',
#            main = 'Test Data Simulation for Each Class')
#       U <- seq(InfInt, SupInt, length.out = Kdrift[i] + M)
#       abline(h = U, col = color_i, lty = 2)
#       for (j in seq(2, ncol(class_data))) {
#         lines(time_grid, class_data[, j], col = color_i)
#       }
#     }
#     legend("topright", legend = paste("Testing Class", labels), col = rainbow(NbClass), lty = 1)
#
#     # Plot predictions
#     plot(time_grid, Xtest[, 1], type = 'l', col = rainbow(NbClass)[PredClass[1]],
#          ylim = range(Xtest), xlab = 'Time Steps', ylab = 'Value', main = 'Predictions')
#     for (j in seq(2, ncol(Xtest))) {
#       lines(time_grid, Xtest[, j], col = rainbow(NbClass)[PredClass[j]])
#     }
#     legend("topright", legend = paste("Predicted Class", labels), col = rainbow(NbClass), lty = 1)
#
#     if (overfit_check) {
#       par(mfrow = c(1, 2))
#
#       # Plot training data
#       for (i in seq_len(NbClass)) {
#         class_data <- Xtrain[, ytrain == labels[i]]
#         color_i <- rainbow(NbClass)[i]
#         plot(time_grid, class_data[, 1], type = 'l', col = color_i, ylim = range(Xtrain),
#              xlab = 'Time Steps', ylab = 'Values of the Mixed process',
#              main = 'Train Data Simulation for Each Class')
#         U <- seq(InfInt, SupInt, length.out = Kdiff + M)
#         abline(h = U, col = color_i, lty = 2)
#         for (j in seq(2, ncol(class_data))) {
#           lines(time_grid, class_data[, j], col = color_i)
#         }
#       }
#       legend("topright", legend = paste("Training Class", labels), col = rainbow(NbClass), lty = 1, cex = 0.5, bty = "n")
#
#       # Plot training predictions
#       plot(time_grid, Xtrain[, 1], type = 'l', col = rainbow(NbClass)[PredClass_train[1]],
#            ylim = range(Xtrain), xlab = 'Time Steps', ylab = 'Value', main = 'Train Predictions')
#       for (j in seq(2, ncol(Xtrain))) {
#         lines(time_grid, Xtrain[, j], col = rainbow(NbClass)[PredClass_train[j]])
#       }
#       legend("topright", legend = paste("Predicted Train Class", labels), col = rainbow(NbClass), lty = 1, cex = 0.5, bty = "n")
#     }
#   }
#
#   if (overfit_check) {
#     return(list(PredClass = PredClass, PredClass_train = PredClass_train))
#   } else {
#     return(PredClass)
#     #return(list(PredClass = PredClass, Drift_a_hat = Drift_a_hat, Diff_a_hat = Diff_a_hat, B_spline = )
#   }
# }


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
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


sde_classif_authors <- function(Xtrain, ytrain, Xtest, ytest, SetK, labels, M = 3,
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
  cat('\n-----Beginning Drift estimation-----\n')

  for (i in seq_len(NbClass)) {
    cat('\n-----Class number', i, '-----\n')
    iXtrain <- Xtrain[, which(ytrain == labels[i])]
    classified_data_list[[i]] <- iXtrain
    N_i <- ncol(iXtrain)
    SupInt <- log(N_i)
    InfInt <- -SupInt
    iZ <- zfun_authors(iXtrain, TimeStep)
    cat('\n-----K value selection-----\n')
    Kdrift[i] <- selectdimdrift_authors(iXtrain, 0.1, SetK,
                                        M, InfInt, SupInt, TimeStep, L)
    cat('\n-----Selected value K : ', Kdrift[i], '-----\n')
    #iB <- bfun(iXtrain[1:(n-1), ], InfInt, SupInt, Kdrift[i], M)
    iB <- bfun_authors(iXtrain, InfInt, SupInt, Kdrift[i], M)
    Drift_a_hat[[i]] <- optimfun(iB, iZ, Kdrift[i], M, L)
  }

  # Estimation of the diffusion coefficient
  cat('\n-----Beginning Diffusion estimation-----\n')
  U <- ufun_authors(Xtrain, TimeStep)
  SupInt_ <- log(N)
  InfInt_ <- -SupInt_
  cat('\n-----K value selection-----\n')
  Kdiff <- selectdimdiff_authors(iXtrain, 5, SetK,
                                 M, InfInt, SupInt, TimeStep, L)
  cat('\n-----Selected value K : ', Kdiff, '-----\n')
  B <- bfun_authors(Xtrain, InfInt_, SupInt_, Kdiff, M)
  Diff_a_hat <- optimfun(B, U, Kdiff, M, L)

  # Computation of probabilities
  cat('\n-----Computation of probabilities (test)-----\n')
  n <- nrow(Xtest)
  ExpFspline <- matrix(0, nrow = ncol(Xtest), ncol = NbClass)
  VectorPi <- matrix(0, nrow = ncol(Xtest), ncol = NbClass)
  VectorProba <- sapply(labels, function(k) sum(ytrain == k) / length(ytrain))

  for (j in pblapply(seq_len(ncol(Xtest)), identity)) {
    for (k in seq_len(NbClass)) {
      kXtrain <- classified_data_list[[k]]
      drift_estim <- driftspline_authors(Xtest[1:(n-2), j], Drift_a_hat[[k]],
                                         InfInt, SupInt, Kdrift[k], M, L)
      diff_estim <- diffspline_authors(Xtest[1:(n-2), j], Diff_a_hat,
                                       InfInt_, SupInt_, Kdiff, M, L)
      result <- ((drift_estim / diff_estim) * diff(Xtest[, j])[-1] - (TimeStep / 2) * (drift_estim^2 / diff_estim))
      ExpFspline[j, k] <- VectorProba[k] * exp(sum(result))
    }
    VectorPi[j, ] <- ExpFspline[j, ] / sum(ExpFspline[j, ])
  }


  # Prediction of labels
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
        drift_estim <- driftspline_authors(Xtrain[1:(n-2), j], Drift_a_hat[[k]],
                                           InfInt, SupInt, Kdrift[k], M, L)
        diff_estim <- diffspline_authors(Xtrain[1:(n-2), j], Diff_a_hat,
                                         InfInt_, SupInt_, Kdiff, M, L)

        real_drift <- sapply(Xtrain[1:(n - 2), j], function(x) b(x, theta[k]))
        real_diff <- sapply(Xtrain[1:(n - 2), j], function(x) sigma(x))

        drift_norm2 <- sqrt(sum((drift_estim - real_drift)^2))
        diff_norm2 <- sqrt(sum((diff_estim - real_diff)^2))

        norm_list_drift[j,k] <- drift_norm2
        norm_list_diff <- c(norm_list_diff, diff_norm2)

        result <- ((drift_estim / diff_estim) * diff(Xtrain[, j])[-1] - (TimeStep / 2) * (drift_estim^2 / diff_estim))
        ExpFspline[j, k] <- VectorProba[k] * exp(sum(result))
      }
      VectorPi[j, ] <- ExpFspline[j, ] / sum(ExpFspline[j, ])
    }


    cat('\n-----Prediction of labels-----\n')
    argmax_vector_train <- max.col(VectorPi)
    PredClass_train <- labels[argmax_vector_train]
  }

  # # Plotting
  # if (plot) {
  #   cat('\n-----Plotting results-----\n')
  #   time_grid <- seq(0, 1, length.out = n)
  #
  #   par(mfrow = c(1, 2))
  #
  #   # Plot testing data
  #   for (i in seq_len(NbClass)) {
  #     class_data <- Xtest[, ytest == labels[i]]
  #     color_i <- rainbow(NbClass)[i]
  #     plot(time_grid, class_data[, 1], type = 'l', col = color_i, ylim = range(Xtest),
  #          xlab = 'Time Steps', ylab = 'Values of the Mixed process',
  #          main = 'Test Data Simulation for Each Class')
  #     U <- seq(InfInt, SupInt, length.out = Kdrift[i] + M)
  #     abline(h = U, col = color_i, lty = 2)
  #     for (j in seq(2, ncol(class_data))) {
  #       lines(time_grid, class_data[, j], col = color_i)
  #     }
  #   }
  #   legend("topright", legend = paste("Testing Class", labels), col = rainbow(NbClass), lty = 1)
  #
  #   # Plot predictions
  #   plot(time_grid, Xtest[, 1], type = 'l', col = rainbow(NbClass)[PredClass[1]],
  #        ylim = range(Xtest), xlab = 'Time Steps', ylab = 'Value', main = 'Predictions')
  #   for (j in seq(2, ncol(Xtest))) {
  #     lines(time_grid, Xtest[, j], col = rainbow(NbClass)[PredClass[j]])
  #   }
  #   legend("topright", legend = paste("Predicted Class", labels), col = rainbow(NbClass), lty = 1)
  #
  #   if (overfit_check) {
  #     par(mfrow = c(1, 2))
  #
  #     # Plot training data
  #     for (i in seq_len(NbClass)) {
  #       class_data <- Xtrain[, ytrain == labels[i]]
  #       color_i <- rainbow(NbClass)[i]
  #       plot(time_grid, class_data[, 1], type = 'l', col = color_i, ylim = range(Xtrain),
  #            xlab = 'Time Steps', ylab = 'Values of the Mixed process',
  #            main = 'Train Data Simulation for Each Class')
  #       U <- seq(InfInt, SupInt, length.out = Kdiff + M)
  #       abline(h = U, col = color_i, lty = 2)
  #       for (j in seq(2, ncol(class_data))) {
  #         lines(time_grid, class_data[, j], col = color_i)
  #       }
  #     }
  #     legend("topright", legend = paste("Training Class", labels), col = rainbow(NbClass), lty = 1, cex = 0.5, bty = "n")
  #
  #     # Plot training predictions
  #     plot(time_grid, Xtrain[, 1], type = 'l', col = rainbow(NbClass)[PredClass_train[1]],
  #          ylim = range(Xtrain), xlab = 'Time Steps', ylab = 'Value', main = 'Train Predictions')
  #     for (j in seq(2, ncol(Xtrain))) {
  #       lines(time_grid, Xtrain[, j], col = rainbow(NbClass)[PredClass_train[j]])
  #     }
  #     legend("topright", legend = paste("Predicted Train Class", labels), col = rainbow(NbClass), lty = 1, cex = 0.5, bty = "n")
  #   }
  # }

  # Plotting
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

  if (overfit_check) {
    if (norm_fun) {
      return(list(PredClass = PredClass, PredClass_train = PredClass_train,
                  NormsDrift = norm_list_drift, NormsDiff = norm_list_diff)
             )
    }
    return(list(PredClass = PredClass, PredClass_train = PredClass_train))
  } else {
      return(PredClass)
    }

}
