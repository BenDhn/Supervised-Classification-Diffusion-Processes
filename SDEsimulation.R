library(ggplot2)
library(reshape2)

sde_sim <- function(t0, T, X0, n, Ni, b, sigma, theta) {
  dt <- (T - t0) / n
  x <- matrix(0, nrow = n, ncol = Ni)
  x[1, ] <- X0
  for (j in 2:n) {
    x[j, ] <- x[j-1, ] + b(x[j-1, ], theta) * dt + sigma(x[j-1, ]) * rnorm(Ni) * sqrt(dt)
  }
  return(x)
}

# Drift and diffusion functions
b <- function(x, theta) {
  theta * (1/4 + (3/4) * (cos(x)^2))
}
sigma <- function(x) {
  0.1 + 0.9 / sqrt(1 + x^2)
}

b_prime <- function(x, theta) {
  theta * (1/4 + (3/4) * (cos(x + (pi/4) * theta)^2))
}
sigma_simple <- function(x) {
  0.1
}

################################################################################
############################ Simulation of Data ################################
################################################################################

simulation_sde <- function(t0, T, x0, n, Ni_train_list, Ni_test_list, labels, b, sigma, theta, plot = TRUE, save_plot = FALSE) {
  ################################################################################
  ################################################################################
  #### Given drift and diffusion coefficients, simulates from the corresponding
  #### SDE
  ####
  #### T, t0 : Starting and ending time points
  #### x0 : Space starting point, i.e deterministic value of X0
  #### n : Discretization stamp
  #### Ni_train_list,labels_train : Arrays such that Ni_train_list[i] = corresponding Ni of the class label[i]
  #### Ni_test_list,labels_test : Same
  #### b,sigma : Respectively drift and diffusion coefficients
  #### theta : Array such that theta[i] is the corresponding theta for class i
  ################################################################################
  ################################################################################

  TimeStep <- 1 / n
  NbClass <- length(labels)

  Xtrain <- sde_sim(t0, T, x0, n, Ni_train_list[1], b, sigma, theta[1])
  ytrain <- rep(labels[1], Ni_train_list[1])
  Xtest <- sde_sim(t0, T, x0, n, Ni_test_list[1], b, sigma, theta[1])
  ytest <- rep(labels[1], Ni_test_list[1])

  for (i in 2:NbClass) {
    Xtrain <- cbind(Xtrain, sde_sim(t0, T, x0, n, Ni_train_list[i], b, sigma, theta[i]))
    ytrain <- c(ytrain, rep(labels[i], Ni_train_list[i]))
    Xtest <- cbind(Xtest, sde_sim(t0, T, x0, n, Ni_test_list[i], b, sigma, theta[i]))
    ytest <- c(ytest, rep(labels[i], Ni_test_list[i]))
  }

  if (plot) {
    time_grid <- seq(t0, T, length.out = n)
    # Prepare the training data for ggplot

    # Prepare the training data for ggplot
    train_data <- data.frame(Time = rep(time_grid, ncol(Xtrain)),
                             Value = as.vector(Xtrain),
                             Class = factor(rep(ytrain, each = length(time_grid))),
                             Group = rep(1:ncol(Xtrain), each = length(time_grid))) # Group to differentiate lines

    # Plot training data using ggplot2
    p1 <- ggplot(train_data, aes(x = Time, y = Value, color = Class, group = interaction(Class, Group))) +
      geom_line() +
      scale_color_manual(values = rainbow(length(labels))) +
      labs(title = "Training Data Simulation for Each Class", x = "Time", y = "Values of the Mixed process") +
      theme_minimal() +
      theme(legend.position = "right",
            panel.background = element_rect(fill = "gray90", color = NA),  # Change panel background to gray
            plot.background = element_rect(fill = "gray95", color = NA)) +
      guides(color = guide_legend(title = "Training Class"))
    print(p1)
    if (save_plot){
      ggsave(filename = "./Plots/train_data.png", plot = p1, width = 8, height = 5, dpi = 300)
    }
    # Prepare the testing data for ggplot
    test_data <- data.frame(Time = rep(time_grid, ncol(Xtest)),
                            Value = as.vector(Xtest),
                            Class = factor(rep(ytest, each = length(time_grid))),
                            Group = rep(1:ncol(Xtest), each = length(time_grid))) # Group to differentiate lines

    # Plot testing data using ggplot2
    p2 <- ggplot(test_data, aes(x = Time, y = Value, color = Class, group = interaction(Class, Group))) +
      geom_line() +
      scale_color_manual(values = rainbow(length(labels))) +
      labs(title = "Testing Data Simulation for Each Class", x = "Time", y = "Value") +
      theme_minimal() +
      theme(legend.position = "right",
            panel.background = element_rect(fill = "gray90", color = NA),  # Change panel background to gray
            plot.background = element_rect(fill = "gray95", color = NA)) +
      guides(color = guide_legend(title = "Testing Class"))
    print(p2)
    if (save_plot){
      ggsave(filename = "./Plots/train_data.png", plot = p2, width = 8, height = 5, dpi = 300)
    }
  }

  return(list(Xtrain = Xtrain, ytrain = ytrain, Xtest = Xtest, ytest = ytest))
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
############### Generalization to non homogeneous diffusions ###################
################################################################################
################################################################################

# Define the SDE simulation function with time dependence
sde_sim_time <- function(t0, T, X0, n, Ni, b, sigma, theta) {
  dt <- (T - t0) / n
  x <- matrix(0, nrow = n, ncol = Ni)
  x[1, ] <- X0

  for (j in 2:n) {
    # Compute the current time step
    current_time <- j * dt
    # Update the process using the drift and diffusion functions
    x[j, ] <- x[j - 1, ] + b(x[j - 1, ], theta, current_time) * dt +
      sigma(x[j - 1, ], current_time) * rnorm(Ni) * sqrt(dt)
  }

  return(x)
}

# Define the functions for drift and diffusion with time dependence
b_time <- function(x, theta, t) {
  t * theta * (1/4 + (3/4) * (cos(x)^2))
}

sigma_time <- function(x, t) {
  0.1 + t*(0.9 / sqrt(1 + x^2))
}

simulation_sde_time <- function(t0, T, x0, n, Ni_train_list, Ni_test_list, labels, b, sigma, theta, plot = TRUE) {
  ################################################################################
  ################################################################################
  #### Given drift and diffusion coefficients, simulates from the corresponding
  #### SDE
  ####
  #### T, t0 : Starting and ending time points
  #### x0 : Space starting point, i.e deterministic value of X0
  #### n : Discretization stamp
  #### Ni_train_list,labels_train : Arrays such that Ni_train_list[i] = corresponding Ni of the class label[i]
  #### Ni_test_list,labels_test : Same
  #### b,sigma : Respectively drift and diffusion coefficients
  #### theta : Array such that theta[i] is the corresponding theta for class i
  ################################################################################
  ################################################################################

  TimeStep <- 1 / n
  NbClass <- length(labels)

  Xtrain <- sde_sim_time(t0, T, x0, n, Ni_train_list[1], b_time, sigma_time, theta[1])
  ytrain <- rep(labels[1], Ni_train_list[1])
  Xtest <- sde_sim_time(t0, T, x0, n, Ni_test_list[1], b_time, sigma_time, theta[1])
  ytest <- rep(labels[1], Ni_test_list[1])

  for (i in 2:NbClass) {
    Xtrain <- cbind(Xtrain, sde_sim_time(t0, T, x0, n, Ni_train_list[i], b_time, sigma_time, theta[i]))
    ytrain <- c(ytrain, rep(labels[i], Ni_train_list[i]))
    Xtest <- cbind(Xtest, sde_sim_time(t0, T, x0, n, Ni_test_list[i], b_time, sigma_time, theta[i]))
    ytest <- c(ytest, rep(labels[i], Ni_test_list[i]))
  }

  if (plot) {
    time_grid <- seq(t0, T, length.out = n)

    # Prepare the training data for ggplot
    train_data <- data.frame(Time = rep(time_grid, ncol(Xtrain)),
                             Value = as.vector(Xtrain),
                             Class = factor(rep(ytrain, each = length(time_grid))),
                             Group = rep(1:ncol(Xtrain), each = length(time_grid))) # Group to differentiate lines

    # Plot training data using ggplot2
    p1 <- ggplot(train_data, aes(x = Time, y = Value, color = Class, group = interaction(Class, Group))) +
      geom_line() +
      scale_color_manual(values = rainbow(length(labels))) +
      labs(title = "Training Data Simulation for Each Class", x = "Time", y = "Values of the Mixed process") +
      theme_minimal() +
      theme(legend.position = "right",
            panel.background = element_rect(fill = "gray90", color = NA),  # Change panel background to gray
            plot.background = element_rect(fill = "gray95", color = NA)) +
      guides(color = guide_legend(title = "Training Class"))
    print(p1)
    # Prepare the testing data for ggplot
    test_data <- data.frame(Time = rep(time_grid, ncol(Xtest)),
                            Value = as.vector(Xtest),
                            Class = factor(rep(ytest, each = length(time_grid))),
                            Group = rep(1:ncol(Xtest), each = length(time_grid))) # Group to differentiate lines

    # Plot testing data using ggplot2
    p2 <- ggplot(test_data, aes(x = Time, y = Value, color = Class, group = interaction(Class, Group))) +
      geom_line() +
      scale_color_manual(values = rainbow(length(labels))) +
      labs(title = "Testing Data Simulation for Each Class", x = "Time", y = "Value") +
      theme_minimal() +
      theme(legend.position = "right",
            panel.background = element_rect(fill = "gray90", color = NA),  # Change panel background to gray
            plot.background = element_rect(fill = "gray95", color = NA)) +
      guides(color = guide_legend(title = "Testing Class"))
    print(p2)
  }

  return(list(Xtrain = Xtrain, ytrain = ytrain, Xtest = Xtest, ytest = ytest))
}


library(ggplot2)

# # Define the complex drift function
# b_complex <- function(x, theta) {
#   theta * (1/4 + (3/4) * (cos(30*x)^2))
# }
# sigma_complex <- function(x) {
#   0.1 + 0.9 / sqrt(1 + x^2)
# }
#
# # Create a sequence of x values from 0 to 1
# x_values <- seq(0, 1, length.out = 100)
#
# # Set theta and t for plotting
# theta <- 1  # Adjust this value as needed
# t <- 0.5    # Time value for the plots
#
# # Compute values for drift and diffusion functions
# b_values <- sapply(x_values, b_complex, theta = theta)
# sigma_values <- sapply(x_values, sigma_complex)
#
# # Create a data frame for plotting
# plot_data <- data.frame(
#   x = x_values,
#   Drift = b_values,
#   Diffusion = sigma_values
# )
#
# # Plot the Drift function
# p1 <- ggplot(plot_data, aes(x = x)) +
#   geom_line(aes(y = Drift), color = 'blue') +
#   labs(title = expression(paste("Drift Function ", b(x, theta, t))),
#        x = "x", y = "Drift Value") +
#   theme_minimal()
#
# # Plot the Diffusion function
# p2 <- ggplot(plot_data, aes(x = x)) +
#   geom_line(aes(y = Diffusion), color = 'red') +
#   labs(title = expression(paste("Diffusion Function ", sigma(x, t))),
#        x = "x", y = "Diffusion Value") +
#   theme_minimal()
#
# # Print both plots
# print(p1)
# #print(p2)
#


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
############### Orstein-Ulhenbeck ###################
################################################################################
################################################################################


sde_sim_orstein <- function(t0, T, X0, n, Ni, b, sigma, label, labels) {
  dt <- (T - t0) / n
  x <- matrix(0, nrow = n, ncol = Ni)
  x[1, ] <- X0
  for (j in 2:n) {
    x[j, ] <- x[j-1, ] + b(x[j-1, ], label, labels) * dt + sigma * rnorm(Ni) * sqrt(dt)
  }
  return(x)
}

# b_orstein <- function(x, class){
#   return((1-x)*(class==1)
#          + (-1-x)*(class==2)
#          + (-x)*(class==3))
# }


b_orstein <- function(x, label, labels) {
  return ((1 - x) * (label == labels[1]) +
            (-1 - x) * (label == labels[2]) +
            (-x) * (label == labels[3]))
}

################################################################################
############################ Simulation of Data ################################
################################################################################

simulation_sde_orstein <- function(t0, T, x0, n, Ni_train_list, Ni_test_list, labels, b, sigma, theta, plot = TRUE, save_plot = FALSE) {
  ################################################################################
  ################################################################################
  #### Given drift and diffusion coefficients, simulates from the corresponding
  #### SDE
  ####
  #### T, t0 : Starting and ending time points
  #### x0 : Space starting point, i.e deterministic value of X0
  #### n : Discretization stamp
  #### Ni_train_list,labels_train : Arrays such that Ni_train_list[i] = corresponding Ni of the class label[i]
  #### Ni_test_list,labels_test : Same
  #### b,sigma : Respectively drift and diffusion coefficients
  #### theta : Array such that theta[i] is the corresponding theta for class i
  ################################################################################
  ################################################################################

  TimeStep <- 1 / n
  NbClass <- length(labels)

  Xtrain <- sde_sim_orstein(t0, T, x0, n, Ni_train_list[1], b, sigma, labels[1], labels)
  ytrain <- rep(labels[1], Ni_train_list[1])
  Xtest <- sde_sim_orstein(t0, T, x0, n, Ni_test_list[1], b, sigma, labels[1], labels)
  ytest <- rep(labels[1], Ni_test_list[1])

  for (i in 2:NbClass) {
    Xtrain <- cbind(Xtrain, sde_sim_orstein(t0, T, x0, n, Ni_train_list[i],
                                            b, sigma, labels[i], labels))
    ytrain <- c(ytrain, rep(labels[i], Ni_train_list[i]))
    Xtest <- cbind(Xtest, sde_sim_orstein(t0, T, x0, n, Ni_test_list[i],
                                          b, sigma, labels[i], labels))
    ytest <- c(ytest, rep(labels[i], Ni_test_list[i]))
  }

  if (plot) {
    time_grid <- seq(t0, T, length.out = n)
    # Prepare the training data for ggplot

    # Prepare the training data for ggplot
    train_data <- data.frame(Time = rep(time_grid, ncol(Xtrain)),
                             Value = as.vector(Xtrain),
                             Class = factor(rep(ytrain, each = length(time_grid))),
                             Group = rep(1:ncol(Xtrain), each = length(time_grid))) # Group to differentiate lines

    # Plot training data using ggplot2
    p1 <- ggplot(train_data, aes(x = Time, y = Value, color = Class, group = interaction(Class, Group))) +
      geom_line() +
      scale_color_manual(values = rainbow(length(labels))) +
      labs(title = "Training Data Simulation for Each Class", x = "Time", y = "Values of the Mixed process") +
      theme_minimal() +
      theme(legend.position = "right",
            panel.background = element_rect(fill = "gray90", color = NA),  # Change panel background to gray
            plot.background = element_rect(fill = "gray95", color = NA)) +
      guides(color = guide_legend(title = "Training Class"))
    print(p1)
    if (save_plot){
      ggsave(filename = "./Plots/train_data.png", plot = p1, width = 8, height = 5, dpi = 300)
    }
    # Prepare the testing data for ggplot
    test_data <- data.frame(Time = rep(time_grid, ncol(Xtest)),
                            Value = as.vector(Xtest),
                            Class = factor(rep(ytest, each = length(time_grid))),
                            Group = rep(1:ncol(Xtest), each = length(time_grid))) # Group to differentiate lines

    # Plot testing data using ggplot2
    p2 <- ggplot(test_data, aes(x = Time, y = Value, color = Class, group = interaction(Class, Group))) +
      geom_line() +
      scale_color_manual(values = rainbow(length(labels))) +
      labs(title = "Testing Data Simulation for Each Class", x = "Time", y = "Value") +
      theme_minimal() +
      theme(legend.position = "right",
            panel.background = element_rect(fill = "gray90", color = NA),  # Change panel background to gray
            plot.background = element_rect(fill = "gray95", color = NA)) +
      guides(color = guide_legend(title = "Testing Class"))
    print(p2)
    if (save_plot){
      ggsave(filename = "./Plots/train_data.png", plot = p2, width = 8, height = 5, dpi = 300)
    }
  }

  return(list(Xtrain = Xtrain, ytrain = ytrain, Xtest = Xtest, ytest = ytest))
}






simulation_sde_mix <- function(t0, T, x0, n, Ni_train_list, Ni_test_list, labels, b, b_or, sigma, theta, plot = TRUE, save_plot = FALSE) {
  ################################################################################
  ################################################################################
  #### Given drift and diffusion coefficients, simulates from the corresponding
  #### SDE
  ####
  #### T, t0 : Starting and ending time points
  #### x0 : Space starting point, i.e deterministic value of X0
  #### n : Discretization stamp
  #### Ni_train_list,labels_train : Arrays such that Ni_train_list[i] = corresponding Ni of the class label[i]
  #### Ni_test_list,labels_test : Same
  #### b,sigma : Respectively drift and diffusion coefficients
  #### theta : Array such that theta[i] is the corresponding theta for class i
  ################################################################################
  ################################################################################

  TimeStep <- 1 / n
  NbClass <- length(labels)

  Xtrain <- sde_sim_orstein(t0, T, x0, n, Ni_train_list[1], b, sigma, theta[1])
  ytrain <- rep(labels[1], Ni_train_list[1])
  Xtest <- sde_sim_orstein(t0, T, x0, n, Ni_test_list[1], b, sigma, theta[1])
  ytest <- rep(labels[1], Ni_test_list[1])

  for (i in 2:NbClass) {
    Xtrain <- cbind(Xtrain, sde_sim_orstein(t0, T, x0, n, Ni_train_list[i],
                                            b_or, sigma, labels[i]))
    ytrain <- c(ytrain, rep(labels[i], Ni_train_list[i]))
    Xtest <- cbind(Xtest, sde_sim_orstein(t0, T, x0, n, Ni_test_list[i],
                                          b_or, sigma, labels[i]))
    ytest <- c(ytest, rep(labels[i], Ni_test_list[i]))
  }

  if (plot) {
    time_grid <- seq(t0, T, length.out = n)
    # Prepare the training data for ggplot

    # Prepare the training data for ggplot
    train_data <- data.frame(Time = rep(time_grid, ncol(Xtrain)),
                             Value = as.vector(Xtrain),
                             Class = factor(rep(ytrain, each = length(time_grid))),
                             Group = rep(1:ncol(Xtrain), each = length(time_grid))) # Group to differentiate lines

    # Plot training data using ggplot2
    p1 <- ggplot(train_data, aes(x = Time, y = Value, color = Class, group = interaction(Class, Group))) +
      geom_line() +
      scale_color_manual(values = rainbow(length(labels))) +
      labs(title = "Training Data Simulation for Each Class", x = "Time", y = "Values of the Mixed process") +
      theme_minimal() +
      theme(legend.position = "right",
            panel.background = element_rect(fill = "gray90", color = NA),  # Change panel background to gray
            plot.background = element_rect(fill = "gray95", color = NA)) +
      guides(color = guide_legend(title = "Training Class"))
    print(p1)
    if (save_plot){
      ggsave(filename = "./Plots/train_data.png", plot = p1, width = 8, height = 5, dpi = 300)
    }
    # Prepare the testing data for ggplot
    test_data <- data.frame(Time = rep(time_grid, ncol(Xtest)),
                            Value = as.vector(Xtest),
                            Class = factor(rep(ytest, each = length(time_grid))),
                            Group = rep(1:ncol(Xtest), each = length(time_grid))) # Group to differentiate lines

    # Plot testing data using ggplot2
    p2 <- ggplot(test_data, aes(x = Time, y = Value, color = Class, group = interaction(Class, Group))) +
      geom_line() +
      scale_color_manual(values = rainbow(length(labels))) +
      labs(title = "Testing Data Simulation for Each Class", x = "Time", y = "Value") +
      theme_minimal() +
      theme(legend.position = "right",
            panel.background = element_rect(fill = "gray90", color = NA),  # Change panel background to gray
            plot.background = element_rect(fill = "gray95", color = NA)) +
      guides(color = guide_legend(title = "Testing Class"))
    print(p2)
    if (save_plot){
      ggsave(filename = "./Plots/train_data.png", plot = p2, width = 8, height = 5, dpi = 300)
    }
  }

  return(list(Xtrain = Xtrain, ytrain = ytrain, Xtest = Xtest, ytest = ytest))
}

