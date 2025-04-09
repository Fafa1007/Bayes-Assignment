#### Required Libraries ####
library(ggplot2)
library(MASS)

#### Data Generating Functions ####
generate_lost <- function(grid_size, nsims){
  # Function to generate the prior distribution for the 
  # location of the lost fisherman.
  # Args: 
  #       grid_size: the dimensions of the square search grid.
  #       nsims: number of samples to base the prior distribution on.
  
  mu_vec  <- c(grid_size/2, grid_size/2)
  sig_mat <- matrix(c(2, 1, 5, 5), 2, 2)
  
  dat <- mvrnorm(nsims, mu_vec, sig_mat)
  dat <- round(abs(dat))
  
  prior <- matrix(rep(0, grid_size^2), grid_size, grid_size)
  for (i in 1:nrow(dat)){
    # Ensure the simulated location is within the grid
    if (dat[i,1] < grid_size & dat[i,2] < grid_size){
      prior[dat[i,1], dat[i,2]] <- prior[dat[i,1], dat[i,2]] + 1
    }
  }
  prior <- prior/sum(prior)
  return(prior)
}

generate_fisherman <- function(grid_size){
  # Function to generate the true location of the lost fisherman.
  # This does not affect the search decision as it is unknown
  # to the search crew.
  # Args: 
  #       grid_size: the dimensions of the square search grid.
  
  mu_vec  <- c(grid_size/2, grid_size/2)
  sig_mat <- matrix(c(2, 1, 5, 5), 2, 2)
  
  location  <- round(mvrnorm(1, mu_vec, sig_mat))
  # Ensure the location is within bounds
  location[location >= grid_size] <- grid_size - 1
  true_grid <- matrix(rep(0, grid_size^2), grid_size, grid_size)
  true_grid[location[1], location[2]] <- 1
  
  return(true_grid)
}

#### Simulation Setup ####
search_size <- 20
# Generate the prior distribution using 1000 simulations
prior_grid <- generate_lost(search_size, nsims = 1000)
# Flatten the prior grid to a vector (cell ordering: row-wise)
theta <- as.vector(prior_grid)

# Generate detection probabilities for each cell:
# Using a uniform random value between 0.6 and 0.9 for each cell.
unifs <- runif(search_size^2, min = 0.6, max = 0.9)
detect_pr <- matrix(unifs, ncol = search_size)
p <- as.vector(detect_pr)

# Generate the true location of the fisherman (unknown to the search crew)
true_grid <- generate_fisherman(search_size)
# Identify the true cell (row-wise indexing)
true_cell <- which(as.vector(true_grid) == 1)
cat("True cell where the fisherman is located:", true_cell, "\n")

#### Helper Function: Heatmap Plotting ####
plot_heatmap <- function(theta, grid_size, title = "Heatmap of Occurrence Probabilities"){
  prob_matrix <- matrix(theta, nrow = grid_size, ncol = grid_size, byrow = TRUE)
  ggplot(melt(prob_matrix), aes(Var2, Var1, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "yellow", high = "red") +
    labs(x = "Grid Column", y = "Grid Row", title = title) +
    scale_y_reverse() +
    theme_minimal()
}

#### Function to Update Probabilities ####
update_probabilities <- function(theta, p, searched_cell) {
  # For the searched cell: update using Equation 3
  theta_new <- theta
  theta_new[searched_cell] <- (1 - p[searched_cell]) * theta[searched_cell] /
    (1 - p[searched_cell] * theta[searched_cell])
  
  # For all other cells: update using Equation 4
  theta_new[-searched_cell] <- theta[-searched_cell] / (1 - p[searched_cell] * theta[searched_cell])
  
  # Normalize to ensure probabilities sum to 1 (numerical stability)
  theta_new <- theta_new / sum(theta_new)
  return(theta_new)
}

#### Simulation of the Search Process ####
n_hours <- 48  # maximum number of search steps (1 search per hour)
true_cell_prob <- numeric()  # Track the posterior probability for the true cell

# For saving snapshots of the heatmap
theta_first <- NULL
theta_last <- NULL

# Search process simulation loop
for (hour in 1:n_hours) {
  # Calculate the chance of successful detection in each cell:
  # product of occurrence probability and detection probability.
  detection_chance <- theta * p
  
  # Select the cell with the highest chance of detection
  search_cell <- which.max(detection_chance)
  
  # Simulate the search in the chosen cell with a Bernoulli trial:
  detected <- rbinom(1, 1, p[search_cell])
  cat("Hour", hour, ": Searching cell", search_cell, "- Detected:", detected, "\n")
  
  # Record the current probability for the true cell
  true_cell_prob <- c(true_cell_prob, theta[true_cell])
  
  # Save heatmap snapshot at the first search
  if (hour == 1) {
    theta_first <- theta
  }
  
  # If the searched cell is the true cell and the fisherman is detected, stop the search.
  if (search_cell == true_cell && detected == 1) {
    cat("Fisherman found in cell", search_cell, "at hour", hour, "\n")
    break
  }
  
  # Update probabilities based on the search outcome (no detection)
  theta <- update_probabilities(theta, p, search_cell)
}

# Save the final heatmap snapshot
theta_last <- theta

#### Visualization ####
# Plot heatmaps: at the first search step and at the final search step
# (if running interactively, use print() to display ggplot objects)
library(reshape2)  # for melting matrices for ggplot

p1 <- plot_heatmap(theta_first, search_size, title = "Heatmap at First Search Step")
p2 <- plot_heatmap(theta_last, search_size, title = "Heatmap at Final Search Step")

print(p1)
print(p2)

# Plot the evolution of the posterior probability for the true cell over time
df_true <- data.frame(Hour = 1:length(true_cell_prob), 
                      True_Cell_Probability = true_cell_prob)
ggplot(df_true, aes(x = Hour, y = True_Cell_Probability)) +
  geom_line(color = "blue") +
  geom_point(color = "blue") +
  labs(title = "Posterior Probability for the True Cell Over Time",
       x = "Search Step (Hour)", y = "Probability") +
  theme_minimal()

