#### Required Functions ####

library(ggplot2)
library(MASS)

#### Data Generating Functions ####

generate_lost <- function(grid_size, nsims){
  
  # Function to generate the prior distribution for the 
  # location of the lost fisherman
  # Args: 
  #       grid_size: the dimensions of the square search grid
  #       nsims: number of samples to base the prior distribution on
  
  mu_vec  <- c(grid_size/2, grid_size/2)
  sig_mat <- matrix(c(2, 1, 5, 5), 2,2)
  
  dat <- mvrnorm(nsims, mu_vec, sig_mat)
  dat <- round(abs(dat))
  
  prior <- matrix(rep(0,grid_size^2), grid_size, grid_size)
  for (i in 1:NROW(dat)){
    
    if (dat[i,1] < grid_size & dat[i,2] < grid_size){
      prior[dat[i,1], dat[i,2]] <- prior[dat[i,1], dat[i,2]] + 1
    }
    
  }
  prior <- prior/sum(prior)
  return(prior)
}

generate_fisherman <- function(grid_size){
  
  # Function to generate the true location of the lost fisherman.
  # This should not effect the search decision in any way!! It is unkown
  # to the search crew.
  # Args: 
  #       grid_size: the dimensions of the square search grid

  
  mu_vec  <- c(grid_size/2, grid_size/2)
  sig_mat <- matrix(c(2, 1, 5, 5), 2,2)
  
  location  <- round(mvrnorm(1, mu_vec, sig_mat))
  true_grid <- matrix(rep(0, grid_size^2), grid_size, grid_size)
  true_grid[location[1], location[2]] <- 1
  
  return(true_grid)
}

#### Simulation ####

search_size <- 20
unifs <- runif(search_size^2, min = 0.6, max = 0.9)
detect_pr <- matrix(unifs, ncol = search_size)

#...








