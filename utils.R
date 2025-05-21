library(pheatmap)
library(RColorBrewer)

library(ggplot2)
library(ggbeeswarm)
library(tidyr)
library(ggpubr)


window_average <- function(x, window_size=10) {
  n_rows_per = length(x)/5
  
  return(c(
    window_average_helper(x[1:n_rows_per], window_size),
    window_average_helper(x[(n_rows_per+1):(2*n_rows_per)], window_size),
    window_average_helper(x[(2*n_rows_per+1):(3*n_rows_per)], window_size),
    window_average_helper(x[(3*n_rows_per+1):(4*n_rows_per)], window_size),
    window_average_helper(x[(4*n_rows_per+1):(5*n_rows_per)], window_size)
  ))
}

window_average_helper <- function(x, window_size) {
  seqs <- seq(1, length(x), by=window_size)
  avg_vals <- rep(0, length(seqs)-1)
  for (i in c(2:length(seqs))) {
    avg_vals[i-1] <- mean(x[seqs[i-1]:seqs[i]])
  }
  return(avg_vals)
}

rolling_average <- function(x, window_size=10) {
  n_rows_per = length(x)/5
  
  return(c(
    rolling_average_helper(x[1:n_rows_per], window_size),
    rolling_average_helper(x[(n_rows_per+1):(2*n_rows_per)], window_size),
    rolling_average_helper(x[(2*n_rows_per+1):(3*n_rows_per)], window_size),
    rolling_average_helper(x[(3*n_rows_per+1):(4*n_rows_per)], window_size),
    rolling_average_helper(x[(4*n_rows_per+1):(5*n_rows_per)], window_size)
  ))
}

rolling_average_helper <- function(x, window_size) {
  avg_vals <- rep(0, length(x) - window_size)
  for (i in c(1:(length(x) - window_size))) {
    avg_vals[i] <- mean(x[i:(i+window_size-1)])
  }
  return(avg_vals)
}

rolling_difference <- function(x, window_size = 2) {
  diff_vals <- rep(0, length(x) - window_size)
  for (i in c(1:(length(x) - window_size))) {
    diff_vals[i] <- abs(x[(i+window_size-1)] - x[i])
  }
  return(diff_vals)
}

create_truncated_indices <- function(remove_first=100) {
  return(c(remove_first:10000, 
           (10000+remove_first):20000,
           (20000+remove_first):30000,
           (30000+remove_first):40000,
           (40000+remove_first):50000))
}