# Simulations for Comparing Generative Network Algorithms
# Author: Srikar Katta
# Date: 04/04/2021
rm(list = ls())
source('generative_algorithms.R')
source('network_characteristic_calculations.R')

library(parallel)
library(data.table)
set.seed(999)

memory.size(max = TRUE)
## change parameters to run simulated networks returning df of summ stats
# n_nodes <- 10000
# n_iter <- 30
# start_prob_a <- 0
# end_prob_a <- 1
# by_prob_a <- 0.01
# prob_a_seq <- seq(start_prob_a, end_prob_a, by_prob_a)
# time_step <- 100
# start_time <- Sys.time()
# pan_graphs_over_time <- rbindlist(mclapply(seq(n_iter), function(iteration) { # simulate network n_iter times
#   rbindlist(mclapply(prob_a_seq, function(prob_a) { # iterate through all prob_a values
#     # returns a single data frame so need to row bind w/ rbindlist
#     cat(iteration, prob_a, '\n')
#     graph_char_time(
#       n = n_nodes,
#       prob_a = prob_a,
#       stop_n = n_nodes,
#       iteration = iteration,
#       time_step = time_step
#     )
#   }, mc.cores = getOption('mc.cores', 20L)))
# }, mc.cores = getOption('mc.cores', 3L)))
# stop_time <- Sys.time()
# time_dif <- stop_time - start_time
# cat('time taken', time_dif)
# # write out results to csv
# fwrite(pan_graphs_over_time, paste0('simulation_data/pan_graphs_iter_', n_nodes, '_', n_iter, '_start', start_prob_a, '_end', end_prob_a, '_by', by_prob_a, '.csv'))

# change parameters to run simulated networks to return networks
# n_nodes <- 5000
# n_iter <- 30
# start_prob_a <- 0
# end_prob_a <- 1
# by_prob_a <- 0.01
# prob_a_seq <- seq(start_prob_a, end_prob_a, by_prob_a)
# # time_step <- 100
# start_time <- Sys.time()
# pan_graphs_over_time <- mclapply(seq(n_iter), function(iteration) { # simulate network n_iter times
#   mclapply(prob_a_seq, function(prob_a) { # iterate through all prob_a values
#     # returns a single data frame so need to row bind w/ rbindlist
#     cat(iteration, prob_a, '\n')
#     if(!(paste0('./simulation_data/pad_networks/', prob_a) %in% list.dirs('./simulation_data/pad_networks', recursive = F))) {
#       dir.create(paste0('./simulation_data/pad_networks/', prob_a))
#     }
#     pan_algo(n = n_nodes, prob_a = prob_a, stop_n = n_nodes)$adjacency_matrix %>%
#       #   # as.vector() %>%
#       #   graph_from_adjacency_matrix(adjmatrix = ., weighted = FALSE) %>%
#       # rbinom(n = n_nodes * n_nodes, size = 2, prob = prob_a) %>%
#       # matrix(data = ., nrow = n_nodes) %>%
#       write.table(x = ., file = paste0('./simulation_data/pad_networks/',
#                                        prob_a, '/nodes_', n_nodes, '_iter_', iteration, '.txt'))
#   }, mc.cores = getOption('mc.cores', 10L))
# }, mc.cores = getOption('mc.cores', 4L))
# stop_time <- Sys.time()
# time_dif <- stop_time - start_time
# cat('time taken', time_dif)

## grow already created networks
old_nodes <- 0
n_nodes <- 1000000
n_iter <- 30
start_prob_a <- 0.2
end_prob_a <- 0.2
by_prob_a <- 0.0
prob_a_seq <- seq(start_prob_a, end_prob_a, by_prob_a)
start_time <- Sys.time()
mclapply(seq(n_iter), function(iteration) { # simulate network n_iter times
  
    mclapply(prob_a_seq, function(prob_a) { # iterate through all prob_a values
      cat(iteration, prob_a, '\n')
      
      # read in old adjacency matrix
      if(old_nodes > 0) {
        new_adj <- read.table(paste0('./simulation_data/pad_networks/',
                          prob_a, '/nodes_', old_nodes, '_iter_', iteration, '.txt'), header = TRUE) %>%
          as.matrix() %>%
          # add new nodes
          pan_algo_cont(n = n_nodes, prob_a = prob_a, old_adj_mat = .) %>%
          .$adjacency_matrix
      } else {
        new_adj <- pan_algo(n = n_nodes, prob_a = prob_a) %>%
          .$adjacency_matrix
      }
      # cat(nrow(new_adj))
      # write out adjacency matrix
      write.table(x = new_adj, paste0('./simulation_data/pad_networks/',
                         prob_a, '/nodes_', n_nodes, '_iter_', iteration, '.txt'))
      rm('new_adj')
      gc()
      
      }, mc.cores = getOption('mc.cores', 1L))
  gc()
}, mc.cores = getOption('mc.cores', 30L))
stop_time <- Sys.time()
time_dif <- stop_time - start_time
cat('time taken', time_dif)



