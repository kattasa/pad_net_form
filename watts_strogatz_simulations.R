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

# change parameters to run simulated networks to return networks
n_nodes <- 10000
n_iter <- 30
rewire_prob_seq <- seq(0, 1, 0.01)

file_pathway <- './simulation_data/ws_networks/'
# time_step <- 100
start_time <- Sys.time()
ws_graphs_over_time <- mclapply(1:n_iter, function(iteration) { # simulate network n_iter times
  mclapply(rewire_prob_seq, function(prob_a) { # iterate through all anchor number values
    cat(iteration, prob_a, '\n')
    if(!(paste0(file_pathway, prob_a) %in% list.dirs(file_pathway, recursive = F))) {
      dir.create(paste0(file_pathway, prob_a))
    }
    network_sim <- create_smallworld(prob_rewire = prob_a, net_size = n_nodes) %>%
      igraph::get.adjacency() %>%
      as.matrix()
    write.table(x = network_sim, file = paste0(file_pathway,
                                               prob_a, '/nodes_', n_nodes, '_iter_', iteration, '.txt'))
    rm('network_sim')
    gc()
  
    
  }, mc.cores = getOption('mc.cores', 1L))
}, mc.cores = getOption('mc.cores', 30L))
stop_time <- Sys.time()
time_dif <- stop_time - start_time
cat('time taken', time_dif, '\n')
