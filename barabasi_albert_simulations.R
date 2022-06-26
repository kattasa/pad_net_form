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
file_pathway <- './simulation_data/ba_networks/'

# time_step <- 100
start_time <- Sys.time()
pan_graphs_over_time <- mclapply(1:n_iter, function(iteration) { # simulate network n_iter times
  mclapply(1:20, function(m_nodes) { # iterate through all anchor number values
    cat(iteration, m_nodes, '\n')
    if(!(paste0(file_pathway, m_nodes) %in% list.dirs(file_pathway, recursive = F))) {
      dir.create(paste0(file_pathway, m_nodes))
    }
    
    network_sim <- barabasi_albert(m = m_nodes, n = n_nodes)
    write.table(x = network_sim, file = paste0(file_pathway,
                                               m_nodes, '/nodes_', n_nodes, '_iter_', iteration, '.txt'))
    rm('network_sim')
    gc()
  
    
  }, mc.cores = getOption('mc.cores', 20L))
}, mc.cores = getOption('mc.cores', 1L))
stop_time <- Sys.time()
time_dif <- stop_time - start_time
cat('time taken', time_dif, '\n')
