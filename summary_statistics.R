# Generative Network Formation Algorithms
# Author: Srikar Katta
# Date: 12/29/2021
rm(list = ls())
source('generative_algorithms.R')
source('network_characteristic_calculations.R')

library(parallel)
library(data.table)
set.seed(999)

memory.size(max = TRUE)

# set parameters to read in files
n_nodes <- 10000
n_iter <- 30

start_prob_a <- 0
end_prob_a <- 1
by_prob_a <- 0.01
prob_a_seq <- seq(start_prob_a, end_prob_a, by_prob_a)
n_node_seq <- c(seq(100, 2000, 100), seq(3000, 10000, 1000))

read_dir <- './simulation_data/dd_networks/'
write_dir <- './simulation_data/dd_networks_summary/'

if(!(write_dir %in% list.dirs('./simulation_data/'))) {
  dir.create(paste0(write_dir))
}


# read in file
mclapply(1:n_iter, function(iteration) {
  mclapply(prob_a_seq, function(prob_a) {
    if(!(paste0(write_dir, prob_a) %in% list.dirs(write_dir, recursive = F))) {
      dir.create(paste0(write_dir,  prob_a))
    }
    # go through each iteration for each prob_a
    cat(iteration, prob_a, '\n')
    start_time <- Sys.time()
    # read in old adjacency matrix
    adj_mat <- read.table(paste0(read_dir, prob_a, '/nodes_', n_nodes, '_iter_', iteration, '.txt'), header = TRUE) %>%
      as.matrix()
    node_level_summs <- lapply(n_node_seq, function(net_size) {
      if(net_size == 3000) cat(iteration, prob_a, ' reached 3000 nodes \n')
      # make network from subset of adjacency matrix
      network_store <- graph_from_adjacency_matrix(adjmatrix = adj_mat[1:net_size, 1:net_size], mode = 'undirected')
      # calculate summary statistics
      node_level_summ_df <- node_level_characteristics(network_store = network_store) %>%
        mutate(n_node = net_size,
               node_order = 1:net_size,
               prob_a = prob_a,
               iteration = iteration)
      # remove network from memory
      rm('network_store')
      gc()
      
      # return summary stats
      return(node_level_summ_df)
    }) %>%
      rbindlist() # join summary statistics from various iterations together
    end_time <- Sys.time()
    cat(prob_a, iteration, 'Time taken: ', end_time - start_time, '\n')
    
    # write out adjacency matrix
    write.table(x = node_level_summs, paste0(write_dir, prob_a, '/nodes_', n_nodes, '_iteration_', iteration, '.txt'))
    rm('adj_mat')
    rm('node_level_summs')
    rm('start_time')
    rm('end_time')
    gc()
    
  }, mc.cores = getOption('mc.cores', 1L))
}, mc.cores = getOption('mc.cores', 30L))
