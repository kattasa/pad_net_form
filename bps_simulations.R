# Srikar Katta
# Bidirectional Preferential Attachment Algorithm
# 4/27/2021
rm(list = ls())
source('generative_algorithms.R')
source('network_characteristic_calculations.R')

library(poweRlaw)
library(data.table)
library(parallel)
set.seed(999)

n_nodes_list <- c(seq(100, 2000, 100), seq(3000, 10000, 1000))
n_iter <- 30
start_pl_exp <- 2
end_pl_exp <- 8
by_pl_exp <- 0.5
time_step <- 100

file_pathway <- './simulation_data/bps_networks/'

bps_df <- lapply(seq(n_iter), function(iteration) {
  cat('new', n_iter, '\n')
  mclapply(seq(start_pl_exp, end_pl_exp, by = by_pl_exp), function(pl_exp) {
    if(!(paste0(file_pathway, pl_exp) %in% list.dirs(file_pathway, recursive = F))) {
      dir.create(paste0(file_pathway, pl_exp))
    }
    mclapply(n_nodes_list, function(n_nodes) {
      cat('iteration, pl_exp, n_nodes: ', iteration, pl_exp, n_nodes, '\n')
      mclapply(seq(n_nodes/2, 4 * n_nodes, by = n_nodes/2), function(n_edges) {
        # number of nodes and number of edges
        network_sim <- create_bps(n_nodes = n_nodes,
                   n_edges = n_edges,
                   power_law_exp = pl_exp)
        write.table(x = network_sim, file = paste0(file_pathway,
                                                   pl_exp, '/nodes_', n_nodes, '_edges_', n_edges, '_iter_', iteration, '.txt'))
        
        rm('network_sim')
        gc()
      }, mc.cores = getOption('mc.cores', 4L))
    }, mc.cores = getOption('mc.cores', 4L))
  }, mc.cores = getOption('mc.cores', 4L))
})
# 
# bps_df %>%
#   pivot_wider(names_from = node_seq, values_from = degree_distribution) %>% 
#   fwrite(., paste0('simulation_data/bps_iter_', n_nodes, '_', n_iter, '_start', start_pl_exp, '_end', end_pl_exp, '_by', by_pl_exp, '.csv'))
