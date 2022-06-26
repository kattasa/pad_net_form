# Generative Network Formation Algorithms
# Author: Srikar Katta
# Date: 04/03/2021
library(tidyverse)
library(igraph)
# library(netcom)
## Preferential Attachment with Neighbors Algorithm ----
pan_algo <- function(n, prob_a, stop_n){
  # parameters:
  #   n: number of nodes in network
  #   prob_a: probability of adopting neighbors of anchor node
  #   stop_n: number of nodes to add in this sequence
  # returns:
  #   adjacency matrix and sequence of anchor nodes
  
  ## parameter checks:
  if(n < 3) stop('`n` must be greater than 2.') 
  if(!(prob_a >= 0 && prob_a <= 1)) stop('`prob_a` must be between 0 and 1') 
  if(missing(stop_n)) stop_n <- n # if no stop number is specified, make entire matrix
  else if(stop_n > n){
    warning('`stop_n` is greater than number of nodes in the matrix. Overriding stop_n.')
    stop_n <- n
  }
  ## run algorithm
  
  # initialize adjacency matrix and list of anchor nodes
  adjacency_matrix <- matrix(data = 0, nrow = n, ncol = n)
  anchor_seq <- c(NA, 1) # NA for first node to enter graph
  # start network with two connected nodes
  adjacency_matrix[1, 2] <- 1
  adjacency_matrix[2, 1] <- 1
  
  # add nodes until stop_n
  for(i in seq(3, stop_n)){
    ## preferential attachment step: choose anchor node based on node degree ----
    graph_degree <- sum(adjacency_matrix)
    deg_prop <- (adjacency_matrix %*% rep(1, n))/graph_degree
    anchor_node <- sample(n, size = 1, replace = FALSE, prob = deg_prop)
    adjacency_matrix[, i] <- adjacency_matrix[, anchor_node] # copy anchor node neighbors to new node
    
    ## neighbor selection step ----
    # decide which neighbors to adopt; prob of gaining neighbor is prob_a
    adjacency_matrix[, i] <- adjacency_matrix[, anchor_node] * rbinom(n, 1, prob_a) 
    # make matrix symmetric
    adjacency_matrix[i, ] <- adjacency_matrix[, i] # make matrix symmetric
    
    ## add edge to anchor nodes ----
    adjacency_matrix[anchor_node, i] <- 1
    adjacency_matrix[i, anchor_node] <- 1
    
    ## add anchor to sequence of anchor nodes
    anchor_seq <- c(anchor_seq, anchor_node)
  }
  # return adjacency matrix and sequence of anchor nodes
  list(adjacency_matrix = adjacency_matrix, anchor_seq = anchor_seq) %>% return()
}


## Continue Preferential Attachment with Neighbors Algorithm Given Sequence of Anchors and Adjacency Matrix
pan_algo_cont <- function(old_adj_mat, prob_a, n){
  # parameters:
  #   old_adj_mat: an adjacency matrix returned from a previously run PAwN algo
  #   prob_a: probability of adopting neighbors of anchor node
  #   n: number of nodes in new matrix
  # returns:
  #   adjacency matrix
  
  ## parameter checks:
  if(!is.matrix(old_adj_mat)) stop('`old_adj_mat` is not a matrix')
  if(!(prob_a >= 0 && prob_a <= 1)) stop('`prob_a` must be between 0 and 1')
  if(n < nrow(old_adj_mat)){
    warning('`n` is less than adjacency matrix. Overriding n')
    quit()
    n <- nrow(old_adj_mat)
  }
  ## find where to start network and find the number of nodes in the network
  start_node <- 1 + nrow(old_adj_mat) # 1 anchor for each node that joined the network
  
  adjacency_matrix <- matrix(data = 0, nrow = n, ncol = n)
  adjacency_matrix[1:nrow(old_adj_mat), 1:nrow(old_adj_mat)] <- old_adj_mat
  
  graph_degree <- sum(adjacency_matrix)
  degree_seq <- adjacency_matrix %*% rep(1, n)
  
  # sample neighbors from outset
  
  for(i in seq(start_node, n)){
    ## preferential attachment step: choose anchor node based on node degree ----
    # start_time <- Sys.time()
    deg_prop <- degree_seq / graph_degree
    anchor_node <- sample(n, size = 1, replace = FALSE, prob = deg_prop)
    
    ## neighbor selection step ----
    # decide which neighbors to adopt; prob of gaining neighbor is prob_a
    adjacency_matrix[, i] <- adjacency_matrix[, anchor_node] * rbinom(n, 1, prob_a) 
    # make matrix symmetric
    adjacency_matrix[i, 1:i] <- adjacency_matrix[1:i, i] # make matrix symmetric
    
    ## add edge between anchor and new nodes ----
    adjacency_matrix[anchor_node, i] <- 1
    adjacency_matrix[i, anchor_node] <- 1
    
    ## update degree sequences ----
    graph_degree <- graph_degree + 2 * sum(adjacency_matrix[, i]) # add 2 to graph degree for each of anchor node and new node
    degree_seq <- degree_seq + adjacency_matrix[, i] # add 1 to degree of each node that is now also connected to new
    degree_seq[i] <- degree_seq[i] + sum(adjacency_matrix[, i]) # add degree of new node to degree sequence
    # end_time <- Sys.time()
  }
  # return adjacency matrix and sequence of anchor nodes
  list(adjacency_matrix = adjacency_matrix) %>% return()
}

## Barabasi Albert Preferential Attachment Algorithm ----
barabasi_albert <- function(n, m) {
  # parameters:
  #   n: number of nodes in network
  #   m: number of nodes to add in each step
  # returns:
  #   adjacency matrix
  
  # initialize adjacency matrix and list of anchor nodes
  adjacency_matrix <- matrix(data = 0, nrow = n, ncol = n)
  # start network with two connected nodes
  # adjacency_matrix[1, 2] <- 1
  # adjacency_matrix[2, 1] <- 1
  adjacency_matrix[(1:(m + 1)), (1:(m + 1))] <- 1
  diag(adjacency_matrix) <- 0
  
  ## preferential attachment step: choose anchor node based on node degree ----
  graph_degree <- sum(adjacency_matrix)
  degree_seq <- (adjacency_matrix %*% rep(1, n))
  
  for(i in seq(m + 1, n)) {
    ## preferential attachment step: choose anchor node based on node degree -----
    deg_prop <- degree_seq / graph_degree
    anchor_nodes <- sample(n, size = m, replace = FALSE, prob = deg_prop)
    # adjacency_matrix[, i] <- adjacency_matrix[, anchor_node] # copy anchor node neighbors to new node
    adjacency_matrix[anchor_nodes, i] <- 1
    adjacency_matrix[i, anchor_nodes] <- 1
    ## update graph degree and degree sequence
    graph_degree <- graph_degree + 2 * m
    degree_seq <- degree_seq + adjacency_matrix[, i] # add 1 to degree of each node that is now also connected to new
    degree_seq[i] <- degree_seq[i] + sum(adjacency_matrix[, i]) # add degree of new node to degree sequence
    
  }
  # return adjacency matrix and sequence of anchor nodes
  return(adjacency_matrix)
}

## Barabasi Albert Preferential Attachment Algorithm Given Sequence of Anchors and Adjacency Matrix
barabasi_albert_cont <- function(adjacency_matrix, anchor_seq, stop_n){
  # parameters:
  #   adjacency_matrix: an adjacency matrix returned from a previously run Barabasi Albert algo
  #   anchor_sequence: a list of nodes denoting anchor nodes
  #   stop_n: number of nodes to add in this sequence
  # returns:
  #   adjacency matrix and sequence of anchor nodes
  
  ## parameter checks:
  if(!is.matrix(adjacency_matrix)) stop('`adjacency_matrix` is not a matrix')
  if(length(anchor_seq) > nrow(adjacency_matrix)) stop('`anchor_seq` invalid: more anchors than nodes')
  if(stop_n > nrow(adjacency_matrix)){
    warning('`stop_n` is greater than adjacency matrix. Overriding stop_n.')
    stop_n <- nrow(adjacency_matrix)
  }
  ## find where to start network and find the number of nodes in the network
  start_node <- 1 + length(anchor_seq) # 1 anchor for each node that joined the network
  n <- nrow(adjacency_matrix)
  
  for(i in seq(start_node, stop_n)){
    ## preferential attachment step: choose anchor node based on node degree ----
    graph_degree <- sum(adjacency_matrix)
    deg_prop <- (adjacency_matrix %*% rep(1, n))/graph_degree
    anchor_node <- sample(n, size = 1, replace = FALSE, prob = deg_prop)
    adjacency_matrix[, i] <- adjacency_matrix[, anchor_node] # copy anchor node neighbors to new node
    
    ## add edge between anchor and new nodes ----
    adjacency_matrix[anchor_node, i] <- 1
    adjacency_matrix[i, anchor_node] <- 1
    
    ## add anchor to sequence of anchor nodes
    anchor_seq <- c(anchor_seq, anchor_node)
  }
  # return adjacency matrix and sequence of anchor nodes
  list(adjacency_matrix = adjacency_matrix, anchor_seq = anchor_seq) %>% return()
}

### Small world methods ----
create_smallworld <- function(prob_rewire, net_size) {
  small_world_sample <- sample_smallworld(dim = 1, size = net_size, nei = 2, p = prob_rewire) %>%
    simplify() %>%
    return()
  # clustering_sw <- transitivity(small_world_sample, type = 'global')
  # diameter_sw <- diameter(small_world_sample)
  # degree_dist <- degree(small_world_sample, mode = 'total', loops = FALSE)
  # return(
  #   list(
  #     degree_distribution = degree_dist,
  #     clustering_coefficient = clustering_sw,
  #     diameter = diameter_sw,
  #     average_degree = mean(degree_dist)
  #   )
  # )
}

### Preferential Attachment methods ----
create_pa <- function(m_nodes, net_size) {
  # parameters:
  #   m_nodes: number of nodes new node connects to at each time step
  #   net_size: number of total nodes in a network
  # cat('a')
  pref_attach_sample <- barabasi_albert(n = net_size, m = m_nodes) %>%
    igraph::graph_from_adjacency_matrix(mode = 'undirected') %>%
    # igraph::make_full_graph(n = m_nodes, directed = FALSE, loops = FALSE) %>%
    # sample_pa(
    #   n = net_size,
    #   power = 1,
    #   m = m_nodes,
    #   directed = FALSE
    # ) %>%
    simplify() %>%
    return()
  # clustering_pref_attach <- transitivity(pref_attach_sample, type = 'global')
  # diameter_pref_attach <- diameter(pref_attach_sample)
  # degree_dist <- degree(pref_attach_sample, mode = 'total', loops = FALSE)
  # return(
  #   list(
  #     degree_distribution = degree_dist,
  #     clustering_coefficient = clustering_pref_attach,
  #     diameter = diameter_pref_attach,
  #     average_degree = mean(degree_dist)
  #   )
  # )
}

### Duplication Divergence methods ----
# make_DD from https://github.com/langendorfr/netcom/blob/master/R/make_DD.R 04/25
dd_algo <- function(n, prob_a, stop_n){
  # parameters:
  #   n: number of nodes in network
  #   prob_a: probability of dropping duplicated neighbors of anchor node
  #   stop_n: number of nodes to add in this sequence
  # returns:
  #   adjacency matrix and sequence of anchor nodes
  
  ## parameter checks:
  if(n < 3) stop('`n` must be greater than 2.') 
  if(!(prob_a >= 0 && prob_a <= 1)) stop('`prob_a` must be between 0 and 1') 
  if(missing(stop_n)) stop_n <- n # if no stop number is specified, make entire matrix
  else if(stop_n > n){
    warning('`stop_n` is greater than number of nodes in the matrix. Overriding stop_n.')
    stop_n <- n
  }
  
  # flip prob_a so that we look at prob of adding neighbors rather than dropping
  prob_a <- 1 - prob_a
  ## run algorithm
  
  # initialize adjacency matrix and list of anchor nodes
  adjacency_matrix <- matrix(data = 0, nrow = n, ncol = n)
  anchor_seq <- c(NA, 1) # NA for first node to enter graph
  # start network with two connected nodes
  adjacency_matrix[1, 2] <- 1
  adjacency_matrix[2, 1] <- 1
  
  # add nodes until stop_n
  for(i in seq(3, stop_n)){
    ## choose anchor node randomly from pre-existing nodes
    anchor_node <- sample(n, size = 1, replace = FALSE)
    adjacency_matrix[, i] <- adjacency_matrix[, anchor_node] # copy anchor node neighbors to new node
    
    ## neighbor selection step ----
    # decide which neighbors to adopt; prob of gaining neighbor is prob_a
    adjacency_matrix[, i] <- adjacency_matrix[, anchor_node] * rbinom(n, 1, prob_a) 
    # make matrix symmetric
    adjacency_matrix[i, ] <- adjacency_matrix[, i] # make matrix symmetric
    
    ## add edge to anchor nodes ----
    adjacency_matrix[anchor_node, i] <- 1
    adjacency_matrix[i, anchor_node] <- 1

  }
  # return adjacency matrix and sequence of anchor nodes
  return(adjacency_matrix)
}

create_dd <- function(div_prob, net_size) {
  # parameters:
  #   m_nodes: number of nodes new node connects to at each time step
  #   div_prob: probability of losing edge of duplicated network
  dd_algo(size = net_size, prob_a = div_prob) %>%
    igraph::graph_from_adjacency_matrix() %>%
    simplify() %>%
    return()
}


## k2 Algorithm ----
k2_algo <- function(n, m) {
  # parameters:
  #   n: number of nodes in network
  #   m: number of nodes to add in each step
  # returns:
  #   adjacency matrix
  
  # initialize adjacency matrix and list of anchor nodes
  adjacency_matrix <- matrix(data = 0, nrow = n, ncol = n)
  # start network with two connected nodes
  # adjacency_matrix[1, 2] <- 1
  # adjacency_matrix[2, 1] <- 1
  adjacency_matrix[(1:(m + 1)), (1:(m + 1))] <- 1
  diag(adjacency_matrix) <- 0
  
  for(i in seq(m + 1, n)){
    ## preferential attachment step: choose anchor node based on number of nodes 1 or 2 steps away ----
    # a_ij in A^2 is the number of paths of length two from node i to j; if greater than 1, path exists
    d2_matrix <- ifelse(((adjacency_matrix %*% adjacency_matrix) >= 1) | (adjacency_matrix), 1, 0)
    graph_d2 <- sum(d2_matrix)
    d2_prop <- (d2_matrix %*% rep(1, n))/graph_d2
    # print(deg_prop)
    anchor_nodes <- sample(n, size = m, replace = FALSE, prob = d2_prop)
    # adjacency_matrix[, i] <- adjacency_matrix[, anchor_node] # copy anchor node neighbors to new node
    adjacency_matrix[anchor_nodes, i] <- 1
    adjacency_matrix[i, anchor_nodes] <- 1
    ## add anchor to sequence of anchor nodes
    # anchor_seq <- c(anchor_seq, anchor_node)
  }
  # return adjacency matrix and sequence of anchor nodes
  return(adjacency_matrix)
}

### k2 methods ----
create_k2 <- function(m_nodes, net_size) {
  # parameters:
  #   m_nodes: number of nodes new node connects to at each time step
  #   net_size: number of total nodes in a network
  k2_sample <- k2_algo(n = net_size, m = m_nodes) %>%
    igraph::graph_from_adjacency_matrix(mode = 'undirected') %>%
    # igraph::make_full_graph(n = m_nodes, directed = FALSE, loops = FALSE) %>%
    # sample_pa(
    #   n = net_size,
    #   power = 1,
    #   m = m_nodes,
    #   directed = FALSE
    # ) %>%
    simplify() %>%
    return()
  # clustering_pref_attach <- transitivity(pref_attach_sample, type = 'global')
  # diameter_pref_attach <- diameter(pref_attach_sample)
  # degree_dist <- degree(pref_attach_sample, mode = 'total', loops = FALSE)
  # return(
  #   list(
  #     degree_distribution = degree_dist,
  #     clustering_coefficient = clustering_pref_attach,
  #     diameter = diameter_pref_attach,
  #     average_degree = mean(degree_dist)
  #   )
  # )
}

## Holme and Kim algorithm ------
holme_kim_algo <- function(n, m) {
  # parameters:
  #   n: number of total nodes in the network
  #   m: 1/2 the number of nodes each new node connects to at each time step
  # returns:
  #   adjacency matrix
  
  # initialize adjacency matrix
  adjacency_matrix <- matrix(data = 0, nrow = n, ncol = n)
  # start network with 2*m connected nodes
  adjacency_matrix[(1:(2*m)), (1:(2*m))] <- 1
  diag(adjacency_matrix) <- 0
  
  # initialize degree sequence
  graph_degree <- sum(adjacency_matrix)
  degree_seq <- rep(0, n)
  degree_seq[1:(2*m)] <- (2*m) - 1

  for(i in seq(2*m + 1, n)){
    ## preferential attachment step: choose m anchor nodes based on node degree ----
    deg_prop <- degree_seq/graph_degree
    
    anchor_nodes <- sample(n, size = m, replace = FALSE, prob = deg_prop)
    adjacency_matrix[anchor_nodes, i] <- 1
    adjacency_matrix[i, anchor_nodes] <- 1
    # replace deg_prop of anchor nodes with 0s so we don't consider them
    deg_prop[anchor_nodes] <- 0
    # find neighbors of the anchor nodes
    for(anchor in anchor_nodes) {
      
      # find anchor's neighbors that are not also neighbors of the new node
      anchor_neighbors <- which((adjacency_matrix[anchor, 1:n] - adjacency_matrix[i, 1:n]) == 1)
      # if there are elements in anchor_neighbors, randomly select one and add edge between it and new node
      if(length(anchor_neighbors) > 0) {
        random_neighbor <- sample(anchor_neighbors, 1)
      } else {
        # if there are no elements, randomly select node via preferential attachment
        random_neighbor <- sample(1:n, size = 1, replace = FALSE, prob = deg_prop)
      }
      # add edge between new and random neighbor
      adjacency_matrix[i, random_neighbor] <- 1
      adjacency_matrix[random_neighbor, i] <- 1
      # replace deg_prop for random_neighbor with 0
      deg_prop[random_neighbor] <- 0
      
      ## update degree sequences ----
      graph_degree <- graph_degree + 4 # add 4 to graph degree for each of anchor node, new node, and random
    }
    
    degree_seq <- degree_seq + adjacency_matrix[, i] # add 1 to degree of each node that is now also connected to new
    degree_seq[i] <- degree_seq[i] + sum(adjacency_matrix[, i]) # add degree of new node to degree sequence
  }
  
  # return adjacency matrix and sequence of anchor nodes
  return(adjacency_matrix)
}

### k2 methods ----
create_holme_kim <- function(m_nodes, net_size) {
  # parameters:
  #   m_nodes: number of nodes new node connects to at each time step
  #   net_size: number of total nodes in a network
  holme_kim_sample <- holme_kim_algo(n = net_size, m = m_nodes) %>%
    igraph::graph_from_adjacency_matrix(mode = 'undirected') %>%
    # igraph::make_full_graph(n = m_nodes, directed = FALSE, loops = FALSE) %>%
    # sample_pa(
    #   n = net_size,
    #   power = 1,
    #   m = m_nodes,
    #   directed = FALSE
    # ) %>%
    simplify() %>%
    return()
  # clustering_pref_attach <- transitivity(pref_attach_sample, type = 'global')
  # diameter_pref_attach <- diameter(pref_attach_sample)
  # degree_dist <- degree(pref_attach_sample, mode = 'total', loops = FALSE)
  # return(
  #   list(
  #     degree_distribution = degree_dist,
  #     clustering_coefficient = clustering_pref_attach,
  #     diameter = diameter_pref_attach,
  #     average_degree = mean(degree_dist)
  #   )
  # )
}



##
padm_algo <- function(prob_a, prob_m, n){
  # parameters:
  #   prob_a: probability of adopting neighbors of anchor node
  #   prob_m: probability of losing neighbor of anchor node
  #   n: number of nodes in new matrix
  # returns:
  #   adjacency matrix
  
  ## parameter checks:
  if(!(prob_a >= 0 && prob_a <= 1)) stop('`prob_a` must be between 0 and 1')
  if(!(prob_m >= 0 && prob_m <= 1)) stop('`prob_a` must be between 0 and 1')
  if(n < 3) stop('`n` must be greater than 2.') 
  
  # initialize adjacency matrix and list of anchor nodes
  adjacency_matrix <- matrix(data = 0, nrow = n, ncol = n)
  anchor_seq <- c(NA, 1) # NA for first node to enter graph
  # start network with two connected nodes
  adjacency_matrix[1, 2] <- 1
  adjacency_matrix[2, 1] <- 1
  
  graph_degree <- sum(adjacency_matrix)
  degree_seq <- adjacency_matrix %*% rep(1, n)
  
  # sample neighbors from outset
  
  for(i in seq(1, n)){
    ## preferential attachment step: choose anchor node based on node degree ----
    # start_time <- Sys.time()
    deg_prop <- degree_seq / graph_degree
    anchor_node <- sample(n, size = 1, replace = FALSE, prob = deg_prop)
    
    ## neighbor selection step ----
    # decide which neighbors to adopt; prob of gaining neighbor is prob_a
    adjacency_matrix[, i] <- adjacency_matrix[, anchor_node] * rbinom(n, 1, prob_a) 
    # make matrix symmetric
    adjacency_matrix[i, 1:i] <- adjacency_matrix[1:i, i] # make matrix symmetric
    
    ## neighbor deletion step ----
    # prob of keeping neighbor is 1 - prob of deleting neighbor
    anchor_neighbors <- adjacency_matrix[, anchor_node]
    adjacency_matrix[, anchor_node] <- adjacency_matrix[, anchor_node] * rbinom(n, 1, 1 - prob_m)
    # make matrix symmetric
    adjacency_matrix[anchor_node, ] <- adjacency_matrix[, anchor_node]
    ## add edge between anchor and new nodes ----
    adjacency_matrix[anchor_node, i] <- 1
    adjacency_matrix[i, anchor_node] <- 1
    
    ## update degree sequences ----
    graph_degree <- graph_degree + 2 * sum(adjacency_matrix[, i]) # add 2 to graph degree for each of anchor node and new node
    graph_degree <- graph_degree - 2 * degree_seq[anchor_node] + 2 * sum(adjacency_matrix[, anchor_node])
    degree_seq <- degree_seq + adjacency_matrix[, i] # add 1 to degree of each node that is now also connected to new
    degree_seq <- degree_seq - anchor_neighbors + adjacency_matrix[, anchor_node] # del 1 from degree of each node prev connected to anchor but not anymore
    degree_seq[i] <- degree_seq[i] + sum(adjacency_matrix[, i]) # add degree of new node to degree sequence
    degree_seq[anchor_node] <- sum(adjacency_matrix[, anchor_node]) # add degree of anchor node to degree seq
    
    rm(anchor_neighbors)
    gc()
  }
  # return adjacency matrix and sequence of anchor nodes
  adjacency_matrix %>% return()
}



## BPS algorithm ------
create_bps <- function(n_nodes, n_edges, power_law_exp) {
  # parameters:
  #   n_nodes: number of nodes in the network
  #   n_edges: number of edges in the network
  #   power_law_exp: exponent of power law to assign node importance from
  
  # initialize a matrix with N nodes
  adjacency_matrix <- matrix(0, nrow = n_nodes, ncol = n_nodes)
  # print(adjacency_matrix)
  # assign importance weights to each node of the network
  importance_weights <- poweRlaw::rpldis(n_nodes, 1, alpha = power_law_exp)
  weights_sum <- sum(importance_weights)
  importance_weights_prop <- importance_weights/weights_sum
  names(importance_weights_prop) <- seq(n_nodes)
  
  # keep adding edges until there are n_edges
  total_edges <- 0
  while(total_edges < n_edges) {
    # first choose nodes to add edge to
    nodes_time_t <- sample(n_nodes, 2, replace = FALSE, prob = importance_weights_prop)
    # if an edge exists between selected nodes, resample
    while((adjacency_matrix[nodes_time_t[1], nodes_time_t[2]]) == 1) {
      nodes_time_t <- sample(n_nodes, 2, replace = FALSE, prob = importance_weights_prop)
    }
    # add an edge between two nodes
    adjacency_matrix[nodes_time_t[1], nodes_time_t[2]] <- 1
    adjacency_matrix[nodes_time_t[2], nodes_time_t[1]] <- 1
    # increase total number of edges by 1
    total_edges <- total_edges + 1
  }
  # return the adjacency matrix
  return(adjacency_matrix)
}