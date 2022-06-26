# Wrangling functions for Generative Network Formation Algorithms
# Author: Srikar Katta
# Date: 04/03/2021
library(tidyverse)
library(igraph)
library(influenceR)
source('./generative_algorithms.R')

## create function for degree-degree-distance distribution
deg_deg_dist_distr <- function(adjacency_matrix) {
  # parameters:
  #   adjacency_matrix: NxN matrix where 1 means edge exists btw nodes i and j; 0 otherwise
  # returns:
  #   degree degree distribution ==> max(deg(i), deg(j))/min(deg(i), deg(j))
  # calculate degree distribution
  degree_dist <- rowSums(adjacency_matrix)
  
  # get only upper triangle
  adjacency_matrix[lower.tri(adjacency_matrix)] <- 0
  
  # find connected nodes (i.e, adj mat == 1)
  edges <- which(adjacency_matrix == 1, arr.ind = T)

  # get degree degree distance for connected nodes
  deg_deg_dist <- function(node1, node2, d_dist) {
    max(d_dist[node1], d_dist[node2])/min(d_dist[node1], d_dist[node2]) %>%
      return()
  } 
  
  mapply(FUN = function(x, y) deg_deg_dist(x, y, degree_dist),
         x = edges[, 'row'],
         y = edges[, 'col']) %>%
    return()
}


## calcualte network characteristics for a single network
classification_characteristics <- function(network_store, degree_dist) {
  # parameters:
  #   network_store: igraph graph
  #   degree_dist: a sequence of degree per node in the network
  # returns:
  #   a list of network characteristics
  if(missing(degree_dist))
    degree_dist <- igraph::degree(network_store, mode = 'total', loops = FALSE)
 
  # degree centrality
  # degree_dist <- degree(network_store, mode = 'total', loops = FALSE)
  degree_dist <- degree_dist[which(degree_dist > 0)]
  min_deg <- min(degree_dist)[1]
  max_deg <- max(degree_dist)[1]
  avg_deg <- mean(degree_dist)[1]
  std_deg <- sd(degree_dist)[1]

  # betweenness centrality
  betweenness_dist <- igraph::betweenness(network_store, directed = F)
  is.na(betweenness_dist) <- sapply(betweenness_dist, is.infinite)
  min_btw <- min(betweenness_dist, na.rm = TRUE)[1]
  max_btw <- max(betweenness_dist, na.rm = TRUE)[1]
  avg_btw <- mean(betweenness_dist, na.rm = TRUE)[1]
  std_btw <- sd(betweenness_dist, na.rm = TRUE)[1]
  
  # closeness centrality
  closeness_dist <- igraph::closeness(network_store, mode = 'all')
  is.na(closeness_dist) <- sapply(closeness_dist, is.infinite)
  min_cls <- min(closeness_dist, na.rm = TRUE)[1]
  max_cls <- max(closeness_dist, na.rm = TRUE)[1]
  avg_cls <- mean(closeness_dist, na.rm = TRUE)[1]
  std_cls <- sd(closeness_dist, na.rm = TRUE)[1]

  # eigenvector centrality
  eigenvector_dist <- igraph::eigen_centrality(network_store, directed = F)$vector
  is.na(eigenvector_dist) <- sapply(eigenvector_dist, is.infinite)
  min_eig <- min(eigenvector_dist, na.rm = TRUE)[1]
  max_eig <- max(eigenvector_dist, na.rm = TRUE)[1]
  avg_eig <- mean(eigenvector_dist, na.rm = TRUE)[1]
  std_eig <- sd(eigenvector_dist, na.rm = TRUE)[1]

  # clustering coefficient
  local_cc_dist <- igraph::transitivity(network_store, type = 'localundirected', isolates = 'zero')
  is.na(local_cc_dist) <- sapply(local_cc_dist, is.infinite)
  min_lcc <- min(local_cc_dist, na.rm = TRUE)[1]
  max_lcc <- max(local_cc_dist, na.rm = TRUE)[1]
  avg_lcc <- mean(local_cc_dist, na.rm = TRUE)[1]
  std_lcc <- sd(local_cc_dist, na.rm = TRUE)[1]
  
  # effective network size: size around a node
  effective_net_size_dist <- influenceR::ens(network_store)
  is.na(effective_net_size_dist) <- sapply(effective_net_size_dist, is.infinite)
  min_ens <- min(effective_net_size_dist, na.rm = TRUE)[1]
  max_ens <- max(effective_net_size_dist, na.rm = TRUE)[1]
  avg_ens <- mean(effective_net_size_dist, na.rm = TRUE)[1]
  std_ens <- sd(effective_net_size_dist, na.rm = TRUE)[1]
  
  # Burt's constraint
  constraint_dist <- igraph::constraint(network_store)
  is.na(constraint_dist) <- sapply(constraint_dist, is.infinite)
  min_cns <- min(constraint_dist, na.rm = TRUE)[1]
  max_cns <- max(constraint_dist, na.rm = TRUE)[1]
  avg_cns <- mean(constraint_dist, na.rm = TRUE)[1]
  std_cns <- sd(constraint_dist, na.rm = TRUE)[1]

  return_list <- list(
    min_deg,
    max_deg,
    avg_deg,
    std_deg,
    min_btw,
    max_btw,
    avg_btw,
    std_btw,
    min_cls,
    max_cls,
    avg_cls,
    std_cls,
    min_eig,
    max_eig,
    avg_eig,
    std_eig,
    min_lcc,
    max_lcc,
    avg_lcc,
    std_lcc,
    min_ens,
    max_ens,
    avg_ens,
    std_ens,
    min_cns,
    max_cns,
    avg_cns,
    std_cns
  )
  names(return_list) <- c(
    'min_deg',
    'max_deg',
    'avg_deg',
    'std_deg',
    'min_btw',
    'max_btw',
    'avg_btw',
    'std_btw',
    'min_cls',
    'max_cls',
    'avg_cls',
    'std_cls',
    'min_eig',
    'max_eig',
    'avg_eig',
    'std_eig',
    'min_lcc',
    'max_lcc',
    'avg_lcc',
    'std_lcc',
    'min_ens',
    'max_ens',
    'avg_ens',
    'std_ens',
    'min_cns',
    'max_cns',
    'avg_cns',
    'std_cns'
  )
  # return_list <- list(
  #   # degree_distribution = degree_dist,
  #   betweenness_distribution = betweenness_dist,
  #   closeness_distribution = closeness_dist,
  #   eigenvector_distribution = eigenvector_dist,
  #   local_cc_distribution = local_cc_dist,
  #   effective_net_size_distribution = effective_net_size_dist,
  #   constraint_distribution = constraint_dist
  # )
  return(return_list)
}

## calculate network characteristics
network_characteristics <- function(network_store, adjacency_matrix) {
  # parameters:
  #   network_store: igraph graph
  #   adjacency_matrix: an adjacency matrix
  # returns:
  #   degree_distribution, clustering coefficient, diameter
  ## parameter checks:
  if(missing(network_store)) {
    if(!is.matrix(adjacency_matrix)) stop('`adjacency_matrix` is not a matrix')
    ## make network from adjacency matrix
    network_store <- igraph::graph_from_adjacency_matrix(adjacency_matrix, mode = 'undirected')
  } else {
    if(!is_igraph(network_store)) stop('`network_store` is not an igraph graph.')
  }
  
  ## find degree distribution
  degree_dist <- degree(network_store, mode = 'total', loops = FALSE)
  ## find global clustering coefficient
  global_cc <- transitivity(network_store, type = 'globalundirected')
  ## find diameter
  network_diameter <- diameter(network_store, directed = FALSE)

  return(
    list(
      node_seq = paste0('node_', seq(length(degree_dist))),
      degree_distribution = degree_dist,
      # btw_distribution = igraph::betweenness(network_store, directed = FALSE),
      global_cc = global_cc,
      diameter = network_diameter,
      average_degree = mean(degree_dist)
    ) %>%
      c(., classification_characteristics(network_store, degree_dist)) %>%
      as.data.frame()
  )
}

## calculate network characteristics without the degree distribution
network_characteristics_no_deg <- function(network_store, adjacency_matrix) {
  # parameters:
  #   network_store: igraph graph
  #   adjacency_matrix: an adjacency matrix
  # returns:
  #   degree_distribution, clustering coefficient, diameter
  ## parameter checks:
  if(missing(network_store)) {
    if(!is.matrix(adjacency_matrix)) stop('`adjacency_matrix` is not a matrix')
    ## make network from adjacency matrix
    network_store <- igraph::graph_from_adjacency_matrix(adjacency_matrix, mode = 'undirected')
  } else {
    if(!is_igraph(network_store)) stop('`network_store` is not an igraph graph.')
  }
  
  ## find degree distribution
  degree_dist <- degree(network_store, mode = 'total', loops = FALSE)
  ## find global clustering coefficient
  global_cc <- transitivity(network_store, type = 'globalundirected')
  ## find diameter
  network_diameter <- diameter(network_store, directed = FALSE)
  
  return(
    list(
      # node_seq = paste0('node_', seq(length(degree_dist))),
      # degree_distribution = degree_dist,
      global_cc = global_cc,
      diameter = network_diameter,
      average_degree = mean(degree_dist)
    ) %>%
      c(., classification_characteristics(network_store, degree_dist)) %>%
      as.data.frame()
  )
}

## calculate characteristics over time
network_characteristics_time_series <- function(adjacency_matrix, time_step) {
  # parameters:
  #   adjacency_matrix: an adjacency matrix
  #   time_step: the number of steps after which network characteristics should be calculated
  # returns:
  #   degree_distribution, clustering coefficient, diameter over time
  ## parameter checks:
  if(!is.matrix(adjacency_matrix)) stop('`adjacency_matrix` is not a matrix')
  
  ## initialize network
  network_store <- igraph::graph_from_adjacency_matrix(adjacency_matrix[1:2, 1:2], mode = 'undirected') %>%
    set_vertex_attr(., name = 'label', index = V(.), value = c(1, 2))

  ## find initial network characteristics
  net_char_df <- network_characteristics(network_store) %>%
    mutate(time = 2)
  
  ## add new vertices and edges to network sequentially
  for(i in seq(3, nrow(adjacency_matrix))) {
    new_edges_to_add <- c()
    ## find edges to add for new node 
    for(neighbor in which(adjacency_matrix[i, 1:i] == 1)) new_edges_to_add <- c(new_edges_to_add, c(i, neighbor))
    
    network_store <- network_store %>%
      add_vertices(1, label = i) %>%
      add_edges(new_edges_to_add)
    ## find network characteristics for network as is after every time_step nodes added
    if((i %% time_step) == 0 || i <= min(time_step, 25)) {
      net_char_df <- network_characteristics(network_store) %>%
        mutate(time = i) %>%
        bind_rows(., net_char_df)
    }
  }
  # calculate network statistics if it hasn't been calculated when network is fully grown
  if((i %% time_step) != 0) {
    net_char_df <- network_characteristics(network_store) %>%
      mutate(time = i) %>%
      bind_rows(., net_char_df)
  }
  # net_char_df[which(is.na(net_char_df)), grep('node_', colnames(net_char_df))] <- 0
  return(net_char_df)
}

## calculate characteristics over time without the degree distribution
network_characteristics_time_series_no_deg <- function(adjacency_matrix, time_step) {
  # parameters:
  #   adjacency_matrix: an adjacency matrix
  #   time_step: the number of steps after which network characteristics should be calculated
  # returns:
  #   degree_distribution, clustering coefficient, diameter over time
  ## parameter checks:
  if(!is.matrix(adjacency_matrix)) stop('`adjacency_matrix` is not a matrix')
  
  ## initialize network
  network_store <- igraph::graph_from_adjacency_matrix(adjacency_matrix[1:2, 1:2], mode = 'undirected') %>%
    set_vertex_attr(., name = 'label', index = V(.), value = c(1, 2))
  
  ## find initial network characteristics
  net_char_df <- network_characteristics_no_deg(network_store) %>%
    mutate(time = 2)
  
  ## add new vertices and edges to network sequentially
  for(i in seq(3, nrow(adjacency_matrix))) {
    new_edges_to_add <- c()
    ## find edges to add for new node 
    for(neighbor in which(adjacency_matrix[i, 1:i] == 1)) new_edges_to_add <- c(new_edges_to_add, c(i, neighbor))
    
    network_store <- network_store %>%
      add_vertices(1, label = i) %>%
      add_edges(new_edges_to_add)
    ## find network characteristics for network as is after every time_step nodes added
    if((i %% time_step) == 0 || i <= min(time_step, 25)) {
      net_char_df <- network_characteristics_no_deg(network_store) %>%
        mutate(time = i) %>%
        bind_rows(., net_char_df)
    }
  }
  # calculate network statistics if it hasn't been calculated when network is fully grown
  if((i %% time_step) != 0) {
    net_char_df <- network_characteristics_no_deg(network_store) %>%
      mutate(time = i) %>%
      bind_rows(., net_char_df)
  }
  # net_char_df[which(is.na(net_char_df)), grep('node_', colnames(net_char_df))] <- 0
  return(net_char_df)
}


## network characteristics dataframe
network_char_to_df <- function(network_char_results, iteration = 1, prob_a) {
  # parameters:
  #   network_char_results: the data frame from network characteristics
  #   iteration: the iteration in which this dataframe is being constructed; default is 1
  #   prob_a: probability with which node connects to neighbors
  # returns:
  #   dataframe version of the network characteristics results
  network_char_results %>%
    # pivot_wider(names_from = node_seq, values_from = degree_distribution) %>%
    mutate(iteration = iteration,
           prob_a = prob_a) %>%
  return()
}

## calculate network characterstics per node, not as a summary
node_level_characteristics <- function(network_store, degree_dist) {
  # parameters:
  #   network_store: igraph graph
  #   degree_dist: a sequence of degree per node in the network
  # returns:
  #   a list of network characteristics
  
  # degree centrality
  if(missing(degree_dist))
    degree_dist <- igraph::degree(network_store, mode = 'total', loops = FALSE)
  
  # betweenness centrality
  betweenness_dist <- igraph::betweenness(network_store, directed = F)
  
  # closeness centrality
  closeness_dist <- igraph::closeness(network_store, mode = 'all')
  
  # eigenvector centrality
  # eigenvector_dist <- igraph::eigen_centrality(network_store, directed = F)$vector
  
  # clustering coefficient
  local_cc_dist <- igraph::transitivity(network_store, type = 'localundirected', isolates = 'zero')
  
  # effective network size: size around a node
  effective_net_size_dist <- influenceR::ens(network_store)
  
  # Burt's constraint
  constraint_dist <- igraph::constraint(network_store)
  
  ## find global clustering coefficient
  global_cc <- transitivity(network_store, type = 'globalundirected')
  ## find diameter
  network_diameter <- diameter(network_store, directed = FALSE)
  
  return_df <- tibble(
    node_seq = paste0('node_', 1:length(degree_dist)),
    degree_distribution = degree_dist,
    betweenness_distribution = betweenness_dist,
    closeness_distribution = closeness_dist,
    # eigenvector_distribution = eigenvector_dist,
    local_cc_distribution = local_cc_dist,
    effective_net_size_distribution = effective_net_size_dist,
    constraint_distribution = constraint_dist
  ) %>%
    mutate(global_cc = global_cc, diameter = network_diameter)
  return(return_df)
}

## calculate node level characteristics over time
node_level_characteristics_time_series <- function(adjacency_matrix, time_step) {
  # parameters:
  #   adjacency_matrix: an adjacency matrix
  #   time_step: the number of steps after which network characteristics should be calculated
  # returns:
  #   degree_distribution, clustering coefficient, diameter over time
  ## parameter checks:
  if(!is.matrix(adjacency_matrix)) stop('`adjacency_matrix` is not a matrix')
  
  ## initialize network
  network_store <- igraph::graph_from_adjacency_matrix(adjacency_matrix[1:2, 1:2], mode = 'undirected') %>%
    set_vertex_attr(., name = 'label', index = V(.), value = c(1, 2))
  
  ## find initial network characteristics
  net_char_df <- node_level_characteristics(network_store) %>%
    mutate(time = 2)
  
  ## add new vertices and edges to network sequentially
  for(i in seq(3, nrow(adjacency_matrix))) {
    new_edges_to_add <- c()
    ## find edges to add for new node 
    for(neighbor in which(adjacency_matrix[i, 1:i] == 1)) new_edges_to_add <- c(new_edges_to_add, c(i, neighbor))
    
    network_store <- network_store %>%
      add_vertices(1, label = i) %>%
      add_edges(new_edges_to_add)
    ## find network characteristics for network as is after every time_step nodes added
    if((i %% time_step) == 0 || i <= min(time_step, 25)) {
      net_char_df <- node_level_characteristics(network_store) %>%
        mutate(time = i) %>%
        bind_rows(., net_char_df)
    }
  }
  # calculate network statistics if it hasn't been calculated when network is fully grown
  if((i %% time_step) != 0) {
    net_char_df <- node_level_characteristics(network_store) %>%
      mutate(time = i) %>%
      bind_rows(., net_char_df)
  }
  return(net_char_df)
}

## node level characteristics dataframe
node_level_char_to_df <- function(network_char_results, iteration = 1, prob_a) {
  # parameters:
  #   network_char_results: the data frame from network characteristics
  #   iteration: the iteration in which this dataframe is being constructed; default is 1
  #   prob_a: probability with which node connects to neighbors
  # returns:
  #   dataframe version of the network characteristics results
  network_char_results %>%
    # pivot_wider(names_from = node_seq, values_from = degree_distribution) %>%
    mutate(iteration = iteration,
           prob_a = prob_a) %>%
    return()
}

pan_single_graph_steps <- function(n, prob_a) {
  # parameters:
  #   n: number of nodes in the network
  #   prob_a: probability of attaching nodes to network
  # returns:
  #   dataframe with summary results in dataframe
  pan_results <- pan_algo(n = n, prob_a = prob_a)
  pan_results_summary <- node_level_characteristics(pan_results$adjacency_matrix)
  node_level_char_to_df(pan_results_summary, prob_a = prob_a) %>% return()
  # pan_results_summary <- network_characteristics(pan_results$adjacency_matrix)
  # network_char_to_df(pan_results_summary, prob_a = prob_a) %>% return()
}

pan_multiple_graph_summary <- function(n_iter, n_nodes, prob_a_seq) {
  # parameters:
  #   n_iter: number of iterations to run for each graph
  #   n_nodes: number of nodes in each graph
  #   prob_a_seq: list of values that describes probability of adopting neighbors
  # results:
  #   dataframe with summary results in dataframe for all networks
  bind_rows(lapply(seq(n_iter), function(x)
    bind_rows(
      lapply(prob_a_seq, function(x)
        pan_single_graph_steps(n = n_nodes, prob_a = x))
    ))) %>%
    return()
}

pan_mult_graphs_deg_dist <- function(mult_graph_summary, n_iter, n_nodes) {
  # parameters:
  #   mult_graph_summary: dataframe returned from pan_multiple_graph_summary
  #   n_iter: number of iterations to run for each graph
  #   n_nodes: number of nodes in each graph
  # results:
  #   plot of degree distribution
  
  return_plot_df <- mult_graph_summary %>%
    select('prob_a', 'degree_distribution') %>%
    # pivot_longer(!prob_a, names_to = 'degree_id', values_to = 'degree') %>%
    # filter(degree_id %in% paste0('node_', seq(n_nodes))) %>%
    group_by(prob_a) %>%
    mutate(n_nodes_iter = n()) %>%
    group_by(prob_a, degree_distribution) %>%
    summarize(pr_degree = n() / n_nodes_iter)
  return_plot <- ggplot() +
    geom_density(data = return_plot_df,
      aes(
        x = degree_distribution,
        fill = factor(prob_a)
        # y = pr_degree
      ),
      position = 'identity',
      alpha = 0.6,
      # color = NA,
      binwidth = 1,
      trim = TRUE
      # se = FALSE
    ) +
    labs(
      x = 'Degree',
      y = 'Pr(Degree = k)',
      fill = 'Prob. of Adopting\nNeighbors',
      title = 'Degree Distributions',
      subtitle = 'Preferential Attachment with Neighbor Adoption',
      caption = paste0(n_iter, ' networks made for each graph with ', n_nodes, ' nodes each.')
    ) +
    # guides(fill = guide_legend(ncol = 2)) +
    theme_bw()
  return(return_plot)
}

pan_mult_graphs_deg_dist_loglog <- function(mult_graph_summary, n_iter, n_nodes) {
  # parameters:
  #   mult_graph_summary: dataframe returned from pan_multiple_graph_summary
  #   n_iter: number of iterations to run for each graph
  #   n_nodes: number of nodes in each graph
  # results:
  #   plot of degree distribution with log-log scales
  return_plot_df <- mult_graph_summary %>%
    # select(contains('node_'), 'prob_a') %>%
    # pivot_longer(!prob_a, names_to = 'degree_id', values_to = 'degree') %>%
    group_by(prob_a) %>%
    mutate(n = n()) %>%
    group_by(prob_a, degree_distribution, n) %>%
    summarize(pr_degree = n()/n)
  
  return_plot <- ggplot() +
    geom_point(data = return_plot_df,
      aes(
        x = degree_distribution,
        color = factor(prob_a),
        y = pr_degree
      ),
      # binwidth = 1,
      position = 'identity',
      alpha = 0.5,
      size = 0.5
      # stat = 'count'
    ) +
    labs(
      x = 'log(Degree)',
      y = 'log(Pr(Degree = k))',
      color = 'Prob. of Adopting\nNeighbors',
      title = 'Degree Distributions',
      subtitle = 'Preferential Attachment with Neighbor Adoption',
      caption = paste0(n_iter, ' networks made for each graph with ', n_nodes, ' nodes each.')
    ) +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw()
  return(return_plot)
}

real_world_deg_dist_loglog <- function(mult_graph_summary, real_world_df, n_iter, n_nodes) {
  # parameters:
  #   mult_graph_summary: dataframe returned from pan_multiple_graph_summary
  #   real_world_df: dataframe with columns degree and pr_degree: k and Pr(deg = k)
  #   n_iter: number of iterations to run for each graph
  #   n_nodes: number of nodes in each graph
  # results:
  #   plot of degree distribution with log-log scales
  return_plot_df <- mult_graph_summary %>%
    select(contains('node_'), 'prob_a') %>%
    pivot_longer(!prob_a, names_to = 'degree_id', values_to = 'degree') %>%
    mutate(n = n()) %>%
    group_by(prob_a, degree, n) %>%
    summarize(pr_degree = n()/n)
  
  return_plot <- ggplot() +
    geom_point(data = return_plot_df,
               aes(
                 x = degree,
                 color = factor(prob_a),
                 y = pr_degree
               ),
               # binwidth = 1,
               position = 'identity',
               alpha = 0.5,
               size = 0.5
               # stat = 'count'
    ) +
    geom_point(data = real_world_df,
               aes(x = degree, y = pr_degree), position = 'identity', alpha = 0.5, size = 0.5) +
    labs(
      x = 'log(Degree)',
      y = 'log(Pr(Degree = k))',
      color = 'Prob. of Adopting\nNeighbors',
      title = 'Degree Distributions',
      subtitle = 'Preferential Attachment with Neighbor Adoption',
      caption = paste0(n_iter, ' networks made for each graph with ', n_nodes, ' nodes each.')
    ) +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw()
  return(return_plot)
}

pan_mult_graphs_deg_cdf <- function(mult_graph_summary, n_iter, n_nodes, linetype = 'solid') {
  # parameters:
  #   mult_graph_summary: dataframe returned from pan_multiple_graph_summary
  #   n_iter: number of iterations to run for each graph
  #   n_nodes: number of nodes in each graph
  # results:
  #   plot of CDF of degree distribution
  return_data <- mult_graph_summary %>%
    select(contains('node_'), 'prob_a') %>%
    pivot_longer(!prob_a, names_to = 'degree_id', values_to = 'degree') %>%
    mutate(linetype = linetype)
  return_plot <- ggplot() +
    stat_ecdf(data = return_data,
      aes(
        x = degree,
        color = factor(prob_a),
        linetype = linetype,
        # y = ..density..
      ),
      # binwidth = 1,
      position = 'identity',
      alpha = 0.5
    ) +
    labs(
      x = 'Degree',
      y = 'Pr(Deg ≤ k)',
      color = 'Prob. of Adopting\nNeighbors',
      title = 'CDF of the Degree Distributions',
      subtitle = 'Preferential Attachment with Neighbor Adoption',
      caption = paste0(n_iter, ' networks made for each graph with ', n_nodes, ' nodes each.')
    ) +
    theme_bw()
  return(return_plot)
}


pan_mult_graphs_clustering <- function(mult_graph_summary, n_iter, n_nodes) {
  # parameters:
  #   mult_graph_summary: dataframe returned from pan_multiple_graph_summary
  #   n_iter: number of iterations to run for each graph
  #   n_nodes: number of nodes in each graph
  # results:
  #   plot of clustering coefficient across different prob_a values
  return_graph <- (
    mult_graph_summary %>%
      select(global_cc, prob_a) %>%
      group_by(prob_a) %>%
      summarize(
        mean = mean(global_cc),
        sd = sd(global_cc),
        n = n(),
        ci_upper = mean + 1.96 * sd / sqrt(n),
        ci_lower = mean - 1.96 * sd / sqrt(n)
      )
  ) %>%
    ggplot() +
    geom_point(aes(x = prob_a, y = mean)) +
    geom_line(aes(x = prob_a, y = mean), linetype = 'dashed') +
    # geom_errorbar(aes(x = prob_a, ymin = ci_lower, ymax = ci_upper), width = 0) +
    geom_ribbon(aes(x = prob_a, ymin = ci_lower, ymax = ci_upper), alpha = 0.5) +
    labs(x = 'Prob. of Adopting Neighbors',
         y = 'Clustering Coefficient',
         title = 'Clustering Coefficients',
         subtitle = 'Preferential Attachment with Neighbor Adoption',
         caption = paste0(n_iter, ' networks made for each graph with ', n_nodes, ' nodes each.')) +
    theme_bw()
  return(return_graph)
}

pan_mult_graphs_diameter <- function(mult_graph_summary, n_iter, n_nodes) {
  # parameters:
  #   mult_graph_summary: dataframe returned from pan_multiple_graph_summary
  #   n_iter: number of iterations to run for each graph
  #   n_nodes: number of nodes in each graph
  # results:
  #   plot of diameter across different prob_a values
  return_graph <- (
    mult_graph_summary %>%
      select(diameter, prob_a) %>%
      group_by(prob_a) %>%
      summarize(
        mean = mean(diameter),
        sd = sd(diameter),
        n = n(),
        ci_upper = mean + 1.96 * sd / sqrt(n),
        ci_lower = mean - 1.96 * sd / sqrt(n)
      )
  ) %>%
    ggplot() +
    geom_point(aes(x = prob_a, y = mean)) +
    geom_line(aes(x = prob_a, y = mean), linetype = 'dashed') +
    # geom_errorbar(aes(x = prob_a, ymin = ci_lower, ymax = ci_upper), width = 0) +
    geom_ribbon(aes(x = prob_a, ymin = ci_lower, ymax = ci_upper), alpha = 0.5) +
    labs(x = 'Probability of Adopting Neighbors',
         y = 'Diameter',
         title = 'Diameter',
         subtitle = 'Preferential Attachment with Neighbor Adoption',
         caption = paste0(n_iter, ' networks made for each graph with ', n_nodes, ' nodes each.')) +
    theme_bw()
  return(return_graph)
}

graph_char_time <- function(n, prob_a, stop_n, iteration, time_step) {
  # parameters:
  #   n: number of nodes in the network
  #   prob_a: probability of attaching nodes to network
  #   stop_n: number of nodes at which to stop
  #   iteration: number of times this network has been simulated
  #   time_step: how many time steps in between each calculation
  # returns:
  #   dataframe with summary results over different time points
  
  pan_results_adj <- pan_algo(n = n, prob_a = prob_a, stop_n = stop_n)$adjacency_matrix
  pan_results_time <-
    node_level_characteristics_time_series(pan_results_adj, time_step = time_step) %>%
    node_level_char_to_df(iteration = iteration, prob_a = prob_a) 
  pan_results_time <- pan_results_time %>%
    pivot_longer(., cols = contains('_distribution'), #c(degree_distribution, btw_distribution),
                 names_to = 'distribution', 
                 values_to = 'value') %>%
    pivot_wider(., names_from = node_seq, values_from = value)
  
  
  # pan_results_time <-
  #   network_characteristics_time_series(pan_results_adj, time_step = time_step) %>%
  #   network_char_to_df(iteration = iteration, prob_a = prob_a) 
  # pan_results_time <- pan_results_time %>%
  #   pivot_longer(., cols = c(degree_distribution, btw_distribution),
  #                names_to = 'distribution', 
  #                values_to = 'value') %>%
  #   pivot_wider(., names_from = node_seq, values_from = value)
  # pan_results_time$degree_distribution %>% 
  #   as.data.frame() %>%
  #   `colnames<-`(paste0('node_', seq(n))) %>%
  #   mutate(global_cc = pan_results_time$global_cc,
  #          diameter = pan_results_time$diameter,
  #          average_degree = pan_results_time$average_degree,
  #          time = seq(n),
  #          prob_a = prob_a) %>%
    return(pan_results_time)
}

dd_graph_char_time <- function(n, prob_a, stop_n, iteration, time_step) {
  # parameters:
  #   n: number of nodes in the network
  #   prob_a: probability of dropping duplicated neighbors of anchor
  #   stop_n: number of nodes at which to stop
  #   iteration: number of times this network has been simulated
  #   time_step: how many time steps in between each calculation
  # returns:
  #   dataframe with summary results over different time points
  dd_results_adj <- dd_algo(n = n, prob_a = prob_a, stop_n = stop_n)
  dd_results_time <-
    network_characteristics_time_series_no_deg(dd_results_adj, time_step = time_step) %>%
    network_char_to_df(iteration = iteration, prob_a = prob_a)
  # pan_results_time$degree_distribution %>% 
  #   as.data.frame() %>%
  #   `colnames<-`(paste0('node_', seq(n))) %>%
  #   mutate(global_cc = pan_results_time$global_cc,
  #          diameter = pan_results_time$diameter,
  #          average_degree = pan_results_time$average_degree,
  #          time = seq(n),
  #          prob_a = prob_a) %>%
  return(dd_results_time)
}

clustering_over_time <- function(graph_char_time_summary, n_iter, n_nodes, time_step) {
  # parameters:
  #   mult_graph_summary: dataframe returned from graph_char_time
  #   n_iter: number of iterations to run for each graph
  #   n_nodes: number of nodes in each graph
  #   time_step: number of time periods between each calculation
  # results:
  #   plot of clustering over time
  
  ## find average nad 95% CF for clustering coefficient for each time point and prob_a value
  pan_time_sum <- graph_char_time_summary %>%
    filter(n_node %in% seq(0, n_nodes, by = time_step)) %>%
    group_by(n_node, prob_a) %>%
    summarize(
      mean = mean(global_cc, na.rm = T),
      sd = sd(global_cc, na.rm = T),
      n = n(),
      ci_upper = mean + 1.96 * sd / sqrt(n),
      ci_lower = mean - 1.96 * sd / sqrt(n)
    )
  
  
  return_plot <- ggplot(pan_time_sum) +
    geom_line(aes( # plot time series split by prob_a
      x = n_node,
      y = mean,
      color = factor(prob_a)
    )) +
    geom_ribbon(aes( # plot 95% CI for each time series
      x = n_node,
      ymin = ci_lower,
      ymax = ci_upper,
      fill = factor(prob_a)
    ),
    alpha = 0.2) +
    labs(
      x = 'Time',
      y = 'Clustering Coefficient',
      fill = 'Prob. of\nAdopting\nNeighbors',
      title = 'Clustering Coefficient as the Network Grows',
      subtitle = 'Preferential Attachment with Neighbor Adoption',
      caption = paste0(
        n_iter,
        ' networks made for each graph with ',
        n_nodes,
        ' nodes each.\nRibbons represent 95% confidence interval.'
      )
    ) +
    guides(color = FALSE) +
    theme_bw()
  return(return_plot)
}

diameter_over_time <- function(graph_char_time_summary, n_iter, n_nodes, time_step) {
  # parameters:
  #   mult_graph_summary: dataframe returned from graph_char_time
  #   n_iter: number of iterations to run for each graph
  #   n_nodes: number of nodes in each graph
  #   time_step: number of time periods between each calculation
  # results:
  #   plot of diameter over time
  
  ## find average nad 95% CF for diameter for each time point and prob_a value
  pan_time_sum <- graph_char_time_summary %>%
    filter(n_node %in% c(seq(0, time_step), seq(time_step, n_nodes, by = time_step))) %>%
    group_by(n_node, prob_a) %>%
    summarize(
      mean = mean(diameter, na.rm = T),
      sd = sd(diameter, na.rm = T),
      n = n(),
      ci_upper = mean + 1.96 * sd / sqrt(n),
      ci_lower = mean - 1.96 * sd / sqrt(n)
    )
  
  
  return_plot <- ggplot(pan_time_sum) +
    geom_line(aes( # plot time series split by prob_a
      x = n_node,
      y = mean,
      color = factor(prob_a)
    )) +
    geom_ribbon(aes( # plot 95% CI for each time series
      x = n_node,
      ymin = ci_lower,
      ymax = ci_upper,
      fill = factor(prob_a)
    ),
    alpha = 0.2) +
    labs(
      x = 'Time',
      y = 'Diameter',
      fill = 'Prob. of\nAdopting\nNeighbors',
      title = 'Diameter as the Network Grows',
      subtitle = 'Preferential Attachment with Neighbor Adoption',
      caption = paste0(
        n_iter,
        ' networks made for each graph with ',
        n_nodes,
        ' nodes each.\nRibbons represent 95% confidence interval.'
      )
    ) +
    guides(color = FALSE) +
    theme_bw()
  return(return_plot)
}

### create small world dataframe from n_iter of simulation
smallworld_df <- function(prob_rewire, net_size, iteration = 1) {
  # parameters:
  #   prob_rewire: probability of rewiring edges
  #   net_size: number of nodes in the network
  #   iteration: how many of this class of graphs were already made?
  # returns:
  #   dataframe with the...
  #     degree of each node of network w/ net_size nodes sim from Small World
  #     clustering coefficient
  #     distance
  #     average degree
  #     iteration
  #     prob_a that describes number of out edges for simulated network
  return(
    create_smallworld(prob_rewire = prob_rewire, net_size = net_size) %>%
      network_characteristics(network_store = .) %>%
      mutate(time = net_size,
             prob_a = prob_rewire,
             iteration = iteration)
      # network_char_to_df(
      #   network_char_results = .,
      #   iteration = iteration,
      #   prob_a = prob_rewire
      # )
  )
}

### create barabasi albert dataframe from n_iter of simulation
pref_attach_df <- function(m_nodes, net_size, iteration = 1) {
  # parameters:
  #   m_nodes: number of edges to create in each time step
  #   net_size: number of nodes in the network
  #   iteration: how many of this class of graphs were already made?
  # returns:
  #   dataframe with the...
  #     degree of each node of network w/ net_size nodes sim from Barabasi Albert 
  #     clustering coefficient
  #     distance
  #     average degree
  #     iteration
  #     prob_a that describes number of out edges for simulated network
  return(
    create_pa(m_nodes = m_nodes, net_size = net_size) %>%
      network_characteristics(network_store = .) %>%
      mutate(time = net_size,
             prob_a = m_nodes,
             iteration = iteration)
      # network_char_to_df(
      #   network_char_results = .,
      #   iteration = iteration,
      #   prob_a = m_nodes
      # )
  )
}

### create k2 dataframe from n_iter of simulation
k2_df <- function(m_nodes, net_size, iteration = 1) {
  # parameters:
  #   m_nodes: number of edges to create in each time step
  #   net_size: number of nodes in the network
  #   iteration: how many of this class of graphs were already made?
  # returns:
  #   dataframe with the...
  #     degree of each node of network w/ net_size nodes sim from Barabasi Albert 
  #     clustering coefficient
  #     distance
  #     average degree
  #     iteration
  #     prob_a that describes number of out edges for simulated network
  return(
    create_k2(m_nodes = m_nodes, net_size = net_size) %>%
      network_characteristics_no_deg(network_store = .) %>%
      mutate(time = net_size,
             prob_a = m_nodes,
             iteration = iteration)
    # network_char_to_df(
    #   network_char_results = .,
    #   iteration = iteration,
    #   prob_a = m_nodes
    # )
  )
}

### create Duplication Divergence dataframe from n_iter of simulation
dup_div_df <- function(div_prob, net_size, iteration = 1) {
  # parameters:
  #   div_prob: probability of dropping edge from duplicated node
  #   net_size: number of nodes in the network
  #   iteration: how many of this class of graphs were already made?
  # returns:
  #   dataframe with the...
  #     degree of each node of network w/ net_size nodes sim from Barabasi Albert 
  #     clustering coefficient
  #     distance
  #     average degree
  #     iteration
  #     prob_a that describes number of out edges for simulated network
  return(
    create_dd(div_prob = div_prob, net_size = net_size) %>%
      network_characteristics(network_store = .) %>%
      mutate(time = net_size,
             prob_a = div_prob,
             iteration = iteration)
    # network_char_to_df(
    #   network_char_results = .,
    #   iteration = iteration,
    #   prob_a = m_nodes
    # )
  )
}
# 
# ### create forest fire dataframe from n_iter of simulation
# forest_fire_df <- function(prob_a, net_size, iteration = 1) {
#   # parameters:
#   #   prob_a: the forward burning probability
#   #   net_size: number of nodes in the network
#   #   iteration: how many of this class of graphs were already made?
#   # returns:
#   #   dataframe with the...
#   #     degree of each node of network w/ net_size nodes sim from Barabasi Albert 
#   #     clustering coefficient
#   #     distance
#   #     average degree
#   #     iteration
#   #     prob_a that describes number of out edges for simulated network
#   return(
#     create_ff(prob_a = prob_a, net_size = net_size) %>%
#       network_char_to_df(
#         network_char_results = .,
#         iteration = iteration,
#         prob_a = prob_a
#       )
#   )
# }

### create holme kim dataframe from n_iter of simulation
holme_kim_df <- function(m_nodes, net_size, iteration = 1) {
  # parameters:
  #   m_nodes: 1/2 number of edges to create in each time step
  #   net_size: number of nodes in the network
  #   iteration: how many of this class of graphs were already made?
  # returns:
  #   dataframe with the...
  #     degree of each node of network w/ net_size nodes sim from Barabasi Albert 
  #     clustering coefficient
  #     distance
  #     average degree
  #     iteration
  #     prob_a that describes number of out edges for simulated network
  return(
    create_holme_kim(m_nodes = m_nodes, net_size = net_size) %>%
      network_characteristics(network_store = .) %>%
      mutate(time = net_size,
             prob_a = m_nodes,
             iteration = iteration)
  )
}

### calculate Earth Mover's Distance between different distributions
emdw_distAB <- function(A, B) {
  # parameters:
  #   A: degree distribution of first network as a dataframe with 2 cols
  #     col1 of A: location that describes the degree
  #     col2 of A: weight that describes Pr(deg = k)
  #   B: degree distribution of second network as a dataframe with 2 cols
  #     col1 of B: location that describes the degree
  #     col2 of B: weight that describes Pr(deg = k)
  # returns:
  #   a numeric distance metric from the Earth Mover's Distance Formula
  full_AB <- full_join(A, B, by = 'location')
  full_AB[is.na(full_AB)] <- 0
  emdist::emdw(
    A = full_AB[1] %>% as.matrix(),
    wA = full_AB[2] %>% as.matrix(),
    B = full_AB[1] %>% as.matrix(),
    wB = full_AB[3] %>% as.matrix()
  ) %>% return()
}


emdw_dist_sim_real <- function(sim, real) {
  # parameters:
  #   sim: the degrees from all simulated networks of interest
  #   real: the degrees of the real world network of interest
  # returns:
  #   a numeric distance metric from the Earth Mover's Distance Formula
  sim <- sim %>%
    pivot_longer(cols = contains('node_'),
                 names_to = 'degree_id',
                 values_to = 'location') %>%
    arrange(location) %>%
    count(location) %>%
    mutate(weight = n / sum(n)) %>%
    select(location, weight)
  real <- real %>%
    count(degree) %>%
    mutate(n = n / sum(n)) %>%
    rename(location = degree, weight = n)
  emdw_distAB(sim, real) %>% return()
}

find_emd_dist <- function(pan_sim, sw_sim, pa_sim, real_deg_dist) {
  # parameters:
  #   pan_sim: simulated networks to test for Pref Attach w/ Neighborhood adoption algorithm
  #   sw_sim: simulated networks to test for Small World algorithm
  #   pa_sim: simulated networks to test for Pref Attach algorithm
  #   real_deg_dist: real world degree distribution
  #   real_world_name: name of the real world networks
  # returns:
  #   list of grid searched values for the earth mover's distance
  emd_pan_real_result <-
    sapply(unique(pan_sim$prob_a), function(x)
      emdw_dist_sim_real(pan_sim %>% filter(prob_a == x), real_deg_dist))
  emd_sw_real_result <-
    sapply(unique(sw_sim$prob_a), function(x)
      emdw_dist_sim_real(
        sw_sim %>% filter(prob_a == x),
        real_deg_dist
      ))
  emd_pa_real_result <-
    sapply(unique(pa_sim$prob_a), function(x)
      emdw_dist_sim_real(pa_sim %>% filter(prob_a == x),
                         real_deg_dist))
  list(pana = emd_pan_real_result,
       sw = emd_sw_real_result,
       pa = emd_pa_real_result) %>%
    return()
}


color_palette_table <- c(
  "0" = "#3182BD",
  "0.01" = "#E6550D",
  "0.02" = "#31A354",
  "0.03" = "#756BB1",
  "0.04" = "#636363",
  "0.05" = "#6BAED6",
  "0.06" = "#FD8D3C",
  "0.07" = "#74C476",
  "0.08" = "#9E9AC8",
  "0.09" = "#969696",
  "0.1" = "#9ECAE1",
  "0.2" = "#FDAE6B",
  "0.3" = "#A1D99B",
  "0.4" = "#BCBDDC",
  "0.5" = "#BDBDBD",
  "0.6" = "#C6DBEF",
  "0.7" = "#FDD0A2",
  "0.8" = "#C7E9C0",
  "0.9" = "#DADAEB",
  "1" = "#D9D9D9"
)