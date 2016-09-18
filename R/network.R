#' Create a data frame of edges from a GO term matrix
#' 
#' @param m GeneMatrix object representing GO terms
#' @param sets Gene sets to use in the network. If NULL, use
#' \code{\link{ThresholdSets}} to find them
#' @param ancestors Whether to include all ancestors of significant
#' in the matrix
#' @param ontology Which of BP, MF, and CC to use (default all three)
#' @param ... Extra arguments passed on to threshold_sets
#' 
#' @return An igraph
#' 
#' @export
GetEdgesTable <- function(m, sets = NULL, ancestors = TRUE,
                          ontology = c("BP", "MF", "CC"), ...) {
  # include ancestors
  if (is.null(sets)) {
    sets <- ThresholdSets(m, ...)
  }
  if (length(sets) == 0) {
    stop("No thresholded sets to include in network")
  }
  
  if (ancestors) {
    sets <- GetGORelatives(sets, m@colData$ID) 
  }
  
  go_terms <- m@colData %>%
    dplyr::filter(ID %in% sets) %>%
    dplyr::filter(Ontology %in% ontology)
  
  # get direct edges going downward to create the graph, and combine with column data
  edges <- get_ancestry_matrix(go_terms$ID, tbl = TRUE, type = "CHILDREN", upward = FALSE) %>%
    mutate(go_id1 = as.character(go_id1)) %>%
    mutate(go_id2 = as.character(go_id2)) %>%
    left_join(dplyr::select(go_terms, ID, Ontology), by = c(go_id1 = "ID")) %>%
    tbl_df()
  
  edges
}


#' Create an igraph from a GeneMatrix representing GO terms
#' 
#' @param m A GeneMatrix object
#' @param edges A table of edges to plot. If not given, compute with
#' \code{\link{GetEdgesTable}}
#' @param ... Extra arguments passed on to GetEdgesTable
#' 
#' @return An igraph object
#' 
#' @export
GenerateNetwork <- function(m, edges = NULL, ...) {
  if (is.null(edges)) {
    edges <- GetEdgesTable(m, ...)
  }
  
  g <- igraph::graph.data.frame(edges)
  
  go_terms <- m@colData %>%
    filter(ID %in% c(edges$go_id1, edges$go_id2))
  
  go_node_data <- data.frame(ID = names(igraph::V(g)), stringsAsFactors = FALSE) %>%
    inner_join(go_terms, by = "ID")
  
  # transfer column information from nodes and edges
  for (col in colnames(go_node_data)) {
    # a hack, but necessary since V(g)[[col]] doesn't work
    eval(substitute(igraph::V(g)$replace <- go_node_data[[col]], list(replace = col)))
  }

  for (col in colnames(edges)) {
    eval(substitute(igraph::E(g)$replace <- edges[[col]], list(replace = col)))
  }
  
  g
}


#' Plot a network from an igraph representation
#' 
#' Return a ggplot2 object representing a network of sets.
#' This is meant to be a sensible default for plotting
#' GO networks, it will not cover all cases. We recommend learning
#' to use ggraph or another igraph plotting functionality.
#' 
#' @param g An igraph object, generally computed by \link{\code{GenerateNetwork}}
#' @param algorithm Algorithm to use for layout
#' @param arrow Arrow 
#' 
#' @import ggplot2
#' @import ggraph
#' 
#' @export
PlotNetwork <- function(g, algorithm = 'kk', arrow = grid::arrow(length = grid::unit(.1, "inches"))) {
  ggraph(g, 'igraph', algorithm = algorithm) +
    geom_edge_link(arrow = arrow) +
    geom_node_point(aes(size = Size)) +
    geom_node_point(aes(color = MeanDifference, size = Size)) +
    geom_node_text(aes(label = Term), check_overlap = TRUE, size = 3) +
    ggforce::theme_no_axes() +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar")
}