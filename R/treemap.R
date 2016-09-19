#' Generate Treemap
#' 
#' Visualize genes nested within specific and general GO sets using treemaps
#' 
#' Each gene should maximize the absolute sum of relevant node weights but every GO category will
#' not be included since a gene can only be assigned to a single ancestry. A small GO category may
#' fall into multiple ancestral paths some which are parents of other meaningful categories and others
#' that are dead-ends.  Should favor attachment to pathways with the most informative children.
#' To allow this to occur, three edge weights are used
#' 1) the path with the highest overall absolute score * number of genes affected is chosen
#' 2) ancestors of LASSO nodes all get a smaller score: [min(abs(beta))]*0.01*N_lasso_children
#' - applied to edge where ancestors are children
#' 3) all remaining edges that are not predictive recieve a tiny weight so that they are not filtered: [min(abs(beta))]*1e-4
#' 
#' @inheritParams GetEdgesTable
#' 
#' @export
#' @return Generates a treemap hierarchy
GenerateTreemap <- function(m, edges = NULL, ...) {
  if (is.null(edges)) {
    signif_sets <- ThresholdSets(m, ...)
    edges <- GetEdgesTable(m, sets = m@colData$ID, ...)
    # orphaned parents connect to root
    edge_roots <- data_frame(go_id1 = "root", go_id2 = setdiff(edges$go_id1, edges$go_id2)) %>%
      dplyr::left_join(edges %>%
                         dplyr::select(go_id2 = go_id1, Ontology) %>%
                         dplyr::distinct(), by = "go_id2")
    edges <- rbind(edges, edge_roots)
    
    ancestors <- get_ancestry_matrix(m@colData$ID, tbl = TRUE, type = "ANCESTOR", upward = FALSE) %>%
      dplyr::mutate(go_id1 = as.character(go_id1)) %>%
      dplyr::mutate(go_id2 = as.character(go_id2)) %>%
      dplyr::filter(go_id1 %in% m@colData$ID, go_id2 %in% m@colData$ID)
    ancestors <- bind_rows(ancestors, data_frame(go_id1 = "root", go_id2 = m@colData$ID))
  }
  
  if(length(signif_edges) == 0){
    stop("No thresholded sets to include in treemap")
  }
  
  signif_sets_weight <- m@colData %>%
    dplyr::filter(ID %in% signif_sets) %>%
    dplyr::select_("ID", "Size", Weight = m@plottingMetric) %>%
    dplyr::mutate(Weight = abs(Weight),
                  nWeight = Size * Weight) %>%
    dplyr::select(-Size)
  
  nodes <- data_frame(ID = unique(c(edges$go_id1, edges$go_id2))) %>%
    dplyr::left_join(m@colData %>% dplyr::select(ID, Size), by = "ID") %>%
    dplyr::left_join(signif_sets_weight, by = "ID") %>%
    dplyr::mutate(Weight = ifelse(is.na(Weight), 0, Weight),
                  nWeight = ifelse(is.na(nWeight), 0, nWeight))
  
  # ancestor weights encourage the same ancestor nodes to be used when significant descendents exist
  ancestor_weights <- ancestors %>%
    dplyr::left_join(nodes %>% dplyr::select(go_id2 = ID, weight = nWeight), by = "go_id2") %>%
    dplyr::mutate(weight = ifelse(is.na(weight), 0, weight)) %>%
    dplyr::group_by(go_id1) %>%
    dplyr::summarize(weight = sum(weight)*0.01) %>%
    dplyr::arrange(desc(weight))
  
  nodes <- nodes %>%
    dplyr::left_join(ancestor_weights, by = c("ID" = "go_id1")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(nWeight = sum(nWeight, weight, na.rm = T)) %>%
    dplyr::select(-weight) %>%
    dplyr::ungroup()
  
  # find a unique GO hierarchy that 
  minimal_edge_set <- find_minimal_rooted_tree(nodes, edges)
  
  # find the optimal path to assign each gene
  gene_paths <- assign_genes_to_paths(m, minimal_edge_set)
  
  # prune gene sets that are not used
  all_used_go_terms <- unique(gene_paths$full_path$ID)
  utilized_go_network <- igraph::delete_vertices(minimal_edge_set, setdiff(igraph::V(minimal_edge_set)$name, all_used_go_terms))
  
    
  #' Treemap GO Plot
  #' 
  #' Generate a treemap plot using the minimal spanning GO tree with genes assigned to individual GO terms.
  #' 
  #' Default color based on GO significance and effect sizes of individual genes
  treemap_GO_plot <- function(m, gene_paths, utilized_go_network, treemap_thresh = 4, ...){
    
    # gene effects
    gene_treemap_data <- gene_paths$assigned_genes %>%
      dplyr::left_join(m@geneData %>% dplyr::select(gene = ID, gene_effect = y), by = "gene") %>%
      dplyr::mutate(gene_effect = pmax(pmin(gene_effect, treemap_thresh), -1*treemap_thresh))
    
    gene_treemap_data <- gene_treemap_data %>%
      dplyr::select(ID = gene, gene_effect) %>%
      dplyr::mutate(go_name = NA_character_,
                    go_effect = NA_real_,
                    data_type = "gene")
    
    # go effects
    go_treemap_data <- m@colData %>%
      dplyr::filter(ID %in% intersect(signif_sets, igraph::V(utilized_go_network)$name)) %>%
      dplyr::select_("ID", go_name = "Term", go_effect = m@plottingMetric)
    
    go_treemap_data <- tibble::data_frame(ID = igraph::V(utilized_go_network)$name) %>%
      dplyr::left_join(go_treemap_data, by = "ID") %>%
      dplyr::mutate(gene_effect = NA_real_,
                    data_type = "go")
    
    all_node_data <- dplyr::bind_rows(gene_treemap_data, go_treemap_data)
    
    # update network
    updated_edges <- dplyr::bind_rows(
      igraph::get.edgelist(utilized_go_network) %>% as.data.frame(stringsAsFactors = F),
      unname(gene_paths$assigned_genes) %>% as.matrix() %>% as.data.frame(stringsAsFactors = F))
    
    # use utilized_go_network with genes assigned to categories
    full_treemap_data <- igraph::graph_from_data_frame(updated_edges, directed = TRUE, vertices = all_node_data)
    
    # update network 
    
    #library(ggraph)
    #library(ggplot2)
    
    #ggraph(full_treemap_data, 'igraph', algorithm = 'tree', circular = TRUE) + 
    # geom_edge_diagonal(aes(alpha = ..index..)) +
    #  coord_fixed() + 
    #  scale_edge_alpha('Direction', guide = 'edge_direction') +
    #  geom_node_point(aes(color = gene_effect, filter = igraph::degree(full_treemap_data, mode = 'out') == 0), size = 1) +
    #  ggforce::theme_no_axes() +
    #  scale_color_gradient2("MeanDifference", low = "blue", high = "red")
    
    # tree
    
    full_treemap_data <- ggraph::treeApply(full_treemap_data, function(node, parent, depth, tree) {
      tree <- igraph::set_vertex_attr(tree, 'depth', node, depth)
      if (depth == 1) {
        tree <- igraph::set_vertex_attr(tree, 'class', node, igraph::V(tree)$shortName[node])
      } else if (depth > 1) {
        tree <- igraph::set_vertex_attr(tree, 'class', node, igraph::V(tree)$class[parent])
      }
      tree
    })
  igraph::V(full_treemap_data)$leaf <- igraph::degree(full_treemap_data, mode = 'out') == 0

ggraph(full_treemap_data, 'treemap') + 
  geom_treemap(aes(fill = gene_effect, filter = leaf), colour = NA) + 
  geom_treemap(aes(size = depth * ifelse(data_type == "gene", 1, 0), alpha = depth,
                   colour = ifelse(data_type == "gene", "white", ifelse(go_effect != 0, "yellow", "black"))), fill = NA) + 
  geom_node_text(aes(label = go_name), size = 3, check_overlap = T, repel = T) +
  scale_fill_gradient2("MeanDifference", low = "green3", high = "firebrick1") +
  scale_color_identity() +
  scale_size(range = c(1, 0), guide = 'none') +
  scale_alpha(range = c(1, 0.2), guide = 'none') +
  ggforce::theme_no_axes()
    
  }
}


get_minimal_gene_path <- function(a_gene, m, nodes, minimal_edge_set){
  
  gene_GO <- m@matrix[a_gene, m@matrix[a_gene,] == 1, drop = F] %>% colnames()
  # add ancestors that may have been missed for the gene
  gene_GO <- union(gene_GO, ancestors$go_id1[ancestors$go_id2 %in% gene_GO])
  gene_GO <- gene_GO[gene_GO %in% nodes$ID]
  
  gene_subgraph <- igraph::induced_subgraph(minimal_edge_set, gene_GO)
  
  gene_paths <- igraph::get.shortest.paths(gene_subgraph, from = "root", to = igraph::V(gene_subgraph))$vpath
  gene_paths <- lapply(seq_along(gene_paths), function(i){
    tibble::data_frame(ID = gene_paths[[i]]$name) %>%
      dplyr::mutate(gene = a_gene, terminus = igraph::V(gene_subgraph)$name[i], step = 1:n())
  }) %>%
    dplyr::bind_rows()
  
  gene_paths
}

#' Assign Genes to Paths
#' 
#' Assign genes to a GO hierarchy which has the greatest signal
#' 
#' @inheritParams GenerateTreemap
#' @param minimal_edge_set a network which is a weighted directed minimal spanning tree of gene ontologies
#' 
#' @return a data_frame which contains which terminal gene set each gene is assigned to
assign_genes_to_paths <- function(m, minimal_edge_set, ...){
  
  # find all paths for each gene
  all_gene_paths <- mclapply(m@geneData$ID, function(a_gene){
    print(a_gene)
    get_minimal_gene_path(a_gene, m, nodes, minimal_edge_set)
  }, mc.cores = parallel::detectCores()) %>%
    dplyr::bind_rows()
  
  # all shortest paths found, now find the path that maximizes the path weight
  
  assigned_paths <- all_gene_paths %>%
    dplyr::left_join(nodes, by = "ID") %>%
    dplyr::group_by(gene, terminus) %>%
    dplyr::summarize(Weight = sum(Weight),
                     nWeight = sum(nWeight),
                     nSteps = n()) %>%
    dplyr::group_by(gene) %>%
    dplyr::arrange(desc(Weight), desc(nWeight), nSteps) %>%
    dplyr::slice(1)
  
  gene_paths <- list()
  gene_paths$full_path <- all_gene_paths %>%
    dplyr::semi_join(assigned_paths, by = c("terminus", "gene"))
  gene_paths$assigned_genes <- assigned_paths %>%
    dplyr::select(terminus, gene)
  gene_paths
}

#' Find Minimal Rooted Tree
#' 
#' Removes all loops so that gene sets can be nested within one another: creates an optimal directed acyclic graph
#' pointing from general GO terms down to specific ones.
#' 
#' @param nodes data_frame significance of gene sets
#' @param edges data_frame go_id1 are parents of go_id2
#' 
#' @return an igraph object of the gene set minimal spanning tree
find_minimal_rooted_tree <- function(nodes, edges){
  
  if(clusters(graph_from_data_frame(edges))$no != 1){
    stop("edges must be a single connected network") 
  }
  
  edge_weights <- edges %>%
    dplyr::left_join(nodes %>% dplyr::select(go_id1 = ID, nWeight_1 = nWeight), by = "go_id1") %>%
    dplyr::left_join(nodes %>% dplyr::select(go_id2 = ID, nWeight_2 = nWeight), by = "go_id2") %>%
    dplyr::mutate(weight = nWeight_1 * 0.1 + nWeight_2,
                  weight = ifelse(weight == 0, min(weight[weight != 0])*1e-4, weight))
  
  GO_graph_NEL <- new("graphNEL", nodes = nodes$ID, edgemode="directed")
  GO_graph_NEL <- graph::addEdge(from = edge_weights$go_id1, to = edge_weights$go_id2, graph = GO_graph_NEL, weights = edge_weights$weight)
  
  GO_optim_branching <- RBGL::edmondsOptimumBranching(GO_graph_NEL)
  
  # recreate graphNEL object from edmonds output 
  edmonds_edgeList <- as.data.frame(t(GO_optim_branching$edgeList))
  minimum_igraph_network <- graph_from_data_frame(edmonds_edgeList)
  
  # 
  GO_graph_parsimony <- graph_from_data_frame(edges, vertices = nodes)
  pruned_edges <- attr(igraph::E(GO_graph_parsimony), "vnames")[!(attr(igraph::E(GO_graph_parsimony), "vnames") %in% attr(igraph::E(minimum_igraph_network), "vnames"))]
  GO_graph_parsimony <- igraph::delete_edges(GO_graph_parsimony, pruned_edges)
  
  if(clusters(minimum_igraph_network)$no != 1 | clusters(GO_graph_parsimony)$no != 1){
    stop("a single connected network was not found in the minimal network")
  }
  if(setdiff(igraph::get.edgelist(minimum_igraph_network)[,1],igraph::get.edgelist(minimum_igraph_network)[,2]) != "root"){
    stop("a single root was not found in the minimal network") 
  }
  
  GO_graph_parsimony
}

update_treemap_root <- function(GO_graph_parsimony){
  
  GO_graph_clusters <- igraph::clusters(GO_graph_parsimony)
  
  new_roots <- c()
  for(i in seq_along(GO_graph_clusters$csize)){
    if(GO_graph_clusters$membership['root'] == i){next}
    
    subgraph_nodes <- names(which(GO_graph_clusters$membership == i))
    
    subgraph_graph <- igraph::induced_subgraph(GO_graph_parsimony, subgraph_nodes)
    
    subgraph_root <- setdiff(igraph::get.edgelist(subgraph_graph)[,1], igraph::get.edgelist(subgraph_graph)[,2])
    if(length(subgraph_root) == 0){
      # if there is no parent then either there is one node or the graph is not a DAG
      if(length(subgraph_nodes) == 1){
        subgraph_root <- subgraph_nodes 
      }else{
        # solve the mst problem
        mst_graph <- igraph::minimum.spanning.tree(igraph::induced_subgraph(GO_graph_parsimony, subgraph_nodes))
        
        # figure out which edges were removed
        GO_graph_parsimony <- igraph::delete_vertices(GO_graph_parsimony, E(subgraph_graph)[setdiff(E(subgraph_graph), E(mst_graph))])
        
        subgraph_root <- setdiff(igraph::get.edgelist(mst_graph)[,1], igraph::get.edgelist(mst_graph)[,2])
        
        if(length(subgraph_root) == 0){
          subgraph_root <- subgraph_nodes
        }
      }
    }
    if(length(subgraph_root) == 0){stop("no root found")}
    new_roots <- c(new_roots, subgraph_root)
  }
  
  if(length(new_roots) != 0){
  GO_graph_parsimony <- igraph::add_edges(GO_graph_parsimony, c(t(data_frame("root", new_roots))))
  }
  
  GO_graph_parsimony
}
