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
  
  function(m, minimal_edge_set){
    
    for(a_gene in m@geneData$ID){
      
      gene_GO <- m@matrix[a_gene, m@matrix[a_gene,] == 1] %>% names()
      gene_GO <- c(gene_GO[gene_GO %in% nodes$ID], "root")
      
      gene_subgraph <- igraph::induced_subgraph(minimal_edge_set, gene_GO)
      gene_subgraph_clusters <- igraph::clusters(gene_subgraph)
      gene_subgraph <- igraph::induced_subgraph(minimal_edge_set, names(gene_subgraph_clusters$membership)[gene_subgraph_clusters$membership == gene_subgraph_clusters$membership['root']])
      
      gene_paths <- get.shortest.paths(gene_subgraph, from = "root", to = V(gene_subgraph))$vpath
      lapply(seq_along(gene_paths), function(i){
        gene_paths[[i]]
      })
      
      
      gene_subset_path <- reshape2::melt(get.shortest.paths(gene_subgraph, from = "root", to = V(gene_subgraph))[[1]], level = 1)
      colnames(gene_subset_path) <- c("path", "leaf")
      gene_subset_path$path <- V(gene_graph_subset)$name[gene_subset_path$path]
      gene_subset_path$leaf <- V(gene_graph_subset)$name[gene_subset_path$leaf]
      
      
      
    }
        
        
    
    
    
    
    
    
  }
  
  gene_graph_subset <- delete.vertices(gene_graph_subset, V(gene_graph_subset)$name[!(V(gene_graph_subset)$name %in% systematic_to_GOIDs$go_id[systematic_to_GOIDs$systematic_name == a_gene])])
    
    
  
  
  
  
  
  
ggraph(GO_graph_parsimony, 'igraph', algorithm = 'tree') + 
  geom_edge_link() +
  ggforce::theme_no_axes()


}

#' Find Minimal Rooted Tree
#' 
#' Removes all loops so that gene sets can be nested within one another: creates an optimal directed acyclic graph
#' pointing from general GO terms down to specific ones.
#' 
#' @param nodes data_frame significance of gene sets
#' @param edges data_frame go_id1 are parents of go_id2
#' @param ancestors data_frame go_id1 are ancestors of go_id2
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
  NEL_edges <- list()
  for(vertex in GO_optim_branching$nodes){
    # a named list pointing to children, including associated weights
    NEL_edges[[vertex]]$edges <- edmonds_edgeList$to[edmonds_edgeList$from == vertex]
    NEL_edges[[vertex]]$weigths <- GO_optim_branching$weights[edmonds_edgeList$from == vertex]
  }
  
  GO_graph_parsimony <- igraph::igraph.from.graphNEL(new("graphNEL", nodes = GO_optim_branching$nodes, edgeL = NEL_edges, edgemode="directed"))
  GO_graph_edgelist <- igraph::get.edgelist(GO_graph_parsimony)
  
  # remove self-links that were added by graphNel
  GO_graph_parsimony <- delete.edges(GO_graph_parsimony, which(GO_graph_edgelist[,1] == GO_graph_edgelist[,2]))
  
  GO_graph_clusters <- igraph::clusters(GO_graph_parsimony)
  
  if(length(unique(GO_graph_edgelist[,2])) != nrow(GO_graph_edgelist)){
    stop("Some GO categories have multiple parents; edmonds algorithm failed")
  }
  if(GO_graph_clusters$no != 1){
    message("Some GO categories were seperated from the root; these branches will be rerooted")
    
    new_roots <- c()
    for(i in seq_along(GO_graph_clusters$csize)){
      if(GO_graph_clusters$membership['root'] == i){next}
      
      subgraph_nodes <- names(which(GO_graph_clusters$membership == i))
      
      subgraph_edgelist <- igraph::induced_subgraph(GO_graph_parsimony, subgraph_nodes) %>%
        igraph::get.edgelist()
      
      subgraph_root <- setdiff(subgraph_edgelist[,1], subgraph_edgelist[,2])
      if(length(subgraph_root) == 0){
       subgraph_root <-  subgraph_nodes
      }
      
      new_roots <- c(new_roots, subgraph_root)
    }
    
    $
    
    }
  
  GO_graph_parsimony  
  }


LASSO_guided_GO_tree <- function(GO_inheritance_subset, GOweights, involvedGenes){
  
  # tkplot(GO_graph_parsimony)
  
  #### Use the simplified tree to assign each gene to the most predicitve subgraph (Highest sum of LASSO) ####
  # We now know that each child node has one parents
  
  unique_gene_ancestry <- NULL
  
  for(a_gene in unique(systematic_to_GOIDs$systematic_name)){
    
    gene_graph_subset <- GO_graph_parsimony
    # delete all categories nodes not associated with the gene of interest
    gene_graph_subset <- delete.vertices(gene_graph_subset, V(gene_graph_subset)$name[!(V(gene_graph_subset)$name %in% systematic_to_GOIDs$go_id[systematic_to_GOIDs$systematic_name == a_gene])])
    
    if(length(V(gene_graph_subset)$name) == 1){
      
      unique_gene_ancestry <- rbind(unique_gene_ancestry, data.frame(gene = a_gene, sequence = V(gene_graph_subset)$name, level = 1, L1 = 0))
      
    }else{
      
      edgeL_subset <- as.data.frame(get.edgelist(gene_graph_subset))
      # find the weights of the chosen edges
      edgeL_subset <- edgeL_subset %>% dplyr::select(parent = V1, child = V2) %>% left_join(GO_inheritance_subset, by = c("parent", "child"))
      edgeL_subset$LASSO[edgeL_subset$Category == "Uninformative"] <- 0
      
      root <- unique(edgeL_subset$parent[!(edgeL_subset$parent %in% edgeL_subset$child)]) # identify root by it only serving as a source
      
      gene_subset_path <- melt(get.shortest.paths(gene_graph_subset, from = root, to = V(gene_graph_subset))[[1]], level = 1)
      colnames(gene_subset_path) <- c("path", "leaf")
      gene_subset_path$path <- V(gene_graph_subset)$name[gene_subset_path$path]
      gene_subset_path$leaf <- V(gene_graph_subset)$name[gene_subset_path$leaf]
      
      path_score <- function(path){
        data.frame(parent = path[1:(length(path)-1)], child = path[2:length(path)]) %>%
          left_join(edgeL_subset, by = c("parent", "child")) %>% summarize(LASSO = sum(LASSO)) %>% unlist() %>% unname()  
      }
      
      path_scores <- gene_subset_path %>% group_by(leaf) %>% dplyr::summarize(
        LASSO = path_score(path), steps = n()
      )
      
      if(any(edgeL_subset$Category == "LASSO")){
        # At least one path includes a LASSO node - choose the path with the highest score
        GOancestry <- path_scores %>% arrange(-LASSO, steps) %>% dplyr::slice(1) %>% left_join(gene_subset_path, by = "leaf") 
        
        GOancestry <- GOancestry %>% left_join(
          GOancestry %>% dplyr::select(path) %>% left_join(GOweights %>% dplyr::select(path = ID, L1 = beta_1se), by = "path") %>%
            dplyr::mutate(L1 = ifelse(!is.na(L1), abs(L1), 0)),
          by = c("path"))
        
      }else{
        # No path includes a LASSO node - choose the path with the highest ancestor weight and then prune at the 4th level
        GOancestry <- path_scores %>% arrange(-LASSO, -steps) %>% slice(1) %>% left_join(gene_subset_path, by = "leaf") %>%
          dplyr::mutate(L1 = 0)
        if(nrow(GOancestry) > 4){
          GOancestry <- GOancestry %>% slice(1:4) 
        }
      }
      
      GOancestry <- GOancestry %>% dplyr::mutate(gene = a_gene, level = 1:length(path)) %>% dplyr::select(gene, sequence = path, level, L1)
      
      unique_gene_ancestry <- rbind(unique_gene_ancestry, GOancestry)
      
    }
  }
}