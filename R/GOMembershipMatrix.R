#' Construct a new gene set membership matrix
#' 
#' Construct a gene set membership matrix (one row per gene, one column per
#' gene set, with 1 indicating membership) from a GO annotation map.
#' 
#' @param annotations A Go3AnnDbBimap object, typically from a Bioconductor Annotation
#' package. For example, org.Hs.egGO from org.Hs.eg.db for human genes.
#' @param ontology A vector of ontologies to include; should be a subset
#' of c("BP", "MF", "CC")
#' @param evidence A vector of GO evidence codes to include
#' @param min_size Minimum size of a gene set to be included
#' @param max_size Maximum size of a gene set to be included
#' @param ancestors Whether a gene included in a gene set should also be included
#' as being in all ancestors of that gene set (default TRUE)
#' @param chosen_genes If provided, restricts genes to those in chosen_genes
#' 
#' The GO annotation maps typically come from the Bioconductor AnnotationDB packages.
#' A couple of notable examples are:
#' 
#' Homo sapiens: org.Hs.egGO
#' S. cerevisiae: org.Sc.sgdGO
#' E coli K12: org.EcK12.egGO
#' 
#' You can restrict the sets to the BP (Biological Process), MF (Molecular Function),
#' or CC (Cellular Compartment) ontologies (by default all are included).
#' 
#' @import Matrix AnnotationDbi
#' @importFrom dplyr %>% data_frame tbl_df filter mutate group_by arrange distinct_
#' 
#' @examples
#' 
#' # yeast membership matrix
#' library(org.Sc.sgd.db)
#' mm <- GOMembershipMatrix(org.Sc.sgdGO, min_size = 5, max_size = 250)
#' 
#' # restrict to Biological Process ontology:
#' mm <- GOMembershipMatrix(org.Sc.sgdGO, ontology = "BP", min_size = 5, max_size = 250)
#' 
#' # human membership matrix
#' library(org.Hs.eg.db)
#' mm <- GOMembershipMatrix(org.Hs.sgdGO, min_size = 5, max_size = 250)
#' 
#' @export
GOMembershipMatrix = function(annotations,
                                ontology = c("BP", "MF", "CC"),
                                evidence = NULL,
                                min_size = 1,
                                max_size = Inf,
                                ancestors = TRUE,
                                chosen_genes = NULL) {
  # create membership matrix
  frame_df <- annotations %>%
    AnnotationDbi::toTable() %>%
    dplyr::tbl_df() %>%
    dplyr::filter(Ontology %in% ontology)
    #filter(go_id %in% available_go_ids)
  
  ontology_tbl <- frame_df %>%
    dplyr::select(ID = go_id, Ontology) %>%
    dplyr::distinct(ID, .keep_all = TRUE)

  if (!is.null(evidence)) {
    frame_df <- frame_df %>%
      dplyr::filter(evidence %in% Evidence)
  }
  
  id_field <- colnames(frame_df)[1]
  frame_df <- frame_df %>%
    dplyr::distinct_(id_field, "go_id", .keep_all = TRUE)
  
  if (!is.null(chosen_genes)) {
    frame_df <- frame_df[frame_df[[id_field]] %in% genes, ]
  }
  
  # turn into a sparse membership matrix
  m <- sparse_cast_(frame_df, id_field, "go_id", 1)

  if (ancestors) {
    # include relationships recursively: if a gene is included in a set,
    # also include it in sets that are ancestors of that set
    ancestry_matrix <- get_ancestry_matrix(colnames(m))
    
    # combine with descendants, turn back into binary matrix
    m <- m + m %*% ancestry_matrix
    m = (m > 0) + 0
  }
  
  # eliminate sets with too few genes, and then genes with no sets
  set_size = colSums(m)
  m <- m[, set_size >= min_size & set_size <= max_size]
  m <- m[rowSums(m) > 0, ]
  sets <- colnames(m)
  genes <- rownames(m)
  
  # create table of set data, matching the membership matrix
  set_data <- data_frame(ID = sets,
                        Term = Term(sets),
                        Definition = Definition(sets),
                        Size = colSums(m)) %>%
    dplyr::inner_join(ontology_tbl, by = "ID")
  
  # just make sure nothing was reordered
  stopifnot(!all(sets == set_data))
  
  # create table of gene data
  gene_data <- data_frame(ID = rownames(m), Size = rowSums(m))
  
  ret <- new("GeneMatrix",
             matrix = m,
             colData = set_data,
             geneData = gene_data)

  ret
}


#' Build an offspring matrix of GO terms
#' 
#' A sparse binary Matrix object with one row and column for each pair
#' of GO terms provided, where each row represents an ancestor and each column
#' represents a descendant, with 1 marking ancestor/descendant pairs.
#' 
#' @param terms IDs of GO terms that should be included in the ancestry matrix
#' @param ontology Ontologies to use
#' @param type Either OFFSPRING (default), CHILDREN, ANCESTOR, or PARENT
#' 
#' @import GO.db
#' 
#' @export
get_ancestry_matrix <- function(terms, ontology = c("BP", "MF", "CC"),
                                type = "OFFSPRING", upward = TRUE, tbl = FALSE) {
  offspring_df <- do.call(rbind, (lapply(ontology, function(o) {
    toTable(get(paste0("GO", o, type)))
  })))[, 1:2]
  
  ret <- offspring_df %>%
    setNames(c("go_id1", "go_id2")) %>%
    filter(go_id1 %in% terms, go_id2 %in% terms) %>%
    mutate(go_id1 = factor(go_id1, levels = terms),
           go_id2 = factor(go_id2, levels = terms)) %>%
    tbl_df()
  
  if (!upward) {
    tmp <- ret$go_id1
    ret$go_id1 <- ret$go_id2
    ret$go_id2 <- tmp
  }
  
  if (!tbl) {
    ret <- sparse_cast(ret, go_id1, go_id2)
  }
  
  ret
}


#' Get all GO IDs in all_IDs that are ancestors/descendants of those in IDs
#' 
#' @param x A vector of GO IDs
#' @param all_IDs A vector of GO IDs to look for relatives of x
#' @param ancestors If true look among ancestors, if false look for children
#' @param combine Whether to combine it with the sets in the IDs parameter
GetGORelatives <- function(IDs, all_IDs,
                             ancestors = TRUE,
                             direct = FALSE,
                             combine = TRUE) {
  type <- if (direct) { "CHILDREN" } else { "OFFSPRING" }
  am <- get_ancestry_matrix(all_IDs, type = type, upward = ancestors)
  ret <- names(which(colSums(am[IDs, , drop = FALSE]) > 0))
  
  if (combine) {
    ret <- unique(c(IDs, ret))
  }
  ret
}
