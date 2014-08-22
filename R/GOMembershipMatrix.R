#' Construct a new gene set membership matrix to use for performing enrichment
#' analysis
#' 
#' Construct a gene set membership matrix (one row per gene, one column per
#' gene set) from a GO annotation map.
#' 
#' @param annotations A Go3AnnDbBimap object, typically from a Bioconductor Annotation
#' package: for example, org.Hs.egGO from org.Hs.eg.db
#' @param ontology A vector of ontologies to include; should be a subset
#' of c("BP", "MF", "CC")
#' @param evidence A vector of GO evidence codes to include
#' @param min.size Minimum size of a gene set to be included
#' @param max.size Maximum size of a gene set to be included
#' @param ancestors Whether a gene included in a gene set should also be included
#' as being in all ancestors of that gene set
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
#' @import data.table Matrix AnnotationDbi
#' @importFrom GO.db GOTERM
#' 
#' @examples
#' 
#' # yeast membership matrix
#' library(org.Sc.sgd.db)
#' mm = GOMembershipMatrix(org.Sc.sgdGO, min.size=5, max.size=250)
#' 
#' # restrict to Biological Process ontology:
#' mm = GOMembershipMatrix(org.Sc.sgdGO, ontology="BP", min.size=5, max.size=250)
#' 
#' # human membership matrix
#' library(org.Hs.eg.db)
#' mm = GOMembershipMatrix(org.Sc.sgdGO, min.size=5, max.size=250)
#' 
#' @export
GOMembershipMatrix = function(annotations, ontology=NULL, evidence=NULL,
                            min.size=0, max.size=1e6, ancestors=TRUE) {
    # create membership matrix
    frame.dt = data.table(toTable(annotations))
    if (is.null(ontology)) {
        ontology = c("BP", "MF", "CC")
    }
    frame.dt = frame.dt[Ontology %in% ontology, ]
    
    if (!is.null(evidence)) {
        frame.dt = frame.dt[Evidence %in% evidence, ]
    }
    
    id.field = colnames(frame.dt)[1]
    frame.dt = frame.dt[!duplicated(frame.dt[, list(get(id.field), go_id)]), ]
    
    genes = sort(unique(frame.dt[[id.field]]))
    sets = sort(unique(frame.dt$go_id))
    
    m = Matrix(0, nrow=length(genes), ncol=length(sets),
               dimnames=list(genes, sets))
    m[cbind(match(frame.dt[[id.field]], genes),
            match(frame.dt$go_id, sets))] = 1
    
    if (ancestors) {
        # include relationships recursively: if a gene is included in a set,
        # also include it in that set's ancestors
        offspring.tab = do.call(rbind, lapply(ontology, function(o) {
            toTable(get(paste0("GO", o, "OFFSPRING")))
        }))
        offspring.indices = cbind(match(offspring.tab[, 1], sets),
                                  match(offspring.tab[, 2], sets))
        # remove those that don't match anything already in our set
        offspring.indices = offspring.indices[complete.cases(offspring.indices), ]
        offspring.m = Matrix(0, nrow=length(sets), ncol=length(sets),
                             dimnames=list(sets, sets))
        offspring.m[offspring.indices] = 1
        # combine with descendants
        m = m + m %*% offspring.m
        # turn back into binary matrix
        m = (m > 0) + 0
    }

    cs = colSums(m)
    m = m[, cs >= min.size & cs <= max.size]
    
    # create table of set data, matching the membership matrix
    goterms = GOTERM[colnames(m)]
    setData = data.table(ID=colnames(m), Term=Term(goterms),
                         Definition=Definition(goterms), Count=colSums(m))
    
    # create table of gene data
    geneData = data.table(ID=rownames(m), Count=rowSums(m))
    
    return(new("GeneMatrix", matrix=m, colData=setData,
               geneData=geneData))
}
