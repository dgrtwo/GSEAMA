\name{GOMembershipMatrix}
\alias{GOMembershipMatrix}
\title{Construct a new gene set membership matrix to use for performing enrichment
analysis}
\usage{
GOMembershipMatrix(annotations, ontology = NULL, evidence = NULL,
  min.size = 0, max.size = 1e+06, ancestors = TRUE)
}
\arguments{
  \item{annotations}{A Go3AnnDbBimap object, typically from
  a Bioconductor Annotation package: for example,
  org.Hs.egGO from org.Hs.eg.db}

  \item{ontology}{A vector of ontologies to include; should
  be a subset of c("BP", "MF", "CC")}

  \item{evidence}{A vector of GO evidence codes to include}

  \item{min.size}{Minimum size of a gene set to be
  included}

  \item{max.size}{Maximum size of a gene set to be
  included}

  \item{ancestors}{Whether a gene included in a gene set
  should also be included as being in all ancestors of that
  gene set

  The GO annotation maps typically come from the
  Bioconductor AnnotationDB packages. A couple of notable
  examples are:

  Homo sapiens: org.Hs.egGO S. cerevisiae: org.Sc.sgdGO E
  coli K12: org.EcK12.egGO

  You can restrict the sets to the BP (Biological Process),
  MF (Molecular Function), or CC (Cellular Compartment)
  ontologies (by default all are included).}
}
\description{
Construct a gene set membership matrix (one row per gene,
one column per gene set) from a GO annotation map.
}
\examples{
# yeast membership matrix
library(org.Sc.sgd.db)
mm = GOMembershipMatrix(org.Sc.sgdGO, min.size=5, max.size=250)

# restrict to Biological Process ontology:
mm = GOMembershipMatrix(org.Sc.sgdGO, ontology="BP", min.size=5, max.size=250)

# human membership matrix
library(org.Hs.eg.db)
mm = GOMembershipMatrix(org.Sc.sgdGO, min.size=5, max.size=250)
}

