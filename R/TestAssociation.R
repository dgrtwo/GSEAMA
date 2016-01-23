#' Test association of a per-gene attribute with each column in a GeneMatrix
#' 
#' Test the association of a per-gene attribute y with each column, where a column
#' can represent a gene set, a motif, transcription factor targets, or other genes
#' that are functionally related
#' 
#' @param m A GeneMatrix object
#' @param genes A vector of gene names (must match the gene names in the GeneMatrix)
#' @param y A per-gene metric of the same length as the vector of gene names
#' @param method Method used to test association between the per-gene metric and
#' each column - either c("lasso", "t.test", "wilcoxon" or "hypergeometric")
#' @param ... Additional arguments that will be given to the test in question
#' 
#' @details A t-test or Wilcoxon test can take an additional argument of \code{alternative},
#' which allows a one-sided test with an alternative hypothesis of "less" or "greater".
#' A t=test can also take the \code{var.equal} argument, indicating whether to treat the
#' variance in and outside each group as identical.
#' 
#' @export
TestAssociation = function(m, genes, y, method = "lasso", ...) {
    if (!inherits(m, "GeneMatrix")) {
        stop("TestAssociation should be given a GeneMatrix object")
    }

    if (!is.null(m@fit)) {
        warning(paste0("Performing a test on a model that has already been ",
                       "tested, previous results may be overwritten"))
    }
    
    # prefiltering of membership matrix
    v <- genes %in% m@geneData$ID
    y <- y[v]
    genes <- genes[v]
    m <- m[genes, ]
    m <- m[, colSums(m@matrix != 0) > 0]  # filter columns *after* rows
    m@geneData$y <- y

    # apply desired method
    method = match.arg(method, c("lasso", "wilcoxon", "hypergeometric", "t.test"))
    if (method == "lasso") {
        m@fit <- glmnet::cv.glmnet(m@matrix, y, ...)
        m@colData$beta <- apply(m@fit$glmnet.fit$beta, 1,
                                function(r) r[r != 0][1])
        m@colData$step <- apply(m@fit$glmnet.fit$beta, 1,
                                function(r) which(r != 0)[1])
        m@rankingMetric <- "step"
    }
    else if (method == "t.test") {
        m@colData$p.value <- vectorized_t_test(m@matrix > 0, y, ...)
        m@rankingMetric <- "pvalue"
    }
    else if (method == "wilcoxon") {
        m@colData$p.value <- vectorized_wilcoxon(m@matrix > 0, y, ...)
        m@rankingMetric <- "p.value"
    }
    else if (method == "hypergeometric") {
        if (!is.logical(y) & !(is.numeric(y) & all(y == 0 | y == 1))) {
            stop("For hypergeometric test, y must be either logical or binary")
        }
        mat <- m@matrix > 0
        overlaps <- colSums(mat[as.logical(y), ])
        totals <- colSums(mat)
        m@colData$p.value <- phyper(overlaps - 1,      # number of white balls in each set
                                    sum(y),            # white balls in urn
                                    nrow(m) - sum(y),  # black balls in urn
                                    totals,            # number drawn in each set
                                    lower.tail = FALSE)
        m@rankingMetric <- "p.value"
    }
    m
}
