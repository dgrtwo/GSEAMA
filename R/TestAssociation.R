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
        beta_ind <- which(m@fit$lambda.1se == m@fit$lambda)
        m@colData$beta <- apply(m@fit$glmnet.fit$beta, 1,
                                function(r) r[r != 0][1])
        m@colData$beta1se <- as.numeric(m@fit$glmnet.fit$beta[, as.numeric(beta_ind)])
        m@colData$step <- apply(m@fit$glmnet.fit$beta, 1,
                                function(r) which(r != 0)[1])
        m@colData$absbeta1se <- abs(m@colData$step)
        
        m@rankingMetric <- "absbeta1se"
        m@effectMetric <- "beta"
    }
    else if (method == "t.test") {
        tt <- vectorized_t_test(m@matrix > 0, y, tbl = TRUE, ...)
        m@colData$estimate <- tt$estimate
        m@colData$p.value <- tt$p.value
        m@rankingMetric <- "p.value"
        m@effectMetric <- "estimate"
    }
    else if (method == "wilcoxon") {
        w <- vectorized_wilcoxon(m@matrix > 0, y, tbl = TRUE, ...)
        m@colData$auc <- w$auc
        m@colData$p.value <- w$p.value
        m@rankingMetric <- "p.value"
        m@effectMetric <- "auc"

    }
    else if (method == "hypergeometric") {
        if (!is.logical(y) & !(is.numeric(y) & all(y == 0 | y == 1))) {
            stop("For hypergeometric test, y must be either logical or binary")
        }
        stop("Not yet implemented")
        mat <- m@matrix > 0
        
        overlaps <- colSums(mat[as.logical(y), ])
        
        total_genes <- nrow(m)
        total_y <- sum(y)
        totals <- colSums(mat)
        
        #m@colData$phi <- 
        m@colData$p.value <- phyper(overlaps - 1,      # number of white balls in each set
                                    total_y,            # white balls in urn
                                    total_genes - total_genes,  # black balls in urn
                                    totals,            # number drawn in each set
                                    lower.tail = FALSE)
        m@rankingMetric <- "p.value"
        m@effectMetric <- "phi"
    }
    
    # finally, compute the difference within and between the sets
    mat <- m@matrix
    not_mat <- 1 - mat
    mean_within <- colSums(y * mat) / colSums(mat)
    mean_outside <- colSums(y * not_mat) / colSums(not_mat)
    m@colData$MeanDifference <- mean_within - mean_outside

    m
}


#' Compute mean differences between "in a set" and "not in a set"
#' 
#' For each set in a GeneMatrix, compute the difference in y between
#' "in the set"
#' and "outside the set". It will be set in the MeanDifference column
#' 
#' @param m GeneMatrix object
#' 
#' @export
mean_difference <- function(m) {
  y <- m@geneData$y
  
  if (is.null(y)) {
    stop("y not yet set, has TestAssociation been performed?")
  }
  
  sub_m <- m@matrix
  sub_m_not <- 1 - sub_m
  mean_within <- colSums(y * sub_m) / colSums(sub_m)
  mean_outside <- colSums(y * sub_m_not) / colSums(sub_m_not)
  m@colData$MeanDifference <- mean_within - mean_outside
  
  m
}


#' Return a vector of genes that pass a threshold for includions
#' 
#' Return a vector of genes that pass a threshold of "significance"
#' or other kind of inclusion. For tests with a p-value (t-test,
#' Wilcoxon, hypergeometric), this is FDR-controlled p-values.
#' For LASSO, this is all cases where beta1sd != 0.
#' 
#' @param m GeneMatrix object
#' @param alpha Threshold for (corrected) p-values, not used in LASSO
#' @param method Method, passed to \code{\link{p.adjust}}, or "qvalue"
#' to use qvalue. Not used in LASSO
#' @param ... Extra arguments to pass on to correction method
#' 
#' @seealso \link{p.adjust}, \link{p.adjust.methods}
ThresholdSets <- function(m, alpha = .05, method = "fdr", ...) {
  if (!is.null(m@colData$beta1se)) {
    # LASSO, filter by beta1se
    return(m@colData$ID[m@colData$beta1se != 0])
  }
  
  if (method == "qvalue") {
    p.adjusted <- qvalue::qvalue(m@colData$p.value, ...)$qvalues
  } else {
    p.adjusted <- p.adjust(m@colData$p.value, method, ...)
  }
  m@colData$ID[!is.na(p.adjusted) & p.adjusted < alpha]
}
