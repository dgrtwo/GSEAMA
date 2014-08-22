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
#' @details A t.test or wilcoxon test can take an additional argument of \code{alternative},
#' which allows a one-sided test with an alternative hypothesis of "less" or "greater".
#' A t.test can also take the \code{var.equal} argument, indicating whether to treat the
#' variance in and outside each group as identical.
#' 
#' @export
TestAssociation = function(m, genes, y, method="lasso", ...) {
    if (!inherits(m, "GeneMatrix")) {
        stop("TestAssociation should be given a GeneMatrix object")
    }

    if (!is.null(m@fit)) {
        warning(paste0("Performing a test on a model that has already been ",
                       "tested, previous results may be overwritten"))
    }
    
    # prefiltering of membership matrix
    v = genes %in% rownames(m@matrix)
    y = y[v]
    genes = genes[v]
    m = m[genes, ]
    m = m[, colSums(m@matrix) > 0]  # filter columns *after* rows
    m@geneData$y = y
    
    # filter so that every column has at least one nonzero value
    m = m[, colSums(m@matrix != 0) > 0]
    
    # apply desired method
    method = match.arg(method, c("lasso", "wilcoxon", "hypergeometric", "t.test"))
    if (method == "lasso") {
        library(glmnet)
        m@fit = cv.glmnet(m@matrix, y, ...)
        min.lambda.index = m@fit$glmnet.fit$lambda == m@fit$lambda.min
        m@colData$beta= apply(m@fit$glmnet.fit$beta, 1,
                               function(r) r[r != 0][1])
        m@colData$step = apply(m@fit$glmnet.fit$beta, 1,
                                function(r) which(r != 0)[1])
        m@rankingMetric = "step"
        #x = as.matrix(m[genes, ]@matrix)
        #a = lars(x, y)
        #covTest(a, x, y)
    }
    else if (method == "t.test") {
        m@colData$pvalue = vectorized.t.test(m@matrix > 0, y, ...)
        m@rankingMetric = "pvalue"
    }
    else if (method == "wilcoxon") {
        m@colData$pvalue = vectorized.wilcoxon(m@matrix > 0, y, ...)
        m@rankingMetric = "pvalue"
    }
    else if (method == "hypergeometric") {
        if (!is.logical(y) & !(is.numeric(y) & all(y %in% c(0, 1)))) {
            stop("For hypergeometric test, y must be either logical or binary")
        }
        m = m@matrix > 0
        overlaps = colSums(m[as.logical(y), ])
        totals = colSums(m)
        m@colData$pvalue = phyper(overlaps, totals,
                                   NROW(m) - totals, sum(y), lower.tail=FALSE)
        m@rankingMetric = "pvalue"
    }
    m
}

#' Perform a Student's T-test comparing a metric to each column of a matrix
#' 
#' This function is similar to \code{apply(m, 2, function(col) t.test(y ~ col)$p.value)},
#' but is vectorized to make it much faster
vectorized.t.test = function(m, y, var.equal=FALSE, alternative="two.sided") {
    stopifnot(NROW(m) == length(y))
    
    in.m = y * m
    out.m = y * (1 - m)
    n.in = colSums(m)
    n.out = NROW(m) - n.in
    
    in.mu = colSums(in.m) / n.in
    out.mu = colSums(out.m) / n.out
    
    v.in = colSums((in.m - t(t(m) * in.mu))^2) / (n.in - 1)
    v.out =colSums((out.m - t(t(m) * out.mu))^2) / (n.out - 1)
    
    if (var.equal) {
        df = NROW(m) - 2
        v = (n.in - 1) * vx + (n.out - 1) * vx
        v = v / df
        stderr <- sqrt(v*(1/n.in+1/n.out))
    }
    else {
        stderr.in <- sqrt(v.in/n.in)
        stderr.out <- sqrt(v.out/n.out)
        stderr <- sqrt(stderr.in^2 + stderr.out^2)
        df <- stderr^4/(stderr.in^4/(n.in-1) + stderr.out^4/(n.out-1))
    }
    
    tstat <- (in.mu - out.mu) / stderr
    
    switch(alternative, less=pt(tstat, df),
           greater=pt(tstat, df, lower.tail = FALSE),
           two.sided=2 * pt(-abs(tstat), df))
}


#' Perform a Wilcoxon rank-sum test comparing a metric to each column of a matrix
#' 
#' This function is similar to \code{apply(m, 2, function(col) wilcox.test(y ~ col)$p.value)},
#' but is vectorized to make it \emph{much} faster
vectorized.wilcoxon = function(m, y, alternative="two.sided") {
    # given a boolean matrix and a vector y, apply the wilcoxon test to see
    # if y depends on each column of the matrix, returning a vector of
    # p-values
    stopifnot(NROW(m) == length(y))
    
    rk = rank(y)
    n.x = colSums(m)
    n.y = NROW(m) - n.x
    
    STATISTIC = colSums(rk * m) - n.x * (n.x + 1) / 2
    NTIES = table(rk)
    
    z <- STATISTIC - n.x * n.y / 2
    
    CORRECTION <- switch(alternative,
                         "two.sided" = sign(z) * 0.5,
                         "greater" = 0.5,
                         "less" = -0.5)
    
    SIGMA <- sqrt((n.x * n.y / 12) *
                      ((n.x + n.y + 1)
                       - sum(NTIES^3 - NTIES)
                       / ((n.x + n.y) * (n.x + n.y - 1))))
    z <- (z - CORRECTION) / SIGMA
    
    PVAL <- switch(alternative,
                   "less" = pnorm(z),
                   "greater" = pnorm(z, lower.tail=FALSE),
                   "two.sided" = 2 * pmin(pnorm(z),
                                          pnorm(z, lower.tail=FALSE)))
    PVAL
}
