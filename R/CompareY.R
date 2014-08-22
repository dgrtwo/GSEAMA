#' Create a violin plot, boxplot or histogram looking at the metric of interest
#' over one or more columns
#' 
#' This compares the metric of interest (y) across the given genes, as well as
#' comparing it to the overall distribution. This shows a visual representation
#' of how a column (a gene set, a motif, a transcription factor, etc) is associated
#' with the metric.
#' 
#' @param m a GeneMatrix
#' @param columns The ID of one or more columns to plot
#' @param mode Type of plot to create: either "violin", "boxplot" or "histogram"
#' @param charlimit Maximum number of characters to include in column names
#' 
#' @export
CompareY = function(m, columns, mode="violin", charlimit=35) {
    # graph one or more columns 
    # reverse the list (so that it reads top down instead of bottom up)
    columns = rev(columns)
    coldata = m@colData[match(columns, ID), ]

    y = rep(m@geneData$y, length(columns))
    dat = data.frame(y=y, Present=c(as.matrix(m@matrix[, columns]) > 0))
    
    # fix the terms so they are not too long
    terms = strtrim(coldata$Term, charlimit)
    terms[nchar(coldata$Term) > charlimit] = paste0(
        terms[nchar(coldata$Term) > charlimit], "...")
    dat$Set = rep(paste0(terms, " (", coldata$Count, ")"), each=NROW(m))
    dat$Set = factor(dat$Set, levels=dat$Set[!duplicated(dat$Set)])
    
    dat = dat[dat$Present, ]
    dat = rbind(dat, data.frame(y=m@geneData$y, Present=TRUE, Set="Overall"))
    dat$y = as.numeric(as.character(dat$y))
    
    mode = match.arg(mode, c("violin", "boxplot", "histogram"))
    if (mode == "histogram") {
        g = ggplot(dat, aes(y)) + geom_histogram() + facet_wrap(~ Set, ncol=1, scale="free_y")
    }
    else {
        g = (ggplot(dat, aes(x=Set, y=y)) + coord_flip() +
             theme(axis.text.y=element_text(color="black", size=15)))
        if (mode == "violin") {
            g = g + geom_violin()
        }
        else {
            g = g + geom_boxplot()
        }
    }
    g
}

#' Create a graph comparing the top n columns (in terms of association with y)
#' 
#' This compares the metric of interest (y) across the 
#'
#' @param m GeneMatrix object
#' @param n Number of genes to compare
#' @param ... Additional arguments to be passed to CompareY
#' 
#' @details Most methods ("t.test", "wilcoxon", "hypergeometric") use the
#' p-value to find the most significant columns. The "lasso" uses the order
#' in which the terms were added to the regression
#' 
#' @export
CompareTopColumns = function(m, n=9, ...) {
    if (length(m@rankingMetric) == 0) {
        stop("Cannot plot top genes without first running a test")
    }
    CompareY(m, m@colData[order(m@colData[[m@rankingMetric]]), ]$ID[1:n],
             ...)
}
