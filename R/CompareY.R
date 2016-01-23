#' Return all genes, along with the saved metric in each of the given sets
#' 
#' @param m a GeneMatrix object
#' @param columns vector of columns (e.g. gene sets) whose genes should be extracted
#' @param char_limit optionally shorten the Term field to the desired length
#' @param add_count whether the count within each set should be added to the
#' Term field (like "Metabolism (360)")
#' 
#' @importFrom dplyr inner_join do
#' 
#' @export
ExtractGenes <- function(m, columns, char_limit = NULL, add_count = FALSE,
                         overall = FALSE) {
    # repeat gene data, once for each column
    dat <- data_frame(Set = columns) %>%
      group_by(Set) %>%
      do(m@geneData)

    dat$Present <- c(as.matrix(m@matrix[, columns]) > 0)
    
    dat <- dat %>%
      filter(Present == TRUE) %>%
      dplyr::select(-Present, -Size) %>%
      inner_join(m@colData, by = c(Set = "ID")) %>%
      arrange(match(Set, columns))
    
    if (!is.null(char_limit)) {
        # shorten the Term
        dat <- dat %>%
          mutate(Term = trim_ellipses(Term, char_limit))
    }
    if (add_count) {
        # add parenthetical counts to each term
        dat <- dat %>%
          mutate(Term = paste0(Term, " (", Size, ")"))
    }
    if (overall) {
        dat <- bind_rows(dat, cbind(m@geneData, Term = "Overall"))
    }
    
    # place terms in order of appearance
    dat$Term <- factor(dat$Term, levels = dat$Term[!duplicated(dat$Term)])
    dat
}


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
#' @param char_limit Maximum number of characters to include in column names
#' 
#' @import ggplot2
#' 
#' @export
CompareY = function(m, columns, mode = "boxplot", char_limit = 35) {
  # graph one or more columns 
  # reverse the list (so that it reads top down instead of bottom up)
  dat <- ExtractGenes(m, rev(columns),
                      char_limit = char_limit,
                      add_count = TRUE,
                      overall = TRUE)
  
  mode <- match.arg(mode, c("boxplot", "violin", "histogram"))
  if (mode == "histogram") {
    g <- ggplot(dat, aes(y)) +
      geom_histogram() +
      facet_wrap(~Term, ncol = 1, scale = "free_y")
  }
  else {
    g <- ggplot(dat, aes(x = Term, y = y)) +
      coord_flip() +
      theme(axis.text.y = element_text(color = "black", size = 15))
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
#' Plot a bocplot, violin plot, or faceted histogram with the distributions
#' of y within the most associated gene sets.
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
CompareTopColumns = function(m, n = 9, ...) {
  if (length(m@rankingMetric) == 0) {
    stop("Cannot plot top genes without first running a test")
  }

  top_sets <- m@colData %>%
    arrange_(m@rankingMetric) %>%
    head(n)
    
  CompareY(m, top_sets$ID, ...)
}
