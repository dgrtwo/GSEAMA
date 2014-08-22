\name{CompareTopColumns}
\alias{CompareTopColumns}
\title{Create a graph comparing the top n columns (in terms of association with y)}
\usage{
CompareTopColumns(m, n = 9, ...)
}
\arguments{
  \item{m}{GeneMatrix object}

  \item{n}{Number of genes to compare}

  \item{...}{Additional arguments to be passed to CompareY}
}
\description{
This compares the metric of interest (y) across the
}
\details{
Most methods ("t.test", "wilcoxon", "hypergeometric") use
the p-value to find the most significant columns. The
"lasso" uses the order in which the terms were added to the
regression
}

