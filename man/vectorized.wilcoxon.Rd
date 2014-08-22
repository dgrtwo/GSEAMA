\name{vectorized.wilcoxon}
\alias{vectorized.wilcoxon}
\title{Perform a Wilcoxon rank-sum test comparing a metric to each column of a matrix}
\usage{
vectorized.wilcoxon(m, y, alternative = "two.sided")
}
\description{
This function is similar to \code{apply(m, 2, function(col)
wilcox.test(y ~ col)$p.value)}, but is vectorized to make
it \emph{much} faster
}

