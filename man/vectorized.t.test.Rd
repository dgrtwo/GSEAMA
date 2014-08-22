\name{vectorized.t.test}
\alias{vectorized.t.test}
\title{Perform a Student's T-test comparing a metric to each column of a matrix}
\usage{
vectorized.t.test(m, y, var.equal = FALSE, alternative = "two.sided")
}
\description{
This function is similar to \code{apply(m, 2, function(col)
t.test(y ~ col)$p.value)}, but is vectorized to make it
much faster
}

