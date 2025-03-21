% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tables.R
\name{gm_mean}
\alias{gm_mean}
\title{Geometric mean}
\usage{
gm_mean(x, offset = 0, zero.propagate = FALSE, conf.level = 0.95, na.rm = TRUE)
}
\arguments{
\item{x}{numeric vector as input}

\item{offset}{Optional variable: adds a value to the data (1 is often used),
so that no more zeros occur in the data. This is an acceptable way of
handling zeros when, for example, zeros represent undetectable low
concentrations in a biomedical laboratory assay.}

\item{zero.propagate}{if TRUE, the output will be zero if any value is zero.
Overrides the offset parameter is true.}

\item{conf.level}{If not missing, the confidence interval will be given together with the geometric mean}

\item{na.rm}{logical, indicating whether NA values should be excluded}
}
\value{
a vector of length 1 or, if conf.level is not missing, length 3
}
\description{
Calculates the geometric mean with confidence interval.
If conf.level is given, a vector will be output containing the geometric mean,
the lower bound confidence value, and the upper bound confidence value.
Output can be used directly in the function
\code{\link{format_values}} to create summary tables.
Loosely based on \href{https://stackoverflow.com/a/25555105/11856430}{this}
stackoverflow post.
}
\details{
Standard zeros are omitted from the data and mean values are calculated
as follows: \code{exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}.
Side effect can be that, when there are many zeros in the data,
the geometric mean could fall outside the confidence interval, since
zeros are not omitted from length(x). To prevent this, you can use a
different way to handle zeros (for example, \code{offset = 1}).
}
\examples{
gm_mean(c(1,4,2,3))
}
