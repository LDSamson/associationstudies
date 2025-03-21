% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_p_histogram}
\alias{plot_p_histogram}
\title{Function to plot P histogram}
\usage{
plot_p_histogram(
  data,
  p_value = "p.value",
  title.to.plot = "P-values of the association study",
  gamma.hat = 1
)
}
\arguments{
\item{data}{dataframe to use}

\item{p_value}{character string. name of column containing p values}

\item{title.to.plot}{title of the plot to produce}

\item{gamma.hat}{gamma hat to use. Standard set to 1. In future versions, this value could
be changed in order to get more accurate detection of associations.}
}
\value{
histogram
}
\description{
If no associations are found, by definition the pvalues should come from
a uniform distribution. This can be checked with the function below.
}
\examples{
# Create random P values, possible result of an association study
pvals_data <- data.frame(
"p.value" = c(runif(1000, 0, .99999), runif(250, 1.0e-11, 1.0e-8))
)
plot_p_histogram(pvals_data)
}
