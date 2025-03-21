% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{unite_vars}
\alias{unite_vars}
\title{Unite values}
\usage{
unite_vars(x, vars, colname = "block", ...)
}
\arguments{
\item{x}{a data frame}

\item{vars}{character vector, with values to unite}

\item{colname}{name new column}

\item{...}{parsed to underlying function tidyr::unite}
}
\value{
data frame where multiple columns are united
}
\description{
This function unites column variable and makes it a factor,
for easy use in other functions.
}
