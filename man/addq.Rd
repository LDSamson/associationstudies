% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{addq}
\alias{addq}
\title{Add backticks}
\usage{
addq(x)
}
\arguments{
\item{x}{character value that needs backticks}
}
\value{
character value with backticks
}
\description{
Small helper function to add backticks around complicated names
Inspired by https://stackoverflow.com/questions/16674045/as-formula-in-r-doesnt-seem-to-accept-a-name-that-starts-with-a-number-followed
}
\examples{
addq("Long Column Name")

}
