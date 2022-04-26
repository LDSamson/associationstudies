#' Function to plot P histogram
#'
#' If no associations are found, by definition the pvalues should come from
#' a uniform distribution. This can be checked with the function below.
#'
#' @param data dataframe to use
#' @param p_value character string. name of column containing p values
#' @param title.to.plot title of the plot to produce
#' @param gamma.hat
#' gamma hat to use. Standard set to 1. In future versions, this value could
#' be changed in order to get more accurate detection of associations.
#' @importFrom rlang .data
#'
#' @return histogram
#' @export
#'
#' @examples
#' # Create random P values, possible result of an association study
#' pvals_data <- data.frame(
#' "p.value" = c(runif(1000, 0, .99999), runif(250, 1.0e-11, 1.0e-8))
#' )
#' plot_p_histogram(pvals_data)

plot_p_histogram <- function(
    data,
    p_value = "p.value",
    title.to.plot = "P-values of the association study",
    gamma.hat = 1
){
  p.values.hist <- ggplot2::ggplot(data, ggplot2::aes(x = !!rlang::sym(p_value), y = .data$..density..)) +
    ggplot2::geom_histogram(binwidth = 0.1, fill = 'lightblue', col = 'black') +
    ggplot2::geom_hline(yintercept = gamma.hat, col = 'tomato') +
    ggplot2::theme_light() +
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::labs(x = "P value", title = title.to.plot)
  print(p.values.hist)
}
