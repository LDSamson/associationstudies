#' Function to plot P histogram
#'
#' @param data dataframe to use
#' @param p_value character string. name of column containing p values
#' @param title.to.plot title of the plot to produce
#' @param gamma.hat gamma hat to use
#' @importFrom rlang .data
#'
#' @return histogram
#' @export
#'
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
