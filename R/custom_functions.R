

####################### colorblind friendly color palette: #####################

# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# For scientific notation:
# https://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales

#' Scientific notation
#'
#' Create scientific notation in ggplot2 functions
#'
#' @param x any numeric value that should have scientific notation
#'
#' @return value in scientific notation
#' @export
scientific_10 <- function(x) {
  parse(text=gsub("\\+", "", # addition to remove the + sign
                  gsub("e", " %*% 10^", scales::scientific_format()(x))
                  )
        )
}

#' Add pagebreak
#'
#' Function to add pagebreaks in Rmarkdown files.
#' For more info see: https://stackoverflow.com/a/55064070
#'
#' @return either a html or a latex tag that adds a page break
#' @export
#'
pagebreak <- function() {
  if(knitr::is_latex_output())
    return("\\newpage")
  else
    return('<div style="page-break-before: always;" />')
}

#' Weigthed average rho
#'
#' calculates weighted average rho in blocked correlation study, as discussed
#' with Jose Ferreira (RIVM)
#'
#' @param data expects a data frame
#' @param x Character value, first variable. Should be in data
#' @param y Character value, second variable. Should be in data
#' @param stratum stratum (factor) to block the association with
#'
#' @return numeric value (weighted average rho)
#' @export
#'
#' @examples
#' weighted_average_rho(iris, "Petal.Length", "Petal.Width", "Species")
#'
weighted_average_rho <-function(data, x, y, stratum = "block"){
  if(is.null(stratum)){
   # estimates.of.rho <- data.frame()
    global.rho <-  stats::cor(data[[x]],
                       data[[y]],
                       method = "spearman")
  }else{
    # 06-11-2019 LS:
    # Adding rho as weighted average in the results,
    # following the advise from José Ferreira
    # First, calculate rho per block with base R cor function:
    blocks <- unique(data[[stratum]])
    estimates.of.rho <- data.frame()
    for(b in blocks){
      # b <- "70-75yr_CMV-"
      data.set.b <- data[data[[stratum]] == b,]
      rho.b <-  stats::cor(data.set.b[[x]], data.set.b[[y]], method = "spearman")
      estimates.of.rho <- rbind(estimates.of.rho,
                                data.frame(rho.hat=rho.b, block = b,
                                           sample.size=nrow(data.set.b)))
    }
    # Now, calculate weighted average from the individual rho's per block:
    global.rho <- sum(estimates.of.rho$rho.hat*estimates.of.rho$sample.size)/
      sum(estimates.of.rho$sample.size)
  }
  return(global.rho)
}


#' Add backticks
#'
#' Small helper function to add backticks around complicated names
#' Inspired by https://stackoverflow.com/questions/16674045/as-formula-in-r-doesnt-seem-to-accept-a-name-that-starts-with-a-number-followed
#'
#' @param x character value that needs backticks
#'
#' @return character value with backticks
#' @export
#'
#' @examples addq("Long Column Name")
addq <- function(x){paste0("`", x, "`")}



#' Unite values
#'
#' This function unites column variable and makes it a factor,
#' for easy use in other functions.
#'
#' @param x a data frame
#' @param vars character vector, with values to unite
#' @param colname name of the newly created column
#' @param ...  parsed to underlying function tidyr::unite
#'
#' @return  data frame where multiple columns are united
#' @export
#'
unite_vars <- function(x, vars, colname = "block", ...){
  data.to.return <- x %>%
    tidyr::unite(col = {{colname}}, tidyselect::all_of(vars), ...)
  data.to.return[[colname]] <- factor(data.to.return[[colname]])
  return(data.to.return)
}


#' Function to plot P histogram
#'
#' @param data dataframe to use
#' @param p_value character string. name of column containing p values
#' @param title.to.plot title of the plot to produce
#' @param gamma.hat gamma hat to use
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
  p.values.hist <- ggplot2::ggplot(data, ggplot2::aes(x = !!rlang::sym(p_value), y = ..density..)) +
    ggplot2::geom_histogram(binwidth = 0.1, fill = 'lightblue', col = 'black') +
    ggplot2::geom_hline(yintercept = gamma.hat, col = 'tomato') +
    ggplot2::theme_light() +
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::labs(x = "P value", title = title.to.plot)
  print(p.values.hist)
}
