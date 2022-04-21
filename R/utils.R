#' Scientific notation
#'
#' Create scientific notation to use in ggplot2 functions. Inspired by
#' \url{https://stackoverflow.com/a/18530540/11856430}
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
#' For more info see \url{https://stackoverflow.com/a/55064070}
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
#' @examples
#' addq("Long Column Name")
#'
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


#' Create all combinations of variable names
#'
#' Useful in a an apply loop when associations should be tested
#' between many variables in a data frame.
#'
#' @param data a data frame as input
#' @param sel parsed to dplyr::select(). To include/exclude variables.
#'
#' @return a character value matrix with two columns, containing all possible pairs of variable names
#' @export
#'
#' @examples
#' variable_pair_combinations(immune_data, sel = -c(Sex, Frailty.index, Batch))
variable_pair_combinations <- function(data, sel = tidyselect::everything()){
  selected.vars <- names(dplyr::select(data, {{sel}}))
  # Create all possible pairs of columns:
  t(utils::combn(selected.vars, m=2))
}
