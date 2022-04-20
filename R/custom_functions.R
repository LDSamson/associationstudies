

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
