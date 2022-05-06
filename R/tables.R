
#' Format values for table input
#'
#' Formats the results of a summary value so that it can be input into
#' a (latex) table, such as the output of knitr::kable().
#' Aims to give a publication-ready result.
#' Expected input is a numeric vector, with the summary statistic value on
#' the first place, and on the second/third place the standard deviation or
#' confidence interval. It is designed to work with the output of functions
#' such as \code{\link{gm_mean}}, \code{\link{median_iqr}} and
#' \code{\link{mean_sd}}.
#'
#' @param x numeric values that need to be summarized (for now, only geometric mean + confidence interval is supported)
#' @param n.digits digits to preserve in output
#' @param latex_output
#' Logical, Whether dollar signs should be used around the output string,
#' for latex formatting in knitr::kable. Make sure to set escape=FALSE in
#' knitr::kable function when using this.
#' @param na.rm logical, indicating whether NA values should be excluded
#' @param ... other parameters, parsed onto function gm_mean
#' @details
#' Inspired by \href{https://stackoverflow.com/questions/44325464/how-to-control-knitr-kable-scientific-notation}{this post}.
#'
#' @return Character value to input in kable function
#' @export
#'
#' @examples
#' df <- data.frame(
#'   V1 = abs(rnorm(n = 1000, mean = 10000, sd = 10000)),
#'   V2 = exp(rnorm(n = 1000, mean = log(10^6), sd = log(10^4))),
#'   V3 = abs(rnorm(n = 1000, mean = 100, sd = 10))
#' )
#'
#' #V1: median (IQR)
#' format_values(median_iqr(df$V1), n.digits = 1)
#' #V2: geomean (conf)
#' format_values(gm_mean(df$V2), n.digits = 1)
#' #V3: mean (SD)
#' format_values(mean_sd(df$V3), n.digits = 1)
#'
#' ## illustrate how to use summary tables:
#' library(dplyr)
#' library(knitr)
#' library(kableExtra)
#'
#' # create dummy data, with some variables having much higher average values
#' # than others, to illustrate the use of this function:
#' table_data <- immune_data %>%
#'   mutate(across(c(Neutrophils, Monocytes), .fns = ~(.+1)*10^5))
#'
#' summ_table <- immune_data %>%
#'   group_by(Sex) %>%
#'   summarize(
#'     Neutrophils = format_values(gm_mean(Neutrophils)),
#'     Tregs = format_values(gm_mean(Tregs))
#'   )
#' knitr::kable(summ_table) %>%
#'   kableExtra::kable_styling()
#'
#' # The function is also helpful to create many summary values at once:
#' library(tidyr)
#' table_data_long <- table_data %>%
#'   select(-c(Sex, Batch, Frailty.index)) %>%
#'   tidyr::pivot_longer(col = everything())
#'
#' summ_table_long <- table_data_long %>%
#'   group_by(name) %>%
#'   summarize(
#'     "Summary value" = format_values(gm_mean(value, offset = 1), n.digits = 2)
#'   )
#' # offset gives better results when geomean and confidence interval are calculated
#' # and there are many zeros in the data. See details of the function gm_mean()
#' # for more information.
#'
#' knitr::kable(summ_table_long) %>%
#'   kableExtra::kable_styling()
#'
format_values <- function(
    x,
    n.digits = 3,
    latex_output = TRUE,
    na.rm = FALSE,
    ...
    ){
  if(length(x)>3) stop("Only numerical vectors of length <= 3 are supported as input")
  if(all(is.na(x))) stop("Invalid summary input.")
  exponent <- floor(log10(abs(x[1])))
  if(exponent <= n.digits) exponent <- 0
  mean_interval <- format(
    round(x/10^exponent, digits = n.digits),
    nsmall = n.digits # see: https://stackoverflow.com/a/12135122/11856430
  )
  output.string <- paste0(
    mean_interval[1],
    "(",
    paste0(mean_interval[-1], collapse = "-"),
    ")"
    )
  if(exponent != 0) output.string <- paste0(c(output.string, "*10^", exponent),
                                            collapse = "")
  if(latex_output) output.string <- paste0("$", output.string, "$")
  output.string
}

#' Geometric mean
#'
#' Calculates the geometric mean with confidence interval.
#' If conf.level is given, a vector will be output containing the geometric mean,
#' the lower bound confidence value, and the upper bound confidence value.
#' Output can be used directly in the function
#' \code{\link{format_values}} to create summary tables.
#' Loosely based on \href{https://stackoverflow.com/a/25555105/11856430}{this}
#' stackoverflow post.
#'
#' @param x numeric vector as input
#' @param offset
#' Optional variable: adds a value to the data (1 is often used),
#' so that no more zeros occur in the data. This is an acceptable way of
#' handling zeros when, for example, zeros represent undetectable low
#' concentrations in a biomedical laboratory assay.
#' @param zero.propagate
#' if TRUE, the output will be zero if any value is zero.
#' Overrides the offset parameter is true.
#' @param conf.level
#' If not missing, the confidence interval will be given together with the geometric mean
#' @param na.rm logical, indicating whether NA values should be excluded
#'
#' @return a vector of length 1 or, if conf.level is not missing, length 3
#' @details
#' Standard zeros are omitted from the data and mean values are calculated
#' as follows: \code{exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}.
#' Side effect can be that, when there are many zeros in the data,
#' the geometric mean could fall outside the confidence interval, since
#' zeros are not omitted from length(x). To prevent this, you can use a
#' different way to handle zeros (for example, \code{offset = 1}).
#' @export
#'
#' @examples
#' gm_mean(c(1,4,2,3))
gm_mean = function(
    x,
    offset = 0,
    zero.propagate = FALSE,
    conf.level = 0.95,
    na.rm=TRUE
){
  if(any(x < 0, na.rm = TRUE)){warning("Non-positive values in 'x'"); return(NaN)}
  if(zero.propagate && any(x == 0, na.rm = TRUE)) return(0)
  #! Unlike in the example from stackoverflow,
  # this also removes NA's from the length call below, which I think is correct
  if (sum(is.na(x)) > 0){
    if(all(is.na(x))) return(NaN)
    if (na.rm) {x <- stats::na.omit(x)}
    else {return(NA)}
  }
  if(all(x == 0)) return(0)
  if(offset > 0 ) x <- x + offset
  gm <- exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  if(!is.null(conf.level)){
    if(all(x == 0)) gconf <- c(NaN, NaN)
    else gconf <- exp(stats::t.test(log(x[x > 0]), conf.level = conf.level)$conf.int)
    gm <- c(gm, gconf)
  }
  gm
}


#' Mean plus standard deviation
#'
#' A vector will be output containing the arithmetic mean and the standard
#' deviation. The output can be used directly in the function
#' \code{\link{format_values}} to create summary tables.
#'
#' @param var variable name
#' @param na.rm logical, indicating whether NA values should be excluded
#'
#' @return Character string
#' @export
#'
#' @examples
#' mean_sd(rnorm(100, 1000, 10))
mean_sd <- function(
    var,
    na.rm = T
    ){
  c(mean(var, na.rm = na.rm), stats::sd(var, na.rm = na.rm))
}


#' Median + interquartile range
#'
#' Creates a vector with the mean and interquartile range.
#' Output can be used directly in the function
#' \code{\link{format_values}} to create summary tables.
#'
#' @param var numerical vector
#' @param na.rm logical, indicating whether NA values should be excluded
#'
#' @return Character string
#' @export
#'
#' @examples
#' median_iqr(rnorm(100, 1000, 100))
median_iqr <- function(
    var,
    na.rm = TRUE
    ){
  c(stats::median(var, na.rm = na.rm), stats::IQR(var, na.rm = na.rm))
}

#' Percentage and number in ordinal variable.
#'
#' Helper function to create summary statistics of ordinal or categorical
#' variables, for use in a summary table.
#' It will display the occurrence of a condition within an ordinal or
#' categorical variable. Works best with two groups.
#'
#' @param var character of factor vector
#' @param count.condition character string. Variable to group by
#' @param n.digits number of digits to use in output
#' @param na.rm logical, indicating whether NA values should be excluded
#'
#' @return Character string
#' @export
#'
#' @examples
#' groupvar <-  sample(c("Men", "Women"), 80, TRUE)
#' perc_number(groupvar, count.condition = "Men")
perc_number <- function(
    var,
    count.condition,
    n.digits = 1,
    na.rm = TRUE
){
  paste0(
    round(100*sum(var %in% count.condition, na.rm = T)/length(var), digits = n.digits),
    " (",
    sum(var %in% count.condition, na.rm = T),
    ")"
  )
}

#' Mean, SD and range
#'
#' Function for table output giving mean, standard deviation,
#' plus range (minimum value and maximum value). Can be useful
#' in summary tables.
#'
#' @param x numerical value
#' @param digits number of digits to use
#' @param ... other vars, parsed to mean/min/max/sd functions
#'
#' @return character string of mean plus range
#' @export
#'
#' @examples
#' mean_sd_range(rnorm(1000, 1000, 100))
mean_sd_range <- function(
    x,
    digits = 1,
    ...
){
  sumvals <- round(c(
    mean(x, ...),
    stats::sd(x, ...),
    min(x, ...),
    max(x, ...)
  ), digits = digits)
  paste0(
    sumvals[1], " (SD: ", sumvals[2], ", range: ",
    sumvals[3], " - ", sumvals[4], ")"
  )
}
