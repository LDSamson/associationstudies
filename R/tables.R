
#' Geometric mean function
#'
#' Calculates the geometric mean.
#' If conf.level is given, a vector will be output containing the geometric mean,
#' the lower bound confidence value, and the upper bound confidence value.
#' Loosely based on \href{https://stackoverflow.com/a/25555105/11856430}{this}
#' stackoverflow post.
#'
#' @param x numeric vector as input
#' @param na.rm logical, whether missing values are allowed.
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
#'
#' @return a vector of length 1 or, if conf.level is not missing, length 3
#' @export
#' @details
#' Standard zeros are omitted from the data and mean values are calculated
#' as follows: \code{exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}.
#' Side effect can be that, when there are many zeros in the data,
#' the geometric mean could fall outside the confidence interval, since
#' zeros are not omitted from length(x). To prevent this, you can use a
#' different way to handle zeros (for example, \code{offset = 1}).

#'
#' @examples
#' gm_mean(c(1,4,2,3))
gm_mean = function(
    x,
    na.rm=TRUE,
    offset = 0,
    zero.propagate = FALSE,
    conf.level = NA
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
  if(!is.na(conf.level)){
    if(all(x == 0)) gconf <- c(NaN, NaN)
      else gconf <- exp(stats::t.test(log(x[x > 0]), conf.level = conf.level)$conf.int)
    gm <- c(gm, gconf)
  }
  gm
}

#' Format values for table input
#'
#' Formats the results of gm_mean so that it can be input into
#' a (latex) table, such as the output of knitr::kable().
#' Aims to give a publication-ready result
#'
#' Inspired by \href{https://stackoverflow.com/questions/44325464/how-to-control-knitr-kable-scientific-notation}{this post}.
#' @param x numeric values that need to be summarized (for now, only geometric mean + confidence interval is supported)
#' @param n.digits digits to preserve in output
#' @param na.rm whether or not to remove NA's
#' @param conf.level confidence interval, standard 0.95
#' @param latex_output
#' Logical, Whether dollar signs should be used around the output string,
#' for latex formatting in knitr::kable. Make sure to set escape=FALSE in
#' knitr::kable function when using this.
#' @param ... other parameters, parsed onto function gm_mean
#' @details
#' Standard zeros are omitted from the data and mean values are calculated
#' as follows: \code{exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}.
#' Side effect can be that, when there are many zeros in the data,
#' the geometric mean could fall outside the confidence interval, since
#' zeros are not omitted from length(x). To prevent this, you can use a
#' different way to handle zeros (for example, \code{offset = 1}).
#'
#' @return Character value to input in kable function
#' @export
#'
#' @examples
#' data <- abs(rnorm(n = 100, mean = 10000, sd = 10000))
#' format_summary_values(data)
#'
#' library(dplyr)
#' library(knitr)
#' library(kableExtra)
#'
#' # create dummy data, with some variables being, just for illustration of this
#' # function:
#' table_data <- immune_data %>%
#'   mutate(across(c(Neutrophils, Monocytes), .fns = ~(.+1)*10^5))
#'
#' summ_table <- immune_data %>%
#'   group_by(Sex) %>%
#'   summarize(
#'     Neutrophils = format_summary_values(Neutrophils),
#'     Tregs = format_summary_values(Tregs)
#'   )
#' knitr::kable(summ_table) %>%
#'   kableExtra::kable_styling()
#'
#' # The function is also helpful to create many summary values at once:
#' library(tidyr)
#' table_data_long <- immune_data %>%
#'   mutate(across(c(Neutrophils, Monocytes), .fns = ~(.+1)*10^5)) %>%
#'   select(-c(Sex, Batch, Frailty.index)) %>%
#'   tidyr::pivot_longer(col = everything())
#'
#' summ_table_long <- table_data_long %>%
#'   group_by(name) %>%
#'   summarize(
#'     "Summary value" = format_summary_values(value, offset = 1, n.digits = 2)
#'   )
#' # offset gives better results when geomean and confidence interval are calculated
#' # and there are many zeros in the data. See details of the function gm_mean()
#' # for more information.
#'
#' knitr::kable(summ_table_long) %>%
#'   kableExtra::kable_styling()
#'
format_summary_values <- function(
    x,
    n.digits = 3,
    na.rm = FALSE,
    conf.level = 0.95,
    latex_output = TRUE,
    ...
    ){
  sumvals <- gm_mean(x, na.rm = na.rm, conf.level = conf.level, ...)
  if(all(is.na(sumvals))) stop("Invalid summary input. Negative numbers included? Check values with function gm_mean()")
  exponent <- floor(log10(sumvals[1]))
  if(exponent <= n.digits) exponent <- 0
  mean_interval <- format(
    round(sumvals/10^exponent, digits = n.digits),
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

#' Function for table output giving mean plus range:
#'
#' @param x numerical value
#' @param digits number of digits to use
#' @param ... other vars, parsed to mean/min/max/sd functions
#'
#' @return character string of mean plus range
#' @export
#'
mean_plus_range <- function(x, digits = 1, ...){
  mean.value <- round(mean(x, ...), digits = digits)
  sd.value   <- round(stats::sd(x, ...),   digits = digits)
  min.value  <- round(min(x, ...),  digits = digits)
  max.value  <- round(max(x, ...),  digits = digits)
  value.to.return <- paste0(
    mean.value, " (SD: ", sd.value, ", range: ",
    min.value, " - ", max.value, ")"
  )
  return(value.to.return)
}

#' Table output helper
#'
#' older function used for nice table output
#'
#' @param var variable name
#' @param n.digits  number of digits to round number
#' @param na.rm exclude NA's or not
#'
#' @return Character string
#' @export
mean_sd <- function(var, n.digits = 2, na.rm = T){
  mean.sd.calculated <-
    paste0(
      round(mean(var, na.rm = na.rm), digits = n.digits),
      " (",
      round(stats::sd(var, na.rm = na.rm), digits = n.digits),
      ")"
    )
  return(mean.sd.calculated)
}

#' Table output helper
#'
#' older function used for nice table output
#'
#' @param var variable name
#' @param count.condition factor variable to group by
#' @param n.digits number of digits to use in output
#' @param na.rm exclude NA's or not
#'
#' @return Character string
#' @export
perc_number <- function(
    var,
    count.condition,
    n.digits = 1,
    na.rm = T
    ){
  paste0(
    round(100*sum(var %in% count.condition, na.rm = T)/dplyr::n(), digits = n.digits),
    " (",
    sum(var %in% count.condition, na.rm = T),
    ")"
  )
}

#' Table output helper
#'
#' older function used for nice table output
#'
#' @param var variable name
#' @param n.digits number of digits to use in output
#' @param na.rm exclude NA's or not
#'
#' @return Character string
#' @export
median_iqr <- function(var, n.digits = 2, na.rm = T){
  median.iqr.calculated <-
    paste0(
      round(stats::median(var, na.rm = na.rm), digits = n.digits),
      " (",
      round(stats::IQR(var, na.rm = na.rm), digits = n.digits),
      ")"
    )
  return(median.iqr.calculated)
}
