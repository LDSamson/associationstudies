
#' Geometric mean function
#'
#' Calculates the geometric mean.
#' If conf.level is given, a vector will be output containing the geometric mean,
#' the lower bound confidence value, and the upper bound confidence value.
#' Loosely based on https://stackoverflow.com/a/25555105/11856430
#'
#' @param x numeric vector as input
#' @param na.rm logical, whether missing values are allowed.
#' @param offset Optional variable: adds a value to the data (1 is often used), so that no more zeros occur.
#' @param zero.propagate if TRUE, the output will be zero if any value is zero. Overrides offset.
#' @param conf.level If not missing, the confidence interval will be given together with the geometric mean
#'
#' @return a vector of length 1 or, if conf.level is not missing, length 3
#' @export
#'
#' @examples
#' gm_mean(c(1,4,2,3))
gm_mean = function(x, na.rm=TRUE, offset = 0, zero.propagate = FALSE, conf.level = NA){
  if(any(x < 0, na.rm = TRUE)){warning("Non-positive values in 'x'"); return(NaN)}
  if(zero.propagate && any(x == 0, na.rm = TRUE)) return(0)
  if (sum(is.na(x)) > 0){
    if (na.rm) {x <- stats::na.omit(x)}
      else {return(NA)}
  }
  if(offset > 0 ) x <- x + offset
  gm <- exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  if(!is.na(conf.level)){
    if(all(x == 0)) gconf <- c(NaN, NaN)
      else gconf <- exp(stats::t.test(log(x[x > 0]), conf.level = conf.level)$conf.int)
    gm <- c(gm, gconf)
  }
  gm
}

#' Geomean + Conf interval
#'
#' Improved version of geomean_conf for latex tables:
#' Inspired by this post: # https://stackoverflow.com/questions/44325464/how-to-control-knitr-kable-scientific-notation
#' n.b. dollar signs give math output in kable function
#' Make sure to set escape=F in kable function!
#'
#' @param x numeric value
#' @param n.digits digits to preserve in output
#' @param na.rm whether or not to remove NA's
#' @param conf.level confidence interval, standard 0.95
#' @param scientific.above threshold above which scientific notation should be used
#'
#' @return Character value to input in kable function
#' @export
format_summary_values <- function(x, n.digits = 3, na.rm = FALSE, conf.level = 0.95,
                           scientific.above = 10000){
  sumvals <- gm_mean(x, na.rm = na.rm, conf.level = conf.level)
  # 'g' in formatC saves space only when necessary:
  mean.formatted <- formatC(sumvals[1],  format = "g", digits = n.digits)
  if(stringr::str_detect(mean.formatted, ".e")){
    exponent <- as.numeric(stringr::str_remove(mean.formatted, ".*e"))
  }else{exponent <- 0}
  # see link below: force n.digits here:
  # https://stackoverflow.com/a/12135122/11856430
  mean.without.base.ten <- format(
    round(sumvals[1]/10^exponent, digits = n.digits),
    nsmall = n.digits
  )
  interval <- format(
    round(sumvals[2:3]/(10^exponent), digits = n.digits),
    nsmall = n.digits
  )
  interval.formatted <- paste0("(", paste0(interval, collapse = "-"), ")")

  if(exponent == 0){# don't display base 10 when format is 10^1:
    output.string <- paste0("$", mean.without.base.ten, interval.formatted, "$")
  }else{
    output.string <- paste0("$", mean.without.base.ten, interval.formatted, "*10^",
                            exponent, "$")
  }
  # n.b. dollar signs give math output in kable function
  # (set escape=F in kable function!)
  return(output.string)
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
perc_number <- function(var, count.condition, n.digits = 1, na.rm = T){
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
