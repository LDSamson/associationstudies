
#' @export

############ Function to return geomean values with 95% confidence interval: ################

geomean_conf <- function(x, n.digits = 2, na.rm = FALSE, conf.level = 0.95){
  if (sum(is.na(x)) > 0){
    if (na.rm) {x <- na.omit(x)}
    else {return(NA)}
  }
  if (any(x <= 0)) {
    warning("Non-positive values in 'x'")
    return(NA)
  }
  meanvalue <- round(exp(mean(log(x))), digits = n.digits)
  cinterval <- round(exp(t.test(log(x), conf.level = conf.level)$conf.int), digits = n.digits)
  return(paste0(meanvalue, " (", cinterval[1], "-", cinterval[2], ")"))
}

# cleaner display output for geomean_conf:
# https://stackoverflow.com/questions/44325464/how-to-control-knitr-kable-scientific-notation
geomean_conf_2 <- function(x, n.digits = 3, na.rm = FALSE, conf.level = 0.95,
                           scientific.above = 10000){
  if (sum(is.na(x)) > 0){
    if (na.rm) {x <- na.omit(x)}
    else {return(NA)}
  }
  if (any(x <= 0)) {
    warning("Non-positive values in 'x'")
    return(NA)
  }
  meanvalue <- exp(mean(log(x)))
  cinterval <- exp(t.test(log(x), conf.level = conf.level)$conf.int)
  # 'g' in formatC saves space only when necessary:
  mean.formatted <- formatC(meanvalue,  format = "g", digits = n.digits)

  if(str_detect(mean.formatted, ".e")){
    decile.factor <- as.numeric(str_remove(mean.formatted, ".*e"))
  }else{decile.factor <- 1}
  # see link below: force n.digits here:
  # https://stackoverflow.com/a/12135122/11856430
  mean.without.base.ten <- format(
    round(meanvalue/10^decile.factor, digits = n.digits),
    nsmall = n.digits
  )
  interval <- format(
    round(cinterval/(10^decile.factor), digits = n.digits),
    nsmall = n.digits
  )
  interval.formatted <- paste0("(", paste0(interval, collapse = "-"), ")")

  if(decile.factor == 1){# don't display base 10 when format is 10^1:
    output.string <- paste0("$", mean.without.base.ten, interval.formatted, "$")
  }else{
    output.string <- paste0("$", mean.without.base.ten, interval.formatted, "*10^",
                            decile.factor, "$")
  }
  # n.b. dollar signs give math output in kable function
  # (set escape=F in kable function!)
  return(output.string)
}

# other functions for table output:
# mean with range function:
mean.plus.range <- function(x, digits = 1, ...){
  mean.value <- round(mean(x, ...), digits = digits)
  sd.value   <- round(sd(x, ...),   digits = digits)
  min.value  <- round(min(x, ...),  digits = digits)
  max.value  <- round(max(x, ...),  digits = digits)
  value.to.return <- paste0(
    mean.value, " (SD: ", sd.value, ", range: ",
    min.value, " - ", max.value, ")"
  )
  return(value.to.return)
}

# older functions used for tables:
mean.sd <- function(var, n.digits = 2, na.rm = T){
  mean.sd.calculated <-
    paste0(
      round(mean(var, na.rm = na.rm), digits = n.digits),
      " (",
      round(sd(var, na.rm = na.rm), digits = n.digits),
      ")"
    )
  return(mean.sd.calculated)
}

perc.number <- function(var, count.condition, n.digits = 1, na.rm = T){
  paste0(
    round(100*sum(var %in% count.condition, na.rm = T)/n(), digits = n.digits),
    " (",
    sum(var %in% count.condition, na.rm = T),
    ")"
  )
}

median.iqr <- function(var, n.digits = 2, na.rm = T){
  median.iqr.calculated <-
    paste0(
      round(median(var, na.rm = na.rm), digits = n.digits),
      " (",
      round(IQR(var, na.rm = na.rm), digits = n.digits),
      ")"
    )
  return(median.iqr.calculated)
}
