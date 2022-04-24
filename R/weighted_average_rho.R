#' Weighted average rho
#'
#' calculates weighted average rho in a blocked correlation study, as discussed
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
    global.rho <-  stats::cor(data[[x]], data[[y]], method = "spearman")
  }else{
    # First, calculate rho per block with base R cor function:
    blocks <- unique(data[[stratum]])
    estimates.of.rho <- data.frame()
    for(b in blocks){
      data.set.b <- data[data[[stratum]] == b,]
      rho.b <-  stats::cor(data.set.b[[x]], data.set.b[[y]], method = "spearman")
      estimates.of.rho <- rbind(estimates.of.rho,
                                data.frame(rho.hat=rho.b, block = b,
                                           sample.size=nrow(data.set.b)))
    }
    if(any(is.na(estimates.of.rho$rho.hat))){
      warning("in some groups rho cannot be estimated. these groups will be ignored when calculating weighted average rho")
      estimates.of.rho <- stats::na.omit(estimates.of.rho)
    }

    # Now, calculate weighted average from the individual rho's per block:
    global.rho <- sum(estimates.of.rho$rho.hat*estimates.of.rho$sample.size)/
      sum(estimates.of.rho$sample.size)
  }
  return(global.rho)
}
