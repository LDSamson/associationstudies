
#' Test association
#'
#' Function that can be used to test an association between two variables.
#' It uses permutation testing functions from the coin package.
#' It gives a consistent data frame output.
#'
#' Note: Errors can occur when some groups within the blocking variables
#' contain less than two observations.
#' Note2: missing values are (silently) deleted.
#'
#' @param response.var response variable (for example Granulocytes_count)
#' @param dataset dataset (data frame format, original: data.wide)
#' @param explanatory.var explanatory variable (for example Frailty.index)
#' @param stratum blocking variables (character string)
#' @param n.resample number of resamples
#' @param ... other variables that influence the coin::approximate function
#' @importFrom rlang .data
#'
#' @return data frame with results
#' @export
#'
#' @examples
#' test_association(immune_data, "Tregs", "Frailty.index", "Batch")
test_association <- function(
    dataset,
    response.var,
    explanatory.var,
    stratum         = NULL,
    n.resample      = 10^4,
    ...){
  data.to.analyze <- dataset %>%
    dplyr::select(tidyselect::all_of(c(response.var, explanatory.var, stratum))) %>%
    stats::na.omit() %>%
    droplevels()
  check_low_group_numbers(data.to.analyze, stratum)
  string_formula <- paste0(addq(response.var), "~", addq(explanatory.var))
  if(!is.null(stratum)){
    data.to.analyze <- unite_vars(data.to.analyze, vars = stratum, colname = "block")
    string_formula <- paste0(string_formula, "|block")
    }
  formula.to.test <- stats::as.formula(string_formula)

  if(!is.numeric(data.to.analyze[[response.var]])){stop("At the moment, only numerical response variables are supported.")}

  if(is.character(data.to.analyze[[explanatory.var]])){
    message("Warning: The explanatory variable is of type character which cannot be used in analyses.\nIt will be converted to a factor and used in a wilcox_test")
    data.to.analyze[[explanatory.var]] <- factor(data.to.analyze[[explanatory.var]])
  }

  if(is.factor(data.to.analyze[[explanatory.var]])){
    ## test results when comparing a continuous variable between groups/factor:
    test.results <- coin::wilcox_test(formula.to.test, data = data.to.analyze,
                                      distribution=coin::approximate(
                                        nresample=n.resample, ...
                                        )
                                      )
  }else if(is.numeric(data.to.analyze[[explanatory.var]])){
    ## test results when testing two continuous variables:
    test.results <- coin::spearman_test(formula.to.test, data = data.to.analyze,
                                        distribution=coin::approximate(
                                          nresample=n.resample, ...
                                          )
                                        )
    global.rho <- weighted_average_rho(data = data.to.analyze, x = response.var,
                                       y = explanatory.var,
                                       stratum = if(is.null(stratum)) NULL else
                                         "block")
  }
  results.to.return <- data.frame(
    "Response var"    = response.var,
    "Explanatory var" = explanatory.var,
    "sample.size"     = nrow(data.to.analyze),
    "stratified.by"   = paste(stratum, collapse = "_"),
    "n.resample"      = n.resample,
    "Method"          = test.results@method,
    "Direction"       = sign(coin::statistic(test.results)),
    "rho"             = NA_integer_,
    "p.value"         = as.numeric(coin::pvalue(test.results))
    )
  ## Add rho to outcome, if applicable:
  if(results.to.return[["Method"]] == "Spearman Correlation Test"){
    results.to.return[["rho"]] <- global.rho}
  ## Return results:
  return(results.to.return)
}



#' Perform association study
#'
#' This function, in combination with \code{\link{test_association}},
#' can basically replace the function '\code{\link{association_study}}'.
#' I think \code{\link{test_association}} is more intuitive,
#' and this function is just a small wrapper
#' that you can use when data is in long format, mimicking behavior
#' of \code{\link{association_study}}.
#'
#' Note: Response variables should be of the same type
#' (e.g. all numerical, all categorical).
#' Note 2: all blocks should contain enough observations.
#'
#' @param data data frame to use
#' @param expl_var_names character value of column that contains all the names of the response values that need to be investigated separately
#' @param ... other values will be parsed to \code{\link{test_association}}
#'
#' @return data frame with results as output
#' @export
#'
association_study_long <- function(data, expl_var_names, ...){
  purrr::map_dfr(unique(data[[expl_var_names]]), .f = function(x){
    df <- test_association(data[data[[expl_var_names]] == x, ], ...)
    df["Response.var"] <- x
    df
  })
}
