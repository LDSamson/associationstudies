
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
    message(paste0("Warning: The explanatory variable '",
                   explanatory.var,
                   "' is of type character which cannot be used in analyses.",
                   "\nIt will be converted to a factor and used in a wilcox_test"))
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
    "response.var"    = response.var,
    "explanatory.var" = explanatory.var,
    "sample.size"     = nrow(data.to.analyze),
    "stratified.by"   = paste(stratum, collapse = "_"),
    "n.resample"      = n.resample,
    "method"          = test.results@method,
    "direction"       = sign(coin::statistic(test.results)),
    "rho"             = NA_integer_,
    "p.value"         = as.numeric(coin::pvalue(test.results))
    )
  ## Add rho to outcome, if applicable:
  if(results.to.return[["method"]] == "Spearman Correlation Test"){
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
#' @param expl.var.names character value of column that contains all the names of the response values that need to be investigated separately
#' @param ... other values will be parsed to \code{\link{test_association}}
#'
#' @return data frame with results as output
#' @export
#'
association_study_long <- function(data, expl.var.names, ...){
  purrr::map_dfr(unique(data[[expl.var.names]]), .f = function(x){
    df <- test_association(data[data[[expl.var.names]] == x, ], ...)
    df["response.var"] <- x
    df
  })
}

#' Perform association study
#'
#' This function is a wrapper around the function \code{\link{test_association}}.
#' This function tests associations between multiple variables by repeatedly
#' calling the function \code{\link{test_association}}.
#' Standard, all possible single associations between pairs of variables in the
#' given data frame are tested.
#'
#' The response variables are columns in a data frame and their
#' column names should be given in a character string.
#'
#' Note 2: Associations can be stratified (see options for
#' \code{\link{test_association}}). All strata/blocks should contain enough
#' observations.
#'
#' @param data data frame to use
#' @param response.var
#' Character value, optional. Response variable. If provided, associations will
#' only be tested for this single response variable.
#' @param vars.to.select
#' Variables to select. Will be used within dplyr::select.
#' All possible pairs of associations between these variables will be tested.
#' Standard all variables will be selected.
#' @param stratum
#' character vector. When used, the tests will be blocked by these variables.
#' @param ... other parameters that will be parsed to \code{\link{test_association}}.
#'
#' @return data frame with results as output
#' @export
#'
#' @examples
#' library(dplyr)
#' cols_to_analyze <- immune_data %>%
#' select(-c(Batch, Sex, Frailty.index)) %>% names()
#' test_outcome <- association_study(immune_data,
#' "Frailty.index", cols_to_analyze)
#'
association_study <- function(
    data,
    response.var = NULL,
    vars.to.select = tidyselect::everything(),
    stratum = NULL,
    ...
    ){
  if(!is.null(response.var)){
  }
  if(!"data.frame" %in% class(data)) stop("data should be of type data.frame")

  data <- data %>%
    # select 'vars.to.select' first so that also negative selection works:
    dplyr::select({{vars.to.select}},
                  tidyselect::any_of(c(response.var, stratum))) %>%
    # (not essential) now just change variable order:
    dplyr::select(tidyselect::any_of(c(response.var, stratum)),
                  tidyselect::everything())
  expl.var.names <- names(data)
  if(!is.null(stratum)){expl.var.names <- names(data)[!names(data) %in% stratum]}

  if(!is.null(response.var)){
    expl.var.names <- expl.var.names[expl.var.names != response.var]
    test_results <- purrr::map_dfr(expl.var.names, .f = function(x){
      test_association(data, response.var = response.var,  explanatory.var = x,
                       stratum = stratum, ...)
      })
   return(test_results)
  }
  pairs_to_test <- variable_pair_combinations(data, expl.var.names)
  # since at the moment only numerical response variables are supported,
  # below shifts non-numerical variables from response variable
  # to the explanatory variable:
  non_numerical <- unlist(lapply(pairs_to_test[,1],
                                 FUN = function(x){!is.numeric(data[[x]])}))
  pairs_to_test[non_numerical, ] <- pairs_to_test[non_numerical, 2:1]

  test_results <- do.call(
    "rbind",
    apply(
      pairs_to_test,
      MARGIN = 1,
      FUN = function(x){
        test_association(
          dataset = data,
          response.var = x[1],
          explanatory.var = x[2],
          stratum = stratum,
          ...
        )
      }
    )
  )
  test_results
}
