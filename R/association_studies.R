## 08-01-2020 LS:
##    - Removed argument CMV.status != "Missing"
##        (this gives an error when CMV.status is missing)
##    - Added calculation of rho when no blocking variable is used

#' Association studies
#'
#' Function to perform associations between two parameters
#' Associations with p values calculated by permutation, using coin package
#' Note: dataset is expected to be in "long" format for this function.
#' Note2: missing values are (silently) deleted
#' Use the function in an apply loop (or map_dfr loop) to perform it on
#' many different subsets in the data.
#'
#' @param subset.parameter character value that can be found in the subset.column (original: Granulocytes_Count)
#' @param subset.column column containing subset names (original: Cellname)
#' @param dataset data frame to use (original: data.long)
#' @param numerical.value a numerical value for the comparison
#' @param comparison.value the second parameter. If it is a factor, the wilcoxon test will be used. If it is numerical, the spearman_test will be used
#' @param stratum character string. Can be one or several blocking variables.
#' @param n.resample number of random samples for permutation testing
#' @param ... Other parameters that influence the function coin::approximate
#'
#' @return data frame with results
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' cell_vars <- cell_vars <- names(immune_data)[!names(immune_data) %in%
#' c("Sex", "Frailty.index", "Batch")]
#' immune_data_long <- immune_data %>%
#' tidyr::pivot_longer(names_to = "cellnames", values_to = "value", cell_vars)
#'
#' purrr::map_dfr(cell_vars, .f = ~association_study(
#' dataset = immune_data_long,
#' subset.parameter = .,
#' subset.column = "cellnames",
#' numerical.value = "value",
#' comparison.value = "Frailty.index",
#' stratum = c("Sex", "Batch")
#' ))
#'
association_study <- function(subset.parameter,
                              subset.column,
                              dataset,
                              numerical.value,
                              comparison.value,
                              stratum          = NULL,
                              n.resample       = 10^4,
                              ...){
  data.to.analyze <- dataset %>%
    dplyr::select(tidyselect::all_of(c(subset.column, numerical.value,
                                       comparison.value, stratum))) %>%
    dplyr::filter(!!rlang::sym(subset.column) == subset.parameter) %>%
    stats::na.omit() %>%
    droplevels()

  string_formula <- paste0(addq(numerical.value), "~", addq(comparison.value))
  if(!is.null(stratum)){
    data.to.analyze <- unite_vars(data.to.analyze, vars = stratum, colname = "block")
    string_formula <- paste0(string_formula, "|block")
  }
  formula.to.test <- stats::as.formula(string_formula)

  if(is.factor(dataset[[comparison.value]])){
    ## test results when comparing a continuous variable between groups/factor:
    test.results <- coin::wilcox_test(formula.to.test,
                                      data = data.to.analyze,
                                      distribution=coin::approximate(nresample=n.resample,
                                                                     ...)
    )
  }else if(is.numeric(dataset[[comparison.value]])){
    ## test results when testing two continuous variables:
    test.results <- coin::spearman_test(
      formula.to.test,
      data = data.to.analyze,
      distribution=coin::approximate(nresample=n.resample, ...)
    )
    global.rho <- weighted_average_rho(
      data = data.to.analyze,
      x = numerical.value,
      y = comparison.value,
      stratum = if(is.null(stratum)){NULL}else "block"
    )
  }

  results.to.return <- list(
    "Subsets"          = subset.parameter,
    "Comparison value" = comparison.value,
    "sample.size"      = nrow(data.to.analyze),
    "stratified.by"    = paste(stratum, collapse = "_"),
    "n.resample"       = n.resample,
    "Method"           = test.results@method,
    "Direction"        = sign(coin::statistic(test.results)),
    "p.value"          = as.numeric(coin::pvalue(test.results))
  )
  ## Add rho to outcome, if applicable:
  if(results.to.return[["Method"]] == "Spearman Correlation Test"){
    results.to.return[["rho"]] <- global.rho}
  return(results.to.return)
}

#' Test association
#'
#' Function that can be used to test an association between two variables.
#' It uses permutation testing functions from the coin package.
#' It gives a consistent data frame output.
#'
#' Note: Errors can occur when some groups within the blocking variables
#' contain less than two observations.
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
test_association <- function(
    dataset,
    response.var,
    explanatory.var,
    stratum         = NULL,
    n.resample      = 10^4,
    ...){
  data.to.analyze <- dataset %>%
    dplyr::select(response.var, explanatory.var, stratum) %>%
    stats::na.omit() %>%
    droplevels()
  string_formula <- paste0(addq(response.var), "~", addq(explanatory.var))
  if(!is.null(stratum)){
    data.to.analyze <- unite_vars(data.to.analyze, vars = stratum, colname = "block")
    string_formula <- paste0(string_formula, "|block")
  }
  formula.to.test <- stats::as.formula(string_formula)

  if(is.factor(dataset[[explanatory.var]])){
    ## test results when comparing a continuous variable between groups/factor:
    test.results <- coin::wilcox_test(formula.to.test, data = data.to.analyze,
                                      distribution=coin::approximate(
                                        nresample=n.resample, ...
                                      )
    )
  }else if(is.numeric(dataset[[explanatory.var]])){
    ## test results when testing two continuous variables:
    test.results <- coin::spearman_test(formula.to.test, data = data.to.analyze,
                                        distribution=coin::approximate(
                                          nresample=n.resample, ...
                                        )
    )
    global.rho <- weighted_average_rho(data = data.to.analyze, x = response.var,
                                       y = explanatory.var, stratum = "block")
  }
  results.to.return <- list(
    "Response var"    = response.var,
    "Explanatory var" = explanatory.var,
    "sample.size"     = nrow(data.to.analyze),
    "stratified.by"   = paste(stratum, collapse = "_"),
    "n.resample"      = n.resample,
    "Method"          = test.results@method,
    "Direction"       = sign(coin::statistic(test.results)),
    "p.value"         = as.numeric(coin::pvalue(test.results))
    )
  ## Add rho to outcome, if applicable:
  if(results.to.return[["Method"]] == "Spearman Correlation Test"){
    results.to.return[["rho"]] <- global.rho}
  ## Return results:
  return(results.to.return)
}

#' Associations per individual
#'
#' This function will investigate the association between biomarkers over time
#' (longitudinal analyses) with ID/individual as blocking variable.
#'
#' @param x Subsets in which comparison should be made
#' @param dataset a data frame
#' @param comparison.var comparison var (for example Age)
#' @param nresample times to resample for permutation test (see coin::approximate)
#' @importFrom rlang .data
#'
#' @return data frame with results
#' @export
associations_per_ID <- function(
    x,
    dataset,
    comparison.var,
    nresample = 10^4
    ){
  form.to.test <- stats::as.formula(paste0(comparison.var, "~", x, "|ID"))
  data.to.analyze <- dataset %>%
    dplyr::select(c("ID", "DoetRound", x, comparison.var, "Sex")) %>%
    stats::na.omit() %>%
   dplyr::mutate(ID = factor(.data$ID))

  ## Remove people that have less than two complete observations:
  IDs.with.2.or.more.complete.obs <- data.to.analyze %>%
    dplyr::group_by(.data$ID) %>%
    dplyr::summarize(n.obs = dplyr::n()) %>%
    dplyr::filter(.data$n.obs >=2) %>%
    dplyr::pull(.data$ID)
  data.to.analyze <- data.to.analyze %>%
    dplyr::filter(.data$ID %in% IDs.with.2.or.more.complete.obs)

  test.results <- coin::spearman_test(
    form.to.test,
    data = data.to.analyze,
    distribution = coin::approximate(nresample = nresample)
  )
  #### calculate rho: ####
  global.rho <- weighted_average_rho(data = data.to.analyze, x = x,
                                     y = comparison.var, stratum = "ID")
  p.val <- as.numeric(coin::pvalue(test.results))
  Direction <- sign(coin::statistic(test.results))
  results.to.return <- list(
    "Subsets"          = x,
    "Comparison value" = comparison.var,
    "sample.size"      = length(unique(data.to.analyze[["ID"]])),
    "n.resample"       = nresample,
    "Method"           = test.results@method,
    "Direction"        = Direction,
    "rho"              = global.rho,
    "p.value"          = p.val
  )
  return(results.to.return)
}

# single.pair <- c("CXCL10/IP-10", "CXCL11/I-TAC" )
## original dataset: immune.trajectories.long.endpoint
#' perform paired testing of associations.
#'
#' Note: deprecated function. This function can be completely replaced by test_association.
#' Only thing missing is the check of blocks, whether blocks with n <= 2
#' observations exists.
#'
#' @param single.pair character varlue of length 2
#' @param data.set data frame to use
#' @param resampling.number number resamplings (see also coin::approximate)
#' @param stratum character value. blocking vars to use
#'
#' @importFrom rlang .data
#'
perform_single_pair_test <- function(
    single.pair,
    data.set,
    resampling.number = 10^4,
    stratum           = NULL
){
  # single.pair <- c("Granuloytes", "Monocytes") # for debugging
  if(length(single.pair) != 2){warning(
    "single pair should contain exactly two variables to compare"
  )}
  vars.to.select <- c(single.pair, stratum)
  names(single.pair) <- NULL
  # only select complete cases:
  data.for.test.i <- data.set %>%
    dplyr::select(vars.to.select) %>%
    stats::na.omit()

  formula.string <- paste0(addq(single.pair), collapse = "~")
  # Add block to formula if existing:
  if(!is.null(stratum)){
    data.for.test.i <- data.for.test.i %>%
      # if block contains more than one variable, this will concatenate them:
      tidyr::unite(col = "block", tidyselect::all_of(stratum)) %>%
      dplyr::mutate(block = factor(.data$block))
    formula.string <- stats::as.formula(paste0(formula.string, "|block"))
    ## test whether blocks with two or less measurements exist:
    obs.per.block <- data.for.test.i %>%
      dplyr::group_by(.data$block) %>%
      dplyr::summarize('n.per.block' = dplyr::n())
    # omit blocks with not enough measurements (the correlation
    # is perfect/ does not make sense):
    if(any(obs.per.block$n.per.block<=2)){
      blocks.to.print <- obs.per.block %>%
        dplyr::filter(.data$n.per.block<=2) %>%
        tidyr::unite(col ="block", sep = ": ") %>%
        dplyr::pull(.data$block)
      warning(paste0("some blocks with too few observations detected;",
                     "these will be omitted from analysis:\n",
                     paste0(blocks.to.print, collapse = ", ")))
    }
    data.for.test.i <- dplyr::left_join(data.for.test.i, obs.per.block,
                                 by = "block") %>%
      dplyr::filter(.data$n.per.block >2)
  }
  formula.to.test  <- stats::as.formula(formula.string)

  test.results     <- coin::spearman_test(
    formula.to.test,
    data = data.for.test.i,
    distribution = coin::approximate(nresample = resampling.number)
  )
  #### calculate rho: ####
  global.rho <- weighted_average_rho(data = data.for.test.i, x = single.pair[1],
                                     y = single.pair[2], stratum = "block")
  ######### End of calculating rho ###########
  data.to.return.i <- data.frame(
    "first.var"   = single.pair[1],
    "second.var"  = single.pair[2],
    "sample.size" = nrow(data.for.test.i),
    "stratified.by" = "",
    "Direction"   = sign(coin::statistic(test.results)),
    "rho"         = global.rho,
    "p.value"     = as.numeric(coin::pvalue(test.results)),
    stringsAsFactors = F
  )
  if(!is.null(stratum)){
    data.to.return.i$`stratified.by` <- paste(stratum, collapse = "_")
  }
  return(data.to.return.i)
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
#' @param response.names character value of column that contains all the names of the response values that need to be investigated separately
#' @param ... other values will be parsed to \code{\link{test_association}}
#'
#' @return data frame with results as output
#' @export
#'
association_study_long <- function(data, response.names, ...){
  purrr::map_dfr(unique(data[[response.names]]), .f = function(x){
    df <- test_association(data[data[[response.names]] == x, ], ...)
    df["Response var"] <- x
    df
  })
}
