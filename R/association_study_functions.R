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
association_study <- function(subset.parameter = "Granulocytes_Count",
                              subset.column    = "Cellname",
                              dataset,
                              numerical.value  = "Concentration",
                              comparison.value = "Frailty.index",
                              stratum          = c("Age.category", "CMV.status"),
                              n.resample       = 100,
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

## 08-01-2020 LS:
##    - Fixed bug: Now returns warning and a valid table
##      when no associations are found

#' BH selection
#'
#' Function to create a selection based on a set False Discovery Rate
#' I used this function in combination with the association_study function.
#'
#' @param data a data frame
#' @param p_value character value of column containing the p values
#' @param FDR_cutoff numerical value. FDR cutoff
#' @param hist logical. Whether or not to create a histogram
#'
#' @return data frame with added columns with FDR and FDR selection
#' @importFrom rlang .data
#' @export
BH_selection <- function(data, p_value = "p.value",
                                 FDR_cutoff = 0.15, hist=FALSE){
  # calculate FDR
  n.of.tests <- nrow(data)
  data.sorted <- data %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(p_value)) %>%
    dplyr::arrange(!!rlang::sym(p_value)) %>%
    dplyr::mutate(
      'order.of.tests' = dplyr::row_number(),
      'FDR'            = n.of.tests*!!rlang::sym(p_value)/.data$order.of.tests,
      FDR              = round(.data$FDR, digits = 4)
    )
  if(min(data.sorted$FDR) > FDR_cutoff){
    warning(paste0("No associations detected with FDR of ", FDR_cutoff))
    data.with.BH.selection <- data.sorted %>%
     dplyr::mutate(`FDR_selection` = 0) %>%
      dplyr::select(-.data$order.of.tests)
  }else{
    # Define max number of detected associations, based on FDR cutoff
    max.n.of.tests <- max(data.sorted[data.sorted$FDR <= FDR_cutoff,
                                      "order.of.tests"])
    data.with.BH.selection <- data.sorted %>%
      dplyr::mutate(
        FDR_selection = ifelse(.data$order.of.tests <= max.n.of.tests, 1, 0)
      ) %>% dplyr::select(-.data$order.of.tests)
  }
  # (optional): print p value histogram
  if(isTRUE(hist)){plot_p_histogram(data.with.BH.selection)}
  return(data.with.BH.selection)
}

#' Function to plot P histogram
#'
#' @param data dataframe to use
#' @param p_value character string. name of column containing p values
#' @param title.to.plot title of the plot to produce
#' @param gamma.hat gamma hat to use
#'
#' @return histogram
#' @export
#'
plot_p_histogram <- function(
    data,
    p_value = "p.value",
    title.to.plot = "P-values of the association study",
    gamma.hat = 1
){
  p.values.hist <- ggplot2::ggplot(data, ggplot2::aes(x = !!rlang::sym(p_value), y = ..density..)) +
    ggplot2::geom_histogram(binwidth = 0.1, fill = 'lightblue', col = 'black') +
    ggplot2::geom_hline(yintercept = gamma.hat, col = 'tomato') +
    ggplot2::theme_light() +
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::labs(x = "P value", title = title.to.plot)
  print(p.values.hist)
}

#' Association studies
#'
#' Association study function that can be used with data in wide format:
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
association_study_wide <- function(
    dataset,
    response.var,
    explanatory.var,
    stratum         = NULL,
    n.resample      = 100,
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
#' @param comparison.var comparison var (for example flowratio)
#' @param nresample times to resample for permutation test (see coin::approximate)
#' @importFrom rlang .data
#'
#' @return data frame with results
#' @export
associations_per_ID <- function(x, dataset,
                                comparison.var, nresample = 10^5){
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
## Create function for single pair association
## original dataset: immune.trajectories.long.endpoint
#' @importFrom rlang .data
perform_single_pair_test <- function(
    single.pair,
    data.set,
    resampling.number = 10^4,
    stratum         = c("Age.category", "Batch")
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

#' Association studies in long format
#'
#' This function, in combination with association_study_wide,
#' can basically replace the function 'association_study'.
#' While I used the function association_study in my publications,
#' I think this implementation is better. The function association_study_wide
#' is more intuitive, and this function is just a small wrapper
#' that you can use when data is in long format.
#'
#' @param data data frame to use
#' @param response.names character value of column that contains all the names of the response values that need to be investigated separately
#' @param ... other values that will be parsed to the function association_study_wide
#'
#' @return data frame with results as output
#' @export
#'
association_study_long <- function(data, response.names, ...){
  purrr::map_dfr(unique(data[[response.names]]), .f = function(x){
    df <- association_study_wide(
      data %>% dplyr::filter({{response.names}} == x), ...)
    df["Response var"] <- x
    df
  })
}
