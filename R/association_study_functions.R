####### Function to perform associations between two parameters: ############################

## Associations with p values calculated by permutation, using coin package
## Note: dataset is expected to be in "long" format for this function.
## Note2: missing values are (silently) deleted

##  a numerical value and a second parameter (comparison value).
## If the second parameter is a factor, a wilcoxon test will be used.
## If the second parameter is a numerical value, a spearman association will be used.
## Use the function in an apply loop (or map_dfr loop) to perform it on
## many different subsets in the data.

## 08-01-2020 LS:
##    - Removed argument CMV.status != "Missing"
##        (this gives an error when CMV.status is missing)
##    - Added calculation of rho when no blocking variable is used

association.study <- function(subset.parameter = "Granulocytes_Count",
                              subset.column    = "Cellname",
                              dataset          = data.long,
                              numerical.value  = "Concentration",
                              comparison.value = "Frailty.index",
                              Stratum          = c("Age.category", "CMV.status"),
                              n.resample       = 100,
                              ...){
  data.to.analyze <- dataset %>%
    select(all_of(c(subset.column, numerical.value,
                    comparison.value, Stratum))) %>%
    filter(!!sym(subset.column) == subset.parameter) %>%
    na.omit() %>%
    droplevels()
  if(!is.null(Stratum)){
    data.to.analyze <- data.to.analyze %>%
      unite(col = "block", all_of(Stratum)) %>%
      mutate(block = factor(block))
    formula.to.test <- as.formula(paste0(numerical.value, "~", comparison.value, "|block"))
  }else {
    formula.to.test <- as.formula(paste0(numerical.value, "~", comparison.value))
  }
  if(is.factor(dataset[[comparison.value]])){
    ## test results when comparing a continuous variable between groups/factor:
    test.results <- coin::wilcox_test(formula.to.test,
                                      data = data.to.analyze,
                                      distribution=approximate(nresample=n.resample,
                                                               ...)
    )
  }else if(is.numeric(dataset[[comparison.value]])){
    ## test results when testing two continuous variables:
    test.results <- coin::spearman_test(
      formula.to.test,
      data = data.to.analyze,
      distribution=approximate(nresample=n.resample, ...)
    )
    global.rho <- weighted.average.rho(
      data = data.to.analyze,
      x = numerical.value,
      y = comparison.value,
      stratum = if(is.null(Stratum)){NULL}else "block"
    )
  }
  p.val <- as.numeric(pvalue(test.results))
  Direction <- sign(statistic(test.results))
  results.to.return <- list(
    "Subsets"          = subset.parameter,
    "Comparison value" = comparison.value,
    "sample.size"      = nrow(data.to.analyze),
    "stratified.by"    = "",
    "n.resample"       = n.resample,
    "Method"           = test.results@method,
    "Direction"        = Direction,
    "p.value"          = p.val
  )
  ## Add rho to outcome, if applicable:
  if(results.to.return[["Method"]] == "Spearman Correlation Test"){
    results.to.return[["rho"]] <- global.rho}
  if(!is.null(Stratum)){
    results.to.return$`stratified.by` <- paste(Stratum, collapse = "_")
  }
  ## Return results:
  return(results.to.return)
}

############### Function to create a selection based on FDR ####################
## This function can be used after the previous one.

## 08-01-2020 LS:
##    - Fixed bug: Now returns warning and a valid table
##      when no associations are found

perform.BH.selection <- function(data, identifier, p_value = "p.value",
                                 FDR_cutoff = 0.15, hist=FALSE){
  # calculate FDR
  n.of.tests <- nrow(data)
  data.sorted <- data %>%
    ungroup() %>%
    # select(identifier, p_value) %>%
    filter(!is.na(p_value)) %>%
    arrange(!!sym(p_value)) %>%
    mutate(order.of.tests = row_number(),
           FDR            = n.of.tests*!!sym(p_value)/order.of.tests,
           FDR            = round(FDR, digits = 4)
    )
  if(min(data.sorted$FDR) > FDR_cutoff){
    warning(paste0("No associations detected with FDR of ", FDR_cutoff))
    data.with.BH.selection <- data.sorted %>%
      mutate(`FDR_selection` = 0) %>%
      select(-order.of.tests)
  }else{
    # Define max number of detected associations, based on FDR cutoff
    max.n.of.tests <- max(data.sorted[data.sorted$FDR <= FDR_cutoff,
                                      "order.of.tests"])
    data.with.BH.selection <- data.sorted %>%
      mutate(
        FDR_selection = ifelse(order.of.tests <= max.n.of.tests, 1, 0)
      ) %>% select(-order.of.tests)
  }
  # (optional): print p value histogram
  if(isTRUE(hist)){
    gamma.hat <- 1
    p.values.plot <- ggplot(data.with.BH.selection,
                            aes(x = !!sym(p_value), y = ..density..)) +
      geom_histogram(binwidth = 0.1, fill = 'lightblue', col = 'black') +
      geom_hline(yintercept = 1, col = 'tomato') +
      theme_light() +
      scale_y_continuous(expand = c(0,0)) +
      labs(x = "P value",
           title = "P-values of the association study")
    print(p.values.plot)
  }
  return(data.with.BH.selection)
}

plot.p.histogram <- function(
    data = data.with.BH.selection,
    p_value = "p.value",
    title.to.plot = "P-values of the association study",
    gamma.hat = 1
){
  #  gamma.hat <- 1
  p.values.hist <- ggplot(data, aes(x = !!sym(p_value), y = ..density..)) +
    geom_histogram(binwidth = 0.1, fill = 'lightblue', col = 'black') +
    geom_hline(yintercept = gamma.hat, col = 'tomato') +
    theme_light() +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "P value", title = title.to.plot)
  print(p.values.hist)
}


## same association study function but now in wide format:
association.study.wide <- function(
    response.var    = "Granulocytes_Count",
    dataset         = data.wide,
    explanatory.var = "Frailty.index",
    stratum         = NULL,
    n.resample      = 100,
    ...){
  data.to.analyze <- dataset %>%
    select(response.var, explanatory.var, stratum) %>%
    na.omit() %>%
    droplevels()

  if(!is.null(stratum)){
    data.to.analyze <- data.to.analyze %>%
      unite(col = "block", stratum) %>%
      mutate(block = factor(block))
    formula.to.test <- as.formula(paste0(response.var, "~", explanatory.var,
                                         "|block"))
  }else {
    formula.to.test <- as.formula(paste0(response.var, "~", explanatory.var))
  }
  if(is.factor(dataset[[explanatory.var]])){
    ## test results when comparing a continuous variable between groups/factor:
    test.results <- coin::wilcox_test(formula.to.test, data = data.to.analyze,
                                      distribution=approximate(
                                        nresample=n.resample, ...
                                      )
    )
  }else if(is.numeric(dataset[[explanatory.var]])){
    ## test results when testing two continuous variables:
    test.results <- coin::spearman_test(formula.to.test, data = data.to.analyze,
                                        distribution=approximate(
                                          nresample=n.resample, ...
                                        )
    )
    global.rho <- weighted.average.rho(data = data.to.analyze, x = response.var,
                                       y = explanatory.var, stratum = "block")
  }
  p.val <- as.numeric(pvalue(test.results))
  Direction <- sign(statistic(test.results))
  results.to.return <- list("Response var"    = response.var,
                            "Explanatory var" = explanatory.var,
                            "sample.size"     = nrow(data.to.analyze),
                            "n.resample"      = n.resample,
                            "Method"          = test.results@method,
                            "Direction"       = Direction,
                            "p.value"         = p.val)
  ## Add rho to outcome, if applicable:
  if(results.to.return[["Method"]] == "Spearman Correlation Test"){
    results.to.return[["rho"]] <- global.rho}
  ## Return results:
  return(results.to.return)
}


#############
# This function will investigate the association between biomarkers over time
# with ID as blocking variable
# Fev1PercPred, FvcPercPred Fev1FvcPercPred
associations.per.ID <- function(x, dataset = immune.trajectories.long,
                                comparison.var = "flowratio", nresample = 10^5){
  form.to.test <- as.formula(paste0(comparison.var, "~", x, "|ID"))
  data.to.analyze <- dataset %>%
    select(c("ID", "DoetRound", x, comparison.var, "Sex")) %>%
    na.omit() %>%
    mutate(ID = factor(ID))

  ## Remove people that have less than two complete observations:
  IDs.with.2.or.more.complete.obs <- data.to.analyze %>%
    group_by(ID) %>%
    summarize(n.obs = n()) %>%
    filter(n.obs >=2) %>%
    pull(ID)
  data.to.analyze <- data.to.analyze %>%
    filter(ID %in% IDs.with.2.or.more.complete.obs)

  test.results <- spearman_test(
    form.to.test,
    data = data.to.analyze,
    distribution = approximate(nresample = nresample)
  )
  #### calculate rho: ####
  global.rho <- weighted.average.rho(data = data.to.analyze, x = x,
                                     y = comparison.var, stratum = "ID")
  p.val <- as.numeric(pvalue(test.results))
  Direction <- sign(statistic(test.results))
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
perform.single.pair.test <- function(
    single.pair,
    data.set = immune.trajectories.long.endpoint,
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
    na.omit()
  # need backticks in complicated names:
  # https://stackoverflow.com/questions/16674045/as-formula-in-r-doesnt-seem-to-accept-a-name-that-starts-with-a-number-followed
  addq <- function(x) paste0("`", x, "`")

  formula.string <- paste0(addq(single.pair), collapse = "~")
  # Add block to formula if existing:
  if(!is.null(stratum)){
    data.for.test.i <- data.for.test.i %>%
      # if block contains more than one variable, this will concatenate them:
      unite(col = "block", stratum) %>%
      mutate(block = factor(block))
    formula.string <- as.formula(paste0(formula.string, "|block"))
    ## test whether blocks with two or less measurements exist:
    obs.per.block <- data.for.test.i %>%
      group_by(block) %>%
      summarize(n.per.block = n())
    # omit blocks with not enough measurements (the correlation
    # is perfect/ does not make sense):
    if(any(obs.per.block$n.per.block<=2)){
      blocks.to.print <- obs.per.block %>% filter(n.per.block<=2) %>%
        unite( col ="block", sep = ": ") %>%
        pull(block)
      warning(paste0("some blocks with too few observations detected;",
                     "these will be omitted from analysis:\n",
                     paste0(blocks.to.print, collapse = ", ")))
    }
    data.for.test.i <- left_join(data.for.test.i, obs.per.block,
                                 by = "block") %>%
      filter(n.per.block >2)
  }
  formula.to.test  <- as.formula(formula.string)

  test.results     <- spearman_test(
    formula.to.test,
    data = data.for.test.i,
    distribution = approximate(nresample = resampling.number)
  )
  #### calculate rho: ####
  global.rho <- weighted.average.rho(data = data.for.test.i, x = single.pair[1],
                                     y = single.pair[2], stratum = "block")
  ######### End of calculating rho ###########
  data.to.return.i <- data.frame(
    "first.var"   = single.pair[1],
    "second.var"  = single.pair[2],
    "sample.size" = nrow(data.for.test.i),
    "stratified.by" = "",
    "Direction"   = sign(statistic(test.results)),
    "rho"         = global.rho,
    "p.value"     = as.numeric(pvalue(test.results)),
    stringsAsFactors = F
  )
  if(!is.null(stratum)){
    data.to.return.i$`stratified.by` <- paste(stratum, collapse = "_")
  }
  return(data.to.return.i)
}
