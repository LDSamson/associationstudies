association_study()

saveRDS(data.to.test, "data.to.test.rds")

data.to.test <- selected.phosflow.baseline.data %>%
  select(ID, Sex, FI.2017, condition, Unstimulated, Phosflow.batch) %>%
  pivot_wider(names_from = condition, values_from = Unstimulated) %>%
  mutate(ID = row_number())


response.vars <- names(data.to.test)[!names(data.to.test) %in% c(stratum, explanatory.var, "ID", "Sex")]


## maybe use association_study wide and loop this over a column named condition?
# attempt to reduce the code in package
thesisfunctions::association_study_wide(data.to.test,
                       response.var =  "STAT1_B cells",
                       explanatory.var =  "FI.2017",
                       stratum = "Phosflow.batch")

data.to.test.long$condition
data.to.test.long <- selected.phosflow.baseline.data %>%
  select(ID, Sex, FI.2017, condition, Unstimulated, Phosflow.batch)
a <- lapply(unique(data.to.test.long[["condition"]]), FUN = function(x){
  df <- thesisfunctions::association_study_wide(
    data.to.test.long %>% filter(condition == x),
    response.var =  "Unstimulated",
    explanatory.var =  "FI.2017",
    stratum = "Phosflow.batch"
  )
  df["Response var"] <- x
  df
}
)
b <- as.data.frame(do.call("rbind", a))
library(purrr)
a <- map_dfr(unique(data.to.test.long[["condition"]]), .f = function(x){
  df <- thesisfunctions::association_study_wide(
    data.to.test.long %>% filter(condition == x),
    response.var =  "Unstimulated",
    explanatory.var =  "FI.2017",
    stratum = "Phosflow.batch"
  )
  df["Response var"] <- x
  df
}
)

thesisfunctions::association_study_wide(data.to.test,
                                        response.var =  "STAT1_B cells",
                                        explanatory.var =  "FI.2017",
                                        stratum = "Phosflow.batch")
response.vars <- "condition"
data <- data.to.test.long
a <- lapply(unique(data[[response.vars]]),
            FUN = function(x) filter(data, condition == x))
b <- map_dfr(a, thesisfunctions::association_study_wide,
        response.var =  "Unstimulated",
        explanatory.var =  "FI.2017",
        stratum = "Phosflow.batch"
        )

b["Response var"] <-

association_study_long <- function(data, response.names, ...){
  map_dfr(unique(data[[response.names]]), .f = function(x){
    df <- thesisfunctions::association_study_wide(
      data %>% filter(condition == x), ...)
    df["Response var"] <- x
    df
  })
}
b <- association_study_long(
  data.to.test.long,
  response.names = "condition",
  response.var =  "Unstimulated",
  explanatory.var =  "FI.2017",
  stratum = "Phosflow.batch",
  n.resample = 100000
  )

library(purrr)
a <- map_dfr(response.vars,
  .f = ~thesisfunctions::association_study_wide(
  dataset = data.to.test %>% filter(Sex == "Men"),
  response.var =  .,
  explanatory.var =  "FI.2017",
  stratum = "Phosflow.batch",
  n.resample = 10000
  )
) %>%
  BH_selection()
names(a) %in% names(b)

b <- map_dfr(
  all_of(baseline.conditions.to.test),  .f = ~association.study(
    dataset = selected.phosflow.baseline.data %>%
      filter(Sex == "Men"),
    subset.parameter = .,
    response.vars    = "condition",
    numerical.value  = "Unstimulated",
    comparison.value = "FI.2017",
    Stratum = "Phosflow.batch",
    n.resample = 10000
  )
) %>% perform.BH.selection()
phosflow_assoc_baseline_fi_men



association_study_wide()

unite_vars <- function(x, vars, colname = "block", ...){
 # x <- data.to.test
#  vars <- NULL
  data.to.return <- x %>%
    tidyr::unite(col = {{colname}}, all_of(vars))
  data.to.return[[colname]] <- factor(data.to.return[[colname]])
  return(data.to.return)
}

association_study_wide <- function(
    dataset,
    response.var,
    explanatory.var,
    stratum         = NULL,
    n.resample      = 1000,
    ...){
  dataset <- data.to.test
  response.var <- "STAT1_B cells"
  explanatory.var =  "FI.2017"
  data.to.analyze <- dataset %>%
    dplyr::select(c(response.var, explanatory.var, tidyselect::all_of(stratum))) %>%
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
  p.val <- as.numeric(coin::pvalue(test.results))
  Direction <- sign(coin::statistic(test.results))
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
  if(!is.null(stratum)){
    results.to.return$`stratified.by` <- paste(stratum, collapse = "_")}
  ## Return results:
  return(results.to.return)
}

data.to.test.long <- selected.phosflow.baseline.data %>%
  select(ID, Sex, FI.2017, condition, Unstimulated, Phosflow.batch)


association_study_new <- function(dataset,
                              subset.parameter = "",
                              response.vars    = "condition",
                              numerical.value  = "Unstimulated",
                              comparison.value = "FI.2017",
                              stratum          = c("Phosflow.batch"),
                              n.resample       = 1000,
                              ...){
  dataset <- data.to.test.long
  response.vars <- "condition"
  numerical.value  = "Unstimulated"
  comparison.value = "FI.2017"
  stratum          = c("Phosflow.batch")
  n.resample       = 1000

  data.to.analyze <- dataset %>%
    dplyr::select(tidyselect::all_of(c(response.vars, numerical.value,
                                       comparison.value, stratum))) %>%
    stats::na.omit() %>%
    droplevels()

  string_formula <- paste0(addq(numerical.value), "~", addq(comparison.value))
  if(!is.null(stratum)){
    data.to.analyze <- unite_vars(data.to.analyze, vars = stratum, colname = "block")
    string_formula <- paste0(string_formula, "|block")
  }
  formula.to.test <- stats::as.formula(string_formula)

  lapply(unique(dataset[[response.vars]]), FUN = function(i){
   # i <- "STAT1_B cells"
    data.to.analyze.i <- data.to.analyze %>%
      # below is necessary so that only one association is tested:
      dplyr::filter(!!rlang::sym(response.vars) == i) %>%
      droplevels()

      if(is.factor(dataset[[comparison.value]])){
    ## test results when comparing a continuous variable between groups/factor:
    test.results <- coin::wilcox_test(formula.to.test,
                                      data = data.to.analyze.i,
                                      distribution=coin::approximate(nresample=n.resample,
                                                                     ...)
                                      )
  }else if(is.numeric(dataset[[comparison.value]])){
    ## test results when testing two continuous variables:
    test.results <- coin::spearman_test(
      formula.to.test,
      data = data.to.analyze.i,
      distribution=coin::approximate(nresample=n.resample#, ...
                                     )
      )
    global.rho <- weighted_average_rho(
      data = data.to.analyze.i,
      x = numerical.value,
      y = comparison.value,
      stratum = if(is.null(stratum)){NULL}else "block"
      )
  }

## how to return results and combine them properly?=
results.to.return.i <- c(
    "Subsets"          = i,
    "Comparison value" = comparison.value,
    "sample.size"      = nrow(data.to.analyze),
    "stratified.by"    = paste(stratum, collapse = "_"),
    "n.resample"       = n.resample,
    "Method"           = test.results@method,
    "Direction"        = sign(coin::statistic(test.results)),
    "p.value"          = as.numeric(coin::pvalue(test.results))
  )
})

## Add rho to outcome, if applicable:
  if(results.to.return[["Method"]] == "Spearman Correlation Test"){
    results.to.return[["rho"]] <- global.rho}
  ## Return results:
  return(results.to.return)
}
