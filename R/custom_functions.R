####################### colorblind friendly color palette: #####################

# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# For scientific notation:
# https://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales
scientific_10 <- function(x) {
  parse(text=gsub("\\+", "", # addition to remove the + sign
                  gsub("e", " %*% 10^", scales::scientific_format()(x))
                  )
        )
}


## helper function: function to create FDR or rho matrix for input in corrplot.
## Expects input from p values/FDR values or rhos, in format as output from
## association studies (see below for function association.study)
create.matrix <- function(data = correlation.results, x = "first.var",
                          y = "second.var",values = "rho"){
  matrix.1 <- select(.data = data, x, y, values)
  matrix.2 <- select(.data = data, y, x, values)
  names(matrix.2) <- names(matrix.2)[c(2,1,3)]
  matrix.full <- bind_rows(matrix.1, matrix.2) %>%
    pivot_wider(names_from = y, values_from = values)
  matrix.full[[x]] <- factor(matrix.full[[x]], levels = markers.to.test)
  matrix.full <- matrix.full %>%
    arrange(!!sym(x)) %>%
    as.data.frame()
  if(values == "rho"){
    matrix.full <- matrix.full %>%
      # for the middle row (cor=1), replace NA (not tested) with cor=1
      mutate_all(.funs = ~ifelse(is.na(.), 1, .))
  }
  if(values == "FDR_selection"){
    matrix.full <- matrix.full %>%
      # for the middle row, replace NA (not tested) with FDR_selection=0
      mutate_all(.funs = ~ifelse(is.na(.), 0, .)) %>%
      # ! below will change the FDR selection so that it can be used directly
      # in the corrplot function as if they are p values
      mutate_all(.funs = ~ifelse(. == 1, 0.01, 0.99))
  }
  rownames(matrix.full) <- markers.to.test
  matrix.full <- matrix.full %>%
    select(markers.to.test) %>%
    as.matrix()
  return(matrix.full)
}


### helper function: calculates weighted average rho in blocked correlation study
## as discussed with Jose Ferreira.
weighted.average.rho <-function(data = data.for.test.i, x, y, stratum = "block"){
  if(is.null(stratum)){
   # estimates.of.rho <- data.frame()
    global.rho <-  cor(data[[x]],
                       data[[y]],
                       method = "spearman")
  }else{
    # 06-11-2019 LS:
    # Adding rho as weighted average in the results,
    # following the advise from José Ferreira
    # First, calculate rho per block with base R cor function:
    blocks <- unique(data[[stratum]])
    estimates.of.rho <- data.frame()
    for(b in blocks){
      # b <- "70-75yr_CMV-"
      data.set.b <- data[data[[stratum]] == b,]
      rho.b <-  cor(data.set.b[[x]], data.set.b[[y]], method = "spearman")
      estimates.of.rho <- rbind(estimates.of.rho,
                                data.frame(rho.hat=rho.b, block = b,
                                           sample.size=nrow(data.set.b)))
    }
    # Now, calculate weighted average from the individual rho's per block:
    global.rho <- sum(estimates.of.rho$rho.hat*estimates.of.rho$sample.size)/
      sum(estimates.of.rho$sample.size)
  }
  return(global.rho)
}




########## Function to add pagebreaks in Rmarkdown files. ###################
## For more info see: https://stackoverflow.com/a/55064070

pagebreak <- function() {
  if(knitr::is_latex_output())
    return("\\newpage")
  else
    return('<div style="page-break-before: always;" />')
}

####################### Function to perform Random Forest analysis #################
## Note: dataset is expected to be in "wide" format for this function.

RFanalysis <- function(data = ISA.cellular,
                       subset = NULL,
                       predictors = cell.subsets.to.analyse,
                       confounders = confounding.vars,
                       dependent.var = dependent.variable,
                       ntree = 10000,
                       select.best.mtry = TRUE,
                       ...){
  ## set.seed(2019) # Is this necessary?? Not sure
  ### prepare data. Select subset
  if(!is.null(subset)) {data <- data %>% filter(Sex == subset)}
  data <- na.omit(data) # RF cannot handle missing data

  ## Prepare formula for the Random forest function
  ## Remove sex from confounders if the analysis is per sex
  ## (in that case, the variable obviously does not contain information for the model  anymore)
  if(!is.null(confounders)){
    if(!is.null(subset) & ("Sex" %in% confounders) ){
      confounders <- confounders[!confounders ==  "Sex"]
      }
    RF.formula <- paste0(dependent.var, "~", paste0(c(confounders, predictors), collapse="+")) %>% as.formula
  }else{
    RF.formula <- paste0(dependent.var, "~", paste0(predictors, collapse="+")) %>% as.formula
  }

  if(isTRUE(select.best.mtry)){
    ### select the mtry that gives the lowest RMSE:
    set.seed(2019)
    optimize.RF <<- caret::train(form = RF.formula, data = data,
                                 method="rf",
                                 ntree=500,
                                 tuneGrid = data.frame(mtry = 5:20),
                                 control = trainControl(method = "oob")
    )
    print(optimize.RF)
    # mtry with the lowest RMSE found after repeated 'out of bag' bootstrapping:
    best.mtry <- optimize.RF$bestTune$mtry
    set.seed(2019)
    RFOutcome <<- randomForest(formula = RF.formula, data = data,
                               ntree=ntree,importance=TRUE,
                               mtry = best.mtry, ...)
  }else{
    set.seed(2019)
    RFOutcome <<- randomForest(formula = RF.formula, data = data,
                               ntree=ntree,importance=TRUE, ...)
    ## mtry=ncol(data)-1 # other argument to use/originally used?
  }
  print(RFOutcome)

  if(isTRUE(select.best.mtry)){plot(optimize.RF)}

  varImpPlot(RFOutcome,sort=TRUE,type=1,
             main=paste(subset, dependent.var,
                        "associations \nVariable importance", sep = " " ))

  varImpPlot(RFOutcome,sort=TRUE,type=1,
             main=paste(subset, dependent.var,
                        "associations \nVariable importance", sep = " " ),
             n.var = 15)
}

