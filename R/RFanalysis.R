#' Function to perform Random Forest analysis
#' Note: dataset is expected to be in "wide" format for this function.
#'
#' @param data dataframe to use (original: ISA.cellular)
#' @param subset whether or not to use subset of data
#' @param predictors character value, predictors (original: cell.subsets.to.analyse)
#' @param confounders character value, confounding variables (original: confounding.vars)
#' @param dependent.var character value (original: dependent.variable)
#' @param ntree numerical, number of trees to use in random forest
#' @param select.best.mtry logical
#' @param ... other parameters parsed to the randomForest function
#'
#' @return a randomForest object
#' @importFrom rlang .data
#' @export
RFanalysis <- function(data,
                       subset = NULL,
                       predictors,
                       confounders,
                       dependent.var,
                       ntree = 10000,
                       select.best.mtry = TRUE,
                       ...){
  ## set.seed(2019) # Is this necessary?? Not sure
  ### prepare data. Select subset
  if(!is.null(subset)) {data <- data %>% dplyr::filter(.data$Sex == subset)}
  data <- stats::na.omit(data) # RF cannot handle missing data

  ## Prepare formula for the Random forest function
  ## Remove sex from confounders if the analysis is per sex
  ## (in that case, the variable obviously does not contain information for the model  anymore)
  if(!is.null(confounders)){
    if(!is.null(subset) & ("Sex" %in% confounders) ){
      confounders <- confounders[!confounders ==  "Sex"]
    }
    RF.formula <- paste0(dependent.var, "~", paste0(c(confounders, predictors), collapse="+")) %>% stats::as.formula
  }else{
    RF.formula <- paste0(dependent.var, "~", paste0(predictors, collapse="+")) %>% stats::as.formula
  }

  if(isTRUE(select.best.mtry)){
    ### select the mtry that gives the lowest RMSE:
    set.seed(2019)
    optimize.RF <<- caret::train(form = RF.formula, data = data,
                                 method="rf",
                                 ntree=500,
                                 tuneGrid = data.frame(mtry = 5:20),
                                 control = caret::trainControl(method = "oob")
    )
    print(optimize.RF)
    # mtry with the lowest RMSE found after repeated 'out of bag' bootstrapping:
    best.mtry <- optimize.RF$bestTune$mtry
    set.seed(2019)
    RFOutcome <<- randomForest::randomForest(formula = RF.formula, data = data,
                               ntree=ntree,importance=TRUE,
                               mtry = best.mtry, ...)
  }else{
    set.seed(2019)
    RFOutcome <<- randomForest::randomForest(formula = RF.formula, data = data,
                               ntree=ntree,importance=TRUE, ...)
    ## mtry=ncol(data)-1 # other argument to use/originally used?
  }
  print(RFOutcome)

  if(isTRUE(select.best.mtry)){plot(optimize.RF)}

  randomForest::varImpPlot(RFOutcome,sort=TRUE,type=1,
             main=paste(subset, dependent.var,
                        "associations \nVariable importance", sep = " " ))

  randomForest::varImpPlot(RFOutcome,sort=TRUE,type=1,
             main=paste(subset, dependent.var,
                        "associations \nVariable importance", sep = " " ),
             n.var = 15)
}

