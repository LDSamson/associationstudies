#' Function to perform Random Forest analysis
#' Note: dataset is expected to be in "wide" format for this function.
#'
#' @param data
#' @param subset
#' @param predictors
#' @param confounders
#' @param dependent.var
#' @param ntree
#' @param select.best.mtry
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
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

