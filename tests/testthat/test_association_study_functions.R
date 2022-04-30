# library(associationstudies)

tab_names <- c('response.var','explanatory.var','sample.size',
               'stratified.by','n.resample','method','direction', 'rho', 'p.value')

test_that("Ouput is not a data frame", {
  testthat::expect_equal(
    is.data.frame(test_association(immune_data, "Tregs", "Frailty.index")),
    TRUE
    )
})
test_that("Output names incorrect", {

  testthat::expect_equal(names(test_association(immune_data, "Frailty.index", "Tregs")),
                         tab_names)
  testthat::expect_equal(names(test_association(immune_data, "Frailty.index", "Sex")),
                         tab_names)
})

test_that("Error when non-numerical response vars are used", {
  expect_error(test_association(immune_data, "Sex", "Frailty.index"))
  })

test_that("Error when low number of observations in a block", {
  expect_error(test_association(head(immune_data, n=10),
                                "Tregs", "Frailty.index", c("Sex", "Batch")))
})


#### Test association_study_long:
immune_data_long <- immune_data %>%
  tidyr::pivot_longer(-c(Batch, Sex, Frailty.index))

test_that("Unexpected output", {
 test_outcome <- association_study_long(immune_data_long,
                                        "name", "value", "Frailty.index")
  expect_equal(is.data.frame(test_outcome), TRUE)
  expect_equal(nrow(test_outcome), 20)
  expect_equal(names(test_outcome), tab_names)
})

test_that("Error when low number of observations in a block", {
  expect_error(association_study_long(head(immune_data_long, 40),
                                      "name", "value", "Frailty.index", "Batch"))
})


## test association_study:
cols_to_analyze <- unique(immune_data_long$name)
test_that("Unexpected output", {
  test_outcome <- association_study(immune_data,
                                    response.var = "Frailty.index",
                                    vars.to.select = cols_to_analyze
                                    )
  expect_equal(is.data.frame(test_outcome), TRUE)
  expect_equal(nrow(test_outcome), 20)
  expect_equal(names(test_outcome), tab_names)
})


############ Test function association_study
standard_test <- function(data = test_output){
  test_that("Unexpected standard output", {

    expect_true(is.data.frame(data))
    expect_true(all(names(data) == tab_names))
    expect_true(all(
      unlist(lapply(data[, c("sample.size", "n.resample",
                             "direction", "rho", "p.value")],
                    is.numeric))
    ))
    expect_true(all(data$p.value <= 1))
    expect_true(all(na.omit(abs(data$rho)) <= 1))
  })
}

# Test pairwise associations between all variables in a a data frame:
standard_test(association_study(immune_data[1:5]))

#Test associations of a response variable with all other variables in a data frame:
standard_test(association_study(immune_data, "Frailty.index"))

# Test associations with a selection of variables:
test_output <- association_study(immune_data, "Frailty.index", c(Tregs, Neutrophils))
standard_test()
test_that("First var name should be the same", {
  expect_equal(unique(test_output$response.var), "Frailty.index")
})

# Test excluding variables:
test_output <- association_study(immune_data, "Frailty.index", -c(Tregs, Neutrophils))
standard_test()
test_that("variables are not excluded appropriately", {
  expect_false(any(c("Tregs", "Neutrophils") %in% test_output$explanatory.var))
})

# test a blocked design:
test_output <- association_study(immune_data, "Frailty.index",
                                 stratum = c("Batch", "Sex"))
standard_test()
test_that("Blocking unexpected", {
  expect_true(all(test_output$stratified.by == "Batch_Sex"))
})

test_that("not enough data should give an error", {
  expect_error(association_study(head(immune_data), "Frailty.index",
                                 stratum = c("Batch", "Sex")))
})
