# library(associationstudies)

tab_names <- c('Response.var','Explanatory.var','sample.size',
               'stratified.by','n.resample','Method','Direction', 'rho', 'p.value')

test_that("Ouput is not a data frame", {
  testthat::expect_equal(
    is.data.frame(test_association(immune_data, "Tregs", "Frailty.index")),
    TRUE
    )
})
test_that("Output names incorrect", {

  testthat::expect_equal(names(test_association(immune_data, "Tregs", "Frailty.index")),
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
                                         cols_to_analyze, "Frailty.index")
  expect_equal(is.data.frame(test_outcome), TRUE)
  expect_equal(nrow(test_outcome), 20)
  expect_equal(names(test_outcome), tab_names)
})
