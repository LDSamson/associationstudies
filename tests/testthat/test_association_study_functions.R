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
                                    response.var = "Frailty.index",
                                    vars.to.select = cols_to_analyze
                                    )
  expect_equal(is.data.frame(test_outcome), TRUE)
  expect_equal(nrow(test_outcome), 20)
  expect_equal(names(test_outcome), tab_names)
})

# library(dplyr)
# cols_to_analyze <- immune_data %>%
#  select(-c(Batch, Sex, Frailty.index)) %>% names()
# library(corrplot)
#test_outcome <- association_study(immune_data, n.resample = 1000)
# test_outcome <- BH_selection(test_outcome, FDR_cutoff = 0.05)

# plot_p_histogram(test_outcome)

# b <- associationstudies::create_matrix(
#   test_outcome, "response.var", "explanatory.var",
#   markers.to.test = names(immune_data))
# c <- associationstudies::create_matrix(
#   test_outcome, "response.var", "explanatory.var", values = "FDR_selection",
#   markers.to.test = names(immune_data))
# corrplot(b, tl.col = "black",
#          method = 'ellipse', type = 'upper', diag = FALSE,
#          order = "hclust",  addrect = 3,
#          p.mat = c, insig = "label_sig",
#          pch = "*", pch.cex = 1.5 )
#
# a <- association_study(immune_data, exclude = c("Sex"))
# a <- association_study(as.matrix(immune_data), "Frailty.index",
#                        vars.to.select = -c(Tregs, Neutrophils),
#                        stratum = c("Batch", "Sex"))
# association_study(immune_data[1:10])
# a <- immune_data %>%
#   select(-c("Tregs", "Neutrophils"), "Frailty.index")
