

test_that("Normal usage of count_observation", {
  expect_equal(as.data.frame(count_observations(immune_data)), data.frame(n = 100))
  expect_equal(as.data.frame(count_observations(immune_data, "Batch")),
               data.frame(Batch = c("1", "2", "3"), n = c(27, 35, 38)))
  expect_equal(as.data.frame(count_observations(immune_data, "Batch", "Sex")),
               data.frame(Batch_Sex = c("Women_1", "Men_1", "Men_2", "Men_3",
                                        "Women_2", "Women_3"),
                          n = c(12, 15, 16, 17, 19, 21)))
})

test_that("Error when incorrect stratum or response names", {
  expect_error(count_observations(immune_data, "Batchy"))
  expect_error(count_observations(immune_data, "Batch", "Typo_"))
  expect_error(count_observations(immune_data, "Batchy", "Sex"))
})

test_that("error when not enough observations in data", {
expect_error(check_low_group_numbers(head(immune_data),
                                     stratum = c("Batch"), "Sex"))})

