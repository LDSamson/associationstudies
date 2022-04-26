
##### Test gm_ mean functions

test_data <- data.frame(
  x1 = c(1,3,2,5),
  x2 = c(1,3,2,0),
  x3 = c(1,3,NA,5),
  x4 = c(1,3,NA, 0),
  x5 = c(NA, NA, NA, NA),
  x6 = c(0,0,0,0)
)
# gm_mean function shold be tested with test_data above.
# I still need to write unit tests for this function.

##### Test outcome formatting with compare_summ_outcome_format():

# if necessary, loop over the list so that new data is generated:
data_to_test <- list(
  data1 <- rnorm(n = 100, mean = 5, sd = 1),
  data2 <- rnorm(n = 100, mean = 50, sd = 1),
  # check large numbers:
  data3 <- abs(rnorm(n = 100, mean = 10000, sd = 10000)),
  ## check zeroes in dataset:
  data4 <- c(abs(rnorm(n = 100, mean = 10000, sd = 10000)), 0, 0),
  data5 <- c(abs(rnorm(n = 100, mean = 1, sd = 2)), 0, 0),
  # check missing values in dataset:
  data6 <- c(abs(rnorm(n = 100, mean = 1, sd = 2)), 0, 0, NA)
)

compare_summ_outcome_format <- function(i){
  # i <- data6
  gm_outcome <- as.numeric(formatC(gm_mean(i, conf.level = 0.95, na.rm = TRUE),
                                   digits = 5))
  form_outcome <- format_summary_values(i, na.rm = TRUE, n.digits = 4)
  exponent <- stringr::str_extract(form_outcome, "(?<=\\*10\\^)[[1-9]]")
  if(is.na(exponent)) exponent <- "0"
  exponent <- as.numeric(exponent)
  formatted_vals <- strsplit(form_outcome, split = "\\$|\\(|-|\\)")[[1]]
  formatted_vals <- gsub("\\*10\\^[1-9]", "", formatted_vals) # remove exponent
  formatted_vals <- as.numeric(formatted_vals[formatted_vals != ""])*10^exponent
  formatted_vals <- as.numeric(formatC(formatted_vals, digits = 5)) # correct n.digits
  # Rounding errors can occur. Therefore, only test for the first digit:
  list(formatC(gm_outcome, digits = 2),
       formatC(formatted_vals, digits = 2))
}

test_that("formatting gives same values", {
  for(i in data_to_test){
    # i <- data6
    outcomes <- compare_summ_outcome_format(i)
    testthat::expect_equal(outcomes[[1]], outcomes[[2]])
  }
})

# data with negative numbers should fail:
fail_data <- c((rnorm(n = 100, mean = 1, sd = 2)), 0, 0, -1)
test_that("Error with negative values", {
  expect_warning(gm_mean(fail_data, na.rm = TRUE))
})

test_that("expect failure with negative values", {
  #expect_message(format_summary_values(fail_data, na.rm = TRUE, n.digits = 4))
  expect_error(gm_outcome(fail_data, na.rm = TRUE, n.digits = 4))
})
