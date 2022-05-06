
####################### Test gm_mean function ######################

test_data <- data.frame(
  x1 = c(1,3,2,5),
  x2 = c(1,3,2,0),
  x3 = c(1,3,NA,5),
  x4 = c(1,3,NA, 0),
  x5 = c(NA, NA, NA, NA),
  x6 = c(0,0,0,0)
)

####                            X1         x2        x3        x4     x5    x6
expected_outcomes         <- c(2.340347, 1.565085, 2.466212, 1.442250, NaN, 0)
expected_na_outcomes      <- c(2.34035,   1.56508,       NA,       NA,  NA, 0)
expected_NA_zero_outcomes <- c(2.34035,   0.00000,       NA,  0.00000,  NA, 0)
expected_offset_outcomes  <- c(3.46410,   2.21336,  3.63424,  2.00000, NaN, 0)

# lapply(test_data, FUN = function(x){
# list(exp(sum(log(x+1), na.rm = T)/length(na.omit(x))),
#      gm_mean(x, na.rm = T, offset = 1))
# })

check_gm_results <- function(...){
  outcomes <- round(unlist(lapply(test_data, gm_mean, conf.level = NULL, ...)), digits = 5)
  names(outcomes) <- NULL
  outcomes
}

test_that("unexpected gm_mean results", {
  expect_equal(check_gm_results(na.rm = TRUE, zero.propagate = FALSE),
    round(expected_outcomes, digits = 5))
})


test_that("failed gm_mean results with na.rm=FALSE", {
  expect_equal(check_gm_results(na.rm = FALSE, zero.propagate = FALSE),
    round(expected_na_outcomes, digits = 5))
})

test_that("failed gm_mean results with zero propagate", {
  expect_equal(check_gm_results(na.rm = FALSE, zero.propagate = TRUE),
    round(expected_NA_zero_outcomes, digits = 5))
})

test_that("unexpected gm_mean results with offset = 1", {
  expect_equal(check_gm_results(na.rm = TRUE, offset = 1),
    round(expected_offset_outcomes, digits = 5))
})

##### Test outcome formatting with compare_summ_outcome_format():

# if necessary, loop over the list so that new data is generated:
# for(x in 1:100){
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
  form_outcome <- format_values(gm_mean(i), na.rm = TRUE, n.digits = 4)
  exponent <- stringr::str_extract(form_outcome, "(?<=\\*10\\^)[[1-9]]")
  if(is.na(exponent)) exponent <- "0"
  exponent <- as.numeric(exponent)
  formatted_vals <- strsplit(form_outcome, split = "\\$|\\(|-|\\)")[[1]]
  formatted_vals <- gsub("\\*10\\^[1-9]", "", formatted_vals) # remove exponent
  formatted_vals <- as.numeric(formatted_vals[formatted_vals != ""])*10^exponent
  formatted_vals <- as.numeric(formatC(formatted_vals, digits = 5)) # correct n.digits
  # Rounding errors can occur. Therefore, only test for the first digit:
  list(round(gm_outcome, digits = 1),
       round(formatted_vals, digits = 1))
}

test_that("formatting gives same values", {
  for(i in data_to_test){
    # i <- data6
    outcomes <- compare_summ_outcome_format(i)
    testthat::expect_equal(outcomes[[1]], outcomes[[2]])
  }
})
# } #### bracket ends loop if used
# if an error occurs, check it with this function:
# purrr::map(data_to_test, gm_mean, conf.level = 0.95, na.rm = TRUE)
# purrr::map(data_to_test, format_summary_values, n.digits = 4, conf.level = 0.95, na.rm = TRUE)

# Data with negative numbers should fail:
fail_data <- c((rnorm(n = 100, mean = 1, sd = 2)), 0, 0, -1)
test_that("Error with negative values", {
  expect_warning(gm_mean(fail_data, na.rm = TRUE))
})

test_that("expect failure with negative values", {
  expect_error(suppressWarnings(format_values(gm_mean(fail_data, na.rm = TRUE))))
})

# create dummy data, with some variables being, just for illustration of this
# function:
test_that("unexpected format output", {
  test_values <- immune_data$Neutrophils*10^5
  expect_equal(
    format_values(gm_mean(test_values), n.digits = 2),
    "$3.40(4.53-9.77)*10^5$")
  expect_equal(
    format_values(median_iqr(test_values), n.digits = 2),
    "$1.13(3.37)*10^6$")
  expect_equal(
    format_values(mean_sd(test_values), n.digits = 2),
    "$1.75(1.75)*10^6$")
})
