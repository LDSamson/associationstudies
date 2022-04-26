
library(dplyr)
immune_data %>%
  group_by(Sex) %>%
  summarize(Neutrophils = geomean_conf(Neutrophils))
geomean_conf()


test_that("Ouput is not a data frame", {
  testthat::expect_equal(
    is.data.frame(test_association(immune_data, "Tregs", "Frailty.index")),
    TRUE
  )
})

x <- data.frame(
  x1 = c(1,3,2,5),
  x2 = c(1,3,2,0),
  x3 = c(1,3,NA,5),
  x4 = c(1,3,NA, 0),
  x5 = c(NA, NA, NA, NA),
  x6 = c(0,0,0,0)
)

a <- rnorm(n = 100, mean = 5, sd = 1)
b <- rnorm(n = 100, mean = 50, sd = 1)
gm_mean(a)
gm_mean(a, conf.level = 0.95)
sumval_tab_format(a, na.rm = T)
gm_mean(b)
gm_mean(b, conf.level = 0.95)
sumval_tab_format(b, na.rm = T)
