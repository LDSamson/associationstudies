## code to prepare `immune_data` dataset goes here

# I use some data from the following package: https://doi.org/10.3389/fimmu.2021.670070
# Furthermore, I add some variables with random data
set.seed(2022)
immune_data <- head(readRDS(url("https://figshare.com/ndownloader/files/25995722")), n=100)
immune_data$Sex <- sample(c("Men", "Women"), size = 100, replace = TRUE)
immune_data$Frailty.index <- runif(100, min = 0, max = 0.7)
immune_data$Batch <- sample(1:3, size = 100, replace = TRUE)
immune_data$`Macrophages M0` <- NULL

usethis::use_data(immune_data, overwrite = TRUE)
summary(immune_data)
