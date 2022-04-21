#functions to test:
library(tidyr)
library(dplyr)
library(purrr)
library(thesisfunctions)
cell_vars <- names(immune_data)[!names(immune_data) %in%
                                  c("Sex", "Frailty.index", "Batch")]

immune_data_long <- immune_data %>%
  pivot_longer(names_to = "cellnames", values_to = "value",
               -c("Sex", "Frailty.index", "Batch"))

a <- purrr::map_dfr(cell_vars, .f = ~association_study(
  dataset = immune_data_long,
  subset.parameter = .,
  subset.column = "cellnames",
  numerical.value = "value",
  comparison.value = "Frailty.index",
  stratum = c("Sex", "Batch"),
  n.resample = 10^4
)) %>%
  mutate(p.value = round(p.value, digits = 1)) %>%
  BH_selection()

plot_p_histogram(a)

table(data$cellnames)

association_study_long2 <- function(data, response.names, ...){
  data <- immune_data_long
  response.names <- "cellnames"
  x <- "Tregs"
  purrr::map_dfr(unique(data[[response.names]]), .f = function(x){
    df <- test_association(
      data[data[[response.names]] == x, ],
      response.var = "value",
      explanatory.var = "Frailty.index",
      stratum = c("Sex"))
    df["Response var"] <- x
    df
  })
}

# careful: if response.names variable is a factor, the names will not
# be displayed properly in the output data frame!
b <- association_study_long(immune_data_long,
                       response.names = "cellnames",
                       response.var = "value",
                       explanatory.var = "Frailty.index",
                       stratum = c("Sex", "Batch"),
                       n.resample = 10^4) %>%
  mutate(p.value = round(p.value, digits = 1)) %>%
  BH_selection()

c <- map_dfr(cell_vars, .f = ~test_association(
  immune_data,
  response.var = .,
  explanatory.var = "Frailty.index",
  stratum = c("Sex", "Batch"),
  n.resample = 10^4
)) %>%
  BH_selection()
