# library(thesisfunctions)

#View(immune_data)

a <- thesisfunctions::test_association(dataset = immune_data,
                                  response.var = "T cells CD4 naive",
                                  explanatory.var = "Frailty.index")

# data <- immune_data
#  <-
# global.rho <-  stats::cor(data[["T cells CD4 naive"]], data[["Frailty.index"]], method = "spearman")
#
