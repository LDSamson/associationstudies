####################### colorblind friendly color palette: #####################

# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
save(cbbPalette, file= "data/cbbPalette.rda")
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
save(cbPalette, file= "data/cbPalette.rda")
