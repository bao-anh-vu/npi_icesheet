## Combine velocity data
setwd("~/SSA_model/CNN/real_data/")

# library(ncdf4)
# library(CFtime)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(qs)
library(abind)
# library(lattice)
# library(sf) # for shapefiles
# library(sfheaders)
# library(sp)
# library(parallel)

years <- 2000:2020
data_dir <- "./data"

## Missing pattern for surface elevation data
# year <- 2000
surf_elev_mat <- matrix(NA, nrow = J, ncol = length(years))
for (i in 1:length(years)) {
    year <- years[i]
    surf_elev_data <- readRDS(file = paste0(data_dir, "/surface_elev/surf_elev_", year, ".rds"))
    surf_elev_mat[, i] <- unlist(surf_elev_data)[1:J]
}

qsave(surf_elev_mat, "./data/surface_elev/surf_elev_mat.qs")

png("./plots/surface_elev/surf_elev_st_plot.png", width = 800, height = 600)
image(surf_elev_mat)
dev.off()

surf_ev_missing_pattern <- ifelse(is.na(surf_elev_mat), 0, 1)
qsave(surf_ev_missing_pattern, "./data/surface_elev/missing_pattern.qs")
