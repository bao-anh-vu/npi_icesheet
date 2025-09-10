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

years <- 2010:2020
data_dir <- "./data"

## Missing pattern for surface elevation data
# year <- 2000
flowline <- qread(paste0(data_dir, "/flowline_regrid.qs"))
J <- nrow(flowline) # number of grid points

surf_elev_mat <- matrix(NA, nrow = J, ncol = length(years))
for (i in 1:length(years)) {
    year <- years[i]
    # surf_elev_data <- readRDS(file = paste0(data_dir, "/surface_elev/surf_elev_", year, ".rds"))
    surf_elev_data <- qread(file = paste0(data_dir, "/surface_elev/surf_elev_", year, ".qs"))
    surf_elev_mat[, i] <- unlist(surf_elev_data)[1:J]
}

## There may be some mismatch between the surface elevation data and grounding line data
## Discard any observations that are beyond the grounding line
gl_pos <- qread(file = paste0(data_dir, "/grounding_line/gl_pos.qs"))
gl_ind <- gl_pos$ind
surf_elev_mat[(gl_ind+1):nrow(surf_elev_mat), ] <- NA

qsave(surf_elev_mat, "./data/surface_elev/surf_elev_mat.qs")

png("./plots/surface_elev/surf_elev_st_plot.png", width = 800, height = 600)
image(surf_elev_mat)
dev.off()

surf_ev_missing_pattern <- ifelse(is.na(surf_elev_mat), 0, 1)
qsave(surf_ev_missing_pattern, "./data/surface_elev/missing_pattern.qs")
