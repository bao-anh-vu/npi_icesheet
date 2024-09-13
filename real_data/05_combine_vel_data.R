## Combine velocity data
setwd("~/SSA_model/CNN/real_data/")

# library(ncdf4)
# library(CFtime)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
# library(lattice)
# library(sf) # for shapefiles
# library(sfheaders)
# library(sp)
# library(parallel)

years <- 2000:2020
data_dir <- "./data"
velocities <- list()

for (i in 1:length(years)) {
    year <- years[i]
    velocity_i <- readRDS(paste0(data_dir, "/velocity/flowline_vel_mapped_", year, ".rds"))
    velocity_i$year <- year
    velocity_i$gridpt <- 1:nrow(velocity_i)
    velocities[[i]] <- velocity_i
}

## Concatenate all velocity data
velocities_df <- bind_rows(velocities)
tail(velocities_df)

## Plot velocity data over space-time
vel_st_plot <- ggplot(velocities_df) + 
    geom_tile(aes(x = gridpt, y = year, fill = vel_avg)) + 
    scale_fill_distiller(palette = "BuPu", direction = 1) +
    theme_bw() + 
    labs(x = "Grid Point", y = "Year", fill = "Velocity (m/yr)") + 
    ggtitle("Thwaites Glacier Velocity Over Time") + theme(plot.title = element_text(hjust = 0.5))

png("./plots/velocity/vel_st_plot.png", width = 800, height = 600)
print(vel_st_plot)
dev.off()
