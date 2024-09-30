## Map surface elevation data to flowline

# library(ncdf4)
# library(CFtime)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(lattice)
library(sf) # for shapefiles
# library(sfheaders)
library(sp)
library(parallel)

# library(raster)

setwd("~/SSA_model/CNN/real_data/")

reread_data <- T
# remap_vel <- T

data_dir <- "./data"

## Basin shapefilethe7
# unzip("./boundaries/Basins_Antarctica_v02.zip", junkpaths = FALSE)
basin_data <- read_sf(paste0(data_dir, "/boundaries/Basins/Basins_Antarctica_v02.shp"))
thwaites_bound <- basin_data %>% filter(NAME == "Thwaites")

## Grounding line data
gl_thwaites <- readRDS(paste0(data_dir, "/gl_thwaites.rds"))

surf_elev_data <- readRDS(file = paste0(data_dir, "/surface_elev/height_thwaites.rds"))

## Map surface elevation data to flowline

avg_nearest_four <- function(pos, values) {
    near_pts <- values %>% filter(
        x >= (pos[1] - delta) & x <= (pos[1] + delta),
        y >= (pos[2] - delta) & y <= (pos[2] + delta)
    )

    near_height <- near_pts %>%
        mutate(dist = sqrt((x - pos[1])^2 + (y - pos[2])^2)) %>%
        filter(dist > 0) %>%
        arrange(dist) %>%
        slice_min(dist, n = 4) %>%
        # select(all_of(var)) %>%
        summarise(height_avg = mean(.data[[var]])) # , vx_avg = mean(vx), vy_avg = mean(vy)) # %>% # take average vx and vx as vx and vy at current position
    # as.numeric()

    return(near_height)
}

years <- 2001:2020
# year <- 2019
for (year in years) {
    var <- paste0("height", year)
    print("Mapping surface elevation data to flowline...")

    delta <- 1920 # grid size
    flowline <- readRDS(paste0(data_dir, "/flowline_regrid.rds"))
    flowline_pos <- lapply(1:nrow(flowline), function(i) as.numeric(flowline[i, ]))

    gl <- readRDS(paste0(data_dir, "/gl_thwaites.rds"))

    annual_surf_elev <- surf_elev_data %>% select(x, y, all_of(var))

    surf_elev_plot <- annual_surf_elev %>% ggplot() + 
                        geom_point(aes(x = x, y = y, color = .data[[var]])) + 
                        scale_colour_distiller(palette = "Blues", direction = 1, 
                                    name = "Height (m)") + 
                        geom_sf(data = thwaites_bound, color = "black", fill = NA) +
                        geom_point(data = gl, mapping = aes(x = X, y = Y), color = "salmon", size = 0.05) +
                        geom_point(data = flowline, mapping = aes(x = x, y = y), color = "black", size = 0.1) + 
                        theme_bw()

    # png(paste0("./plots/surface_elev/surf_elev_", year, ".png"), width = 800, height = 500)
    # print(surf_elev_plot)
    # dev.off()

    surf_elev <- lapply(flowline_pos, avg_nearest_four, values = annual_surf_elev) #, mc.cores = 12L)

    ## where's the ice shelf??

    png(paste0("./plots/elev_flowline_1d_", year, ".png"), width = 800, height = 500)
    plot(unlist(surf_elev), type = "l")
    dev.off()

    saveRDS(surf_elev, paste0(data_dir, "/surface_elev/surf_elev_", year, ".rds"))
}



# pos <- flowline_pos[[1]]
# near_pts <- surf_elev_data %>% filter(
#     x >= (pos[1] - delta) & x <= (pos[1] + delta),
#     y >= (pos[2] - delta) & y <= (pos[2] + delta)
# )

#     near_height <- near_pts %>%
#         mutate(dist = sqrt((x - pos[1])^2 + (y - pos[2])^2)) %>%
#         filter(dist > 0) %>%
#         arrange(dist) %>%
#         slice_min(dist, n = 4) %>%
#         select(.data[[var]]) %>%
#         summarise(height_avg = mean(.data[[var]])) # , vx_avg = mean(vx), vy_avg = mean(vy)) # %>% # take average vx and vx as vx and vy at current position
#     # as.numeric()

#     return(near_height)