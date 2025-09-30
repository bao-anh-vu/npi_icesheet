## Read bedmachine data

setwd("~/SSA_model/CNN/real_data/")

rm(list = ls())
library(ncdf4)
library(sf)
library(dplyr)
library(qs)
library(parallel)

reread_data <- F

if (reread_data) {
    ## BedMachine data
    bedmachine <- nc_open("./data/bedmachine/BedMachineAntarctica-v3.nc")
    # print(bedmachine)

    x <- ncvar_get(bedmachine, "x")
    y <- ncvar_get(bedmachine, "y")
    bed <- ncvar_get(bedmachine, "bed")

    grid <- expand.grid(x = x, y = y)
    grid$bed <- c(bed)

    ## Basin shapefile
    basin_data <- read_sf(paste0("./data/boundaries/Basins/Basins_Antarctica_v02.shp"))
    thwaites_bound <- basin_data %>% filter(NAME == "Thwaites")

    ## Select only grid points that fall within Thwaites glacier
    thwaites_points <- as.data.frame(st_coordinates(st_geometry(thwaites_bound)))

    xmin <- min(x) # min(thwaites_points$X)
    xmax <- max(thwaites_points$X)
    ymin <- min(thwaites_points$Y)
    ymax <- max(thwaites_points$Y)

    # png(paste0("./plots/vel.png"), width = 800, height = 800)
    # print(plot)
    # dev.off()

    x_thwaites <- x[x >= xmin & x <= xmax]
    y_thwaites <- y[y >= ymin & y <= ymax]
    # grid_thwaites <- expand.grid(x = x_thwaites, y = y_thwaites)

    ## Try filtering grid points in v with coordinates in x_thwaites and y_thwaites
    thwaites_bed <- grid %>% filter(x %in% x_thwaites & y %in% y_thwaites)

    qsave(thwaites_bed, "./data/bedmachine/thwaites_bed.qs")
} else {
    thwaites_bed <- qread("./data/bedmachine/thwaites_bed.qs")
}

## Map bed observations to flowline
cat("Mapping BedMachine data to flowline... \n")
flowline <- qread("./data/flowline_regrid.qs")
flowline <- na.omit(flowline)
flowline_pos <- lapply(1:nrow(flowline), function(i) as.numeric(flowline[i, ]))

avg_nearest_four <- function(pos) {
    near_pts <- thwaites_bed %>% filter(
        x >= (pos[1] - delta) & x <= (pos[1] + delta),
        y >= (pos[2] - delta) & y <= (pos[2] + delta)
    )

    near_bed <- near_pts %>%
        mutate(dist = sqrt((x - pos[1])^2 + (y - pos[2])^2)) %>%
        filter(dist > 0) %>%
        arrange(dist) %>%
        slice_min(dist, n = 4) %>%
        select(bed) %>%
        summarise(bed_nearest = bed[1], bed_avg = mean(bed)) # , vx_avg = mean(vx), vy_avg = mean(vy)) # %>% # take average vx and vx as vx and vy at current position
    # as.numeric()

    return(near_bed)
}

delta <- 500 # grid size

t12 <- system.time({
    bedmachine_obs <- mclapply(flowline_pos, avg_nearest_four, mc.cores = 10L)
})

bedmachine_obs_nearest <- sapply(bedmachine_obs, function(x) x$bed_nearest)
bedmachine_obs_avg <- sapply(bedmachine_obs, function(x) x$bed_avg)
flowline$bed_nearest <- bedmachine_obs_nearest
flowline$bed_avg <- bedmachine_obs_avg

qsave(flowline, "./data/bedmachine/flowline_bedmachine.qs")

## Plot bedmachine data along flowline
png("./plots/temp/bedmachine_flowline.png")
plot(flowline$bed_nearest, type = 'l', lwd = 2, col = 'blue',
    ylab = "Bed elevation (m)", xlab = "Distance along flowline (m)")
lines(flowline$bed_avg, col = 'red', lwd = 2)
legend("topright", legend = c("Nearest", "Average of 4 nearest"), col = c("blue", "red"), lwd = 2)
dev.off()
