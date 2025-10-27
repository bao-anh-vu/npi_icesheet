## Read Bedmap data in shapefiles

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(sp)
library(sf)
library(parallel)
library(qs)

setwd("~/SSA_model/CNN/real_data/")

print("Reading BedMap data...")
data_dir <- "./data/"

## Basin data
basin_data <- read_sf(paste0("./data/boundaries/Basins/Basins_Antarctica_v02.shp"))
thwaites_bound <- basin_data %>% filter(NAME == "Thwaites")

## Grounding line data
gl_thwaites <- qread(paste0(data_dir, "grounding_line/gl_thwaites.qs"))

## Flowline
flowline <- qread(paste0(data_dir, "flowline_regrid.qs"))

## Mark GL position on the flowline
gl_pos <- qread(file = paste0(data_dir, "grounding_line/gl_pos.qs"))

## Bed data
bas_data_2018 <- read_sf(paste0(data_dir, "bedmap/BAS_2018_Thwaites_AIR_BM3/BAS_2018_Thwaites_AIR_BM3_points.shp"))
bas_data_2019 <- read_sf(paste0(data_dir, "bedmap/BAS_2019_Thwaites_AIR_BM3/BAS_2019_Thwaites_AIR_BM3_points.shp"))
# head(bas_data_2018$geometry)

bed_data_2018 <- bas_data_2018 %>% select(Mean_bed, SD_bed, geometry) %>% mutate(year = 2018)
bed_data_2019 <- bas_data_2019 %>% select(Mean_bed, SD_bed, geometry) %>% mutate(year = 2019)

bed_data <- rbind(bed_data_2018, bed_data_2019)

bed_geom <- as.data.frame(st_coordinates(st_geometry(bed_data$geometry)))
# geom_2019 <- as.data.frame(st_coordinates(st_geometry(bed_data_2019$geometry)))

bed_data$X <- bed_geom$X
bed_data$Y <- bed_geom$Y
# bed_data_2019$X <- geom_2019$X
# bed_data_2019$Y <- geom_2019$Y
bed_data$Mean_bed[bed_data$Mean_bed == -9999] <- NA
bed_data$SD_bed[bed_data$SD_bed == -9999] <- NA

## Plot bed geom
bed_plot <- ggplot() +
#   geom_sf(data = bas_data_2018, aes(color = Mean_bed)) +
  geom_sf(data = bed_data, aes(color = Mean_bed)) +
  geom_sf(data = thwaites_bound, color = "black", fill = NA) +
  scale_color_viridis_c() +
  geom_line(data = flowline, mapping = aes(x = x, y = y), colour = "cyan") +
  theme_bw()

png("./plots/bed/bed_geom.png", width = 800, height = 800)
print(bed_plot)
dev.off()

bed_sd_plot <- ggplot() +
#   geom_sf(data = bas_data_2018, aes(color = Mean_bed)) +
  geom_sf(data = bed_data, aes(color = SD_bed)) +
  geom_sf(data = thwaites_bound, color = "black", fill = NA) +
  scale_color_viridis_c() +
  geom_line(data = flowline, mapping = aes(x = x, y = y), colour = "cyan") +
  theme_bw()

png("./plots/bed/bed_sd.png", width = 800, height = 800)
print(bed_sd_plot)
dev.off()

## Map bed data to flowline

delta <- 200 # grid size
# flowline_dist <- sqrt((flowline$x[2:nrow(flowline)] - flowline$x[1:(nrow(flowline) - 1)])^2 +
#     (flowline$y[2:nrow(flowline)] - flowline$y[1:(nrow(flowline) - 1)])^2)
# flowline$ind <- 1:nrow(flowline)
# bed_near_pts <- bed_data %>% filter(
#             x >= (gl_pos[1] - delta) & x <= (gl_pos[1] + delta),
#             y >= (gl_pos[2] - delta) & y <= (gl_pos[2] + delta)) #%>% 
            # mutate(dist = sqrt((x - gl_pos[1])^2 + (y - gl_pos[2])^2)) 
                
cat("Mapping bed elevation data to flowline...")
## Now map velocity to flowline

avg_nearest_four <- function(pos) {
    near_pts <- bed_data %>% filter(
        X >= (pos[1] - delta) & X <= (pos[1] + delta),
        Y >= (pos[2] - delta) & Y <= (pos[2] + delta)
    ) 

    if (nrow(near_pts) > 0) {
        near_vals <- near_pts %>%
        mutate(dist = sqrt((X - pos[1])^2 + (Y - pos[2])^2)) %>%
        filter(dist > 0) %>%
        arrange(dist) %>%
        slice_min(dist, n = 4) %>%
        select(Mean_bed, SD_bed) %>%
        summarise(bed_nearest = Mean_bed[1], bed_avg = mean(Mean_bed), 
                    sd_nearest = SD_bed[1], sd_avg = mean(SD_bed)) # , vx_avg = mean(vx), vy_avg = mean(vy)) # %>% # take average vx and vx as vx and vy at current position
    # as.numeric()
    } else {
        near_vals <- data.frame(bed_nearest = NA, bed_avg = NA,
                                sd_nearest = NA, sd_avg = NA)   
    }

    return(near_vals)
}



# delta <- 500 # grid size
# flowline <- na.omit(flowline)
flowline_pos <- lapply(1:nrow(flowline), function(i) as.numeric(flowline[i, ]))
t12 <- system.time({
    bed_elev <- lapply(flowline_pos, avg_nearest_four) #, mc.cores = 10L)
})

m <- sapply(bed_elev, function(x) nrow(x))

bed_elev_nearest <- sapply(bed_elev, function(x) x$bed_nearest)
bed_elev_avg <- sapply(bed_elev, function(x) x$bed_avg)
bed_sd_nearest <- sapply(bed_elev, function(x) x$sd_nearest)
bed_sd_avg <- sapply(bed_elev, function(x) x$sd_avg)

png(paste0("./plots/bed/bed_elev_flowline_new.png"), width = 800, height = 800)
plot(unlist(bed_elev_avg), type = "l")
dev.off()

qsave(bed_elev_avg, file = "./data/bedmap_obs.qs")
qsave(bed_sd_avg, file = "./data/bedmap_sd.qs")