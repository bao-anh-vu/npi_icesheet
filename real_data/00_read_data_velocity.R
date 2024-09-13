# Read velocity data

## Read surface elev change data from ITS_LIVE

library(ncdf4)
library(CFtime)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(lattice) 
library(sf) # for shapefiles
# library(sfheaders)
library(sp)
# library(raster)

setwd("~/SSA_model/CNN/real_data/")

reread_data <- F

data_dir <- "./data"

## Basin shapefile
# unzip("./boundaries/Basins_Antarctica_v02.zip", junkpaths = FALSE)
basin_data <- read_sf(paste0(data_dir, "/boundaries/Basins/Basins_Antarctica_v02.shp"))

# basin_plot <- ggplot(basin_data) +
#   geom_sf(fill = "#69b3a2", color = "white") +
#   theme_bw()

# png(paste0("./plots/basin_plot.png"), width = 800, height = 800)
# print(basin_plot)
# dev.off()

thwaites_bound <- basin_data %>% filter(NAME == "Thwaites")

# ggplot(thwaites_bound) +
#   geom_sf(color = "black", fill = NA) +
#   theme_bw()


if (reread_data) {
  ## Velocity data
  vel_data <- nc_open(paste0(data_dir, "/velocity/ITS_LIVE_velocity_120m_RGI19A_0000_v02.nc"))
  ## Print variables and their info
  print(vel_data)
  names(vel_data$var)
  # attributes(vel_data)$names

  ## Velocity data at each grid point
  x <- ncvar_get(vel_data, "x")
  y <- ncvar_get(vel_data, "y")
  vx <- ncvar_get(vel_data, "vx")
  vy <- ncvar_get(vel_data, "vy")
  v <- ncvar_get(vel_data, "v")
  grid <- expand.grid(x = x, y = y)
  grid$v <- c(v)
  grid$vx <- c(vx)
  grid$vy <- c(vy)

  # plot_vel_ant <- ggplot(grid) + 
  # geom_point(aes(x = x, y = y, colour = v)) +
  # scale_colour_distiller(palette = "Reds", direction = 1,
  #                       #  limits = c(min_vel, max_vel),
  #                        name = "Velocity (m/a)") +
  # # geom_sf(data = basin_data, color = "black", fill = NA) +
  # theme_bw()

  # print("Saving velocity plot...")
  # png(paste0("./plots/vel_ant.png"), width = 800, height = 800)
  # print(plot_vel_ant)
  # dev.off()

  ## Select only grid points that fall within Thwaites glacier
  thwaites_points <- as.data.frame(st_coordinates(st_geometry(thwaites_bound)))
  # thwaites_wgs <- st_transform(thwaites_bound, crs = 4326)

  # rast <- raster(extent(thwaites_points[, c(X, Y)])) #, ncol = 5, nrow = 5)
  xmin <- min(x) #min(thwaites_points$X)
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
  grid_thwaites <- grid %>% filter(x %in% x_thwaites & y %in% y_thwaites)

  saveRDS(grid_thwaites, file = "./data/vel_thwaites.rds")

} else {
  grid_thwaites <- readRDS("./data/vel_thwaites.rds")
}


## Plot velocity data
plot_thwaites <- ggplot(grid_thwaites) + 
  geom_point(aes(x = x, y = y, colour = v)) +
  scale_colour_distiller(palette = "Reds", direction = 1,
                        #  limits = c(min_vel, max_vel),
                         name = "Velocity (m/a)") +
  geom_point(data = gl_thwaites, aes(x = X, y = Y), color = "salmon", size = 0.5) +
  # geom_sf(data = gl_thwaites_sf, color = "salmon", fill = NA) +
  geom_sf(data = thwaites_bound, color = "black", fill = NA) +
  theme_bw()

png(paste0("./plots/vel_thwaites.png"), width = 800, height = 800)
print(plot_thwaites)
dev.off()

# Select only points within Thwaites glacier

# thwaites_sp <- as(thwaites_bound, "Spatial")
# pts <- SpatialPoints(grid_thwaites[, c("x", "y")], proj4string = CRS(proj4string(thwaites_sp)))
# # ov <- over(pts, as(city_bdry, "SpatialPolygons"))
# ov <- over(pts, thwaites_sp)
# removed_rows <- attr(na.omit(ov), "na.action") # remove NA rows as those are points that do not fall within the shapefile
# retained_rows <- setdiff(1:nrow(grid_thwaites), removed_rows)
# thwaites_vel <- grid_thwaites[retained_rows,]

saveRDS(grid_thwaites, file = "./data/thwaites_vel.rds")





