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
library(qs)
# library(raster)

setwd("~/SSA_model/CNN/real_data/")

reread_data <- T

data_dir <- "./data/"

## Basin shapefile
# unzip("./boundaries/Basins_Antarctica_v02.zip", junkpaths = FALSE)
basin_data <- read_sf(paste0(data_dir, "boundaries/Basins/Basins_Antarctica_v02.shp"))

# png(paste0("./plots/basin_plot.png"), width = 800, height = 800)
# print(basin_plot)
# dev.off()

thwaites_bound <- basin_data %>% filter(NAME == "Thwaites")

# ggplot(thwaites_bound) +
#   geom_sf(color = "black", fill = NA) +
#   theme_bw()

year <- "0000" # "0000" for the composite velocity data

if (reread_data) {
  ## Velocity data
  vel_data <- nc_open(paste0(data_dir, "velocity/ITS_LIVE_velocity_120m_RGI19A_", year, "_v02.nc"))
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
  v_error <- ncvar_get(vel_data, "v_error")
  grid <- expand.grid(x = x, y = y)
  grid$v <- c(v)
  grid$vx <- c(vx)
  grid$vy <- c(vy)
  grid$v_error <- c(v_error)

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

  qsave(grid_thwaites, file = "./data/velocity/vel_thwaites.qs")

} else {
  grid_thwaites <- qread("./data/velocity/vel_thwaites.qs")
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

qsave(grid_thwaites, file = "./data/thwaites_vel.qs")





