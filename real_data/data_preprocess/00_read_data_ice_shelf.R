## Read ice shelf data

library(ncdf4)
library(CFtime)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(lattice) 
library(sf) # for shapefiles
# library(sfheaders)
library(sp)
library(parallel)

setwd("~/SSA_model/CNN/real_data/")

data_dir <- "./data"

## Basin shapefile
basin_data <- read_sf(paste0(data_dir, "/boundaries/Basins/Basins_Antarctica_v02.shp"))
thwaites_bound <- basin_data %>% filter(NAME == "Thwaites")
thwaites_points <- as.data.frame(st_coordinates(st_geometry(thwaites_bound)))

## Grounding line data
gl_thwaites <- readRDS(paste0(data_dir, "/gl_thwaites.rds"))

## Ice shelf elevation change data
shelf_elev_data <- nc_open(paste0(data_dir, "/surface_elev/ice_shelf/NSIDC-0792_19920317-20171216_V01.0.nc"))


sink(paste0(data_dir, "/surface_elev/ice_shelf/ice_shelf_data_desc.txt"))
print(shelf_elev_data)
sink()

## Get x, y coordinates and time units
## Get x, y coordinates and time units
x <- ncvar_get(shelf_elev_data, "x")
y <- ncvar_get(shelf_elev_data, "y")
time <- ncvar_get(shelf_elev_data, "time")
tunits <- ncatt_get(shelf_elev_data, "time", "units") # get time units
# dim(time)

## Need to convert the time to something readable here
# decode time
cf <- CFtime(tunits$value, calendar = "proleptic_gregorian", time) # convert time to CFtime class
# cf
timestamps <- as_timestamp(cf) # get character-string times
# timestamps
time_cf <- CFparse(cf, timestamps) # parse the string into date components
head(time_cf)

time_cf$row <- 1:nrow(time_cf)

## Get height data
dname <- "height_change"
height_change_arr <- ncvar_get(shelf_elev_data, dname)
fillvalue <- ncatt_get(shelf_elev_data, dname,"_FillValue")
height_change_arr[height_change_arr==fillvalue$value] <- NA

time <- 1
shelf_grid <- expand.grid(x = x, y = y)
shelf_grid$height_change <- as.vector(height_change_arr[,, time])

thwaites_points <- as.data.frame(st_coordinates(st_geometry(thwaites_bound)))
# thwaites_wgs <- st_transform(thwaites_bound, crs = 4326)

# rast <- raster(extent(thwaites_points[, c(X, Y)])) #, ncol = 5, nrow = 5)
xmin <- min(x) #min(thwaites_points$X)
xmax <- max(thwaites_points$X)
ymin <- min(thwaites_points$Y)
ymax <- max(thwaites_points$Y)

thwaites_shelf <- shelf_grid %>%
  filter(x >= xmin & x <= xmax & y >= ymin & y <= ymax)

shelf_plot <- thwaites_shelf %>% ggplot() + 
  geom_point(aes(x = x, y = y, colour = height_change)) +
  scale_colour_distiller(palette = "RdBu", direction = 1, limits = c(-20, 20),
                                name = "Surface elevation change (m)") + 
  geom_point(data = gl_thwaites, aes(x = X, y = Y), color = "cyan", size = 0.05) +
  geom_sf(data = thwaites_bound, color = "black", fill = NA) +
  # geom_point(data = thwaites_points, aes(x = X, y = Y), color = "black", size = 0.01) +
  theme_bw() + 
  # coord_fixed() +
  labs(title = "Ice shelf elevation change", x = "Longitude", y = "Latitude", fill = "Elevation change (m)")

png(paste0("./plots/surface_elev/ice_shelf_elev_change.png"), width = 800, height = 800)
print(shelf_plot)
dev.off()

## Thickness data
dname <- "thickness"
thickness_arr <- ncvar_get(shelf_elev_data, dname)
fillvalue <- ncatt_get(shelf_elev_data, dname,"_FillValue")
thickness_arr[height_change_arr==fillvalue$value] <- NA

shelf_grid$thickness <- as.vector(thickness_arr[,, time])

thickness_shelf <- shelf_grid %>%
  filter(x >= xmin & x <= xmax & y >= ymin & y <= ymax)

thickness_shelf_plot <- test %>% ggplot() + 
  geom_point(aes(x = x, y = y, colour = thickness)) +
  scale_colour_distiller(palette = "RdBu", direction = 1, limits = c(-20, 20),
                                name = "Thickness (m)") + 
  # geom_point(data = gl_thwaites, aes(x = X, y = Y), color = "cyan", size = 0.05) +
  geom_sf(data = thwaites_bound, color = "black", fill = NA) +
  # geom_point(data = thwaites_points, aes(x = X, y = Y), color = "black", size = 0.01) +
  theme_bw() + 
  # coord_fixed() +
  labs(title = "Ice shelf thickness", x = "Longitude", y = "Latitude", fill = "Thickness (m)")

png(paste0("./plots/surface_elev/shelf_thickness.png"), width = 800, height = 800)
print(thickness_shelf_plot)
dev.off()
