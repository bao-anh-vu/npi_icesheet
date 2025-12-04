## RACMO surface mass balance data

## Read SMB data

setwd("~/SSA_model/CNN/real_data/")

library(ncdf4)
library(proj4)
library(CFtime)
library(sf)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(parallel)
library(qs)

data_dir <- "./data"

racmo_data <- nc_open(paste0(data_dir, "/SMB/SMB_RACMO2-3p2_yearly_ANT27_1979_2016.nc"))
print(racmo_data)

smb_kg_per_m2 <- ncvar_get(racmo_data, "smb")
lon <- ncvar_get(racmo_data, "lon")
lat <- ncvar_get(racmo_data, "lat")

## Convert from kg m^-2 s^-1 to m a^-1
# secpera <- 31556926 #seconds per annum
ice_density <- 917 # kg m^-3
smb <- smb_kg_per_m2 / ice_density #* secpera

# Extract rotated coordinates
# rot_lats <- ncvar_get(racmo_data, "rlat") # Replace 'rlat' with your variable name
# rot_lons <- ncvar_get(racmo_data, "rlon") # Replace 'rlon' with your variable name

# Extract pole information (adjust names as needed)
# pole_lon <- ncvar_get(racmo_data, "rotated_pole") 
# pole_lon <- ncatt_get(racmo_data, 0, "rotated_pole")$value
# pole_lat <- ncatt_get(racmo_data, 0, "pole_lat")$value

# # Create a projection string for the rotated pole system
# rotated_proj <- "-m 57.295779506 +proj=ob_tran +o_proj=latlon +o_lat_p=-180.0 +lon_0=10.0" # from the data description
# standard_proj <- "+proj=longlat +datum=WGS84"

# # Create a meshgrid of rotated coordinates
# rot_lat_grid <- as.vector(rot_lats)
# rot_lon_grid <- as.vector(rot_lons)

# # Transform coordinates
# coords <- expand.grid(rot_lon_grid, rot_lat_grid)
# transformed_coords <- project(as.matrix(coords), rotated_proj, inverse = TRUE)

# Reshape transformed coordinates back to original grid
# nlat <- length(rot_lats)
# nlon <- length(rot_lons)
# lons_standard <- matrix(transformed_coords[, 1], nrow = nlat, ncol = nlon)
# lats_standard <- matrix(transformed_coords[, 2], nrow = nlat, ncol = nlon)

## Read timestamps
time <- ncvar_get(racmo_data, "time_bnds")
tunits <- ncatt_get(racmo_data, "time_bnds", "units") # get time units

## Need to convert the time to something readable here
# decode time
cf <- CFtime(tunits$value, calendar = "proleptic_gregorian", time[2, ]) # convert time to CFtime class
# cf
timestamps <- as_timestamp(cf) # get character-string times
# timestamps
time_cf <- CFparse(cf, timestamps) # parse the string into date components
head(time_cf)


## Get x, y coordinates
# smb1 <- smb[,,1]

## Reshape smb and lon, lat into a dataframe
# smb_df <- data.frame(lon = as.vector(lon), lat = as.vector(lat), smb = as.vector(smb1) / 1000)

# Convert to sf object in geographic (WGS84) coordinates
lonlat <- data.frame(lon = as.vector(lon), lat = as.vector(lat))
points_sf <- st_as_sf(lonlat, coords = c("lon", "lat"), crs = 4326)  # EPSG:4326 = WGS84

# Define the polar stereographic projection (EPSG:3413 for North Pole)
# Use EPSG:3031 for South Polar Stereographic projection if working with the South Pole
points_projected <- st_transform(points_sf, crs = 3031)  # EPSG:3413 = North Polar Stereographic

# Extract the transformed coordinates
polar_coordinates <- st_coordinates(points_projected)
# head(polar_coordinates)

# smb_df_xy <- data.frame(x = polar_coordinates[, 1], y = polar_coordinates[, 2], smb = smb_df$smb)
# head(smb_df_xy)

## Subset the data to the Thwaites region
basin_data <- read_sf(paste0(data_dir, "/boundaries/Basins/Basins_Antarctica_v02.shp"))
thwaites_bound <- basin_data %>% filter(NAME == "Thwaites")
thwaites_points <- as.data.frame(st_coordinates(st_geometry(thwaites_bound)))

# rast <- raster(extent(thwaites_points[, c(X, Y)])) #, ncol = 5, nrow = 5)
xmin <- min(thwaites_points$X)
xmax <- max(thwaites_points$X)
ymin <- min(thwaites_points$Y)
ymax <- max(thwaites_points$Y)

## Map SMB data to flowline
flowline <- readRDS(paste0(data_dir, "/flowline_regrid.rds"))
# flowline <- na.omit(flowline)
flowline_pos <- lapply(1:nrow(flowline), function(i) as.numeric(flowline[i, ]))

avg_nearest_four <- function(pos, grid_size) {
            
    delta <- grid_size
    
    near_pts <- smb_thwaites %>% filter(
        x >= (pos[1] - delta) & x <= (pos[1] + delta),
        y >= (pos[2] - delta) & y <= (pos[2] + delta)
    )
    near_v <- near_pts %>%
        mutate(dist = sqrt((x - pos[1])^2 + (y - pos[2])^2)) %>%
        filter(dist > 0) %>%
        arrange(dist) %>%
        slice_min(dist, n = 4) %>%
        select(smb) %>%
        summarise(nearest = smb[1], avg = mean(smb)) 
    # as.numeric()

    return(near_v)
}


flowline_smb_annual <- matrix(NA, nrow = dim(smb)[3], ncol = nrow(flowline))
for (t in 1:dim(smb)[3]) {
    cat("t = ", t, "\n")
    smb_t <- smb[,,t]
    smb_df <- data.frame(x = polar_coordinates[, 1], y = polar_coordinates[, 2], smb = as.vector(smb_t))

    ## Try filtering grid points in v with coordinates in x_thwaites and y_thwaites
    smb_thwaites <- smb_df %>% filter(x >= xmin & x <= xmax & y >= ymin & y <= ymax)

    ## Map SMB data to flowline
    smb_mapped <- mclapply(flowline_pos, avg_nearest_four, grid_size = 27000, mc.cores = 10L)

    # smb_nearest <- sapply(smb_mapped, function(x) x$v_nearest)
    smb_avg <- sapply(smb_mapped, function(x) x$avg)

    flowline_smb_annual[t, ] <- smb_avg
    # flowline$smb_nearest <- smb_nearest
    # flowline$smb_avg <- smb_avg

}

## Might need to average SMB data over time?
qsave(flowline_smb_annual, file = paste0(data_dir, "/SMB/flowline_landice_smb.qs"))


## Plot the data
smb_plot <- ggplot() +
  geom_sf(data = thwaites_bound, fill = NA, color = "black") +
  geom_point(data = smb_thwaites, aes(x = x, y = y, color = smb), size = 1) +
  scale_colour_gradient(low = "blue", high = "red") +
  theme_bw()

png(paste0("./plots/smb_plot.png"), width = 800, height = 800)
print(smb_plot)
dev.off()

## Map SMB data to flowline



# # t12 <- system.time({
    