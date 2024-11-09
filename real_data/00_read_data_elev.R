
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
library(parallel)

setwd("~/SSA_model/CNN/real_data/")

data_dir <- "./data"

## Basin shapefile
basin_data <- read_sf(paste0(data_dir, "/boundaries/Basins/Basins_Antarctica_v02.shp"))
thwaites_bound <- basin_data %>% filter(NAME == "Thwaites")
thwaites_points <- as.data.frame(st_coordinates(st_geometry(thwaites_bound)))

## Grounding line data
gl_thwaites <- readRDS(paste0(data_dir, "/gl_thwaites.rds"))

## Surface elevation change data
surf_elev_data <- nc_open(paste0(data_dir, "/surface_elev/ANT_G1920_GroundedIceHeight_v01.nc"))

## Print variables and their info
print(surf_elev_data)

## Get x, y coordinates and time units
x <- ncvar_get(surf_elev_data, "x")
y <- ncvar_get(surf_elev_data, "y")
time <- ncvar_get(surf_elev_data, "time")
tunits <- ncatt_get(surf_elev_data, "time", "units") # get time units
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
dname <- "height"
height_arr <- ncvar_get(surf_elev_data, dname)
# dlname <- ncatt_get(surf_elev_data, dname,"long_name")
# dunits <- ncatt_get(surf_elev_data, dname, "units")
fillvalue <- ncatt_get(surf_elev_data, dname,"_FillValue")
height_arr[height_arr==fillvalue$value] <- NA

dim(height_arr)
height_grid <- expand.grid(x = x, y = y)
height_grid$height <- as.vector(height_arr)

## Filter height data to Thwaites glacier
xmin <- min(x) # min(thwaites_points$X)
xmax <- max(thwaites_points$X)
ymin <- min(thwaites_points$Y)
ymax <- max(thwaites_points$Y)

# png(paste0("./plots/vel.png"), width = 800, height = 800)
# print(plot)
# dev.off()

x_thwaites <- x[x >= xmin & x <= xmax]
y_thwaites <- y[y >= ymin & y <= ymax]
thwaites_height <- height_grid %>% filter(x %in% x_thwaites & y %in% y_thwaites)

thwaites_elev_plot <- thwaites_height %>% ggplot() + 
        geom_point(aes(x = x, y = y, colour = height)) + 
        scale_colour_distiller(palette = "Blues", direction = 1, 
                                name = "Height (m)") + 
        # geom_point(data = thwaites_points, mapping = aes(x = X, y = Y), color = "black", size = 0.1) +
        geom_sf(data = thwaites_bound, color = "black", fill = NA) +
        theme_bw() #+     

png("./plots/surface_elev/height_map_thwaites.png", width = 800, height = 800)
print(thwaites_elev_plot)
dev.off()

# Get height change data
dname <- "height_change"
height_change_arr <- ncvar_get(surf_elev_data, dname)
# dlname <- ncatt_get(surf_elev_data, dname,"long_name")
# dunits <- ncatt_get(surf_elev_data, dname, "units")
fillvalue <- ncatt_get(surf_elev_data, dname,"_FillValue")

dim(height_array)
height_change_arr[height_change_arr==fillvalue$value] <- NA

time_2000 <- time_cf %>% filter(year >= 2000)
data_from_2000 <- height_change_arr[,, time_2000$row]
nyears <- 2020-2000+1
year_ind <- split(time_2000$row, cut(seq_along(time_2000$row), nyears, labels = FALSE))

yearly_height_change <- lapply(year_ind, function(x) {
    height_change_arr[,, x]
})

# dim(yearly_height[[1]])
# test <- apply(yearly_height[[1]], 1:2, mean, na.rm = TRUE)

height_change_list <- mclapply(yearly_height_change, function(x) {
    apply(x, 1:2, mean, na.rm = TRUE)
}, mc.cores = 12L)

# saveRDS(height_list, file = paste0(data_dir, "/surface_elev/height_change_thwaites.rds"))

for (year in 1:length(year_ind)) {
    height_grid[, ncol(height_grid) + 1] <- as.vector(height_change_list[[year]])
}
colnames(height_grid)[4:ncol(height_grid)] <- paste0("change", 2000:2020)

test <- na.omit(height_grid)

height_change_thwaites <- height_grid %>% filter(x %in% x_thwaites & y %in% y_thwaites)

year <- 2020
var <- paste0("change", year)

height_change_plot <- height_change_thwaites %>% 
        ggplot() + geom_point(aes(x = x, y = y, colour = .data[[var]])) + 
        scale_colour_distiller(palette = "RdBu", direction = 1, 
                                limits = c(-40, 40),
                                name = "Height change (m)") + 
        geom_point(data = thwaites_points, mapping = aes(x = X, y = Y), color = "black", size = 0.1) +
        # geom_sf(data = thwaites_bound, color = "black", fill = NA) +
        coord_equal() + 
        theme_bw() #+

png(paste0("./plots/surface_elev/height_change_thwaites", year, ".png"), width = 800, height = 800)
print(height_change_plot)
dev.off()

# saveRDS(height_change_thwaites, file = paste0(data_dir, "/surface_elev/height_change_thwaites.rds"))
add_height <- function(change, height) {
    height + change
}
## Now recover absolute height data
height_thwaites <- height_change_thwaites %>% 
        mutate_at(vars(starts_with("change")), ~ add_height(., height))

height_thwaites <- rename_with(
  height_thwaites,
#   ~ paste0("height", .x),
  ~ sub("change", "height", .x),
  starts_with("change")
)

# saveRDS(height_thwaites, file = paste0(data_dir, "/surface_elev/height_thwaites.rds"))


height_plot <- height_thwaites %>% ggplot() + geom_point(aes(x = x, y = y, colour = .data[[var]])) + 
        scale_colour_distiller(palette = "Blues", direction = 1, 
                                name = "Height (m)") + 
        geom_point(data = thwaites_points, mapping = aes(x = X, y = Y), color = "black", size = 0.1) +
        # geom_sf(data = thwaites_bound, color = "black", fill = NA) +
        coord_equal() + 
        theme_bw() #+


png(paste0("./plots/surface_elev/height_thwaites_", year, ".png"), width = 800, height = 800)
print(height_plot)
dev.off()

