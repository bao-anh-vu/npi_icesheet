## Read surface elev change data from ITS_LIVE

library(ncdf4)
library(CFtime)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(lattice) 


setwd("~/SSA_model/CNN/real_data/")

data_dir <- "./data/"
surf_elev_data <- nc_open(paste0(data_dir, "surface_elev/ANT_G1920_GroundedIceHeight_v01.nc"))

## Print variables and their info
print(surf_elev_data)

## Basin boundary data
basin_data <- read_sf(paste0(data_dir, "/boundaries/Basins/Basins_Antarctica_v02.shp"))
# basin_imbie_data <- read_sf(paste0(data_dir, "/boundaries/Basins_IMBIE/Basins_IMBIE_Antarctica_v02.shp"))

glacier_basin <- ncvar_get(surf_elev_data, "glacier_basin")
# dim(glacier_basin)


## Get x, y coordinates and time units
x <- ncvar_get(surf_elev_data, "x")
y <- ncvar_get(surf_elev_data, "y")
time <- ncvar_get(surf_elev_data, "time")
tunits <- ncatt_get(surf_elev_data, "time", "units") # get time units
dim(time)

## Need to convert the time to something readable here
# decode time
cf <- CFtime(tunits$value, calendar = "proleptic_gregorian", time) # convert time to CFtime class
# cf
timestamps <- as_timestamp(cf) # get character-string times
# timestamps
time_cf <- CFparse(cf, timestamps) # parse the string into date components
head(time_cf)

time_cf$row <- 1:nrow(time_cf)
# ref_timept <- time_cf %>% filter(year == 2013, month == 12) %>% select(row) %>% as.numeric() ## December 16, 2013

# Get height change data
dname <- "height_change"
height_array <- ncvar_get(surf_elev_data, dname)
dlname <- ncatt_get(surf_elev_data, dname,"long_name")
dunits <- ncatt_get(surf_elev_data, dname, "units")
fillvalue <- ncatt_get(surf_elev_data, dname,"_FillValue")

dim(height_array)
height_array[height_array==fillvalue$value] <- NA

## Plot
# get a single slice or layer 
yr <- 2020 #1985
mth <- 12
timept <- time_cf %>% filter(year == yr, month == mth) %>% select(row) %>% as.numeric() #345
height_slice <- height_array[,, timept]

# quick map
# y <- sort(y, decreasing = FALSE)
# png("~/SSA_model/CNN/real_data/plots/surf_elev_change.png", width = 800, height = 800)
# image(x, y, height_slice, col = rev(brewer.pal(10, "RdBu")))
# dev.off()

grid <- expand.grid(x = x, y = y)
grid$height <- c(height_slice)

thwaites_bound <- basin_data %>% filter(NAME == "Thwaites")
thwaites_sp <- as(thwaites_bound, "Spatial")
pts <- SpatialPoints(grid[, c("x", "y")], proj4string = CRS(proj4string(thwaites_sp)))
ov <- over(pts, thwaites_sp)
removed_rows <- attr(na.omit(ov), "na.action") # remove NA rows as those are points that do not fall within the shapefile
retained_rows <- setdiff(1:nrow(grid), removed_rows)
thwaites_elev <- grid[retained_rows, ]

# min_change <- min(grid$height, na.rm = TRUE)
# max_change <- max(grid$height, na.rm = TRUE)

p <- ggplot(thwaites_elev) + 
  geom_point(aes(x = x, y = y, colour = height)) +
  scale_colour_distiller(palette = "RdBu", direction = 1,
                         limits = c(-50, 50), #c(min_change, max_change),
                         name = "Surface elev change (m)") +
  geom_sf(data = thwaites_bound, color = "black", fill = NA) +

  theme_bw() 

png(paste0("./plots/elev_change_thwaites_", yr, formatC(mth, width = 2, flag = "0"), ".png"), width = 800, height = 800)
print(p)
dev.off()


