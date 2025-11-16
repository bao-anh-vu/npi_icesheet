## Read SMB data

setwd("~/SSA_model/CNN/real_data/")

library(sf)
library(ncdf4)
library(CFtime)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(parallel)
library(qs)
# library(fields)
# library(parallel)
# library(Matrix)
# library("qlcMatrix")
# library("fastmatrix")
# library("expm")
# library(R.utils)
# library("sp")
# library("tidyr")
# library(matrixStats) # for the rowMaxs() function
# library(mvtnorm)
# library(FRK)
# library(qs)

data_dir <- "./data/"

reread_height_change_data <- F
reread_smb_data <- F
reread_basal_melt_data <- F
reread_thickness_data <- T

## Basal melt data
basal_melt_data <- nc_open(paste0(data_dir, "basal_melt/NSIDC-0792_19920317-20171216_V01.0.nc"))
print(basal_melt_data)

height_change <- ncvar_get(basal_melt_data, "height_change")
smb <- ncvar_get(basal_melt_data, "smb")
melt <- ncvar_get(basal_melt_data, "melt")
mean_melt <- ncvar_get(basal_melt_data, "melt_mean")
x <- ncvar_get(basal_melt_data, "x")
y <- ncvar_get(basal_melt_data, "y")

thickness_bedmachine <- ncvar_get(basal_melt_data, "thickness")

## Read timestamps
time <- ncvar_get(basal_melt_data, "time")
tunits <- ncatt_get(basal_melt_data, "time", "units") # get time units

## Need to convert the time to something readable here
cf <- CFtime(tunits$value, calendar = "proleptic_gregorian", time) # convert time to CFtime class
timestamps <- as_timestamp(cf) # get character-string times
# timestamps
time_cf <- CFparse(cf, timestamps) # parse the string into date components
head(time_cf)
time_cf$ind <- 1:length(time_cf$year)

n_years <- max(time_cf$year) - min(time_cf$year)
years <- unique(time_cf$year)

annual_height_change_re2014 <- mclapply(0:n_years, function(i) { ## annual height change relative to 2014
    height_change_t <- height_change[, , (i * 4 + 1):((i + 1) * 4)]
    height_change_avg <- apply(height_change_t, c(1, 2), mean, na.rm = TRUE)
}, mc.cores = 10L)

## Annual change (each year relative to the previous year)
annual_height_change <- lapply(2:length(annual_height_change_re2014), function(i) {
    annual_height_change_re2014[[i]] - annual_height_change_re2014[[i - 1]]
})


## Grounding line data
gl_thwaites <- readRDS(paste0(data_dir, "/gl_thwaites.rds"))

## Restrict to Thwaites glacier
## Basin shapefile
basin_data <- read_sf(paste0(data_dir, "/boundaries/Basins/Basins_Antarctica_v02.shp"))
thwaites_bound <- basin_data %>% filter(NAME == "Thwaites")
thwaites_points <- as.data.frame(st_coordinates(st_geometry(thwaites_bound)))
xmin <- min(x) # min(thwaites_points$X)
xmax <- max(thwaites_points$X)
ymin <- min(thwaites_points$Y)
ymax <- max(thwaites_points$Y)

if (reread_height_change_data) {
    ## Create a data frame containing ice height_change data
    # fillvalue <- ncatt_get(basal_melt_data, dname,"_FillValue")
    # height_change[height_change==fillvalue$value] <- NA
    height_change_df <- data.frame(expand.grid(x = x, y = y))
    # thwaites_df <- xy_df %>% filter(x >= xmin & x <= xmax & y >= ymin & y <= ymax)

    for (t in 1:n_years) {
        varname <- paste0("height_change_", years[t])
        height_change_df[, varname] <- as.vector(annual_height_change[[t]])
    }

    # height_change_df <- data.frame(expand.grid(x = x, y = y), height_change = as.vector(height_change1))

    ## Restrict to Thwaites glacier
    height_change_thwaites <- height_change_df %>% filter(x >= xmin & x <= xmax & y >= ymin & y <= ymax)

    height_change_shelf_plot <- ggplot() +
        geom_raster(data = height_change_thwaites, aes(x = x, y = y, fill = height_change_1)) +
        scale_fill_viridis_c() +
        #   geom_sf(data = thwaites_bound, color = "black", fill = NA) +
        geom_point(data = thwaites_points, mapping = aes(x = X, y = Y), color = "black", size = 0.1) +
        #   geom_raster() +
        coord_fixed() +
        theme_bw()

    png("plots/height_change_shelf.png")
    print(height_change_shelf_plot)
    dev.off()

    qsave(height_change_thwaites, file = paste0(data_dir, "/height_change_shelf_thwaites.qs"))
} else {
    height_change_thwaites <- qread(paste0(data_dir, "/height_change_shelf_thwaites.qs"))
}


if (reread_thickness_data) {
    ## Do the same for thickness data
    annual_thickness <- mclapply(0:n_years, function(i) {
        thickness_t <- thickness_bedmachine[, , (i * 4 + 1):((i + 1) * 4)]
        thickness_avg <- apply(thickness_t, c(1, 2), mean, na.rm = TRUE)
    }, mc.cores = 10L)


    ## Restrict ice thickness to Thwaites glacier
    thickness_thwaites <- data.frame(expand.grid(x = x, y = y))
    for (t in 1:(n_years + 1)) {
        varname <- paste0("thickness_", years[t])
        thickness_thwaites[, varname] <- as.vector(annual_thickness[[t]])
    }
    thickness_thwaites <- thickness_thwaites %>% filter(x >= xmin & x <= xmax & y >= ymin & y <= ymax)

    qsave(thickness_thwaites, file = paste0(data_dir, "/thickness/thickness_shelf_thwaites_bedmachine.qs"))
}

browser()

## SMB
if (reread_smb_data) {
    annual_smb <- mclapply(0:n_years, function(i) {
        smb_t <- smb[, , (i * 4 + 1):((i + 1) * 4)]
        smb_avg <- apply(smb_t, c(1, 2), mean, na.rm = TRUE)
    }, mc.cores = 10L)

    smb_df <- data.frame(expand.grid(x = x, y = y))

    for (t in 1:(n_years + 1)) {
        varname <- paste0("smb_", years[t])
        smb_df[, varname] <- as.vector(annual_smb[[t]])
    }

    ## Restrict to Thwaites glacier
    smb_thwaites <- smb_df %>% filter(x >= xmin & x <= xmax & y >= ymin & y <= ymax)

    ## Plot SMB
    smb_shelf_plot <- ggplot() +
        geom_raster(data = smb_thwaites, aes(x = x, y = y, fill = smb_1)) +
        scale_fill_viridis_c() +
        # geom_sf(data = thwaites_bound, color = "black", fill = NA) +
        geom_point(data = thwaites_points, mapping = aes(x = X, y = Y), color = "black", size = 0.1) +
        coord_fixed() +
        theme_bw()

    png("plots/smb_shelf.png")
    print(smb_shelf_plot)
    dev.off()

    qsave(smb_thwaites, file = paste0(data_dir, "/SMB/smb_shelf_thwaites.qs"))
} else {
    smb_thwaites <- qread(paste0(data_dir, "/SMB/smb_shelf_thwaites.qs"))
}


## Melt rate
if (reread_basal_melt_data) {
    annual_melt <- mclapply(0:n_years, function(i) {
        melt_t <- melt[, , (i * 4 + 1):((i + 1) * 4)] # melt rate for each quarter
        melt_avg <- apply(melt_t, c(1, 2), mean, na.rm = TRUE) # average melt rate for each year
    }, mc.cores = 10L)

    melt_df <- data.frame(expand.grid(x = x, y = y))

    for (t in 1:(n_years + 1)) {
        varname <- paste0("melt_", years[t])
        melt_df[, varname] <- as.vector(annual_melt[[t]])
    }

    ## Restrict to Thwaites glacier
    melt_thwaites <- melt_df %>% filter(x >= xmin & x <= xmax & y >= ymin & y <= ymax)

    qsave(melt_thwaites, file = paste0(data_dir, "/SMB/melt_shelf_thwaites.qs"))

    ## Plot SMB
    melt_plot <- ggplot() +
        geom_raster(data = melt_thwaites, aes(x = x, y = y, fill = melt_1)) +
        scale_fill_viridis_c() +
        # geom_sf(data = thwaites_bound, color = "black", fill = NA) +
        geom_point(data = thwaites_points, mapping = aes(x = X, y = Y), color = "black", size = 0.1) +
        coord_fixed() +
        theme_bw()

    png("./plots/temp/melt_shelf.png")
    print(melt_plot)
    dev.off()
} else {
    melt_thwaites <- qread(paste0(data_dir, "/SMB/melt_shelf_thwaites.qs"))
}
