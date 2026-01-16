## Map SMB to flowline

setwd("~/SSA_model/CNN/real_data/")

# library(sf)
# library(ncdf4)
# library(CFtime)
library(dplyr)
library(ggplot2)
# library(RColorBrewer)
library(parallel)
library(qs)


data_dir <- "./data/"

smb_thwaites <- qread(file = paste0(data_dir, "/SMB/smb_shelf_thwaites.qs"))

flowline <- qread(paste0(data_dir, "/flowline_regrid.qs"))
flowline_pos <- lapply(1:nrow(flowline), function(i) as.numeric(flowline[i, ]))

avg_nearest_four <- function(df, pos, grid_size) {
            
    delta <- grid_size
    
    near_pts <- df %>% filter(
        x >= (pos[1] - delta) & x <= (pos[1] + delta),
        y >= (pos[2] - delta) & y <= (pos[2] + delta)
    )
    near_v <- near_pts %>%
        mutate(dist = sqrt((x - pos[1])^2 + (y - pos[2])^2)) %>%
        filter(dist > 0) %>%
        arrange(dist) %>%
        slice_min(dist, n = 4) %>%
        select(val) %>%
        summarise(nearest = val[1], avg = mean(val)) 
    # as.numeric()

    return(near_v)
}

flowline_shelf_smb <- matrix(NA, nrow = ncol(smb_thwaites) - 2, ncol = nrow(flowline))
for (col in 3:ncol(smb_thwaites)) {

    t <- col - 2
    cat("t = ", t, "\n")

    ## Select grid points in v with coordinates in x_thwaites and y_thwaites
    smb_df <- smb_thwaites[, c(1, 2, col)] 
    colnames(smb_df) <- c("x", "y", "val")

    ## Map SMB data to flowline
    smb_mapped <- mclapply(flowline_pos, avg_nearest_four, 
                        df = smb_df, grid_size = 1920, mc.cores = 10L)

    # smb_nearest <- sapply(smb_mapped, function(x) x$v_nearest)
    smb_avg <- sapply(smb_mapped, function(x) x$avg)

    flowline_shelf_smb[t, ] <- smb_avg
    # flowline$smb_nearest <- smb_nearest
    # flowline$smb_avg <- smb_avg

}

qsave(flowline_shelf_smb, file = paste0(data_dir, "/SMB/flowline_shelf_smb.qs"))

##############################
##      Basal melt data     ## 
##############################

melt_data <- qread(file = paste0(data_dir, "/SMB/melt_shelf_thwaites.qs"))

flowline_shelf_melt <- matrix(NA, nrow = ncol(melt_data) - 2, ncol = nrow(flowline))
for (col in 3:ncol(melt_data)) {
# for (col in 3:4) {

    t <- col - 2
    cat("t = ", t, "\n")

    ## Try filtering grid points in v with coordinates in x_thwaites and y_thwaites
    melt_df <- melt_data[, c(1, 2, col)] 
    colnames(melt_df) <- c("x", "y", "val")

    ## Map SMB data to flowline
    melt_mapped <- mclapply(flowline_pos, avg_nearest_four, 
                            df = melt_df, grid_size = 1920, mc.cores = 10L)

    # smb_nearest <- sapply(smb_mapped, function(x) x$v_nearest)
    melt_avg <- sapply(melt_mapped, function(x) x$avg)

    flowline_shelf_melt[t, ] <- melt_avg

}

qsave(flowline_shelf_melt, file = paste0(data_dir, "/SMB/flowline_shelf_melt.qs"))

##################################
##      Height change data      ##
##################################

height_change_data <- qread(file = paste0(data_dir, "/height_change_shelf_thwaites.qs"))

flowline_shelf_height_change <- matrix(NA, nrow = ncol(height_change_data) - 2, ncol = nrow(flowline))

for (col in 3:ncol(height_change_data)) {
# for (col in 3:4) {

    t <- col - 2
    cat("t = ", t, "\n")

    ## Try filtering grid points in v with coordinates in x_thwaites and y_thwaites
    height_change_df <- height_change_data[, c(1, 2, col)] 
    colnames(height_change_df) <- c("x", "y", "val")

    ## Map SMB data to flowline
    height_change_mapped <- mclapply(flowline_pos, avg_nearest_four, 
                            df = height_change_df, grid_size = 1920, mc.cores = 10L)

    # smb_nearest <- sapply(smb_mapped, function(x) x$v_nearest)
    height_change_avg <- sapply(height_change_mapped, function(x) x$avg)

    flowline_shelf_height_change[t, ] <- height_change_avg

}

qsave(flowline_shelf_height_change, file = paste0(data_dir, "/SMB/flowline_shelf_height_change.qs"))
