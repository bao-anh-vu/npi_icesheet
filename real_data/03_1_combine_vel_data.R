## Combine velocity data
setwd("~/SSA_model/CNN/real_data/")

# library(ncdf4)
# library(CFtime)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(qs)
library(abind)
# library(lattice)
# library(sf) # for shapefiles
# library(sfheaders)
# library(sp)
# library(parallel)

years <- 2010:2020
data_dir <- "./data"
velocities <- list()
recombine_data <- T
censor_data <- T ## censor velocity data from GL onwards

# flowline <- readRDS(paste0(data_dir, "/flowline_regrid.rds"))
flowline <- qread(paste0(data_dir, "/flowline_regrid.qs"))
J <- nrow(flowline) # number of grid points

## Grounding line data
gl_pos <- qread(paste0(data_dir, "/grounding_line/gl_pos.qs")) ## grounding line position
gl_ind <- gl_pos$ind

# if (recombine_data) {
for (i in 1:length(years)) {
    year <- years[i]
    # velocity_i <- readRDS(paste0(data_dir, "/velocity/flowline_vel_mapped_", year, ".rds"))
    velocity_i <- qread(paste0(data_dir, "/velocity/flowline_vel_mapped_", year, ".qs"))
    velocity_i$year <- year
    velocity_i$gridpt <- 1:nrow(velocity_i)
    velocities[[i]] <- velocity_i
}

velocities_only <- lapply(velocities, function(x) x %>% select(vel_avg))
velocities_truncated <- lapply(velocities_only, function(x) x[1:J, ])
velocity_arr <- abind(velocities_truncated, along = 2)

qsave(velocity_arr, "./data/velocity/all_velocity_arr.qs")

# png("./plots/velocity/velocity_arr.png", width = 750, height = 500)
# matplot(velocity_arr, col = "grey", type = "l")
# abline(v = gl_ind, col = "black", lty = 2)
# dev.off()

###################################################
##  Smooth velocity data using median polishing  ##
###################################################
flowline_dist <- sqrt((flowline$x[2:J] - flowline$x[1:(J - 1)])^2 + (flowline$y[2:J] - flowline$y[1:(J - 1)])^2)
flowline_dist <- c(0, cumsum(na.omit(flowline_dist)))

vel_smoothed <- matrix(NA, nrow = J, ncol = length(years))

# gl_pos <- readRDS(paste0(data_dir, "/grounding_line/gl_pos.rds")) ## grounding line position
# gl_pos <- qread(paste0(data_dir, "/grounding_line/gl_pos.qs")) ## grounding line position
# gl_ind <- gl_pos$ind
# delta <- 120 # grid size
# flowline$ind <- 1:nrow(flowline)
# gl_near_pts <- flowline %>% filter(
#                 x >= (gl_pos[1] - delta) & x <= (gl_pos[1] + delta),
#                 y >= (gl_pos[2] - delta) & y <= (gl_pos[2] + delta))
# gl_ind <- gl_near_pts$ind

## If before GL, use 20 intervals; if after GL, use 50 intervals
# n_intervals <- 20
interval_size_before_gl <- 100 # gl_ind %/% 5
interval_size_after_gl <- 200 # (J - gl_ind) %/% 3

# i <- 1
smooth_years <- years[1:3] # only smooth first 3 years (2010-2012) as the other years look fine
for (i in 1:length(smooth_years)) {
    year <- years[i]
    vel_yr <- velocity_arr[, i]

    ## Divide into 10 roughly equal intervals
    before_gl <- vel_yr[1:gl_ind]
    after_gl <- vel_yr[(gl_ind + 1):J]
    intervals_before_gl <- split(before_gl, ceiling(seq_along(before_gl) / interval_size_before_gl))
    intervals_after_gl <- split(after_gl, ceiling(seq_along(after_gl) / interval_size_after_gl))
    intervals <- c(intervals_before_gl, intervals_after_gl)
    # intervals <- split(vel_yr, cut(1:length(vel_yr), breaks = seq(0, length(vel_yr), by = interval_size)))

    # interval_medians_before_gl <- sapply(intervals_before_gl, function(x) median(x, na.rm = T))
    # interval_medians_after_gl <- sapply(intervals_after_gl, function(x) median(x, na.rm = T))
    # interval_medians <- c(interval_medians_before_gl, interval_medians_after_gl)
    interval_medians <- sapply(intervals, function(x) median(x, na.rm = T))

    interval_outliers <- sapply(intervals, function(x) {
        which(x > (median(x, na.rm = T) + 1.5 * IQR(x, na.rm = T)) |
            x < (median(x, na.rm = T) - 1.5 * IQR(x, na.rm = T)))
    })
    ## Replace outliers with the median
    for (j in 1:length(intervals)) {
        intervals[[j]][interval_outliers[[j]]] <- interval_medians[j]
    }
    vel_smoothed[, i] <- unlist(intervals)

    png(paste0("./plots/temp/vel_smoothed_", year, ".png"), width = 750, height = 500)
    plot_range <- 1:J
    plot(velocity_arr[plot_range, i], type = "l", col = "grey")
    lines(vel_smoothed[plot_range, i], col = "red")
    abline(v = gl_ind, col = "black", lty = 2)
    dev.off()
}

vel_smoothed2 <- cbind(vel_smoothed[, 1:3], velocity_arr[, 4:ncol(velocity_arr)])

## Plot velocity for individual years
plots <- list()
n_years <- length(years)
for (i in 1:n_years) {
    year <- years[i]
    vel_df <- data.frame(x = flowline_dist / 1000, 
                        vel = velocity_arr[, i], 
                        vel_smooth = vel_smoothed2[, i])
    p <- ggplot(vel_df, aes(x = x, y = vel)) +
        geom_line(color = "black") +
        # geom_line(aes(y = vel_smooth), color = "salmon") +
        # xlim(c(0, 150)) +
        theme_bw() +
        labs(x = "Distance along flowline (km)", y = "Velocity (m/yr)") +
        ggtitle(paste("Velocity in", year)) 
    plots[[i]] <- p

    png(paste0("./plots/velocity/vel_yr_", year, ".png"), width = 1000, height = 600, res = 200)
    print(p)
    dev.off()

}

vel_smoothed2[, 1:3] <- NA # discard 2010-2012 data as well due to noisy observations
    
qsave(vel_smoothed2, "./data/velocity/vel_smoothed.qs")

######################################
##          Data masking            ## 
######################################

if (censor_data) { # mask data after GL
    vel_smoothed2[(gl_ind + 1):J, ] <- NA
    # velocity_arr[(gl_ind + 1):J, ] <- NA
    ## Manually "mask" some unreliable velocity values in 2010 and 2011
    # vel_smoothed2[1:500, 1] <- NA # discard first 500 grid points in 2010 as the values seem unreliable
    
}


png("./plots/velocity/vel_smoothed.png", width = 750, height = 500)
matplot(velocity_arr, col = "grey", type = "l")
matlines(vel_smoothed2, col = "salmon")
abline(v = gl_ind, col = "black", lty = 2)
dev.off()

## Concatenate all velocity data
# velocities_df <- bind_rows(vel_smoothed2)
# tail(velocities_df)

# Turn vel_smoothed2 into a data frame for plotting
velocities_df <- data.frame(
    gridpt = rep(1:J, times = length(years)),
    year = rep(years, each = J),
    vel_avg = as.vector(vel_smoothed2)
)

## Plot velocity data over space-time
vel_st_plot <- ggplot(velocities_df) +
    geom_tile(aes(x = gridpt, y = year, fill = vel_avg)) +
    scale_fill_distiller(palette = "Reds", direction = 1) +
    theme_bw() +
    labs(x = "Grid Point", y = "Year", fill = "Velocity (m/yr)") +
    # scale_y_discrete(limits=rev) +
    ggtitle("Thwaites Glacier Velocity Over Time") +
    theme(plot.title = element_text(hjust = 0.5))

png("./plots/velocity/vel_st_plot.png", width = 800, height = 600)
print(vel_st_plot)
dev.off()


###############################
##  Extract missing pattern  ##
###############################

vel_missing_pattern <- ifelse(is.na(vel_smoothed2), 0, 1)
qsave(vel_missing_pattern, "./data/velocity/missing_pattern.qs")


velocities_df <- data.frame(
    gridpt = rep(1:J, times = length(years)),
    year = rep(years, each = J),
    vel_avg = as.vector(vel_smoothed2)
)
# vel_mat <- as.matrix(velocities_df %>% select(year, vel_avg))

velocities_df <- velocities_df %>% mutate(nonmissing = ifelse(is.na(vel_avg), 0, 1))

vel_mp_plot <- ggplot(velocities_df) +
    geom_tile(aes(x = gridpt, y = year, fill = factor(nonmissing))) +
    # scale_fill_distiller(palette = "BuPu", direction = 1) +
    theme_bw() +
    # scale_y_discrete(limits=rev) +
    labs(x = "Grid Point", y = "Year", fill = "Non-missing") +
    ggtitle("Thwaites Glacier Velocity Over Time") +
    theme(plot.title = element_text(hjust = 0.5))

png("./plots/missing_pattern/vel_missing_ptn.png", width = 800, height = 600)
print(vel_mp_plot)
dev.off()

## Plot actual velocity data with missing pattern overlay
vel_masked_plot <- ggplot(velocities_df) +
    geom_tile(aes(x = gridpt, y = year, fill = vel_avg)) +
    scale_fill_distiller(palette = "Reds", direction = 1) +
    theme_bw() +
    # scale_y_discrete(limits=rev) +
    labs(x = "Grid Point", y = "Year", fill = "Non-missing") +
    ggtitle("Thwaites Glacier Velocity Over Time") +
    theme(plot.title = element_text(hjust = 0.5))

png("./plots/missing_pattern/vel_masked.png", width = 800, height = 600)
print(vel_masked_plot)
dev.off()


# ## Turn missing pattern into a matrix
# vel_missing_pattern <- velocities_df %>%
#     filter(gridpt <= J) %>%
#     select(nonmissing) %>%
#     as.matrix() %>%
#     matrix(nrow = J, ncol = length(years))

# } else {
#     velocity_arr <- qread("./data/velocity/all_velocity_arr.qs")
#     vel_missing_pattern <- qread("./data/velocity/missing_pattern.qs")
# }





## Smooth velocity data
# vel_smoothed <- matrix(NA, nrow = nrow(velocity_arr), ncol = ncol(velocity_arr))

# for (i in 1:length(years)) {
#     year <- years[i]
#     vel_df <- data.frame(x = flowline_dist, vel = velocity_arr[, i])
#     test <- loess(vel ~ x, data = vel_df, span = 0.05)#$fitted
#     imputed <- predict(test, newdata = vel_df)

#     png(paste0("./plots/temp/vel_imputed_", year, ".png"), width = 750, height = 500)
#     plot(vel_df$x, vel_df$vel, type = "l", col = "grey", lwd = 2,
#             ylim = c(0, max(imputed, na.rm = T)),
#             main = paste("Velocity in ", year))
#     lines(vel_df$x, imputed, col = "red")
#     dev.off()

# }



# ## Missing pattern for surface elevation data
# # year <- 2000
# surf_elev_mat <- matrix(NA, nrow = J, ncol = length(years))
# for (i in 1:length(years)) {
#     year <- years[i]
#     surf_elev_data <- readRDS(file = paste0(data_dir, "/surface_elev/surf_elev_", year, ".rds"))
#     surf_elev_mat[, i] <- unlist(surf_elev_data)[1:J]
# }

# qsave(surf_elev_mat, "./data/surface_elev/surf_elev_mat.qs")

# png("./plots/surface_elev/surf_elev_st_plot.png", width = 800, height = 600)
# image(surf_elev_mat)
# dev.off()

# surf_ev_missing_pattern <- ifelse(is.na(surf_elev_mat), 0, 1)
# qsave(surf_ev_missing_pattern, "./data/surface_elev/missing_pattern.qs")
