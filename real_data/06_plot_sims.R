setwd("~/SSA_model/CNN/real_data/")

rm(list = ls())

# library(parallel)
# library(Matrix)
# library(R.utils)
# library(fields)
library(dplyr)
# library(matrixStats) # for the rowMaxs() function
# library(mvtnorm)
library(abind)
library(ggplot2)
library(gridExtra)
# library(FRK)
library(qs)

## Flags
use_missing_pattern <- TRUE

## Plot simulations and also plot real data
set <- 51
setf <- formatC(set, width = 2, flag = "0")

## Directory for training data
train_data_dir <- "./data/training_data"
data_date <- "20241111"

### True thickness, friction, bed, grounding line
surface_obs_arr <- qread(file = paste0(train_data_dir, "/surface_obs_arr_", setf, "_", data_date, ".qs"))
gl_arr <- qread(file = paste0(train_data_dir, "/gl_arr_", setf, "_", data_date, ".qs"))

friction_arr <- qread(file = paste0(train_data_dir, "/friction_arr_", setf, "_", data_date, ".qs"))
bed_arr <- qread(file = paste0(train_data_dir, "/bed_arr_", setf, "_", data_date, ".qs"))

friction_basis <- qread(file = paste0("data/training_data/friction_basis_", setf, "_", data_date, ".qs"))
bed_basis <- qread(file = paste0("data/training_data/bed_basis_", setf, "_", data_date, ".qs"))

if (use_missing_pattern) {
    print("Using missing patterns to mask the observations")

    surf_elev_missing_pattern <- qread("./data/surface_elev/missing_pattern.qs")
    vel_missing_pattern <- qread("./data/velocity/missing_pattern.qs")

    ## Replace 0/1 by NA/1
    surf_elev_missing_pattern <- ifelse(surf_elev_missing_pattern == 0, NA, 1)
    vel_missing_pattern <- ifelse(vel_missing_pattern == 0, NA, 1)

    surface_obs_list <- lapply(1:dim(surface_obs_arr)[1], function(i) { surface_obs_arr[i,,,]})

    ## Multiply the missing patterns by the surface_obs array
    surface_obs_list_missing <- lapply(surface_obs_list, function(arr) {
    se <- arr[,,1] * surf_elev_missing_pattern
    vel <- arr[,,2] * vel_missing_pattern
    abind(se, vel, along = 3)
    })

    surface_obs_arr <- abind(surface_obs_list_missing, along = 0)
    rm(surface_obs_list)
}


## Load data on flowline
ssa_steady <- qread(file = paste0(train_data_dir, "/steady_state/steady_state_", data_date, ".qs"))
domain <- ssa_steady$domain

## Number of grid points and time points
J <- length(domain)
years <- dim(surface_obs_arr)[3]

## BedMap observations
bed_obs_df <- qread(file = paste0("./data/bedmap/bed_obs_df_all.qs"))

## Plot surface observations for a simulation
png("plots/cnn/input/surface_obs_sim.png", width = 6, height = 8, units = "in", res = 300)
par(mfrow = c(2, 1))
matplot(surface_obs_arr[1, , , 1], type = "l", main = "Surface elevation (m)", ylab = "Surface elev. (m)", xlab = "Domain (grid points)")
matplot(surface_obs_arr[1, , , 2], type = "l", main = "Velocity (m/yr)", ylab = "Velocity (m/yr)", xlab = "Domain (grid points)")
dev.off()

## Check if the generated data is close to real data here
surf_elev_mat <- qread("./data/surface_elev/surf_elev_mat.qs")
# vel_mat <- qread(file = "./data/velocity/all_velocity_arr.qs")
vel_mat <- qread(file = "./data/velocity/vel_smoothed.qs")

sim <- 1
png("./plots/cnn/input/sim_vs_real.png", width = 500, height = 600)
par(mfrow = c(2, 1))

matplot(surface_obs_arr[sim, , , 1],
    type = "l", col = "grey",
    xlab = "Grid point", ylab = "Surface elevation (m)"
)
matlines(surf_elev_mat, col = "salmon", lty = 2)
lines(ssa_steady$current_top_surface, col = "blue")

matplot(vel_mat,
    ylim = c(0, 6000), type = "l", col = "salmon",
    xlab = "Grid point", ylab = "Velocity (m/a)"
)
# matlines(vel_mat2, col = "black", lty = 2)
matlines(surface_obs_arr[sim, , , 2], col = "grey", lty = 2)
lines(ssa_steady$current_velocity, col = "blue")
dev.off()

## Hovmoller plots of some simulations
plots <- list()

nsamples <- 4
sims <- sample(1:dim(surface_obs_arr)[1], size = nsamples)

space <- domain / 1000
time <- 1:years # 1:dim(thickness_velocity_arr)[3]
grid_test <- expand.grid(space, time)
head(grid_test)
names(grid_test) <- c("space", "time")

inds <- matrix(1:(nsamples * 4), nsamples, 4, byrow = T)

gl <- ceiling(gl_arr[1, 1] / (domain[length(domain)] / 1000) * length(domain))

# s <- 1
for (s in 1:nsamples) {
    sim <- sims[[s]]
    ind <- inds[s, ]
    surface_elev <- surface_obs_arr[sim, , , 1]
    velocity <- surface_obs_arr[sim, , , 2]
    grid_test$surface_elev <- as.vector(surface_elev)
    grid_test$velocity <- as.vector(velocity)

    surface_elev_plot <- ggplot(grid_test) +
        geom_tile(aes(space, time, fill = surface_elev)) +
        scale_y_reverse() +
        scale_fill_distiller(palette = "Blues", direction = 1) +
        theme_bw() +
        theme(text = element_text(size = 24)) +
        xlab("Domain (km)") +
        ylab("Year") +
        # labs(fill="Thickness (m)")
        labs(fill = "Surface elev. (m)")

    velocity_plot <- ggplot(grid_test) +
        geom_tile(aes(space, time, fill = velocity)) +
        scale_y_reverse() +
        theme_bw() +
        xlab("Domain (km)") +
        ylab("Year") +
        theme(text = element_text(size = 24)) +
        scale_fill_distiller(palette = "Reds", direction = 1) +
        labs(fill = bquote("Velocity (m" ~ a^-1 ~ ")"))


    # if (log_transform) {
    #   # fitted_fric_sim <- exp(fitted_friction[sim, 1:gl])
    #   friction_sim <- exp(friction_arr[sim, 1:gl])
    # } else {
    # fitted_fric_sim <- fitted_friction[sim, 1:gl]
    friction_sim <- friction_arr[sim, 1:gl]
    # }

    df <- data.frame(
        domain = ssa_steady$domain[1:gl] / 1000, friction = friction_sim # ,
        # fitted_fric = fitted_fric_sim
    )
    friction_plot <- ggplot(df, aes(x = domain, y = friction)) +
        geom_line() +
        # geom_line(aes(x = domain, y = fitted_fric), col = "red") +
        theme_bw() +
        xlim(c(0, domain[gl] / 1e3)) +
        xlab("Domain (km)") +
        ylab(bquote("Friction (M Pa m"^"-1/3" ~ "a"^"1/3" ~ ")")) +
        theme(text = element_text(size = 24))

    bed_sim <- bed_basis$fitted_values[sim, ] + bed_basis$mean

    #   fitted_bed_sim <- fitted_bed[sim, ] + bed_mean
    bed_df <- data.frame(domain = ssa_steady$domain / 1000, bed = bed_sim) # , fitted_bed = fitted_bed_sim)

    bed_plot <- ggplot(bed_df, aes(x = domain, y = bed)) +
        geom_line() +
        geom_point(data = bed_obs_df, aes(x = loc/1000, y = bed_elev), col = "red", size = 2) +
        # geom_line(aes(x = domain, y = fitted_bed), col = "red") +
        theme_bw() +
        ylim(c(-1500, -500)) +
        xlim(c(0, domain[gl] / 1e3)) +
        xlab("Domain (km)") +
        ylab("Bed (m)") +
        theme(text = element_text(size = 24))

    # plots[[ind[1]]] <- thickness_plot
    plots[[ind[1]]] <- surface_elev_plot
    plots[[ind[2]]] <- velocity_plot
    plots[[ind[3]]] <- friction_plot
    plots[[ind[4]]] <- bed_plot
}

png(
    file = paste0("./plots/cnn/input/simulations_", setf, "_", data_date, ".png"),
    width = 2800, height = 400 * nsamples, res = 100
)
grid.arrange(grobs = plots, nrow = nsamples, ncol = 4)
dev.off()


# png(file = paste0("./plots/cnn/test.png"))
# matplot(surface_obs_arr[2,,,1], type = "l", col = "grey")
# matlines(surf_elev_mat, col = "red")
# dev.off()

## Hovmoller plots of real data

vel_mat <- vel_mat * vel_missing_pattern
velocities_df <- data.frame(
    gridpt = rep(1:J, times = years),
    dist = rep(domain/1000, times = length(years)),
    year = rep(1:years, each = J),
    vel = as.vector(vel_mat)
)

se_df <- data.frame(
    gridpt = rep(1:J, times = years),
    dist = rep(domain/1000, times = length(years)),
    year = rep(1:years, each = J),
    se = as.vector(surf_elev_mat)
)

# velocities_df <- velocities_df %>% mutate(nonmissing = ifelse(is.na(vel_mat), 0, 1))
# se_df <- se_df %>% mutate(nonmissing = ifelse(is.na(surf_elev_mat), 0, 1))

# vel_mp_plot <- ggplot(velocities_df) +
#     geom_tile(aes(x = gridpt, y = year, fill = factor(nonmissing))) +
#     # scale_fill_distiller(palette = "BuPu", direction = 1) +
#     theme_bw() +
#     # scale_y_discrete(limits=rev) +
#     labs(x = "Grid Point", y = "Year", fill = "Non-missing") +
#     ggtitle("Thwaites Glacier Velocity Over Time") +
#     theme(plot.title = element_text(hjust = 0.5))

# png("./plots/cnn/input/real_vel_hovmoller.png", width = 800, height = 600)
# print(vel_mp_plot)
# dev.off()

## Plot actual velocity data with missing pattern overlay
vel_masked_plot <- ggplot(velocities_df) +
    geom_tile(aes(x = dist, y = year, fill = vel)) +
    scale_fill_distiller(palette = "Reds", direction = 1) +
    theme_bw() +
    scale_y_reverse() +
    # scale_y_discrete(limits=rev) +
    labs(x = "Dist (km)", y = "Year", fill = "Velocity (m/yr)") +
    # ggtitle("Thwaites Glacier Velocity Over Time") +
    theme(text = element_text(size = 24))

se_masked_plot <- ggplot(se_df) +
    geom_tile(aes(x = dist, y = year, fill = se)) +
    scale_fill_distiller(palette = "Blues", direction = 1) +
    theme_bw() +
    scale_y_reverse() +
    # scale_y_discrete(limits=rev) +
    labs(x = "Dist (km)", y = "Year", fill = "Surf. elev. (m)") +
    # ggtitle("Thwaites Glacier Surface Elevation Over Time") +
    theme(text = element_text(size = 24))

png("./plots/cnn/input/real_data_hovmoller.png", width = 2000, height = 600, res = 150)
grid.arrange(grobs = list(se_masked_plot, vel_masked_plot), nrow = 1)
dev.off()
