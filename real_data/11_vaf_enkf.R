## Approximate the area of ice using quadrature
setwd("~/SSA_model/CNN/real_data/")

library(qs)
library(ggplot2)
library(abind)
library(dplyr)

source("./source/surface_elev.R")

## Calculate VAF based on bed and EnKF-adjusted elevation data

data_date <- "20241111" # "20241103"
output_date <- "20241111" # "20240518"

sets <- 1:50 #51:100 # 6:20
correct_model_discrepancy <- T
avg_over_time <- T

setsf <- paste0("sets", sets[1], "-", sets[length(sets)])

data_dir <- paste0("./data/training_data/", setsf, "/")
output_dir <- paste0("./output/cnn/", setsf, "/")

## Directories
if (correct_model_discrepancy) {
    if (avg_over_time) {
        pred_output_dir <- paste0(output_dir, "pred/discr_avg/")
        plot_dir <- paste0("./plots/cnn/", setsf, "/pred/discr_avg/")
    } else {
        pred_output_dir <- paste0(output_dir, "pred/discr/")
        plot_dir <- paste0("./plots/cnn/", setsf, "/pred/discr/")
    }
    enkf_output_dir <- paste0(pred_output_dir, "enkf/")
} else {
    pred_output_dir <- paste0(output_dir, "pred/")
    plot_dir <- paste0("./plots/cnn/", setsf, "/pred/")
}

## Load flowline data here so we have the grid spacing
ssa_steady <- qread(file = paste0("./data/training_data/steady_state/steady_state_", data_date, ".qs"))
domain <- ssa_steady$domain  # in m
# J <- length(domain)
dx <- diff(domain)[1]

## Load real data
surf_elev_data <- qread(file = "./data/surface_elev/surf_elev_mat.qs")
velocity_data <- qread(file = "./data/velocity/vel_smoothed.qs")
n_years <- ncol(surf_elev_data)

## Model discrepancy
vel_discr_mat <- qread(file = paste0("./data/discrepancy/", setsf, "/vel_discr_smooth_", data_date, ".qs"))
se_discr_mat <- qread(file = paste0("./data/discrepancy/", setsf, "/se_discr_smooth_", data_date, ".qs"))

# Extrapolate to 11 years if needed
se_discr_mat <- cbind(se_discr_mat, se_discr_mat[, ncol(se_discr_mat)]) # repeat last column to make 11 years
vel_discr_mat <- cbind(vel_discr_mat, vel_discr_mat[, ncol(vel_discr_mat)]) # repeat last column to make 11 years

## Posterior samples
post_bed_samples <- qread(file = paste0(pred_output_dir, "bed_samples_real_", data_date, ".qs"))
prior_bed_samples <- qread(file = paste0(pred_output_dir, "prior_bed_samples_", data_date, ".qs"))

## EnKF ensembles 
samples <- 1:50
post_se_ens <- post_vel_ens <- list()
for (p in 1:length(samples)) {
cat("Posterior sample ", samples[p], "\n")
    # thick_file <- paste0(enkf_output_dir, "/enkf_thickness_p", p, "_", output_date, ".qs")
    post_se_file <- paste0(enkf_output_dir, "/enkf_surf_elev_posterior_p", samples[p], "_", output_date, ".qs")
    post_vel_file   <- paste0(enkf_output_dir, "/enkf_velocities_posterior_p", samples[p], "_", output_date, ".qs")

    # Thickness file
    if (file.exists(post_se_file)) {
        post_se_ens[[p]] <- qread(post_se_file)
    } else {
        post_se_ens[[p]] <- "MISSING"
    }

    # Velocity file
    if (file.exists(post_vel_file)) {
        post_vel_ens[[p]] <- qread(post_vel_file)
    } else {
        post_vel_ens[[p]] <- "MISSING"
    }
    
}

## Do the same for prior samples
# samples <- 1:10
prior_se_ens <- prior_vel_ens <- list()
for (p in 1:length(samples)) {
    prior_se_file <- paste0(enkf_output_dir, "enkf_surf_elev_prior_p", samples[p], "_", data_date, ".qs")
    prior_vel_file   <- paste0(enkf_output_dir, "enkf_velocities_prior_p", samples[p], "_", data_date, ".qs")

    # Thickness file
    if (file.exists(prior_se_file)) {
        prior_se_ens[[p]] <- qread(prior_se_file)
    } else {
        prior_se_ens[[p]] <- "MISSING"
    }

    # Velocity file
    if (file.exists(prior_vel_file)) {
        prior_vel_ens[[p]] <- qread(prior_vel_file)
    } else {
        prior_vel_ens[[p]] <- "MISSING"
    }
}

quadrature <- function(dx, y1, y2) {
    # Compute the difference between curves
    diff_y <- y2 - y1
    
    if (any(diff_y < 0)) {
        warning("Warning: y2 is less than y1 at some points. Area will be computed as absolute value.")
        diff_y <- abs(diff_y)
    }
    # Compute area using the trapezoidal rule
    # trapezoidal rule: sum((y[i] + y[i+1])/2 * dx)
    # dx <- diff(x)[1]  # uniform grid spacing
    area <- sum((diff_y[-1] + diff_y[-length(diff_y)]) / 2) * dx
    return(area)
}

calc_vaf <- function(dx, bed, se, rho_i = 917, rho_w = 1028) {
    # Compute the difference between curves
    h <- se - bed
    
    total_vol <- sum((h[-1] + h[-length(h)]) / 2) * dx
    h_under_water <- pmin(bed, 0)
    vol_under_water <- sum((h_under_water[-1] + h_under_water[-length(h_under_water)]) / 2) * dx

    vaf <- total_vol + (rho_w / rho_i) * vol_under_water

    h_per_cell <- (h[-1] + h[-length(h)]) / 2 # take height in the middle of each grid cell
    bed_per_cell <- (bed[-1] + bed[-length(bed)]) / 2

    vaf2 <- sum(pmax(h_per_cell + pmin(bed_per_cell, 0) * (rho_w / rho_i)), 0) * dx

    # Compute area using the trapezoidal rule
    # trapezoidal rule: sum((y[i] + y[i+1])/2 * dx)
    # dx <- diff(x)[1]  # uniform grid spacing
    # area <- sum((diff_y[-1] + diff_y[-length(diff_y)]) / 2) * dx
    return(list(vaf = vaf, vaf2 = vaf2))
}


## Filter out any runs where EnKF failed (i.e., no ensemble members)
nonmissing <- lapply(post_se_ens, function(x) !identical(x, "MISSING")) # Or sapply(my_list, length)
non_empty_indices <- unlist(nonmissing)

post_se_ens <- post_se_ens[non_empty_indices]
post_vel_ens <- post_vel_ens[non_empty_indices]
post_bed_samples <- post_bed_samples[, non_empty_indices]

prior_vaf_ls <- post_vaf_ls <- list()
for (p in 1:length(post_se_ens)) {
    # cat("Posterior sample ", p, "\n")

    ## Calculate surface elevation
    post_vaf_p <- list()
    for (yr in 1:n_years) {
        # cat("  Year ", yr, "\n")
        # se <- apply(thickness_ens[[p]][[yr]], 2, get_surface_elev,
        #                         b = post_bed_samples[, p])
        # se <- se + matrix(rep(se_discr_mat[, yr], ncol(se)), 
        #                        nrow = nrow(se), ncol = ncol(se))
        # post_se_ens[[yr]] <- se

        ## Calculate the VAF 
        post_vaf_p[[yr]] <- apply(post_se_ens[[p]][[yr+1]], 2, function(se) 
                            calc_vaf(dx = dx, bed = post_bed_samples[, p], se = se)$vaf)

    }
    post_vaf_ls[[p]] <- post_vaf_p
}

nonmissing <- lapply(prior_se_ens, function(x) !identical(x, "MISSING")) # Or sapply(my_list, length)
non_empty_indices <- unlist(nonmissing)

prior_se_ens <- prior_se_ens[non_empty_indices]
prior_vel_ens <- prior_vel_ens[non_empty_indices]
prior_bed_samples <- prior_bed_samples[, non_empty_indices]

for (p in 1:length(prior_se_ens)) {

    ## Calculate surface elevation
    prior_vaf_p <- list()
    for (yr in 1:n_years) {
        # cat("  Year ", yr, "\n")
        # se <- apply(thickness_ens[[p]][[yr]], 2, get_surface_elev,
        #                         b = post_bed_samples[, p])
        # se <- se + matrix(rep(se_discr_mat[, yr], ncol(se)), 
        #                        nrow = nrow(se), ncol = ncol(se))
        # post_se_ens[[yr]] <- se

        ## Calculate the VAF 
        prior_vaf_p[[yr]] <- apply(prior_se_ens[[p]][[yr+1]], 2, function(se) 
                            calc_vaf(dx = dx, bed = prior_bed_samples[, p], se = se)$vaf)
    }
    prior_vaf_ls[[p]] <- prior_vaf_p
}

## Concatenate ensemble members across all posterior samples
post_se_concat <- post_vel_concat <- post_vaf_concat <- list()
for (yr in 1:n_years) {

    se_year <- lapply(post_se_ens, function(x) x[[yr+1]]) # skip initial year
    post_se_concat[[yr]] <- abind(se_year, along = 2)
    
    velocity_year <- lapply(post_vel_ens, function(x) x[[yr+1]])
    post_vel_concat[[yr]] <- abind(velocity_year, along = 2)

    vaf_year <- lapply(post_vaf_ls, function(x) x[[yr]])
    post_vaf_concat[[yr]] <- unlist(vaf_year)
}

## Do the same for the prior samples
prior_se_concat <- prior_vel_concat <- prior_vaf_concat <- list()
for (yr in 1:n_years) {
    se_year <- lapply(prior_se_ens, function(x) x[[yr+1]]) # skip initial year
    prior_se_concat[[yr]] <- abind(se_year, along = 2)
    
    velocity_year <- lapply(prior_vel_ens, function(x) x[[yr+1]])
    prior_vel_concat[[yr]] <- abind(velocity_year, along = 2)

    vaf_year <- lapply(prior_vaf_ls, function(x) x[[yr]])
    prior_vaf_concat[[yr]] <- unlist(vaf_year)
}


## Plot the ensemble surface elevation samples for each year
## And compare against actual observations
plot_range <- which(!is.na(surf_elev_data[, 1]))
year_num <- 2010:2020
png(paste0(plot_dir, "/enkf/post_se_ensemble_", data_date, ".png"), width = 1500, height = 2500, res = 150)
par(mfrow = c(6, 2))    
for (yr in 1:n_years) {
    # cols <- colorRampPalette(c("lightblue", "blue"))(10)
    matplot(domain[plot_range]/1e3, prior_se_concat[[yr]][plot_range, ], type = "l", lty = 1, col = "lightblue",
            xlab = "Distance along flowline (km)", ylab = "Surface elevation (m)",
            main = paste0("Surface elevation ensemble for ", year_num[yr]))
    matlines(domain[plot_range]/1e3, post_se_concat[[yr]][plot_range, ], type = "l", lty = 1, col = "lightpink")
    # Add observed surface elevation data
    lines(domain[plot_range]/1e3, surf_elev_data[plot_range, yr], pch = 16, col = "black")
    legend("topright", legend = c("Prior", "Posterior", "Observed"), col = c("lightblue", "lightpink", "black"), lty = c(1, 1))
    
}
dev.off()

## Plot the ensemble velocity samples for each year
png(paste0(plot_dir, "/enkf/post_velocity_ensemble_", data_date, ".png"), width = 1500, height = 2500, res = 150)
par(mfrow = c(6, 2))    
for (yr in 1:n_years) {
    # cols <- colorRampPalette(c("lightblue", "blue"))(10)
    matplot(domain[plot_range]/1e3, prior_vel_concat[[yr]][plot_range, ], type = "l", lty = 1, col = "lightblue",
            xlab = "Distance along flowline (km)", ylab = "Velocity (m/yr)",
            main = paste0("Velocity ensemble for year ", year_num[yr]))
    matlines(domain[plot_range]/1e3, post_vel_concat[[yr]][plot_range, ], type = "l", lty = 1, col = "lightpink")
    # Add observed velocity data
    lines(domain[plot_range]/1e3, velocity_data[plot_range, yr], pch = 16, col = "black")
    legend("topright", legend = c("Prior", "Posterior", "Observed"), col = c("lightblue", "lightpink", "black"), lty = c(1, 1))
}
dev.off()

## Calculate the surface elevation based on thickness and bed samples


## Just subset the bed samples to the number of samples in the predictive dist
# post_se_samples <- qread(file = paste0(pred_output_dir, "post_pred_obs_", data_date, ".qs"))
# n_post_samples <- dim(post_se_samples)[1]

## For the first sample in post_se_samples,
## compute the change in elevation per year to check
# cols <- colorRampPalette(c("turquoise", "blue"))(10)
# sp <- 6
# se_change <- post_se_samples[sp, , 2:11, 1] - post_se_samples[sp, , 1:10, 1]
# png(paste0(plot_dir, "se_change_sample1_", data_date, ".png"), width = 1000, height = 500, res = 150)
# matplot(domain/1e3, se_change, type = "l", lty = 1, col = cols,
#         xlab = "Distance along flowline (km)", ylab = "Change in surface elevation (m)",
#         main = paste0("Surface elevation change per year for posterior sample ", sp))
# dev.off()



# ## Calculate volume above floatation for each posterior sample
# vol_ls <- post_vaf_ls <- post_vaf_ls2 <- list()
# for (s in 1:n_post_samples) {
#     post_bed <- post_bed_samples[, s]
#     post_se <- post_se_samples[s, , , 1]
#     n_grounded <- which(!is.na(post_se[, 1]))
    
#     ## Subset the bed to only where surface elevation data is available
#     post_bed <- post_bed[n_grounded]
#     # post_bed_mat <- matrix(rep(post_bed, dim(post_se)[2]), nrow = length(post_bed), ncol = dim(post_se)[2])
#     post_se <- post_se[n_grounded, ]
#     # test <- quadrature(dx = dx, y1 = post_bed, y2 = post_se)
#     vol <- sapply(1:dim(post_se)[2], function(i) quadrature(dx = dx, y1 = post_bed, y2 = post_se[, i]))
#     vaf <- sapply(1:dim(post_se)[2], function(i) calc_vaf(dx = dx, bed = post_bed, se = post_se[, i])$vaf )
#     vaf2 <- sapply(1:dim(post_se)[2], function(i) calc_vaf(dx = dx, bed = post_bed, se = post_se[, i])$vaf2)
    
#     # if (vol[length(vol)] > vol[1]) {
#     #     warning("Warning: Ice volume is increasing over time in sample ", s, ". Check the model predictions.")
#     #     browser()

#     #     png("./plots/temp/ice_volume_increasing_warning.png", width = 800, height = 600)
#     #     plot(domain[n_grounded], post_se[, 1], type = "l", lty = 1, col = "lightblue",
#     #          xlab = "Distance along flowline (km)", ylab = "Surface elevation (m)",
#     #          main = paste("Posterior surface elevation samples for sample", s))
#     #     lines(domain[n_grounded], post_se[, dim(post_se)[2]], col = "blue", lwd = 2)
#     #     # matplot(domain[n_grounded], post_se, type = "l", lty = 1, col = rgb(0, 0, 1, 0.3),
#     #     #         xlab = "Distance along flowline (km)", ylab = "Surface elevation (m)",
#     #     #         main = paste("Posterior surface elevation samples for sample", s))
#     #     lines(domain[n_grounded], post_bed, col = "red", lwd = 2)
#     #     legend("topright", legend = c("Surface elevation in year 1", "Surface elevation in year 11", "Bed elevation"), col = c("lightblue", "blue", "red"), lwd = c(1, 2))
#     #     dev.off()   
#     # }
#     vol_ls[[s]] <- vol
#     post_vaf_ls[[s]] <- vaf  
#     post_vaf_ls2[[s]] <- vaf2  
    
#     # cat("Sample:", s, "Mean ice volume (m^2):", mean(test2) * 1e9, "\n")
# }

## Create box plot of ice volume
post_vaf_mat <- do.call(cbind, post_vaf_concat)

## Convert to km^2
post_vaf_mat <- post_vaf_mat * 1e-6

## Convert to data frame for ggplot
N <- length(post_vaf_concat[[1]])
vaf_df <- data.frame(year = rep(2010:2020, each = N),
                     volume = as.vector(post_vaf_mat)) 

vaf_plot <- ggplot(data = vaf_df, aes(x = as.factor(year), y = volume)) +
      geom_boxplot() +
      theme_bw() +
    #   ylim (c(95e3, 120e3)) +
      labs(x = "Year", y = expression("Ice area (km"^2*")")) +
      ggtitle("Plot of \"volume\" above flotation (VAF)") + 
        theme(text = element_text(size = 24))

png(paste0(plot_dir, "/enkf/ice_vaf_enkf_boxplot_", data_date, ".png"), width = 1500, height = 900, res = 150)
print(vaf_plot)
dev.off()

## Also calculate change in ice volume year to year
vaf_diff <- t(apply(post_vaf_mat, 1, diff))

vaf_diff_df <- data.frame(year = rep(2011:2020, each = N),
                          volume_change = as.vector(vaf_diff))

vaf_change_plot <- ggplot(data = vaf_diff_df, aes(x = as.factor(year), y = volume_change)) +
      geom_boxplot() +
      theme_bw() +
      labs(x = "Year", y = expression("Change in ice area (km"^2*")")) +
      ggtitle("Change in area (volume) above flotation") +
      theme(text = element_text(size = 24))

png(paste0(plot_dir, "/enkf/ice_vaf_change_enkf_boxplot_", data_date, ".png"), width = 800, height = 600, res = 100)
print(vaf_change_plot)
dev.off()

# ## Calculate actual VAF using observed surface elevation data
surf_elev_data <- qread(file = "./data/surface_elev/surf_elev_mat.qs")

# ## but we don't have the true bed...
# ## so maybe compare the change in height to observed change in height?
years <- ncol(surf_elev_data)
obs_se_change <- surf_elev_data[, 2:years] - surf_elev_data[, 1:(years - 1)]
avg_obs_se_change <- rowMeans(obs_se_change, na.rm = T)

# total_obs_se_change <- surf_elev_data[, years] - surf_elev_data[, 1]
# ## Calculate mean predicted change in surface elevation from posterior samples
# # mean_post_se <- apply(post_se_samples, c(2, 3), mean)
# # pred_se_change <- mean_post_se[, 2:years] - mean_post_se[,

# ## Plot avg observed change in surface elevation
# png(paste0(plot_dir, "avg_obs_se_change_", data_date, ".png"), width = 800, height = 600, res = 100)
# plot(domain/1e3, avg_obs_se_change, type = "b", pch = 16,
#      xlab = "Year", ylab = "Average observed change in surface elevation (m)",
#      main = "Average observed change in surface elevation over time")
# dev.off()
