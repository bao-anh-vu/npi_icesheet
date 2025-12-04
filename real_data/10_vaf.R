## Approximate the area of ice using quadrature
setwd("~/SSA_model/CNN/real_data/")

library(qs)
library(ggplot2)
## Toy example to test quadrature methods

# Example: area between two curves using trapezoidal integration

# Generate a regular grid
x <- seq(0, 2 * pi, length.out = 200)

# Define two curves
y1 <- sin(x)
y2 <- cos(x)

# Compute the difference between curves
diff_y <- abs(y2 - y1)

# Compute area using the trapezoidal rule
# trapezoidal rule: sum((y[i] + y[i+1])/2 * dx)
dx <- diff(x)[1]  # uniform grid spacing
area <- sum((diff_y[-1] + diff_y[-length(diff_y)]) / 2) * dx

cat("Approximate area between curves:", area, "\n")

# Optional: visualize
png("./plots/temp/area_between_curves.png", width = 800, height = 600)
plot(x, y1, type = "l", col = "blue", lwd = 2, ylim = c(-1, 1))
lines(x, y2, col = "red", lwd = 2)
polygon(c(x, rev(x)), c(y1, rev(y2)), col = rgb(0.2, 0.7, 0.9, 0.3), border = NA)
legend("topright", legend = c("y1 = sin(x)", "y2 = cos(x)"), col = c("blue", "red"), lwd = 2)

dev.off()

## Now do the same but for bed and posterior predicted elevation data

data_date <- "20241111" # "20241103"
sets <- 1:50 #51:100 #1:50 #51:100 # 6:20
# use_missing_pattern <- Tth
# use_basal_melt_data <- T
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
} else {
    pred_output_dir <- paste0(output_dir, "pred/")
    plot_dir <- paste0("./plots/cnn/", setsf, "/pred/")
}

## Load flowline data here so we have the grid spacing
ssa_steady <- qread(file = paste0("./data/training_data/steady_state/steady_state_", data_date, ".qs"))
domain <- ssa_steady$domain  # in m
# J <- length(domain)
dx <- diff(domain)[1]

## Posterior samples
post_bed_samples <- qread(file = paste0(pred_output_dir, "bed_samples_real_", data_date, ".qs"))

## Just subset the bed samples to the number of samples in the predictive dist
post_se_samples <- qread(file = paste0(pred_output_dir, "post_pred_obs_", data_date, ".qs"))
n_post_samples <- dim(post_se_samples)[1]

## For the first sample in post_se_samples,
## compute the change in elevation per year to check
cols <- colorRampPalette(c("turquoise", "blue"))(10)
sp <- 1
se_change <- post_se_samples[sp, , 2:11, 1] - post_se_samples[sp, , 1:10, 1]
png(paste0(plot_dir, "se_change_sample1_", data_date, ".png"), width = 1000, height = 500, res = 150)
matplot(domain/1e3, se_change, type = "l", lty = 1, col = cols,
        xlab = "Distance along flowline (km)", ylab = "Change in surface elevation (m)",
        main = paste0("Surface elevation change per year for posterior sample ", sp))
dev.off()


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
    
    # total_vol <- sum((h[-1] + h[-length(h)]) / 2) * dx
    # h_under_water <- pmin(bed, 0)
    # vol_under_water <- sum((h_under_water[-1] + h_under_water[-length(h_under_water)]) / 2) * dx

    # vaf <- total_vol + (rho_w / rho_i) * vol_under_water

    h_per_cell <- (h[-1] + h[-length(h)]) / 2 # take height in the middle of each grid cell
    bed_per_cell <- (bed[-1] + bed[-length(bed)]) / 2

    vaf <- sum(pmax(h_per_cell + pmin(bed_per_cell, 0) * (rho_w / rho_i)), 0) * dx
    return(vaf)

    # Compute area using the trapezoidal rule
    # trapezoidal rule: sum((y[i] + y[i+1])/2 * dx)
    # dx <- diff(x)[1]  # uniform grid spacing
    # area <- sum((diff_y[-1] + diff_y[-length(diff_y)]) / 2) * dx
    # return(list(vaf = vaf, vaf2 = vaf2))
}

## Calculate volume above floatation for each posterior sample
vol_ls <- vaf_ls <- vaf_ls2 <- list()
for (s in 1:n_post_samples) {
    post_bed <- post_bed_samples[, s]
    post_se <- post_se_samples[s, , , 1]
    n_grounded <- which(!is.na(post_se[, 1]))
    
    ## Subset the bed to only where surface elevation data is available
    post_bed <- post_bed[n_grounded]
    # post_bed_mat <- matrix(rep(post_bed, dim(post_se)[2]), nrow = length(post_bed), ncol = dim(post_se)[2])
    post_se <- post_se[n_grounded, ]
    # test <- quadrature(dx = dx, y1 = post_bed, y2 = post_se)
    vol <- sapply(1:dim(post_se)[2], function(i) quadrature(dx = dx, y1 = post_bed, y2 = post_se[, i]))
    vaf <- sapply(1:dim(post_se)[2], function(i) calc_vaf(dx = dx, bed = post_bed, se = post_se[, i]))
    # vaf2 <- sapply(1:dim(post_se)[2], function(i) calc_vaf(dx = dx, bed = post_bed, se = post_se[, i])$vaf2)
    
    # if (vol[length(vol)] > vol[1]) {
    #     warning("Warning: Ice volume is increasing over time in sample ", s, ". Check the model predictions.")
    #     browser()

    #     png("./plots/temp/ice_volume_increasing_warning.png", width = 800, height = 600)
    #     plot(domain[n_grounded], post_se[, 1], type = "l", lty = 1, col = "lightblue",
    #          xlab = "Distance along flowline (km)", ylab = "Surface elevation (m)",
    #          main = paste("Posterior surface elevation samples for sample", s))
    #     lines(domain[n_grounded], post_se[, dim(post_se)[2]], col = "blue", lwd = 2)
    #     # matplot(domain[n_grounded], post_se, type = "l", lty = 1, col = rgb(0, 0, 1, 0.3),
    #     #         xlab = "Distance along flowline (km)", ylab = "Surface elevation (m)",
    #     #         main = paste("Posterior surface elevation samples for sample", s))
    #     lines(domain[n_grounded], post_bed, col = "red", lwd = 2)
    #     legend("topright", legend = c("Surface elevation in year 1", "Surface elevation in year 11", "Bed elevation"), col = c("lightblue", "blue", "red"), lwd = c(1, 2))
    #     dev.off()   
    # }
    vol_ls[[s]] <- vol
    vaf_ls[[s]] <- vaf  
    # vaf_ls2[[s]] <- vaf2  
    
    # cat("Sample:", s, "Mean ice volume (m^2):", mean(test2) * 1e9, "\n")
}

## Calculate actual VAF using observed surface elevation data
surf_elev_data <- qread(file = "./data/surface_elev/surf_elev_mat.qs")
obs_vaf <- sapply(1:dim(post_se)[2], function(i) calc_vaf(dx = dx, bed = post_bed, se = na.omit(surf_elev_data[, i])))
obs_vaf <- obs_vaf * 1e-6 ## Convert to km^2

## Create box plot of ice volume
vaf_mat <- do.call(rbind, vaf_ls)
vaf_mat <- vaf_mat * 1e-6 ## Convert to km^2

## Convert to data frame for ggplot
vaf_df <- data.frame(year = rep(2010:2020, each = n_post_samples),
                     volume = as.vector(vaf_mat))  

obs_vaf_df <- data.frame(year = 2010:2020,
                         volume = obs_vaf)

vaf_plot <- ggplot(data = vaf_df, aes(x = as.factor(year), y = volume)) +
  geom_boxplot() +
  geom_point(data = obs_vaf_df,
             aes(x = as.factor(year), y = volume),
             color = "red", size = 3) +
  theme_bw() +
  labs(x = "Year", y = expression("Ice area above flotation (km"^2*")")) +
  theme(text = element_text(size = 19))


# vaf_plot <- ggplot(data = vaf_df, aes(x = as.factor(year), y = volume)) +
#       geom_boxplot() +
#       theme_bw() +
#     #   ylim (c(95e3, 120e3)) +
#       labs(x = "Year", y = expression("Ice area above flotation (km"^2*")")) +
#     #   ggtitle("Plot of \"volume\" above flotation (VAF)") + 
#         theme(text = element_text(size = 19))

png(paste0(plot_dir, "ice_vaf_boxplot_", data_date, ".png"), width = 1500, height = 600, res = 150)
print(vaf_plot)
dev.off()


## Also calculate change in ice volume year to year
vaf_diff <- t(apply(vaf_mat, 1, diff))
vaf_diff_df <- data.frame(year = rep(2011:2020, each = n_post_samples),
                          volume_change = as.vector(vaf_diff))

obs_vaf_diff <- obs_vaf[2:length(obs_vaf)] - obs_vaf[1:(length(obs_vaf) - 1)]
obs_vaf_diff_df <- data.frame(year = 2011:2020,
                              volume_change = obs_vaf_diff)

vaf_change_plot <- ggplot(data = vaf_diff_df, aes(x = as.factor(year), y = volume_change)) +
      geom_boxplot() +
        geom_point(data = obs_vaf_diff_df,
                     aes(x = as.factor(year), y = volume_change),
                     color = "red", size = 3) +
        theme_bw() +
      labs(x = "Year", y = expression("Change in ice area (km"^2*")")) +
    #   ggtitle("Change in area (volume) above flotation") +
    theme(text = element_text(size = 20))

png(paste0(plot_dir, "ice_vaf_change_boxplot_", data_date, ".png"), width = 1500, height = 600, res = 150)
print(vaf_change_plot)
dev.off()

## Median change in VAF per year
# median_vaf <- apply(vaf_diff, 2, median)

    
## so maybe compare the change in height to observed change in height?
years <- ncol(surf_elev_data)
obs_se_change <- surf_elev_data[, 2:years] - surf_elev_data[, 1:(years - 1)]
avg_obs_se_change <- colMeans(obs_se_change, na.rm = T) #rowMeans(obs_se_change, na.rm = T)

total_obs_se_change <- surf_elev_data[, years] - surf_elev_data[, 1]

## Calculate mean predicted change in surface elevation from posterior samples
# mean_post_se <- apply(post_se_samples, c(2, 3), mean)
# pred_se_change <- mean_post_se[, 2:years] - mean_post_se[,

## Plot avg observed change in surface elevation
# png(paste0(plot_dir, "avg_obs_se_change_", data_date, ".png"), width = 800, height = 600, res = 100)
# plot(2011:2020, avg_obs_se_change, type = "b", pch = 16,
#      xlab = "Year", ylab = "Average observed change in surface elevation (m)",
#      main = "Average observed change in surface elevation over time")
# dev.off()
