setwd("~/SSA_model/CNN/real_data/")

library(qs)

# data_dir <- "./data/"
data_date <- "20241111" # "20241103"
sets <- 1:50 #51:100 #51:100 #6:20
setsf <- paste0("sets", sets[1], "-", sets[length(sets)])

data_dir <- paste0("./data/training_data/", setsf, "/")
output_dir <- paste0("./output/cnn/", setsf, "/")

pred_output_dir <- paste0(output_dir, "pred/discr/")
plot_dir <- paste0("./plots/cnn/", setsf, "/pred/discr/")

avg_over_time <- T

# Load surface elevation data
surf_elev_mat <- qread(file = "./data/surface_elev/surf_elev_mat.qs")
vel_mat <- qread(file = "./data/velocity/vel_smoothed.qs")

## Model discrepancy
vel_discr <- qread(file = paste0("./data/discrepancy/", setsf, "/vel_discr_", data_date, ".qs"))
se_discr <- qread(file = paste0("./data/discrepancy/", setsf, "/se_discr_", data_date, ".qs"))

avg_vel_discr <- rowMeans(vel_discr, na.rm = T)
avg_se_discr <- rowMeans(se_discr, na.rm = T)

nyears <- 10
  if (avg_over_time) { # just use average discrepancy for all years
      vel_discr_mat <- matrix(rep(avg_vel_discr, nyears), nrow = length(avg_vel_discr), ncol = nyears)
      se_discr_mat <- matrix(rep(avg_se_discr, nyears), nrow = length(avg_se_discr), ncol = nyears)
  } else { 
        vel_discr_mat <- vel_discr
        se_discr_mat <- se_discr
    }

adj_se_mat <- surf_elev_mat[, 1:10] - se_discr_mat
adj_vel_mat <- vel_mat[, 1:10] - vel_discr_mat

## Load prior simulations
prior_pred_obs <- qread(file = paste0(pred_output_dir, "prior_pred_obs_unadj_", data_date, ".qs"))
# prior_pred_obs <- prior_pred_obs[1:100, , ,]

## Take average over prior simulations
prior_mean <- apply(prior_pred_obs, c(2,3,4), mean)

# ## Compare prior sims against adjusted data
# png(paste0(plot_dir, "prior_vs_adj_obs.png"), width = 1200, height = 1200, res = 100)
# par(mfrow = c(2,1))
# matplot(prior_mean[,,1], type = "l", lty = 1, col = "lightblue",
#         xlab = "Adjusted observed surface elevation (m)",
#         ylab = "Prior predicted surface elevation (m)",
#         main = "Prior predictions vs adjusted observed surface elevation")
# matlines(adj_se_mat, col = "lightpink")

# ## Same for the velocity
# matplot(prior_mean[,,2], type = "l", lty = 1, col = "lightblue",
#         xlab = "Adjusted observed velocity (m/yr)",
#         ylab = "Prior predicted velocity (m/yr)",
#         main = "Prior predictions vs adjusted observed velocity")
# matlines(adj_vel_mat, col = "lightpink")
# dev.off()

## For each year, compare prior mean against adjusted data
years <- 1:10
png(paste0(plot_dir, "prior_vs_adj_se_", data_date, ".png"), width = 800, height = 1200, res = 100)
par(mfrow = c(5, 2))
    
for (yr in years) {
    # Surface elevation
    plot(prior_mean[, yr, 1], type = "l", lwd = 1.5, 
         xlab = "Domain (m)",
         ylab = "Surface elevation (m)",
         main = paste0("Year ", yr),
         col = "grey", pch = 16)
    lines(adj_se_mat[, yr], col = "red", lwd = 1.5, lty = 3)
    legend("topright", legend = c("Prior mean", "Adjusted obs"), col = c("grey", "red"), lwd = 1.5, lty = c(1,3))
    # abline(0, 1, col = "red", lty = 2)
    
}
dev.off()

## Same for velocity
png(paste0(plot_dir, "prior_vs_adj_vel_", data_date, ".png"), width = 800, height = 1000, res = 100)
par(mfrow = c(5, 2))
for (yr in years) {
    # Velocity
    plot(prior_mean[, yr, 2], type = "l", lwd = 1.5, 
         xlab = "Domain (m)",
         ylab = "Velocity (m/yr)",
         main = paste0("Year ", yr),
         col = "grey", pch = 16)
    lines(adj_vel_mat[, yr], col = "red", lwd = 1.5, lty = 3)
    legend("topleft", legend = c("Prior mean", "Adjusted obs"), col = c("grey", "red"), lwd = 1.5, lty = c(1,3))
    # abline(0, 1, col = "red", lty = 2)
}
dev.off()