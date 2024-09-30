## Compare results from EnKF conditional on CNN posterior and state-augmented EnKF

library(ggplot2)
library(mvtnorm)

setwd("/home/babv971/SSA_model/CNN/simbed/")
rm(list = ls())

source("./source/enkf/surface_elev.R")

## Presets
data_date <- "20220329" # "20230518"
output_date <- "20240320" # "20240518"

use_basis_funs <- T
plot_ice_thickness <- T
plot_velocity <- T
plot_bed <- T
plot_friction <- T

## SSA model info
ssa_steady <- readRDS(file = paste("./training_data/initial_conds/ssa_steady_20220329.rds", sep = ""))
# reference <- readRDS(file = paste("./training_data/initial_conds/reference_20220329.rds", sep = ""))
domain <- ssa_steady$domain
J <- length(domain)

## Read bed and friction from NN output
sets <- 1:50 # 10
setsf <- paste0("sets", sets[1], "-", sets[length(sets)])
s <- 3 # test sample index

years <- 20
Ne <- 1000 # Ensemble size

data_dir <- paste0("./training_data/", setsf)
cnn.output_dir <- paste0("./output/posterior/", setsf)
cnn_enkf.output_dir <- paste0("./output/posterior/", setsf, "/sample", s)

if (use_basis_funs) {
    enkfsa.output_dir <- paste0("./output/stateaug/", setsf, "/basis", "/sample", s)
} else {
    enkfsa.output_dir <- paste0("./output/stateaug/", setsf, "/no_basis", "/sample", s)
}

## Read test data
test_data <- readRDS(file = paste0(data_dir, "/test_data_", output_date, ".rds"))

true_surface_elevs <- test_data$true_surface_elevs_test
true_thicknesses <- test_data$true_thickness_test
true_velocities <- test_data$true_velocity_test

true_bed <- test_data$true_bed
true_fric <- exp(test_data$true_fric)
true_gl <- test_data$grounding_line * test_data$sd_gl + test_data$mean_gl

## Read CNN predictions
pred_fric <- readRDS(file = paste0(cnn.output_dir, "/pred_fric_", output_date, ".rds"))
pred_bed <- readRDS(file = paste0(cnn.output_dir, "/pred_bed_", output_date, ".rds"))
pred_gl <- readRDS(file = paste0(cnn.output_dir, "/pred_gl_", output_date, ".rds"))

print("Reading posterior samples from CNN...")
fric_samples_ls <- readRDS(file = paste0(cnn.output_dir, "/fric_post_samples_", output_date, ".rds"))
bed_samples_ls <- readRDS(file = paste0(cnn.output_dir, "/bed_post_samples_", output_date, ".rds"))
gl_samples_ls <- readRDS(file = paste0(cnn.output_dir, "/gl_post_samples_", output_date, ".rds"))

## Scaling units for friction coefficients
# secpera <- 31556926
# fric_scale <- 1e6 * secpera^(1 / 3)


## Read EnKF-CNN results
cnn.thickness <- readRDS(file = paste0(cnn_enkf.output_dir, "/enkf_thickness_sample", s, "_", output_date, ".rds", sep = ""))
cnn.velocity <- readRDS(file = paste0(cnn_enkf.output_dir, "/enkf_velocities_sample", s, "_", output_date, ".rds", sep = ""))

cnn.bed <- pred_bed[, s]
cnn.fric <- pred_fric[, s]
cnn.gl <- pred_gl[s, ]

cnn.bed_lq <- apply(bed_samples_ls[[s]], 1, quantile, probs = 0.05)
cnn.bed_uq <- apply(bed_samples_ls[[s]], 1, quantile, probs = 0.95)
cnn.fric_lq <- apply(fric_samples_ls[[s]], 1, quantile, probs = 0.05)
cnn.fric_uq <- apply(fric_samples_ls[[s]], 1, quantile, probs = 0.95)


## Read EnKF-StateAug results
enkfsa.thickness <- readRDS(file = paste0(enkfsa.output_dir, "/enkf_thickness_sample", s, "_Ne", Ne, "_", output_date,  ".rds", sep = ""))
enkfsa.bed <- readRDS(file = paste0(enkfsa.output_dir, "/enkf_bed_sample", s, "_Ne", Ne, "_", output_date,  ".rds", sep = ""))
enkfsa.friction <- readRDS(file = paste0(enkfsa.output_dir, "/enkf_friction_sample", s, "_Ne", Ne, "_", output_date,  ".rds", sep = ""))
enkfsa.velocity <- readRDS(file = paste0(enkfsa.output_dir, "/enkf_velocities_sample", s, "_Ne", Ne, "_", output_date,  ".rds", sep = ""))

enkfsa.friction_fin <- exp(enkfsa.friction[[years + 1]])
enkfsa.bed_fin <- enkfsa.bed[[years + 1]]

if (use_basis_funs) {
    plot_dir <- paste0("./plots/combined/Ne", Ne, "/basis/sample", s)
} else {
    plot_dir <- paste0("./plots/combined_Ne", Ne, "/no_basis/sample", s)
}

if (!dir.exists(plot_dir)) {
    dir.create(paste0(plot_dir))
} else { # delete all previously saved plots
    unlink(paste0(plot_dir, "/*"))
}

plot_times <- seq(1, years + 1, 2)

## Ice thickness plot
if (plot_ice_thickness) {
    for (t in plot_times) {
        # t <- plot_times[ind]
        cnn.ens_t <- cnn.thickness[[t]]
        cnn.covmat <- 1 / (ncol(cnn.ens_t) - 1) * tcrossprod(cnn.ens_t - rowMeans(cnn.ens_t)) # diag(enkf_covmats[[t]])
        cnn.lower <- rowMeans(cnn.ens_t[1:J, ]) + qnorm(0.025) * sqrt(diag(cnn.covmat)[1:J])
        cnn.upper <- rowMeans(cnn.ens_t[1:J, ]) + qnorm(0.975) * sqrt(diag(cnn.covmat)[1:J])

        enkfsa.ens_t <- enkfsa.thickness[[t]]
        enkfsa.covmat <- 1 / (ncol(enkfsa.ens_t) - 1) * tcrossprod(enkfsa.ens_t - rowMeans(enkfsa.ens_t)) # diag(enkf_covmats[[t]])
        enkfsa.lower <- rowMeans(enkfsa.ens_t[1:J, ]) + qnorm(0.025) * sqrt(diag(enkfsa.covmat)[1:J])
        enkfsa.upper <- rowMeans(enkfsa.ens_t[1:J, ]) + qnorm(0.975) * sqrt(diag(enkfsa.covmat)[1:J])

        # Grounding line
        # ref.GL <- reference$grounding_line[t]
        # GL <- gl_migrate(H = enkf_means[1:J, t], b = rowMeans(ens[(J+1):(2*J), ]))
        # GL <- GL / J * domain[length(domain)] / 1000

        thickness.df <- data.frame(
            domain = domain / 1000,
            cnn_mean = rowMeans(cnn.ens_t),
            enkfsa_mean = rowMeans(enkfsa.ens_t),
            cnn_lower = cnn.lower, cnn_upper = cnn.upper,
            enkfsa_lower = enkfsa.lower, enkfsa_upper = enkfsa.upper
        )

        true_thickness.df <- data.frame(domain = domain / 1000, thickness = true_thicknesses[s, , t])

        # Plot title
        title <- paste("t = ", t - 1, "a")
        thickness_plot <- ggplot(thickness.df) +
            geom_ribbon(
                # data = thickness.df,
                aes(domain, ymin = enkfsa_lower, ymax = enkfsa_upper),
                fill = "blue", alpha = 0.2
            ) +
            geom_ribbon(
                # data = thickness.df,
                aes(domain, ymin = cnn_lower, ymax = cnn_upper),
                fill = "red", alpha = 0.2
            ) +
            theme_bw() +
            geom_line(data = true_thickness.df, aes(domain, thickness)) +
            geom_line(data = thickness.df, aes(domain, enkfsa_mean), colour = "blue") +
            geom_line(data = thickness.df, aes(domain, cnn_mean), colour = "red") +
            xlab("Domain (km)") +
            ylab("Ice thickness (m)") +
            # xlim(0, 500) +
            ggtitle(title) +
            theme(plot.title = element_text(
                hjust = 0.95, vjust = 0.5, face = "bold",
                margin = margin(t = 20, b = -30)
            ))

        # thickness_plots[[ind]] <- thickness_plot
        png(paste0(plot_dir, "/thickness_enkf_", t - 1, ".png"), width = 2000, height = 1000, res = 300)
        print(thickness_plot)
        # grid.arrange(grobs = thickness_plots, ncol = 2)
        dev.off()
    }
}


## Velocity plot
if (plot_velocity) {
    for (t in plot_times) {
        cnn.ens_t <- cnn.velocity[[t]]
        cnn.covmat <- 1 / (ncol(cnn.ens_t) - 1) * tcrossprod(cnn.ens_t - rowMeans(cnn.ens_t)) # diag(enkf_covmats[[t]])
        cnn.lower <- rowMeans(cnn.ens_t[1:J, ]) + qnorm(0.025) * sqrt(diag(cnn.covmat)[1:J])
        cnn.upper <- rowMeans(cnn.ens_t[1:J, ]) + qnorm(0.975) * sqrt(diag(cnn.covmat)[1:J])

        enkfsa.ens_t <- enkfsa.velocity[[t]]
        enkfsa.covmat <- 1 / (ncol(enkfsa.ens_t) - 1) * tcrossprod(enkfsa.ens_t - rowMeans(enkfsa.ens_t)) # diag(enkf_covmats[[t]])
        enkfsa.lower <- rowMeans(enkfsa.ens_t[1:J, ]) + qnorm(0.025) * sqrt(diag(enkfsa.covmat)[1:J])
        enkfsa.upper <- rowMeans(enkfsa.ens_t[1:J, ]) + qnorm(0.975) * sqrt(diag(enkfsa.covmat)[1:J])

        # Grounding line
        # ref.GL <- reference$grounding_line[t]
        # GL <- gl_migrate(H = enkf_means[1:J, t], b = rowMeans(ens[(J+1):(2*J), ]))
        # GL <- GL / J * domain[length(domain)] / 1000
        # bg.GL <- gl_migrate(H = bg_state_means[1:J, t], b = bg_state_means[(J+1):(2*J), t])
        # bg.GL <- bg.GL / J * domain[length(domain)] / 1000
        #

        velocities.df <- data.frame(
            domain = domain / 1000,
            cnn_mean = rowMeans(cnn.ens_t),
            enkfsa_mean = rowMeans(enkfsa.ens_t),
            cnn_lower = cnn.lower, cnn_upper = cnn.upper,
            enkfsa_lower = enkfsa.lower, enkfsa_upper = enkfsa.upper
        )
        true_velocities.df <- data.frame(domain = domain / 1000, velocity = true_velocities[s, , t])

        # Plot title
        title <- paste("t = ", t - 1, "a")
        velocity_plot <- ggplot(velocities.df) +
            geom_ribbon(
                aes(domain, ymin = enkfsa_lower, ymax = enkfsa_upper),
                fill = "blue", alpha = 0.2
            ) +
            geom_ribbon(
                aes(domain, ymin = cnn_lower, ymax = cnn_upper),
                fill = "red", alpha = 0.2
            ) +
            theme_bw() +
            geom_line(data = true_velocities.df, aes(domain, velocity)) +
            geom_line(data = velocities.df, aes(domain, enkfsa_mean), colour = "blue") +
            geom_line(data = velocities.df, aes(domain, cnn_mean), colour = "red") +
            xlab("Domain (km)") +
            ylab("Velocity (m)") +
            ggtitle(title) +
            theme(plot.title = element_text(
                hjust = 0.95, vjust = 0.5, face = "bold",
                margin = margin(t = 20, b = -30)
            ))

        png(paste0(plot_dir, "/velocity_enkf", t - 1, ".png"), width = 2000, height = 1000, res = 300)
        print(velocity_plot)
        # grid.arrange(grobs = vel_plots, ncol = 2)
        dev.off()
        # print(thickness_plot)
    }
}

## Plot the GL predicted by the CNN, GL using inferred state + flotation condition, and true GL position here
### Need to first calculate GL based on inferred state and bed
enkfsa.mean_thickness <- lapply(enkfsa.thickness, rowMeans)
enkfsa.bed <- rowMeans(enkfsa.bed_fin)
enkfsa.gl <- sapply(enkfsa.mean_thickness, gl_migrate, b = enkfsa.bed) # grid pts
enkfsa.gl <- enkfsa.gl / J * domain[length(domain)] / 1000 # convert to km

cnn.mean_thickness <- lapply(cnn.thickness, rowMeans)
enkf_cnn.gl <- sapply(cnn.mean_thickness, gl_migrate, b = cnn.bed) # grid pts
enkf_cnn.gl <- enkf_cnn.gl / J * domain[length(domain)] / 1000 # convert to km

png(paste0(plot_dir, "/gl_comparison.png"), width = 2000, height = 1000, res = 300)
plot(true_gl[s, ], type = "l", col = "black", lwd = 2, xlab = "Time (years)", ylab = "Grounding line (km)")
lines(cnn.gl, col = "red", lwd = 2) # CNN prediction
lines(enkfsa.gl, col = "blue", lwd = 2) # EnKF-SA prediction
lines(enkf_cnn.gl, col = "goldenrod", lwd = 2) # EnKF-CNN prediction
legend("bottomright", legend = c("True GL", "CNN", "EnKF-SA", "EnKF-CNN"), col = c("black", "red", "blue", "goldenrod"), lwd = 2)
dev.off()

## Bed
if (plot_bed) {
    enkfsa.ens_t <- enkfsa.bed_fin
    enkfsa.covmat <- 1 / (ncol(enkfsa.ens_t) - 1) * tcrossprod(enkfsa.ens_t - rowMeans(enkfsa.ens_t)) # diag(enkf_covmats[[t]])
    enkfsa.lower <- rowMeans(enkfsa.ens_t[1:J, ]) + qnorm(0.025) * sqrt(diag(enkfsa.covmat)[1:J])
    enkfsa.upper <- rowMeans(enkfsa.ens_t[1:J, ]) + qnorm(0.975) * sqrt(diag(enkfsa.covmat)[1:J])

    # Grounding line
    # ref.GL <- reference$grounding_line[t]
    # GL <- gl_migrate(H = enkf_means[1:J, t], b = rowMeans(ens[(J+1):(2*J), ]))
    # GL <- GL / J * domain[length(domain)] / 1000
    # bg.GL <- gl_migrate(H = bg_state_means[1:J, t], b = bg_state_means[(J+1):(2*J), t])
    # bg.GL <- bg.GL / J * domain[length(domain)] / 1000

    bed.df <- data.frame(
        domain = domain / 1000,
        cnn_mean = cnn.bed,
        enkfsa_mean = rowMeans(enkfsa.ens_t),
        cnn_lower = cnn.bed_lq, cnn_upper = cnn.bed_uq,
        enkfsa_lower = enkfsa.lower, enkfsa_upper = enkfsa.upper
    )
    true_bed.df <- data.frame(domain = domain / 1000, bed = true_bed[s, ])

    # Plot title
    # title <- paste("t = ", t - 1, "a")
    bed_plot <- ggplot(bed.df) +
        geom_ribbon(
            aes(domain, ymin = enkfsa_lower, ymax = enkfsa_upper),
            fill = "blue", alpha = 0.2
        ) +
        geom_ribbon(
            aes(domain, ymin = cnn_lower, ymax = cnn_upper),
            fill = "red", alpha = 0.2
        ) +
        theme_bw() +
        geom_line(data = true_bed.df, aes(domain, bed)) +
        geom_line(data = bed.df, aes(domain, enkfsa_mean), colour = "blue") +
        geom_line(data = bed.df, aes(domain, cnn_mean), colour = "red") +
        geom_vline(xintercept = true_gl[s, years], lty = 2, colour = "black") +
        xlab("Domain (km)") +
        ylab("Velocity (m)") +
        ggtitle(title) +
        theme(plot.title = element_text(
            hjust = 0.95, vjust = 0.5, face = "bold",
            margin = margin(t = 20, b = -30)
        ))

    png(paste0(plot_dir, "/bed_enkf.png"), width = 2000, height = 1000, res = 300)
    print(bed_plot)
    # grid.arrange(grobs = vel_plots, ncol = 2)
    dev.off()
}

## Friction
if (plot_friction) {
    enkfsa.ens_t <- enkfsa.friction_fin
    enkfsa.covmat <- 1 / (ncol(enkfsa.ens_t) - 1) * tcrossprod(enkfsa.ens_t - rowMeans(enkfsa.ens_t)) # diag(enkf_covmats[[t]])
    enkfsa.lower <- rowMeans(enkfsa.ens_t[1:J, ]) + qnorm(0.025) * sqrt(diag(enkfsa.covmat)[1:J])
    enkfsa.upper <- rowMeans(enkfsa.ens_t[1:J, ]) + qnorm(0.975) * sqrt(diag(enkfsa.covmat)[1:J])

    enkfsa.log_fric_samples <- rmvnorm(1000, rowMeans(enkfsa.ens_t), enkfsa.covmat)
    enkfsa.fric_samples <- exp(enkfsa.log_fric_samples)
    enkfsa.fric_mean <- colMeans(enkfsa.fric_samples)
    lower <- apply(enkfsa.fric_samples, 2, quantile, probs = 0.05)
    upper <- apply(enkfsa.fric_samples, 2, quantile, probs = 0.95)

    # Grounding line
    # ref.GL <- reference$grounding_line[t]
    # GL <- gl_migrate(H = enkf_means[1:J, t], b = rowMeans(ens[(J+1):(2*J), ]))
    # GL <- GL / J * domain[length(domain)] / 1000
    # bg.GL <- gl_migrate(H = bg_state_means[1:J, t], b = bg_state_means[(J+1):(2*J), t])
    # bg.GL <- bg.GL / J * domain[length(domain)] / 1000
    #

    friction.df <- data.frame(
        domain = domain / 1000,
        cnn_mean = cnn.fric,
        enkfsa_mean = rowMeans(enkfsa.ens_t),
        cnn_lower = cnn.fric_lq, cnn_upper = cnn.fric_uq,
        enkfsa_lower = enkfsa.lower, enkfsa_upper = enkfsa.upper
    )
    true_fric.df <- data.frame(domain = domain / 1000, fric = true_fric[s, ])

    # Plot title
    # title <- paste("t = ", t - 1, "a")
    friction_plot <- ggplot(friction.df) +
        geom_ribbon(
            aes(domain, ymin = enkfsa_lower, ymax = enkfsa_upper),
            fill = "blue", alpha = 0.2
        ) +
        geom_ribbon(
            aes(domain, ymin = cnn_lower, ymax = cnn_upper),
            fill = "red", alpha = 0.2
        ) +
        theme_bw() +
        geom_line(data = true_fric.df, aes(domain, fric)) +
        geom_line(data = friction.df, aes(domain, enkfsa_mean), colour = "blue") +
        geom_line(data = friction.df, aes(domain, cnn_mean), colour = "red") +
        geom_vline(xintercept = true_gl[s, years], lty = 2, colour = "black") +
        xlab("Domain (km)") +
        ylab("Velocity (m)") +
        ggtitle(title) +
        theme(plot.title = element_text(
            hjust = 0.95, vjust = 0.5, face = "bold",
            margin = margin(t = 20, b = -30)
        ))

    png(paste0(plot_dir, "/friction_enkf.png"), width = 2000, height = 1000, res = 300)
    print(friction_plot)
    # grid.arrange(grobs = vel_plots, ncol = 2)
    dev.off()
}

##############################
##    RMSE (time series)    ##
##############################

rmse <- function(estimated, true) {
    stopifnot(length(estimated) == length(true))
    sum(sqrt((estimated - true)^2))
}

## Thickness RMSE
cnn.thickness_rmse <- c()
enkfsa.thickness_rmse <- c()
for (t in 1:years) {
    cnn.thickness_rmse[t] <- rmse(rowMeans(cnn.thickness[[t + 1]]), true_thicknesses[s, , t + 1])
    enkfsa.thickness_rmse[t] <- rmse(rowMeans(enkfsa.thickness[[t + 1]]), true_thicknesses[s, , t + 1])
}
plot(1:years, enkfsa.thickness_rmse,
    type = "o", col = "blue",
    ylim = c(20000, 47000),
    xlab = "Time (years)", ylab = "RMSE",
    main = paste0("RMSE of ice thickness over time for sample ", s)
)
lines(1:years, cnn.thickness_rmse, type = "o", col = "red")

## Velocity RMSE
cnn.vel_rmse <- c()
enkfsa.vel_rmse <- c()
for (t in 1:years) {
    cnn.vel_rmse[t] <- rmse(rowMeans(cnn.velocity[[t + 1]]), true_velocities[s, , t + 1])
    enkfsa.vel_rmse[t] <- rmse(rowMeans(enkfsa.velocity[[t + 1]]), true_velocities[s, , t + 1])
}
plot(1:years, enkfsa.vel_rmse,
    type = "o", col = "blue",
    # ylim = c(30000, 100000),
    xlab = "Time (years)", ylab = "RMSE",
    main = paste0("RMSE of ice velocity over time for sample ", s)
) ## The scale is a bit ridiculous here
lines(1:years, cnn.vel_rmse, type = "o", col = "red")

fin_gl <- true_gl[s, years]
fin_gl_gridpt <- floor(fin_gl * 1000 / domain[length(domain)] * J)
cnn.bed_rmse <- rmse(cnn.bed[1:fin_gl_gridpt], true_bed[s, 1:fin_gl_gridpt])
enkfsa.bed_rmse <- rmse(rowMeans(enkfsa.bed_fin)[1:fin_gl_gridpt], true_bed[s, 1:fin_gl_gridpt])

cnn.fric_rmse <- rmse(cnn.fric[1:fin_gl_gridpt], true_fric[s, 1:fin_gl_gridpt])
enkfsa.fric_rmse <- rmse(rowMeans(enkfsa.friction_fin)[1:fin_gl_gridpt], true_fric[s, 1:fin_gl_gridpt])

rmse_df <- data.frame(cnn_rmse = c(cnn.bed_rmse, cnn.fric_rmse), 
            enkfsa_rmse = c(enkfsa.bed_rmse, enkfsa.fric_rmse))
rmse_df$parameter <- c("bed", "friction")
print(rmse_df)
