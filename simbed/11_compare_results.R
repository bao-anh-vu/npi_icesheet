## Compare results from EnKF conditional on CNN posterior and state-augmented EnKF

library(ggplot2)
library(mvtnorm)
library(qs)

setwd("/home/babv971/SSA_model/CNN/simbed/")
rm(list = ls())

source("./source/surface_elev.R")

## Presets
data_date <- "20220329" # "20230518"
output_date <- "20240320" # "20240518"

use_missing_pattern <- T
use_cov_taper <- F
# use_basis_funs <- T
plot_ice_thickness <- T
plot_velocity <- T
plot_bed <- T
plot_friction <- T
plot_gl <- T

## SSA model info
ssa_steady <- readRDS(file = paste("./training_data/initial_conds/ssa_steady_20220329.rds", sep = ""))
# reference <- readRDS(file = paste("./training_data/initial_conds/reference_20220329.rds", sep = ""))
domain <- ssa_steady$domain
J <- length(domain)

## Read bed observations
bed_obs <- readRDS(file = paste0("./training_data/bed_obs_", output_date, ".rds"))
bed_obs_df <- data.frame(location = domain[bed_obs$locations] / 1000, bed_elev = bed_obs$obs)

## Read bed and friction from NN output
sets <- 1:50 # 10
setsf <- paste0("sets", sets[1], "-", sets[length(sets)])

## Test sample index
set.seed(2024)
chosen_test_samples <- sample(1:500, 50)
set.seed(NULL)
sample_ind <- 1 #commandArgs(trailingOnly = TRUE) # sample index
s <- chosen_test_samples[sample_ind] # the actual number of the sample in the test set

years <- 20
save_points <- c(1, floor(years/2) + 1, years+1) #c(1, 11, 21)

Ne <- 1000 # Ensemble size for EnKFSA
Ne_enkf <- 500 # Ensemble size for EnKF-CNN

## Read test data
if (use_missing_pattern) {
    data_dir <- paste0("./training_data/", setsf, "/missing")
} else {
    data_dir <- paste0("./training_data/", setsf, "/nonmissing")
}

test_data <- qread(file = paste0(data_dir, "/test_data_", output_date, ".qs"))

true_surface_elevs <- test_data$true_surface_elevs_test
true_thicknesses <- test_data$true_thickness_test
true_velocities <- test_data$true_velocity_test

true_bed <- test_data$true_bed
true_fric <- exp(test_data$true_fric)
true_gl <- test_data$grounding_line * test_data$sd_gl + test_data$mean_gl


## Read CNN predictions
if (use_missing_pattern) {
    cnn.output_dir <- paste0("./output/posterior/", setsf, "/missing")
    cnn_enkf.output_dir <- paste0("./output/posterior/", setsf, "/missing/sample", sample_ind)
    enkfsa.output_dir <- paste0("./output/stateaug/", setsf, "/missing/sample", sample_ind)    
} else {
    cnn.output_dir <- paste0("./output/posterior/", setsf, "/nonmissing")
    cnn_enkf.output_dir <- paste0("./output/posterior/", setsf, "/nonmissing/sample", sample_ind)
    enkfsa.output_dir <- paste0("./output/stateaug/", setsf, "/nonmissing/sample", sample_ind)    
}

if (use_cov_taper) {
    cnn_enkf.output_dir <- paste0(cnn_enkf.output_dir, "/taper")
    enkfsa.output_dir <- paste0(enkfsa.output_dir, "/taper")
} else {
    cnn_enkf.output_dir <- paste0(cnn_enkf.output_dir, "/no_taper")
    enkfsa.output_dir <- paste0(enkfsa.output_dir, "/no_taper")
}

pred_fric <- qread(file = paste0(cnn.output_dir, "/pred_fric_", output_date, ".qs"))
pred_bed <- qread(file = paste0(cnn.output_dir, "/pred_bed_", output_date, ".qs"))
pred_gl <- qread(file = paste0(cnn.output_dir, "/pred_gl_", output_date, ".qs"))

# pred_fric <- readRDS(file = paste0(cnn.output_dir, "/pred_fric_", output_date, ".rds"))
# pred_bed <- readRDS(file = paste0(cnn.output_dir, "/pred_bed_", output_date, ".rds"))
# pred_gl <- readRDS(file = paste0(cnn.output_dir, "/pred_gl_", output_date, ".rds"))


print("Reading posterior samples from CNN...")
# fric_samples_ls <- readRDS(file = paste0(cnn.output_dir, "/fric_post_samples_", output_date, ".rds"))
# bed_samples_ls <- readRDS(file = paste0(cnn.output_dir, "/bed_post_samples_", output_date, ".rds"))
# gl_samples_ls <- readRDS(file = paste0(cnn.output_dir, "/gl_post_samples_", output_date, ".rds"))

fric_samples_ls <- qread(file = paste0(cnn.output_dir, "/fric_post_samples_", output_date, ".qs"))
bed_samples_ls <- qread(file = paste0(cnn.output_dir, "/bed_post_samples_", output_date, ".qs"))
gl_samples_ls <- qread(file = paste0(cnn.output_dir, "/gl_post_samples_", output_date, ".qs"))


## Scaling units for friction coefficients
# secpera <- 31556926
# fric_scale <- 1e6 * secpera^(1 / 3)

## Read EnKF-CNN results
cnn.thickness <- qread(file = paste0(cnn_enkf.output_dir, "/enkf_thickness_sample", sample_ind, "_Ne", Ne_enkf, "_", output_date, ".qs", sep = ""))
cnn.velocity <- qread(file = paste0(cnn_enkf.output_dir, "/enkf_velocities_sample", sample_ind, "_Ne", Ne_enkf, "_", output_date, ".qs", sep = ""))

cnn.thickness <- lapply(cnn.thickness, as.matrix)
cnn.velocity <- lapply(cnn.velocity, as.matrix)

# test <- lapply(cnn.thickness, base::rowMeans)

cnn.bed <- pred_bed[, s]
cnn.fric <- pred_fric[, s]
cnn.gl <- pred_gl[s, ]


cnn.bed_q <- apply(bed_samples_ls[[s]], 1, quantile, probs = c(0.025, 0.975))
cnn.fric_q <- apply(fric_samples_ls[[s]], 1, quantile, probs = c(0.025, 0.975))
cnn.bed_lq <- cnn.bed_q[1, ]
cnn.bed_uq <- cnn.bed_q[2, ]
cnn.fric_lq <- cnn.fric_q[1, ]
cnn.fric_uq <- cnn.fric_q[2, ]


## Read EnKF-StateAug results
enkfsa.thickness <- qread(file = paste0(enkfsa.output_dir, "/enkf_thickness_sample", sample_ind, "_Ne", Ne, "_", output_date,  ".qs", sep = ""))
enkfsa.bed <- qread(file = paste0(enkfsa.output_dir, "/enkf_bed_sample", sample_ind, "_Ne", Ne, "_", output_date,  ".qs", sep = ""))
enkfsa.friction <- qread(file = paste0(enkfsa.output_dir, "/enkf_friction_sample", sample_ind, "_Ne", Ne, "_", output_date,  ".qs", sep = ""))
enkfsa.velocity <- qread(file = paste0(enkfsa.output_dir, "/enkf_velocities_sample", sample_ind, "_Ne", Ne, "_", output_date,  ".qs", sep = ""))

enkfsa.thickness <- lapply(enkfsa.thickness, as.matrix)
enkfsa.bed <- lapply(enkfsa.bed, as.matrix)
enkfsa.friction <- lapply(enkfsa.friction, as.matrix)
enkfsa.velocity <- lapply(enkfsa.velocity, as.matrix)

#############################
##    Plotting results     ## 
#############################

if (use_missing_pattern) {
    plot_dir <- paste0("./plots/combined/Ne", Ne, "/missing/sample", sample_ind)
} else {
    plot_dir <- paste0("./plots/combined/Ne", Ne, "/nonmissing/sample", sample_ind)
}

if (!dir.exists(plot_dir)) {
    dir.create(paste0(plot_dir))
}

if (use_cov_taper) {
    plot_dir <- paste0(plot_dir, "/taper")
} else {
    plot_dir <- paste0(plot_dir, "/no_taper")
}

if (!dir.exists(plot_dir)) {
    dir.create(paste0(plot_dir))
} #else { # delete all previously saved plots
#     unlink(paste0(plot_dir, "/*"))
# }

plot_times <- save_points #seq(1, years + 1, 1)

## Ice thickness plot
if (plot_ice_thickness) {
    print("Plotting ice thickness...")
    for (i in 1:length(plot_times)) {
        t <- plot_times[i]
        cnn.ens_t <- cnn.thickness[[i]]
        cnn.covmat <- 1 / (ncol(cnn.ens_t) - 1) * tcrossprod(cnn.ens_t - rowMeans(cnn.ens_t)) # diag(enkf_covmats[[t]])
        cnn.lower <- rowMeans(cnn.ens_t[1:J, ]) + qnorm(0.025) * sqrt(diag(cnn.covmat)[1:J])
        cnn.upper <- rowMeans(cnn.ens_t[1:J, ]) + qnorm(0.975) * sqrt(diag(cnn.covmat)[1:J])

        enkfsa.ens_t <- enkfsa.thickness[[i]]
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
    print("Plotting velocity...")
    for (i in 1:length(plot_times)) {
        t <- plot_times[i]
        cnn.ens_t <- cnn.velocity[[i]]
        cnn.covmat <- 1 / (ncol(cnn.ens_t) - 1) * tcrossprod(cnn.ens_t - rowMeans(cnn.ens_t)) # diag(enkf_covmats[[t]])
        cnn.lower <- rowMeans(cnn.ens_t[1:J, ]) + qnorm(0.025) * sqrt(diag(cnn.covmat)[1:J])
        cnn.upper <- rowMeans(cnn.ens_t[1:J, ]) + qnorm(0.975) * sqrt(diag(cnn.covmat)[1:J])

        enkfsa.ens_t <- enkfsa.velocity[[i]]
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
enkfsa.mean_bed <- rowMeans(enkfsa.bed[[length(save_points)]])
enkfsa.gl <- sapply(enkfsa.mean_thickness, gl_migrate, b = enkfsa.mean_bed) # grid pts
enkfsa.gl <- domain[enkfsa.gl] / 1000 # convert to km

cnn.mean_thickness <- lapply(cnn.thickness, rowMeans)
enkf_cnn.gl <- sapply(cnn.mean_thickness, gl_migrate, b = cnn.bed) # grid pts
enkf_cnn.gl <- domain[enkf_cnn.gl] / 1000 # convert to km

if (plot_gl) {
    png(paste0(plot_dir, "/gl_comparison.png"), width = 2000, height = 1000, res = 300)
    plot(true_gl[s, ], type = "l", col = "black", lwd = 2, ylim = c(360, 363),
        xlab = "Time (years)", ylab = "Grounding line (km)")
    lines(cnn.gl, col = "red", lwd = 2) # CNN prediction
    lines(enkfsa.gl, col = "blue", lwd = 2) # EnKF-SA prediction
    lines(enkf_cnn.gl, col = "goldenrod", lwd = 2) # EnKF-CNN prediction
    legend("bottomright", legend = c("True GL", "CNN", "EnKF-SA", "EnKF-CNN"), col = c("black", "red", "blue", "goldenrod"), lwd = 2)
    dev.off()
}

## Bed
if (plot_bed) {
    print("Plotting bed...")
    # enkfsa.bed[[years+1]] <- enkfsa.bed[[1]]#enkfsa.bed
    t <- save_points[[length(save_points)]] #years + 1
    enkfsa.ens_t <- enkfsa.bed[[length(save_points)]]
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
        enkfsa_lower = enkfsa.lower, enkfsa_upper = enkfsa.upper,
        cnn_s1 = bed_samples_ls[[s]][, 1],
        cnn_s2 = bed_samples_ls[[s]][, 2],
        cnn_s3 = bed_samples_ls[[s]][, 3],
        cnn_s4 = bed_samples_ls[[s]][, 4],
        cnn_s5 = bed_samples_ls[[s]][, 5],
        enkf_s1 = enkfsa.ens_t[, 1],
        enkf_s2 = enkfsa.ens_t[, 2],
        enkf_s3 = enkfsa.ens_t[, 3],
        enkf_s4 = enkfsa.ens_t[, 4],
        enkf_s5 = enkfsa.ens_t[, 5]
    )
    true_bed.df <- data.frame(domain = domain / 1000, bed = true_bed[s, ])

    # Plot title
    # title <- paste("t = ", t - 1, "a")
    bed_plot <- ggplot() +
        # geom_ribbon(
        #     aes(domain, ymin = enkfsa_lower, ymax = enkfsa_upper),
        #     fill = "blue", alpha = 0.2
        # ) +
        # geom_ribbon(
        #     aes(domain, ymin = cnn_lower, ymax = cnn_upper),
        #     fill = "red", alpha = 0.2
        # ) +
        theme_bw() +
        geom_line(data = bed.df, aes(domain, cnn_s1), colour = "red", alpha = 0.3, lwd = 0.5) +
        geom_line(data = bed.df, aes(domain, cnn_s2), colour = "red", alpha = 0.3, lwd = 0.5) +
        geom_line(data = bed.df, aes(domain, cnn_s3), colour = "red", alpha = 0.3, lwd = 0.5) +
        geom_line(data = bed.df, aes(domain, cnn_s4), colour = "red", alpha = 0.3, lwd = 0.5) +
        geom_line(data = bed.df, aes(domain, cnn_s5), colour = "red", alpha = 0.3, lwd = 0.5) +
        geom_line(data = bed.df, aes(domain, enkf_s1), colour = "blue", alpha = 0.3, lwd = 0.5) +
        geom_line(data = bed.df, aes(domain, enkf_s2), colour = "blue", alpha = 0.3, lwd = 0.5) +
        geom_line(data = bed.df, aes(domain, enkf_s3), colour = "blue", alpha = 0.3, lwd = 0.5) +
        geom_line(data = bed.df, aes(domain, enkf_s4), colour = "blue", alpha = 0.3, lwd = 0.5) +
        geom_line(data = bed.df, aes(domain, enkf_s5), colour = "blue", alpha = 0.3, lwd = 0.5) +
        geom_line(data = true_bed.df, aes(domain, bed)) +
        geom_line(data = bed.df, aes(domain, enkfsa_mean), colour = "blue") +
        geom_line(data = bed.df, aes(domain, cnn_mean), colour = "red") +
        geom_vline(xintercept = true_gl[s, years], lty = 2, colour = "black") +
        geom_point(data = bed_obs_df, aes(location, bed_elev), colour = "turquoise") +
        xlab("Domain (km)") +
        ylab("Bed elevation (m)") +
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
# enkfsa.friction_fin <- exp(enkfsa.friction[[length(save_points)]])

if (plot_friction) {
    print("Plotting friction...")
    enkfsa.ens_t <- enkfsa.friction[[length(save_points)]]
    enkfsa.covmat <- 1 / (ncol(enkfsa.ens_t) - 1) * tcrossprod(enkfsa.ens_t - rowMeans(enkfsa.ens_t)) # diag(enkf_covmats[[t]])
    enkfsa.lower <- rowMeans(enkfsa.ens_t[1:J, ]) + qnorm(0.025) * sqrt(diag(enkfsa.covmat)[1:J])
    enkfsa.upper <- rowMeans(enkfsa.ens_t[1:J, ]) + qnorm(0.975) * sqrt(diag(enkfsa.covmat)[1:J])

    enkfsa.log_fric_samples <- rmvnorm(1000, rowMeans(enkfsa.ens_t), as.matrix(enkfsa.covmat))
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
        enkfsa_mean = enkfsa.fric_mean,
        cnn_lower = cnn.fric_lq, cnn_upper = cnn.fric_uq,
        enkfsa_lower = lower, enkfsa_upper = upper,
        cnn_s1 = fric_samples_ls[[s]][, 1],
        cnn_s2 = fric_samples_ls[[s]][, 2],
        cnn_s3 = fric_samples_ls[[s]][, 3],
        cnn_s4 = fric_samples_ls[[s]][, 4],
        cnn_s5 = fric_samples_ls[[s]][, 5],
        enkf_s1 = enkfsa.fric_samples[1, ],
        enkf_s2 = enkfsa.fric_samples[2, ],
        enkf_s3 = enkfsa.fric_samples[3, ],
        enkf_s4 = enkfsa.fric_samples[4, ],
        enkf_s5 = enkfsa.fric_samples[5, ]
    )

    true_fric.df <- data.frame(domain = domain / 1000, fric = true_fric[s, ])

    # Plot title
    # title <- paste("t = ", t - 1, "a")
    friction_plot <- ggplot(friction.df) +
        # geom_ribbon(
        #     aes(domain, ymin = enkfsa_lower, ymax = enkfsa_upper),
        #     fill = "blue", alpha = 0.2
        # ) +
        # geom_ribbon(
        #     aes(domain, ymin = cnn_lower, ymax = cnn_upper),
        #     fill = "red", alpha = 0.2
        # ) +
        theme_bw() +
        geom_line(data = friction.df, aes(domain, cnn_s1), colour = "red", alpha = 0.2, lwd = 0.5) +
        geom_line(data = friction.df, aes(domain, cnn_s2), colour = "red", alpha = 0.2, lwd = 0.5) +
        geom_line(data = friction.df, aes(domain, cnn_s3), colour = "red", alpha = 0.2, lwd = 0.5) +
        geom_line(data = friction.df, aes(domain, cnn_s4), colour = "red", alpha = 0.2, lwd = 0.5) +
        geom_line(data = friction.df, aes(domain, cnn_s5), colour = "red", alpha = 0.2, lwd = 0.5) +
        geom_line(data = friction.df, aes(domain, enkf_s1), colour = "blue", alpha = 0.2, lwd = 0.5) +
        geom_line(data = friction.df, aes(domain, enkf_s2), colour = "blue", alpha = 0.2, lwd = 0.5) +
        geom_line(data = friction.df, aes(domain, enkf_s3), colour = "blue", alpha = 0.2, lwd = 0.5) +
        geom_line(data = friction.df, aes(domain, enkf_s4), colour = "blue", alpha = 0.2, lwd = 0.5) +
        geom_line(data = friction.df, aes(domain, enkf_s5), colour = "blue", alpha = 0.2, lwd = 0.5) +
        geom_line(data = true_fric.df, aes(domain, fric)) +
        geom_line(data = friction.df, aes(domain, enkfsa_mean), colour = "blue") +
        geom_line(data = friction.df, aes(domain, cnn_mean), colour = "red") +
        geom_vline(xintercept = true_gl[s, years], lty = 2, colour = "black") +
        xlim(0, 400) +
        ylim(0, 0.2) +
        xlab("Domain (km)") +
        ylab(bquote("Friction (M Pa m"^{-1/3}~"a"^{1/3}~")")) +
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

## Quick RMSE calc
rmse <- function(estimated, true) {
    stopifnot(length(estimated) == length(true))
    sqrt(mean((estimated - true)^2))
}

# head(friction.df)

gl <- true_gl[s, years]
gl_ind <- floor(gl / 800 * 2001)
cnn.fric_rmse <- rmse(friction.df$cnn_mean[1:gl_ind], true_fric[s, 1:gl_ind])
enkfsa.fric_rmse <- rmse(friction.df$enkfsa_mean[1:gl_ind], true_fric[s, 1:gl_ind])

cnn.bed_rmse <- rmse(bed.df$cnn_mean, true_bed[s, ])
enkfsa.bed_rmse <- rmse(bed.df$enkfsa_mean, true_bed[s, ])


rmse.df <- data.frame(
    method = c("CNN", "EnKF-SA"),
    bed_rmse = c(cnn.bed_rmse, enkfsa.bed_rmse),
    fric_rmse = c(cnn.fric_rmse, enkfsa.fric_rmse)
)
print(rmse.df)
