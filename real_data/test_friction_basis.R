## Test different basis lengthscales
setwd("~/SSA_model/CNN/real_data/")

rm(list = ls())
library(FRK)
library(qs)

source("./source/fit_basis.R")


## Presets
data_date <- "20241111" #"20241103" 
train_data_dir <- "./data/training_data"

## Domain
ssa_steady <- qread(file = paste0(train_data_dir, "/steady_state/steady_state_", data_date, ".qs"))
domain <- ssa_steady$domain
J <- length(domain)

## Friction simulations
fric_sims <- qread(file = paste0(train_data_dir, "/fric_sims_", data_date, ".qs"))

## Fit basis to log(friction)
lengthscales <- c(2.5e3, 5e3, 10e3, 15e3)
nbasis <- 150

output_list <- list()
# for (i in 1:length(lengthscales)) {
for (i in seq_along(lengthscales)) {
    output_list[[i]] <- fit_friction_basis(
      nbasis = nbasis,
      domain = domain,
      fric_arr = fric_sims,
      log_transform = T,
      lengthscale = lengthscales[i]
    )
}



sim <- 10
fitted <- lapply(output_list, function(x) exp(x$fitted_values[sim, ]))
basis_coefs <- lapply(output_list, function(x) x$basis_coefs[sim, ])

## Plot the fitted values for different lengthscales
png(file = paste0("./plots/friction/compare_lengthscales_", data_date, ".png"), width = 1000, height = 750)
plot_domain <- 1:J #1000
par(mfrow = c(length(lengthscales), 2))
cols <- RColorBrewer::brewer.pal(n = length(lengthscales), name = "Set1")
for (i in seq_along(lengthscales)) {
    plot(domain[plot_domain], fric_sims[sim, plot_domain], cex = 1.5,
    type = "l", main = paste0("Friction basis functions with lengthscale = ", lengthscales[i]/1e3, " km"),
    ylab = "Friction", xlab = "Distance (m)")
    lines(domain[plot_domain], fitted[[i]][plot_domain], col  = cols[i], lwd = 2)

    plot(basis_coefs[[i]], type = "l", cex = 1.5,
    main = paste0("Friction basis functions with lengthscale = ", lengthscales[i]/1e3, " km"),
    ylab = "Basis function", xlab = "Distance (m)")

}
legend("topright", legend = paste0("Lengthscale = ", lengthscales/1e3, " km"), col = cols, lty = 1, lwd = 2)
dev.off()

## Plot the basis functions for different lengthscales
# png(file = paste0("./plots/friction/compare_lengthscales_basis_", data_date, ".png"), width = 800, height = 1500)
# plot_domain <- 1:J #1000
# par(mfrow = c(length(lengthscales), 1))
# for (i in seq_along(lengthscales)) {
#     plot(basis_coefs[[i]], type = "l", 
#     main = paste0("Friction basis functions with lengthscale = ", lengthscales[i]/1e3, " km"),
#     ylab = "Basis function", xlab = "Distance (m)")
# }
# # legend("topright", legend = paste0("Lengthscale = ", lengthscales/1e3, " km"), col = "black", lty = 1, lwd = 2)
# dev.off()