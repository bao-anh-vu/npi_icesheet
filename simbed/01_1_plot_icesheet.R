## Plot ice sheet profiles from simulations

setwd("~/SSA_model/CNN/simbed/")

rm(list = ls())

# library(keras)
# reticulate::use_condaenv("myenv", required = TRUE)
# library(tensorflow)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(abind)
library(parallel)
library(qs)

source("./source/surface_elev.R")

data_date <- "20240320" # "20220329"

arg <- commandArgs(trailingOnly = TRUE)
sets <- 1#:10 #1:50 
setf <- formatC(sets, width = 2, flag = "0")
# setsf <- paste0("sets", sets[1], "-", sets[lenhgth(sets)])#formatC(sets, width=2, flag="0

train_data_dir <- "./training_data"

## Steady state
ssa_steady <- readRDS(file = paste0(train_data_dir, "/initial_conds/ssa_steady_20220329.rds", sep = ""))
domain <- ssa_steady$domain

## Read thickness and velocity data
print("Reading surface data...")
surface_elev_arr <- qread(paste0(train_data_dir, "/true_surface_elevs_", setf, "_", data_date, ".qs"))

print("Reading ice thickness data...")
thickness_arr <- qread(paste0(train_data_dir, "/true_thicknesses_", setf, "_", data_date, ".qs"))

print("Reading friction data...")
fric_arr <- qread(paste0(train_data_dir, "/friction_arr_", setf, "_", data_date, ".qs"))

print("Reading bed data...")
bed_arr <- qread(paste0(train_data_dir, "/bed_arr_", setf, "_", data_date, ".qs"))

gl_arr <- qread(file = paste0(train_data_dir, "/gl_arr_", setf, "_", data_date, ".qs"))


s <- 1
H <- thickness_arr[s,,,1]
z <- surface_elev_arr[s,,,]
b <- t(bed_arr[s, ])
fric <- fric_arr[s,1]

years <- dim(surface_elev_arr)[3]

z_b <- z - H

png(paste0("./plots/temp/ice_sheet_prof.png"), width = 800, height = 600, res = 150)

matplot(domain/1000, z, type = "l", ylim = c(-2000, 2000), 
    xlab = "Domain (km)", ylab = "Elevation (m)")
matlines(domain/1000, z_b, col = "black")
#   abline(v = x[GL]/1000, col = "black", lty = 2)
#   abline(v = x[GL_position[1]]/1000, col = "salmon", lty = 2) # initial GL
abline(h = 0, col = "turquoise", lty = 2)
lines(domain/1000, b, col = "grey")
dev.off()

## Same plot but in ggplot
df <- data.frame(x = rep(domain/1000, years),
                year = rep(1:years, each = length(domain)),
                 z = as.vector(z),
                 z_b = as.vector(z_b),
                    b = rep(b, years))

df_small <- df %>% filter(year %% 5 == 0) # only plot every 5 years

ice_geometry_plot <- ggplot(df_small, aes(x = x, y = z, group = year)) +
  geom_line(aes(alpha = year), color = "darkcyan") +
#   scale_alpha(range = c(0.1, 1), guide = 'none') +
  geom_line(aes(x = x, y = z_b,  group = year, alpha = year), color = "darkcyan") +
  scale_alpha(range = c(0.1, 1), guide = 'none') +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "turquoise") +
  geom_line(aes(y = b), color = "black") +
  labs(x = "Domain (km)", y = "Elevation (m)") +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  xlim(c(200, 450)) +
  ylim(c(-1000, 1800)) 
  

# Plot surface elevation with fading lines corresponding to time                    
png(paste0("./plots/temp/ice_sheet_prof_sim", s, ".png"), width = 1000, height = 600, res = 150)
print(ice_geometry_plot)
dev.off()

## Also plot corresponding bed and friction profiles

gl <- ceiling(gl_arr[1, 1] / (domain[length(domain)] / 1000) * length(domain))

## Fitted friction and bed
friction_basis <- qread(file = paste0(train_data_dir, "/friction_basis_", setf, "_", data_date, ".qs"))

bed_basis <- qread(file = paste0(train_data_dir, "/bed_basis_", setf, "_", data_date, ".qs"))

friction_arr <- friction_basis$true_vals
fitted_friction <- friction_basis$fitted_values

bed_arr <- bed_basis$true_vals
fitted_bed <- bed_basis$fitted_values
bed_mean <- bed_basis$mean

fitted_fric_sim <- exp(fitted_friction[s, 1:gl])
friction_sim <- exp(friction_arr[s, 1:gl])
fric_df <- data.frame(
    domain = ssa_steady$domain[1:gl] / 1000, friction = friction_sim,
    fitted_fric = fitted_fric_sim
  )

  friction_plot <- fric_df %>% ggplot() +
    geom_line(aes(x = domain, y = fitted_fric), lwd = 1) +
    # geom_line(aes(x = domain, y = friction), col = "red", lwd = 1) +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    xlab("Domain (km)") +
    ylab(bquote("Friction (M Pa m"^"-1/3" ~ "yr"^"1/3" ~ ")"))

  bed_sim <- bed_arr[s, ] + bed_mean
  fitted_bed_sim <- fitted_bed[s, ] + bed_mean
  bed_df <- data.frame(domain = ssa_steady$domain / 1000, bed = bed_sim, fitted_bed = fitted_bed_sim)
  
  bed_plot <- bed_df %>% ggplot() + 
    geom_line(aes(x = domain, y = fitted_bed), lwd = 1) +
    # geom_line(aes(x = domain, y = bed), col = "red", lwd = 1) +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    xlim(0, 350) + 
    ylim(-1000, -300) +
    xlab("Domain (km)") +
    ylab("Bed (m)")

png(paste0("./plots/temp/param_sim", s, ".png"), width = 2000, height = 600, res = 150)
grid.arrange(grobs = list(friction_plot, bed_plot), nrow = 1, ncol = 2)
dev.off()