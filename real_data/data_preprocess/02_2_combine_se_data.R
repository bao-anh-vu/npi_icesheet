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

## Missing pattern for surface elevation data
flowline <- qread(paste0(data_dir, "/flowline_regrid.qs"))
J <- nrow(flowline) # number of grid points
flowline_dist <- sqrt((flowline$x[2:J] - flowline$x[1:(J-1)])^2 + (flowline$y[2:J] - flowline$y[1:(J-1)])^2)
flowline_dist <- c(0, cumsum(na.omit(flowline_dist)))

surf_elev_mat <- matrix(NA, nrow = J, ncol = length(years))
for (i in 1:length(years)) {
    year <- years[i]
    # surf_elev_data <- readRDS(file = paste0(data_dir, "/surface_elev/surf_elev_", year, ".rds"))
    surf_elev_data <- qread(file = paste0(data_dir, "/surface_elev/surf_elev_", year, ".qs"))
    surf_elev_mat[, i] <- unlist(surf_elev_data)[1:J]
}

## There may be some mismatch between the surface elevation data and grounding line data
## Discard any observations that are beyond the grounding line
gl_pos <- qread(file = paste0(data_dir, "/grounding_line/gl_pos.qs"))
gl_ind <- gl_pos$ind
surf_elev_mat[(gl_ind+1):nrow(surf_elev_mat), ] <- NA
gl_km <- flowline_dist[gl_ind] / 1000

# qsave(surf_elev_mat, "./data/surface_elev/surf_elev_mat.qs")

## Plot surface elevation using ggplot
elev_df <- data.frame(dist = rep(flowline_dist/1000, times = length(years)),
                    year = rep(years, each = nrow(flowline)),
                    elev = c(surf_elev_mat))

ice_geometry_plot <- ggplot(elev_df, aes(x = dist, y = elev, group = year)) +
  geom_line(aes(alpha = year), color = "darkcyan") +
#   scale_alpha(range = c(0.1, 1), guide = 'none') +
#   geom_line(aes(x = x, y = z_b,  group = year, alpha = year), color = "darkcyan") +
  scale_alpha(range = c(0.1, 1), guide = 'none') +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "turquoise") +
#   geom_line(aes(y = b), color = "black") +
  labs(x = "Distance along flowline (km)", y = "Surface elevation (m)") +
  theme_bw() +
  theme(text = element_text(size = 24)) +
  xlim(c(0, 180)) 
#   ylim(c(-1000, 1800)) 

png(paste0("./plots/surface_elev/surf_elev.png"), width = 1500, height = 600, res = 150)
print(ice_geometry_plot)
dev.off()

surf_ev_missing_pattern <- ifelse(is.na(surf_elev_mat), 0, 1)
qsave(surf_ev_missing_pattern, "./data/surface_elev/missing_pattern.qs")

## Plot sủrface elevation missing pattern
surf_elev_df <- as.data.frame(surf_elev_mat)
colnames(surf_elev_df) <- years
surf_elev_df$gridpt <- 1:nrow(surf_elev_df)
surf_elev_df_long <- surf_elev_df %>%
    tidyr::pivot_longer(cols = -gridpt, names_to = "year", values_to = "surf_elev") %>%
    mutate(nonmissing = ifelse(is.na(surf_elev), 0, 1),
           year = as.integer(year))

surf_elev_mp_plot <- ggplot(surf_elev_df_long) +
    geom_tile(aes(x = gridpt, y = year, fill = factor(nonmissing))) +
    # scale_fill_manual(values = c("0" = "white", "1" = "blue"), labels = c("0" = "Missing", "1" = "Observed")) +
    theme_bw() +
    labs(x = "Grid Point", y = "Year", fill = "Data Status") +
    # scale_y_discrete(limits=rev) +
    ggtitle("Thwaites Glacier Surface Elevation Missing Pattern") +
    theme(plot.title = element_text(hjust = 0.5))

png("./plots/missing_pattern/surf_elev_missing_ptn.png", width = 800, height = 600)
print(surf_elev_mp_plot)
dev.off()