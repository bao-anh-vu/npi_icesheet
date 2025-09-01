## Extract data along flowline

library(dplyr)
library(sf)
library(ggplot2)

extract_flowline <- T
regrid_flowline <- T

## Read (averaged) velocity data
setwd("~/SSA_model/CNN/real_data/")
data_dir <- "./data"

basin_data <- read_sf(paste0(data_dir, "/boundaries/Basins/Basins_Antarctica_v02.shp"))
thwaites_bound <- basin_data %>% filter(NAME == "Thwaites")
thwaites_points <- as.data.frame(st_coordinates(st_geometry(thwaites_bound)))
  
## Velocity data on Thwaites glacier
  thwaites_vel <- readRDS(paste0(data_dir, "/vel_thwaites.rds"))
  thwaites_vel <- thwaites_vel %>% filter(!is.na(v)) # remove NA velocities
  
  ## Choose a point in the middle of the glacier
  half_range_x <- min(thwaites_points$X) + 150e3 #400e3 
  half_range_y <- min(thwaites_points$Y) + 450e3 #500e3 

  x_unique <- unique(thwaites_vel$x) # x-coordinates
  delta <- x_unique[2] - x_unique[1] # grid size

  midpts <- thwaites_vel %>%
    filter(x >= (half_range_x - delta) & x <= half_range_x + delta & 
          y >= (half_range_y - delta) & y <= (half_range_y + delta))

  # head(midpts)
  chosen_pt <- midpts[1, ]

  ## Plot data around the chosen point
  margin <- 100e3
  near_pts <- thwaites_vel #%>% filter(x >= (chosen_pt$x - margin) & x <= (chosen_pt$x + margin) &
    # y >= (chosen_pt$y - margin) & y <= (chosen_pt$y + margin))

  png(paste0("./plots/flowline_start_point.png"), width = 800, height = 800)
  ggplot() +
    # geom_point(aes(x = x, y = y, colour = v)) +
    # scale_colour_distiller(
    #   palette = "Reds", direction = 1,
    #   name = "Velocity (m/a)"
    # ) +
    geom_sf(data = thwaites_bound, color = "black", fill = NA) +
    geom_point(data = data.frame(chosen_pt), aes(x = x, y = y), color = "black", size = 5) +
    theme_bw()
  dev.off()


if (extract_flowline) {

  stepsize <- 1

  # Extract data along flowline
  print("Extracting flowline...")
  pos <- as.numeric(chosen_pt[c("x", "y")])

  positions <- list()
  positions[[1]] <- pos
  # velocities <- list()
  # velocities[[1]] <- as.numeric(chosen_pt[c("v", "vx", "vy")])

  while (pos[1] > (min(thwaites_vel$x) + delta) & pos[2] > (min(thwaites_vel$y) + delta)) {
    # for (i in 1:10) {

    ## need to find 4 nearest grid points around the current position
    near_pts <- thwaites_vel %>% filter(
      x >= (pos[1] - delta) & x <= (pos[1] + delta),
      y >= (pos[2] - delta) & y <= (pos[2] + delta)
    )


    near_v <- near_pts %>%
      mutate(dist = sqrt((x - pos[1])^2 + (y - pos[2])^2)) %>%
      filter(dist > 0) %>%
      arrange(dist) %>%
      slice_min(dist, n = 4) %>%
      select(v, vx, vy) %>%
      summarise(v_nearest = v[1], v_avg = mean(v), vx_avg = mean(vx), vy_avg = mean(vy)) # %>% # take average vx and vx as vx and vy at current position
    # as.numeric()

    vx_vy <- as.numeric(near_v[c("vx_avg", "vy_avg")])

    if (is.na(vx_vy[1])) {
      break
    }

    cat("vx_vy = ", vx_vy, "\n")

    pos <- pos + stepsize * vx_vy
    positions <- c(positions, list(pos))
    print(pos)

    # velocities <- c(velocities, list(near_v))
  }

  saveRDS(positions, paste0(data_dir, "/flowline_positions.rds"))
} else {
  positions <- readRDS(paste0(data_dir, "/flowline_positions.rds"))
}

# positions
flowline_x <- sapply(positions, function(x) x[1])
flowline_y <- sapply(positions, function(x) x[2])

gl_thwaites <- readRDS(paste0(data_dir, "/gl_thwaites.rds"))

plot_flowline <- ggplot(thwaites_vel) +
    geom_point(aes(x = x, y = y, colour = log10(v))) +
    scale_colour_distiller(
      palette = "Reds", direction = 1,
      # limits = log10(c(1, 2000)),
      name = "Log-velocity (m/a)"
    ) +
    geom_point(data = gl_thwaites, aes(x = X, y = Y), color = "black", size = 0.2) +
    geom_sf(data = thwaites_bound, color = "black", fill = NA) +
    geom_point(data = data.frame(chosen_pt), aes(x = x, y = y), color = "black", size = 5) +
    geom_line(data = data.frame(x = flowline_x, y = flowline_y),
                aes(x = x, y = y), color = "cyan", linewidth = 2) +
    theme_bw()

# print("Saving flowline plot...")
# png(paste0("./plots/flowline2.png"), width = 800, height = 800)
# print(plot_flowline)
# plot(flowline_x, flowline_y, type = "l", col = "blue", xlab = "x", ylab = "y")
# dev.off()

## Re-grid
if (regrid_flowline) {
  sq_x_diff <- (flowline_x[2:length(flowline_x)] - flowline_x[1:(length(flowline_x) - 1)])^2
  sq_y_diff <- (flowline_y[2:length(flowline_y)] - flowline_y[1:(length(flowline_y) - 1)])^2
  dist <- sqrt(sq_x_diff + sq_y_diff)
  flowline_length <- sum(dist)

  J <- 2100 #2001
  grid_size <- flowline_length/J

  cumul_dist <- c(0, cumsum(dist))
  # cumul_dist <- lapply(1:5, function(i) cumsum(dist[1:i]))

  which.closest <- function(x, vec) {
    sub_vec <- vec[vec <= x]
    which.min(x - sub_vec)
  }

  ini_pt <- c(flowline_x[1], flowline_y[1])

  x_new <- c()
  y_new <- c()
  for (j in 1:J) {
    closest <- which.closest(grid_size * j, cumul_dist)
    leftover <- grid_size * j - cumul_dist[closest]

    l <- leftover / dist[closest]
    m <- 1 - l

    x_new[j] <- (m * flowline_x[closest] + l * flowline_x[closest + 1]) 
    y_new[j] <- (m * flowline_y[closest] + l * flowline_y[closest + 1])

  }

  x_new <- c(flowline_x[1], x_new) ## add initial point
  y_new <- c(flowline_y[1], y_new)

  ## Plot re-gridded flowline
  plot_range <- (length(flowline_x) - 10):length(flowline_x)
  plot(flowline_x[plot_range], flowline_y[plot_range], type = "o", xlab = "x", ylab = "y")
  points(x_new, y_new, col = "red", pch = 16)

  flowline_regrid <- data.frame(x = x_new, y = y_new)
  flowline_regrid <- na.omit(flowline_regrid)

  ## Truncate flowline to 2001 points
  # flowline_regrid <- flowline_regrid[1:2001, ]

  saveRDS(flowline_regrid, paste0(data_dir, "/flowline_regrid.rds"))
} else {
   flowline_regrid <- readRDS(paste0(data_dir, "/flowline_regrid.rds"))
}

plot_flowline_rg <- ggplot(thwaites_vel) +
    geom_point(aes(x = x, y = y, colour = log10(v))) +
    scale_colour_distiller(
      palette = "BuPu", direction = 1,
      # limits = log10(c(1, 2000)),
      name = "Log-velocity (m/a)"
    ) +
    geom_point(data = gl_thwaites, aes(x = X, y = Y), color = "salmon", size = 0.5) +
    geom_sf(data = thwaites_bound, color = "black", fill = NA) +
    geom_point(data = data.frame(chosen_pt), aes(x = x, y = y), color = "black", size = 5) +
    geom_point(data = data.frame(x = flowline_x, y = flowline_y),
                aes(x = x, y = y), color = "black") +
    theme_bw()

png(paste0("./plots/flowline_regridded.png"), width = 800, height = 800)
print(plot_flowline_rg)
dev.off()
