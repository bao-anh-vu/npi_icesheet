# Read velocity data

library(ncdf4)
library(CFtime)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(lattice)
library(sf) # for shapefiles
# library(sfheaders)
library(sp)
library(parallel)

# library(raster)

setwd("~/SSA_model/CNN/real_data/")

reread_data <- F
remap_vel <- F

data_dir <- "./data"

## Basin shapefile
# unzip("./boundaries/Basins_Antarctica_v02.zip", junkpaths = FALSE)
basin_data <- read_sf(paste0(data_dir, "/boundaries/Basins/Basins_Antarctica_v02.shp"))
thwaites_bound <- basin_data %>% filter(NAME == "Thwaites")

## Grounding line data
gl_thwaites <- readRDS(paste0(data_dir, "/gl_thwaites.rds"))
# ggplot(thwaites_bound) +
#   geom_sf(color = "black", fill = NA) +
#   theme_bw()
years <- 2000:2020

for (year in years) {
    cat("Mapping velocity data for year", year, "\n")

    if (reread_data) {

        ## Velocity data
        vel_data <- nc_open(paste0(data_dir, "/velocity/ITS_LIVE_velocity_120m_RGI19A_", year, "_v02.nc"))


        # flowline_x <- sapply(flowline, function(x) x[1])
        # flowline_y <- sapply(flowline, function(x) x[2])

        # plot(flowline_x, flowline_y, type = "l")

        v <- ncvar_get(vel_data, "v")
        x <- ncvar_get(vel_data, "x")
        y <- ncvar_get(vel_data, "y")
        #   vx <- ncvar_get(vel_data, "vx")
        #   vy <- ncvar_get(vel_data, "vy")
        #   v <- ncvar_get(vel_data, "v")
        grid <- expand.grid(x = x, y = y)
        grid$v <- c(v)
        #   grid$vx <- c(vx)
        #   grid$vy <- c(vy)

        ## Select only grid points that fall within Thwaites glacier
        thwaites_points <- as.data.frame(st_coordinates(st_geometry(thwaites_bound)))

        xmin <- min(x) # min(thwaites_points$X)
        xmax <- max(thwaites_points$X)
        ymin <- min(thwaites_points$Y)
        ymax <- max(thwaites_points$Y)

        # png(paste0("./plots/vel.png"), width = 800, height = 800)
        # print(plot)
        # dev.off()

        x_thwaites <- x[x >= xmin & x <= xmax]
        y_thwaites <- y[y >= ymin & y <= ymax]
        # grid_thwaites <- expand.grid(x = x_thwaites, y = y_thwaites)

        ## Try filtering grid points in v with coordinates in x_thwaites and y_thwaites
        thwaites_vel <- grid %>% filter(x %in% x_thwaites & y %in% y_thwaites)
        saveRDS(thwaites_vel, file = paste0(data_dir, "/velocity/vel_thwaites_", year, ".rds"))
    } else {
        thwaites_vel <- readRDS(paste0(data_dir, "/velocity/vel_thwaites_", year, ".rds"))
    }

    if (remap_vel) {
        cat("Mapping velocity to flowline for year ", year, "\n")
        ## Now map velocity to flowline

        avg_nearest_four <- function(pos) {
            near_pts <- thwaites_vel %>% filter(
                x >= (pos[1] - delta) & x <= (pos[1] + delta),
                y >= (pos[2] - delta) & y <= (pos[2] + delta)
            )

            near_v <- near_pts %>%
                mutate(dist = sqrt((x - pos[1])^2 + (y - pos[2])^2)) %>%
                filter(dist > 0) %>%
                arrange(dist) %>%
                slice_min(dist, n = 4) %>%
                select(v) %>%
                summarise(v_nearest = v[1], v_avg = mean(v)) # , vx_avg = mean(vx), vy_avg = mean(vy)) # %>% # take average vx and vx as vx and vy at current position
            # as.numeric()

            return(near_v)
        }


        delta <- 120 # grid size
        flowline <- readRDS(paste0(data_dir, "/flowline_regrid.rds"))
        flowline <- na.omit(flowline)
        flowline_pos <- lapply(1:nrow(flowline), function(i) as.numeric(flowline[i, ]))
        t12 <- system.time({
            velocities <- mclapply(flowline_pos, avg_nearest_four, mc.cores = 10L)
        })

        velocities_nearest <- sapply(velocities, function(x) x$v_nearest)
        velocities_avg <- sapply(velocities, function(x) x$v_avg)
        flowline$vel_nearest <- velocities_nearest
        flowline$vel_avg <- velocities_avg

        saveRDS(flowline, file = paste0(data_dir, "/velocity/flowline_vel_mapped_", year, ".rds"))
    } else {
        flowline <- readRDS(paste0(data_dir, "/velocity/flowline_vel_mapped_", year, ".rds"))
    }

    ## Plot velocity along flowline
    gl_pos <- readRDS(paste0(data_dir, "/grounding_line/gl_pos.rds"))
    
    delta <- 120 # grid size
    # flowline_dist <- sqrt((flowline$x[2:nrow(flowline)] - flowline$x[1:(nrow(flowline) - 1)])^2 +
    #     (flowline$y[2:nrow(flowline)] - flowline$y[1:(nrow(flowline) - 1)])^2)
    flowline$ind <- 1:nrow(flowline)
    gl_near_pts <- flowline %>% filter(
                x >= (gl_pos[1] - delta) & x <= (gl_pos[1] + delta),
                y >= (gl_pos[2] - delta) & y <= (gl_pos[2] + delta)) #%>% 
                # mutate(dist = sqrt((x - gl_pos[1])^2 + (y - gl_pos[2])^2)) 
                

    png(paste0("./plots/velocity/vel_flowline_1d_", year, ".png"), width = 800, height = 500)
    # plot(flowline$vel_nearest, type = "l")
    plot(1:nrow(flowline), flowline$vel_avg, type = "l", 
    main = paste0("Velocity for ", year), 
    xlab = "Point along flowline", ylab = "Velocity (m/a)")
    abline(v = gl_near_pts$ind, col = "salmon", lty = 2, lwd = 2)
    dev.off()

    tail_len <- 500
    tail_flowline <- tail(flowline, tail_len)
    # plot(tail_flowline$vel_avg, type = "l")

    flowline_tail_grid <- thwaites_vel %>% filter(x >= min(tail_flowline$x) & x <= max(tail_flowline$x) &
        y >= min(tail_flowline$y) & y <= max(tail_flowline$y))

    gl_thwaites_fl <- gl_thwaites %>% filter(X >= min(tail_flowline$x) & X <= max(tail_flowline$x) &
        Y >= min(tail_flowline$y) & Y <= max(tail_flowline$y))

    plot_vel_flowline <- ggplot(flowline_tail_grid) +
        geom_point(aes(x = x, y = y, colour = v)) +
        scale_colour_distiller(
            palette = "Reds", direction = 1,
            #  limits = c(min_vel, max_vel),
            name = "Velocity (m/a)"
        ) +
        geom_point(
            data = gl_thwaites_fl,
            aes(x = X, y = Y), color = "black", size = 0.2
        ) +
        geom_line(data = tail_flowline, mapping = aes(x = x, y = y), colour = "cyan") +
        coord_fixed() +
        # geom_line(tail_flowline, mapping = aes(x = x, y = y)) +
        # scale_colour_distiller(palette = "Reds", direction = 1,
        #                         size = 3, stroke = 0.5) +
        theme_bw()

    # png(paste0("./plots/velocity/vel_flowline_mapped_", tail_len, "_", year, ".png"), width = 800, height = 300)
    # print(plot_vel_flowline)
    # dev.off()
}






# vel_nearest <- sapply(velocities, function(x) x$v_nearest)
# vel_avg <- sapply(velocities, function(x) x$v_avg)

# # png(paste0("./plots/vel_flowline_nearest.png"), width = 800, height = 500)
# # plot(1:length(flowline), vel_nearest, type = "l", xlab = "Point along flowline", ylab = "Velocity (m/a)")
# # # lines(1:length(flowline))
# # dev.off()

# # png(paste0("./plots/vel_flowline_avg.png"), width = 800, height = 500)
# # plot(1:length(flowline), vel_avg, type = "l", xlab = "Grid point", ylab = "Velocity (m/a)")
# # dev.off()

# flowline_df <- data.frame(x = flowline_x, y = flowline_y)
# J <- length(flowline_df$x)
# flowline_tail <- flowline_df[(J - 20):J, ]
# # png(paste0("./plots/flowline_grid.png"), width = 800, height = 800)

# tail_gridpts <- thwaites_vel %>% filter(x >= min(flowline_tail$x) & x <= max(flowline_tail$x) &
#                                         y >= min(flowline_tail$y) & y <= max(flowline_tail$y))

# png(paste0("./plots/flowline_tail", year, ".png"), width = 800, height = 400)
# ggplot(tail_gridpts) +
#     geom_point(aes(x = x, y = y, colour = v)) +
#     scale_colour_distiller(palette = "Reds", direction = 1,
#                             #  limits = c(min_vel, max_vel),
#                             name = "Velocity (m/a)") +
#     geom_point(data = flowline_tail,
#     aes(x = x, y = y), color = "black", size = 3) +
#     coord_fixed() +
#     theme_bw()
# dev.off()

