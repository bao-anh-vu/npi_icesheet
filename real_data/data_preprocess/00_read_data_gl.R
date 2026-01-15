## Read grounding line data

# library(ncdf4)
# library(CFtime)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
# library(lattice) 
library(sf) # for shapefiles
library(qs)
# library(sfheaders)
# library(sp)

setwd("~/SSA_model/CNN/real_data/")

reread_data <- T

data_dir <- "./data"

## Basin shapefile
# unzip("./boundaries/Basins_Antarctica_v02.zip", junkpaths = FALSE)
basin_data <- read_sf(paste0(data_dir, "/boundaries/Basins/Basins_Antarctica_v02.shp"))
thwaites_bound <- basin_data %>% filter(NAME == "Thwaites")

## Grounding line data
gl_data <- read_sf(paste0(data_dir, "/grounding_line/InSAR_GL_Antarctica_v02.shp"))
gl_ERS <- gl_data %>% filter(SENSOR == "ERS")
gl_ERS_polar <- st_transform(gl_ERS, crs = st_crs(thwaites_bound))
gl_ERS_geo <- st_geometry(gl_ERS_polar)

gl_SEN <- gl_data %>% filter(SENSOR == "SENTINEL1A")
gl_SEN_polar <- st_transform(gl_SEN, crs = st_crs(thwaites_bound))
gl_SEN_geo <- st_geometry(gl_SEN)

# gl_coords <- gl_ERS_geo %>% st_cast("MULTIPOINT") %>% 
#         st_cast("POINT") %>%
#         as.data.frame()

# test2 <- test %>% st_cast("POINT") 

# orbit1 <- gl_coords[seq(1, nrow(gl_coords), 5), ]       
# orbit2 <- gl_coords[seq(2, nrow(gl_coords), 5), ]       
# orbit3 <- gl_coords[seq(3, nrow(gl_coords), 5), ]
# orbit4 <- gl_coords[seq(4, nrow(gl_coords), 5), ]

# saveRDS(orbit1, file = paste0(data_dir, "/grounding_line/gl_orbit1.rds"))
# saveRDS(orbit2, file = paste0(data_dir, "/grounding_line/gl_orbit2.rds"))
# saveRDS(orbit3, file = paste0(data_dir, "/grounding_line/gl_orbit3.rds"))
# saveRDS(orbit4, file = paste0(data_dir, "/grounding_line/gl_orbit4.rds"))

# gl_plot <- orbit1 %>% ggplot() +
#     geom_sf(color = "black", fill = NA) +
#     geom_sf(data = orbit2, color = "salmon", fill = NA, size = 0.1) +
#     geom_sf(data = orbit3, color = "lightblue", fill = NA, size = 0.1) +
#     geom_sf(data = orbit4, color = "goldenrod", fill = NA, size = 0.1) +
    
#     theme_bw()

# png(paste0("./plots/gl_test2.png"), width = 800, height = 800)
# # plot(orbit1)
# # points(orbit2, col = "red")
# print(gl_plot)
# dev.off()
# # v <- as.vector(test[1,])

# orbit1_df <- as.data.frame(st_coordinates(orbit1))
# orbit2_df <- as.data.frame(st_coordinates(orbit2))
# orbit3_df <- as.data.frame(st_coordinates(orbit3))
# orbit4_df <- as.data.frame(st_coordinates(orbit4))

# saveRDS(orbit1_df, file = paste0(data_dir, "/grounding_line/gl_orbit1_df.rds"))
# saveRDS(orbit2_df, file = paste0(data_dir, "/grounding_line/gl_orbit2_df.rds"))
# saveRDS(orbit3_df, file = paste0(data_dir, "/grounding_line/gl_orbit3_df.rds"))
# saveRDS(orbit4_df, file = paste0(data_dir, "/grounding_line/gl_orbit4_df.rds"))

## Flowline

flowline <- readRDS(paste0(data_dir, "/flowline_regrid.rds"))
# flowline_sf <- st_as_sf(flowline, coords = c("x", "y"), crs = st_crs(gl_ERS_polar))
flowline <- na.omit(flowline)
flowline_str <- st_linestring(as.matrix(flowline))
flowline_str <- st_sfc(flowline_str, crs = st_crs(gl_ERS_geo))
intersect_pts <- st_intersection(flowline_str, gl_ERS_geo)
intersect_pts <- intersect_pts %>% st_coordinates() %>% as.data.frame()# if (reread_data) {
intersect_pts$ind <- as.factor(1:nrow(intersect_pts))

p <- ggplot() +
    geom_sf(data = thwaites_bound, fill = NA, color = "black") +
    geom_sf(data = gl_ERS_geo, color = "red") +
    geom_sf(data = gl_SEN_geo, color = "yellow") +
    geom_sf(data = flowline_str, color = "black", fill = NA) +
    # geom_sf(data = flowline_str, color = "red", size = 0.1) +
    theme_bw()

png(paste0("./plots/gl_ERS.png"), width = 800, height = 800)
print(p)
dev.off()

## Select only grid points that fall within Thwaites glacier
thwaites_points <- as.data.frame(st_coordinates(st_geometry(thwaites_bound)))

# xmin <- min(x) # min(thwaites_points$X)
xmax <- max(thwaites_points$X)
ymin <- min(thwaites_points$Y)
ymax <- max(thwaites_points$Y)

gl_thwaites <- gl_ERS_geo %>% 
                    st_coordinates() %>% 
                    as.data.frame() %>%
                    filter(X <= xmax, Y >= ymin & Y <= ymax)

gl_thwaites_plot <- ggplot() + 
    geom_sf(data = thwaites_bound, fill = NA, color = "black") +
    geom_point(data = gl_thwaites, aes(x = X, y = Y), col = "black", size = 0.2) +
    geom_line(data = flowline, aes(x = x, y = y), col = "royalblue") +
    geom_point(data = intersect_pts, aes(x = X, y = Y, col = factor(ind))) +
    theme_bw()

png(paste0("./plots/gl_thwaites.png"), width = 800, height = 800)
print(gl_thwaites_plot)
dev.off()

tail_len <- 1000
tail_flowline <- tail(flowline, tail_len)

gl_thwaites_tail <- gl_thwaites %>% filter(X >= min(tail_flowline$x) & X <= max(tail_flowline$x) &
        Y >= min(tail_flowline$y) & Y <= max(tail_flowline$y))

gl_intst_plot <- ggplot() + 
    # geom_sf(data = thwaites_bound, fill = NA, color = "black") +
    geom_point(data = gl_thwaites_tail, aes(x = X, y = Y), col = "black", size = 0.2) +
    geom_line(data = tail_flowline, aes(x = x, y = y), col = "blue") +
    geom_point(data = intersect_pts, aes(x = X, y = Y, col = factor(ind))) +
    coord_fixed() +
    theme_bw()

png(paste0("./plots/gl_intersect.png"), width = 800, height = 800)
print(gl_intst_plot)
dev.off()

## Choose 2nd intersection point (more upstream) as gl position
gl_pos <- intersect_pts[2, ]
gl_pos <- as.numeric(gl_pos[, 1:2])
gl_pos

delta <- 120 # grid size
# flowline_dist <- sqrt((flowline$x[2:nrow(flowline)] - flowline$x[1:(nrow(flowline) - 1)])^2 +
#     (flowline$y[2:nrow(flowline)] - flowline$y[1:(nrow(flowline) - 1)])^2)
flowline$ind <- 1:nrow(flowline)
gl_near_pts <- flowline %>% filter(
            x >= (gl_pos[1] - delta) & x <= (gl_pos[1] + delta),
            y >= (gl_pos[2] - delta) & y <= (gl_pos[2] + delta)) %>% 
            mutate(dist = sqrt((x - gl_pos[1])^2 + (y - gl_pos[2])^2)) %>%
            slice_min(dist, n = 1) 
                
gl_df <- data.frame(gl_x = gl_pos[1], gl_y = gl_pos[2],
                    gl_nearest_gridpt_x = gl_near_pts$x, gl_nearest_gridpt_y = gl_near_pts$y,
                    ind = gl_near_pts$ind) 

# saveRDS(gl_pos, file = "./data/grounding_line/gl_pos.rds")
qsave(gl_df, file = "./data/grounding_line/gl_pos.qs")

# } else {
#     gl_thwaites <- readRDS(paste0(data_dir, "/gl_thwaites.rds"))
# }

# gl_plot <- ggplot(gl_thwaites) +
#     geom_point(aes(x = X, y = Y), col = "black", size = 0.2) +
#     geom_sf(thwaites_bound, fill = NA, color = "black") +
#     theme_bw()

# png(paste0("./plots/gl_orbit3_thwaites.png"), width = 800, height = 800)
# print(gl_plot)
# dev.off()
