## BedMap data

# library(readxl)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(sp)
library(sf)
library(parallel)
# library(lattice)

setwd("~/SSA_model/CNN/real_data/")

print("Reading BedMap data...")
data_dir <- "./data/"
bed_data_2018 <- read.csv(paste0(data_dir, "/bedmap/BAS_2018_Thwaites_AIR_BM3.csv"), 
        skip = 18, header = TRUE)
bed_data_2019 <- read.csv(paste0(data_dir, "/bedmap/BAS_2019_Thwaites_AIR_BM3.csv"), 
        skip = 18, header = TRUE)
bed_data_2009 <- read.csv(paste0(data_dir, "/bedmap/CRESIS_2009_Thwaites_AIR_BM3.csv"), 
        skip = 18, header = TRUE)

bedmap2_bas <- read.csv(paste0(data_dir, "/bedmap/BAS_2004_BBAS_AIR_BM2.csv"), 
        skip = 18, header = TRUE)
bedmap2_agasea <- read.csv(paste0(data_dir, "/bedmap/UTIG_2004_AGASEA_AIR_BM2.csv"), 
        skip = 18, header = TRUE)

desc <- read.csv(paste0(data_dir, "/bedmap/UTIG_2004_AGASEA_AIR_BM2.csv"), nrows = 18)

fillValue <- -9999
bed_data_2009[bed_data_2009 == fillValue] <- NA
bed_data_2018[bed_data_2018 == fillValue] <- NA
bed_data_2019[bed_data_2019 == fillValue] <- NA
bedmap2_bas[bedmap2_bas == fillValue] <- NA
bedmap2_agasea[bedmap2_agasea == fillValue] <- NA

rename_cols <- function(df) {
  df <- rename_with(df, ~ sub("..m.", "", .x), ends_with("..m."))
  df <- df %>% rename(lon = longitude..degree_east., lat = latitude..degree_north.)
  return(df)
}

bed_data_2009 <- rename_cols(bed_data_2009)
bed_data_2018 <- rename_cols(bed_data_2018)
bed_data_2019 <- rename_cols(bed_data_2019)
bedmap2_bas <- rename_cols(bedmap2_bas)
bedmap2_agasea <- rename_cols(bedmap2_agasea)

# bed_df_2009 <- bed_data_2009 %>% select(lon, lat, bedrock_altitude) 
# bed_df_2018 <- bed_data_2018 %>% select(lon, lat, bedrock_altitude) 
# bed_df_2019 <- bed_data_2019 %>% select(lon, lat, bedrock_altitude) 

# coordinates(bed_elev_df) <- ~ lon + lat
# proj4string(bed_elev_df) <- CRS("+init=epsg:4326")
# bed_elev_xy <- spTransform(bed_elev_df, CRS("+init=epsg:3031")) 
# bed_elev_xy_df <- as.data.frame(bed_elev_xy)
# names(bed_elev_xy_df) <- c("bedrock_altitude", "x", "y")

print("Converting coordinates...")
lonlat_to_xy <- function(df) {
  coordinates(df) <- ~ lon + lat
  proj4string(df) <- CRS("+init=epsg:4326")
  xy <- spTransform(df, CRS("+init=epsg:3031")) 
  xy_df <- as.data.frame(xy)
  names(xy_df) <- c("bedrock_altitude", "x", "y")
  return(xy_df)
}

bed_df_2009 <- bed_data_2009 %>% select(lon, lat, bedrock_altitude) %>% lonlat_to_xy()
bed_df_2018 <- bed_data_2018 %>% select(lon, lat, bedrock_altitude) %>% lonlat_to_xy()
bed_df_2019 <- bed_data_2019 %>% select(lon, lat, bedrock_altitude) %>% lonlat_to_xy()
bedmap2_bas_df <- bedmap2_bas %>% select(lon, lat, bedrock_altitude) %>% lonlat_to_xy()

bedmap2_agasea_df <- bedmap2_agasea %>% 
                        select(lon, lat, bedrock_altitude) %>% 
                        filter(lon >= -180 & lon <= 180) %>%
                        lonlat_to_xy()

bed_df_2009$year <- 2009
bed_df_2018$year <- 2018
bed_df_2019$year <- 2019
bedmap2_bas_df$year <- 2004
bedmap2_agasea_df$year <- 2004
bed_data <- rbind(bed_df_2009, bed_df_2018, bed_df_2019, bedmap2_bas_df, bedmap2_agasea_df)
bed_data <- bed_data %>% arrange(year)
# bed_data$year <- as.factor(bed_data$year)
# test <- bed_elev_df %>% 
#   sf::st_as_sf( coords = c("lon", "lat") ) %>%
#   sf::st_set_crs( 4326 ) %>%      #current CRS is WSG84 (I think)  <--!! check!!
#   sf::st_transform( 3031 ) #%>%   #transform CRS to 32632
#   # sf::st_coordinates()  #export new coordinates

# # bed_elev_sp <- SpatialPointsDataFrame(coords = bed_elev_df[, c("lon", "lat")], data = bed_elev_df,
# #                                proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
# p <- test %>% ggplot() + 
#     geom_sf() + 
#     # scale_colour_distiller(palette = "RdBu", direction = 1, name = "Bedrock altitude (m)") + 
#     theme_bw()
# print(p)

# surf_elev_data <- readRDS(paste0(data_dir, "/surface_elev/height_thwaites.rds"))


## Basin boundary data
basin_data <- read_sf(paste0(data_dir, "/boundaries/Basins/Basins_Antarctica_v02.shp"))
thwaites_bound <- basin_data %>% filter(NAME == "Thwaites")
# thwaites_points <- as.data.frame(st_coordinates(st_geometry(thwaites_bound)))

## Grounding line data
gl_thwaites <- readRDS(paste0(data_dir, "/gl_thwaites.rds"))

## Flowline
flowline <- readRDS(paste0(data_dir, "/flowline_regrid.rds"))
        
print("Plotting...")
bed_plot <- bed_data %>% ggplot() + 
    geom_point(aes(x = x, y = y, colour = factor(year)), alpha = 0.5) + 
#     scale_colour_distiller(palette = "RdBu", direction = 1, name = "Bedrock altitude (m)") + 
#     geom_sf(data = basin_data, color = "black", fill = NA) +
    geom_sf(data = thwaites_bound, color = "#928484", fill = NA, linewidth = 1) +
    geom_path(data = flowline, aes(x = x, y = y), color = "royalblue") +
    geom_point(data = gl_thwaites, aes(x = X, y = Y), color = "black", size = 0.5) +
    theme_bw()

# png(paste0("./plots/bed/bed_elev_obs_locations.png"), width = 800, height = 800)
# print(bed_plot)
# dev.off()

delta <- 200 # take the closest bed observations within 200 m around each point along the flowline
## Map bed obs to flowline
avg_nearest_four <- function(pos) {
        near_pts <- bed_data %>% filter(
        x >= (pos[1] - delta) & x <= (pos[1] + delta),
        y >= (pos[2] - delta) & y <= (pos[2] + delta)
        )

        near_bed <- near_pts %>%
        mutate(dist = sqrt((x - pos[1])^2 + (y - pos[2])^2)) %>%
        filter(dist > 0) %>%
        arrange(dist) %>%
        slice_min(dist, n = 4) %>%
        select(bedrock_altitude) %>%
        summarise(bed_nearest = bedrock_altitude[1], bed_avg = mean(bedrock_altitude)) # , vx_avg = mean(vx), vy_avg = mean(vy)) # %>% # take average vx and vx as vx and vy at current position
        # as.numeric()

        return(near_bed)
}

# avg_nearest_four(as.numeric(flowline[1,]))
flowline_pos <- lapply(1:nrow(flowline), function(i) as.numeric(flowline[i, ]))
t12 <- system.time({
        bed_elev <- mclapply(flowline_pos, avg_nearest_four, mc.cores = 10L)
})

## Mark GL position
gl_pos <- readRDS(file = "./data/grounding_line/gl_pos.rds")
delta_gl <- 500
pts_near_gl <- flowline %>% filter(
        x >= (gl_pos[1] - delta_gl) & x <= (gl_pos[1] + delta_gl),
        y >= (gl_pos[2] - delta_gl) & y <= (gl_pos[2] + delta_gl)
        ) %>%
        mutate(dist = sqrt((x - gl_pos[1])^2 + (y - gl_pos[2])^2)) %>%
        # arrange(dist) %>%
        slice_min(dist, n = 1) #%>%

gl_ind <- which(flowline$x == pts_near_gl$x)

bed_elev_nearest <- sapply(bed_elev, function(x) x$bed_nearest)
# saveRDS(bed_elev_nearest, file = "./data/bed_elev_nearest.rds")

png(paste0("./plots/bed/bed_elev_nearest.png"), width = 800, height = 800)
plot(bed_elev_nearest, type = "l")
abline(v = gl_ind, lty = 2)
dev.off()


