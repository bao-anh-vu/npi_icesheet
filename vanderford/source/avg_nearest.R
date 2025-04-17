## Map bed obs to flowline
avg_nearest <- function(pos, radius, df, var, neighbours = 4) {
    
        near_pts <- df %>% 
        filter(
            x >= (pos[1] - radius) & x <= (pos[1] + radius),
            y >= (pos[2] - radius) & y <= (pos[2] + radius)
            ) %>%
        mutate(dist = sqrt((x - pos[1])^2 + (y - pos[2])^2)) %>%
        filter(dist > 0) %>%
        arrange(dist) %>%
        slice_min(dist, n = neighbours) %>%
        select(all_of(var)) %>%
        summarise(nearest = all_of(var)[1], avg = mean(all_of(var))) # , vx_avg = mean(vx), vy_avg = mean(vy)) # %>% # take average vx and vx as vx and vy at current position
        # as.numeric()

        return(near_pts)
}
