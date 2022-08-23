# This script draws the map, this is figure 1 in the paper.


korean_archaeological_site_locations <-
  readr::read_csv(here::here("analysis/data/raw_data/sitelocation.csv"))

korean_archaeological_site_locations <-
  tibble::rowid_to_column(korean_archaeological_site_locations, "ID")

library(ggrepel)
# devtools::install_github('3wen/legendMap')
library(legendMap)

# load map tiles we previously saved
map <- readRDS(here::here("analysis/data/raw_data/map_tiles.rds"))

# map with site names next to the points
map_of_sites <-
  ggmap(map)  +
  geom_point(data = korean_archaeological_site_locations,
             aes(long_dd ,
                 lat_dd),
             colour = "red",
             size = 2) +
  # label points with site name
  geom_text_repel(data = korean_archaeological_site_locations,
                  aes(long_dd ,
                      lat_dd,
                      label = full_sitename),
                  size = 2,
                  bg.color = "white",
                  bg.r = 0.15) +
  theme_minimal(base_size = 6) +
  labs(x = "Longitute",
       y = "Latitude") +
  legendMap::scale_bar(
    # edit these numbers to select a suitable location
    # for the scale bar where it does not cover
    # important details on the map
    lon = 124.5,
    lat = 33.5,
    legend_size = 2,
    # distance of one section of scale bar, in km
    distance_lon = 50,
    # height of the scale bar, in km
    distance_lat = 5,
    # distance between scale bar and units, in km
    distance_legend = 20,
    # units of scale bar
    dist_unit = "km",
    # add the north arrow
    orientation = TRUE,
    # length of N arrow, in km
    arrow_length = 60,
    # distance between scale bar & base of N arrow, in km
    arrow_distance = 50,
    # size of letter 'N' on N arrow, in km
    arrow_north_size = 10) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

# inset world map, from
# https://stackoverflow.com/a/62104396/1036500
library(maps)
library(magrittr)
library(maptools)
library(raster)
library(ggplot2)

#Defining a general CRS
mycrs <- "+proj=longlat +datum=WGS84 +no_defs"

#Using the original maps package, then converting map into SpatialPolygons object
world <- maps::map("world", fill=TRUE,  plot = FALSE) %$%
  maptools::map2SpatialPolygons(., IDs=names,proj4string=CRS(mycrs))

#The resulting map has self intersection problems so any further operation reports errors; using buffers of width 0 is a fast fix
while(rgeos::gIsValid(world)==FALSE){
  world <- rgeos::gBuffer(world, byid = TRUE, width = 0, quadsegs = 5, capStyle = "ROUND")
}

#Dissolving polygon's limits
world <- raster::aggregate(world)

#Plotting. I add theme_void to your code to erase any axis, etc
worldplot <-
  ggplot() +
  geom_polygon(data = world,
               aes(x=long, y=lat,
                   group=group),
               fill='NA',
               color='black',
               size=0.1) +
  annotate("rect",
           xmin = 118,
           xmax = 137,
           ymin = 28,
           ymax = 47,
           fill = NA,
           colour = "red",
           size = 0.5
  ) +
  theme_void() +
  coord_fixed(1) +
  labs(x = "", y = "") +
  labs(x = NULL, y = NULL, title = NULL) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# combine site map and inset map
library(cowplot)
ggdraw() +
  draw_plot(map_of_sites) +
  draw_plot(plot = worldplot,
            x = 0.6, # x location of inset placement
            y = -0.09, # y location of inset placement
            width = .4, # Inset width
            height = .35, # Inset height
            scale = 0.5 # Inset scale
  ) +
  theme(panel.background = element_rect(fill='white', colour="white"),
        plot.background = element_rect(fill='white', colour="white"))

ggsave(here::here("analysis/figures/site-map.png"),
       width = 4.45,
       height = 4.5,
       units = "in")



