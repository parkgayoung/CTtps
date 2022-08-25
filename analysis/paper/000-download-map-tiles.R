# get map tiles from the internet, we only need to do this once and save them

library(ggmap)

# bounding box
# 39.534214, 124.159154 ... 39.657415, 129.138604
# 33.694662, 123.491142 ... 34.338454, 130.569658

# this takes a few minutes and needs an internet connection
map <-
  get_stamenmap(bbox = c(left = 124,
                         bottom = 33,
                         right = 	130,
                         top = 39),
                zoom = 10)

saveRDS(map, file = here::here("analysis/data/raw_data/map_tiles.rds"))
