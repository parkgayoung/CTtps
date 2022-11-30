library(magrittr)
library(magick)
library(tesseract)
library(tidyverse)

# preprocessing
img <-
  image_read(here::here("analysis/data/raw_data/b-and-e-table-5.png")) %>%
  image_transparent("white", fuzz=82) %>%
  image_background("white") %>%
  image_negate() %>%
  image_crop(geometry_area(0, 0, 25, 25))

img

# read img and ocr
data <- img %>%
  image_ocr()

# some wrangling
tbl <-
data %>%
  stringi::stri_split(fixed = "\n") %>%
  .[[1]]


Monitor_Valley <- tbl[9:11] # Nevada - indirect bias
Monitor_Valley_tbl <-
tibble(mean = stringi::stri_split(Monitor_Valley[1], fixed = " ")[[1]],
       std =  stringi::stri_split(Monitor_Valley[2], fixed = " ")[[1]],
       n =    stringi::stri_split(Monitor_Valley[3], fixed = " ")[[1]]) %>%
  slice(-1) %>%
  mutate(cv = parse_number(std) / parse_number(mean))

Eastern_California <- tbl[13:15] # Eastern California - guided variation
Eastern_California_tbl <-
  tibble(mean = stringi::stri_split(Eastern_California[1], fixed = " ")[[1]],
         std =  stringi::stri_split(Eastern_California[2], fixed = " ")[[1]],
         n =    stringi::stri_split(Eastern_California[3], fixed = " ")[[1]]) %>%
  slice(-1) %>%
  mutate(cv = parse_number(std) / parse_number(mean))

Rose_Spring_Site <- tbl[17:19] # Eastern California - guided variation
Rose_Spring_Site_tbl <-
  tibble(mean = stringi::stri_split(str_squish(str_remove(Rose_Spring_Site[1], "no data")), fixed = " ")[[1]],
         std =  stringi::stri_split(str_squish(str_remove(Rose_Spring_Site[2], "no data")), fixed = " ")[[1]][1:7],
         n =    stringi::stri_split(str_squish(str_remove(Rose_Spring_Site[3], "no data")), fixed = " ")[[1]]) %>%
  slice(-1) %>%
  mutate(cv = parse_number(std) / parse_number(mean))

bind_rows(list(Monitor_Valley = Monitor_Valley_tbl,
              Eastern_California = Eastern_California_tbl,
              Rose_Spring_Site = Rose_Spring_Site_tbl),
          .id = "site") %>%
  ggplot() +
  aes(site, cv) +
  geom_boxplot() +
  ggbeeswarm::geom_quasirandom(size = 3,
                               alpha = 0.3) +
  theme_minimal()

