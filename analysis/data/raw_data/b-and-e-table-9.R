
library(tidyverse)
library(readxl)

# data from table 9 of
# https://www.researchgate.net/publication/277543080_Point_Typologies_Cultural_Transmission_and_the_Spread_of_Bow-and-Arrow_Technology_in_the_Prehistoric_Great_Basin

# Central Nevada: highly correlated - indirect bias
c_nv <- read_excel("analysis/data/raw_data/b-and-e-table-9.xlsx",
                   sheet = 1) %>%
  pivot_longer(-`C. NV` ) %>%
  dplyr::rename(dim = `C. NV`)

# Eastern California: poorly correlated - guided variation
e_ca <- read_excel("analysis/data/raw_data/b-and-e-table-9.xlsx",
                   sheet = 2) %>%
  pivot_longer(-`E. CA`) %>%
  dplyr::rename(dim = `E. CA`)

list(
c_nv = c_nv ,
e_ca = e_ca
) %>%
  bind_rows( .id = "location") %>%
  ggplot() +
  aes(location,
      value) +
  geom_boxplot() +
  ggbeeswarm::geom_quasirandom(size = 4,
                               alpha = 0.3) +
  theme_minimal()
