
library(tidyverse)
library(ggplot2)
library(ggbeeswarm)

# import the cleaned data ready to work on
load(here::here("analysis/data/derived_data/data_main.RData"))

# Box plot for each attributes

# mutate data to make plot with ggplot
data_box <- gather(DF)

# make box plot for each attributes with better visulazation options with ggplot
ggplot(data_box, aes(key, value)) +
  geom_boxplot() +
  geom_beeswarm(alpha = 0.3,
                size = 0.8,
                cex = 0.7) +
  ylab("Length (cm)") +
  xlab("Variables") +
  theme_bw()

ggsave(here::here("analysis/figures/002-attributes-variation.png"),
       width = 4.45,
       height = 4.5,
       units = "in")

