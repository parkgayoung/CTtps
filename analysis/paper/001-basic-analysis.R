load("data_main.RData")
library(tidyverse)
library(ggplot2)
library(ggbeeswarm)
#Box plot for each attributes

box1 <- boxplot(DF,
                main = "Variation for each attributes",
                col = "gold",
                xlab = "Attributes",
                ylab = "Length (cm)")

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

ggsave(here::here("analysis/figures/003-attributes-variation.png"),
       width = 4.45,
       height = 4.5,
       units = "in")

# take a look at the distribution of each variable
DF %>%
  pivot_longer(everything()) %>%
  ggplot() +
  aes(value) +
  geom_density() +
  facet_wrap( ~ name) +
  theme_minimal()

# fit distributions for each of these, so we can simulate from those distributions
library("fitdistrplus")
library(distributions3)

# what distributions will explore to see if they match with our data?
the_distributions <- c("weibull",   "norm" )
the_distributions3 <- c(Weibull,  Normal)
names(the_distributions3) <- the_distributions

fit_distributions <-
DF %>%
  pivot_longer(everything())  %>%
  nest(-name) %>%
  # compute fit of each variable with a selection of distributions
  mutate(dist_out_w = map(data, ~fitdist(unlist(.x), the_distributions[1]))) %>%
  mutate(dist_out_l = map(data, ~fitdist(unlist(.x), the_distributions[2]))) %>%
 # mutate(dist_out_g = map(data, ~fitdist(unlist(.x), the_distributions[3]))) %>%
 # mutate(dist_out_n = map(data, ~fitdist(unlist(.x), the_distributions[4]))) %>%
  pivot_longer(cols = starts_with("dist"),
               names_to = "dist") %>%
  nest(-name) %>%
  # extract goodness of fit statistics
  mutate(gofstat_out = map(data, ~gofstat(.x$value,
                                          fitnames = the_distributions))) %>%
  # AIC, in particular, which distribution has the lowest value (so, best fit)
  mutate(gofstat_out_aic = map_chr(gofstat_out,
                                   ~names(which.min(.x$aic)))) %>%
  unnest(data) %>%
  distinct(name, .keep_all = TRUE) %>%
  # estimate parameters of distribution to match our data
  mutate(fitdist_out = map2(data, gofstat_out_aic,
                            ~fitdist(deframe(.x), .y))) %>%
  mutate(fit = map(fitdist_out, ~t(as.data.frame(.x$estimate)))) %>%
  unnest_wider(fit) %>%
  # create a function to generate a random distribution, using the distribution
  # name we deduced earlier, and with the parameter values we got earlier
  mutate(rdists_fn = map(gofstat_out_aic, ~the_distributions3[.x])) %>%
  # generate the simulated distribution of the variable
  mutate(rdists_val = map(rdists_fn, ~unname(.x)[[1]](`...1`, `...2`) %>% random(n = 1000)))

# compare with plots
plot_comparisons <-
pmap(list(fit_distributions$data,
          fit_distributions$rdists_val,
          fit_distributions$name),
     function(.x, .y, .z)
     ggplot() +
     geom_density(data = .x, aes(value), colour = "red")  +
     geom_density(data = as_tibble(.y), aes(value), colour = "green") +
     ggtitle(label = .z)
    )

library(cowplot)
plot_grid(plotlist = plot_comparisons)
# not very good match! terrible

#------------------------------------------------------------
# another approach, from https://stackoverflow.com/a/4809744/1036500
library(SuppDists)

fit_distributions2 <-
  DF %>%
  pivot_longer(everything())  %>%
  nest(-name) %>%
  # extract parms
  mutate(parms = map(data, ~JohnsonFit(.x$value, moment="quant"))) %>%
  mutate(rjohn = map(data, ~rJohnson(100, parms)))

# compare with plots
plot_comparisons2 <-
  pmap(list(fit_distributions2$data,
            fit_distributions2$rjohn,
            fit_distributions2$name),
       function(.x, .y, .z)
         ggplot() +
         geom_density(data = .x, aes(value),
                      colour = "red")  +
         geom_density(data = as_tibble(.y),
                      aes(value),
                      colour = "green") +
         ggtitle(label = .z)
  )

library(cowplot)
plot_grid(plotlist = plot_comparisons2)
# also not good.



