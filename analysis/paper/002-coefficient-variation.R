library(tidyr)
library(purrr)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
load("data_main.RData")

# Coefficient of variation (CV): the ratio of the standard deviation divided by mean

cv <- function(x, ... ) sd(x, ...) / mean(x, ...) *100

# we've got DF a data frame with cols as variables
# rows and specimens, and values as dimensions
# need to add a column of site ID

map(DF, cv)

#Same as map but better view for summary
cv_all <- map_df(DF, cv)

#arrage the data properly
cv_data_all <- gather(cv_all)

library(dplyr)
# ummary_table(cv_all)

#Make a table of CVs for all variable grouped by site
cv_by_site_df  <-
  df_full_sitename %>%
  select(-SPstage.Stage, -lat_dd, -long_dd, -sitename) %>%
  group_by(full_sitename) %>%
  add_tally() %>%
  filter( n > 1) %>%
  select(-n) %>%
  nest(-full_sitename) %>%
  mutate(cv_by_site = map(data, ~map_df(.x, cv)))

cv_by_site_df_unnest <-
  cv_by_site_df %>%
  unnest(cv_by_site)

cv_plot_site <- cv_by_site_df_unnest %>%
  select(-data)

# grouped bar plot
cv_plot_site %>%
  pivot_longer(cols = -full_sitename,
               names_to = "group") %>%
  ggplot(aes(x = group, y = value, fill = full_sitename)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(fill = "Site")

ggsave(here::here("analysis/figures/004-cv-sites.png"),
       width = 7,
       height = 4.5,
       units = "in")

# facetted bar plot
cv_plot_site %>%
  pivot_longer(cols = -full_sitename,
               names_to = "group") %>%
  ggplot(aes(x = group, y = value, fill = full_sitename)) +
  geom_col() +
  facet_wrap( ~ full_sitename)
  labs(fill = "Site")

#Make a table of CVs for all variable grouped by stage, according to stemmed point chronology

cv_by_stage_df  <-
  df_sitename %>%
  select(-sitename) %>%
  group_by(SPstage.Stage) %>%
  add_tally() %>%
  filter(SPstage.Stage > 1) %>%
  nest(-SPstage.Stage) %>%
  mutate(cv_by_stage = map(data, ~map_df(.x, cv)))

cv_by_stage_df_unnest <-
  cv_by_stage_df %>%
  unnest(cv_by_stage)

cv_plot_stage <- cv_by_stage_df_unnest %>%
  select(-data, -n)

cv_plot_stage %>%
  pivot_longer(cols = -SPstage.Stage, names_to = "group") %>%
  ggplot(aes(x = group, y = value, fill = as.factor(SPstage.Stage))) +
  geom_bar(stat="identity", position=position_dodge()) + labs(x="Attributes", fill= "Phase")


## CV per each site (n=X)

cv_by_site_df_label <-
  df_sitename %>%
  select(-SPstage.Stage) %>%
  group_by(sitename) %>%
  add_tally() %>%
  filter( n > 1) %>%
  mutate(label = as.factor(paste0(sitename,' (N = ', n, ")"))) %>%   select(-n) %>%
  nest(-sitename, -label) %>%
  mutate(cv_by_site = map(data, ~map_df(.x, cv)))

cv_by_site_df_unnest_label <-
  cv_by_site_df_label %>%
  unnest(cv_by_site)

cv_plot_site_label <- cv_by_site_df_unnest_label %>%
  select(-data)


## Reference: site and number of artifacts
cv_by_site_df_label_1 <-
  df_sitename %>%
  select(-SPstage.Stage) %>%
  group_by(sitename) %>%
  add_tally() %>%
  mutate(label = as.factor(paste0(sitename,'\nN = ',n))) %>%   select(-n) %>%
  nest(-sitename, -label) %>%
  mutate(cv_by_site = map(data, ~map_df(.x, cv)))

cv_by_site_df_unnest_label_1 <-
  cv_by_site_df_label_1 %>%
  unnest(cv_by_site)


p<- cv_plot_site_label %>%
  pivot_longer(cols = -c(sitename,label), names_to = "group") %>%
  ggplot(aes(label, value, fill = as.factor(group))) +
  geom_bar(stat="identity", position=position_dodge())

p + labs(x= "Site", fill = "Attributes")


# plot CVs by site with number of artefacts showing >10
# phase 2 and phase 3 sites only
cv_by_site_df_label_more_than_10 <-
  df_sitename %>%
  unite(col = "site_phase",
        c(sitename,
          SPstage.Stage),
        sep = "_") %>%
  group_by(site_phase) %>%
  add_tally() %>%
  filter( n > 5) %>%
  mutate(label = as.factor(paste0(site_phase,' (N = ', n, ")"))) %>%   select(-n) %>%
  nest(-site_phase, -label) %>%
  mutate(cv_by_site = map(data, ~map_df(.x, cv))) %>%
  unnest(cv_by_site) %>%
  select(-data)

cv_ten_att <-
cv_by_site_df_label_more_than_10 %>%
  pivot_longer(cols = -c(site_phase, label),
               names_to = "lithic_attribute",
               values_to = "CV") %>%
  ggplot(aes(x = lithic_attribute,
             y = CV)) +
  geom_col() +
  facet_wrap(~  label) +
  ggtitle("CVs of point attributes by site (with sites >10 artefacts)")

cv_by_site_df_label_more_than_10_ml_mw_tw <-
  cv_by_site_df_label_more_than_10 %>%
  pivot_longer(cols = -c(site_phase, label),
               names_to = "lithic_attribute",
               values_to = "CV") %>%
  filter(lithic_attribute %in% c("BL", "ML", "MW", "TW"))  %>%
  mutate(phase = parse_number(str_extract(site_phase, "_."))) %>%
  group_by(site_phase) %>%
  mutate(mean_cv = mean(CV))

cv_ten_four_att <-
  ggplot(cv_by_site_df_label_more_than_10_ml_mw_tw,
       aes(x = lithic_attribute,
           y = CV)) +
  geom_col() +
  geom_hline(aes(yintercept = mean_cv),
             colour = "red") +
  facet_wrap(~  label) +
  ggtitle("CVs of point attributes by site (with sites >10 artefacts)") +
  theme_minimal()

# plot CVs by site with number of artefacts showing >=7
cv_by_site_df_label_more_than_10_ml_mw_bl <-
  cv_by_site_df_label_more_than_10 %>%
  pivot_longer(cols = -c(site_phase, label),
               names_to = "lithic_attribute",
               values_to = "CV") %>%
  filter(lithic_attribute %in% c("ML", "BL","MW"))  %>%
  mutate(phase = parse_number(str_extract(site_phase, "_."))) %>%
  group_by(site_phase) %>%
  mutate(mean_cv = mean(CV))

ggplot(cv_by_site_df_label_more_than_10_ml_mw_bl,
       aes(x = lithic_attribute,
           y = CV)) +
  geom_col() +
  geom_hline(aes(yintercept = mean_cv),
             colour = "red") +
  facet_wrap(~  label) +
  ggtitle("CVs of point attributes by site (with sites >=7 artefacts)") +
  theme_minimal()

library(cowplot)
plot_grid(cv_ten_att, cv_ten_four_att)
