library(tidyr)
library(purrr)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
load("data_main.RData")

# Coefficient of variation (CV): the ratio of the standard deviation divided by mean

cv <- function(x, ... ) sd(x, ...) / mean(x, ...) *100

#try corrected cv with different methods
library(cvcqv)
ccv <- function(x, ... ) cv_versatile(x, method = "kelly", ...)
mccv <- function(x, ... ) sd(x, ...) / mean(x, ...) *100 * (1 + 1/(4*length(x)))

# we've got DF a data frame with cols as variables
# rows and specimens, and values as dimensions
# need to add a column of site ID

DFCV<- map_dbl(DF, cv)

#Same as map but better view for summary
cv_all <- map_df(DF, cv)
#arrange the data properly
cv_data_all <- gather(cv_all)

hist(DF$ML)

hist(DFCV)

# to install MKmisc
## Install package BiocManager
# install.packages("BiocManager")
## Use BiocManager to install limma
# BiocManager::install("limma")
# install.packages("MKmisc")

library(MKmisc)
#CV value
sharma_cv <- function(x,...) cvCI(x, conf.level = 0.95, method = "sharma", na.rm = FALSE)$estimate*100

#CV sharma interval values
sharma_int_low <- function(x,...) cvCI(x, conf.level = 0.95, method = "sharma", na.rm = FALSE)$conf.int[1]*100
sharma_int_high <- function(x,...) cvCI(x, conf.level = 0.95, method = "sharma", na.rm = FALSE)$conf.int[2]*100

all_cv_sharma <- map_dbl(DF, sharma_cv)
all_low_sharma <- map_dbl(DF, sharma_int_low)
all_high_sharma <- map_dbl(DF, sharma_int_high)

#create data frame
cv_data_sharma <- data.frame(all_cv_sharma, all_low_sharma, all_high_sharma)

#add column name
library(tibble)
cv_data_sharma_named <- rownames_to_column(cv_data_sharma, var = "attribute")


library(ggplot2)
ggplot(cv_data_sharma_named, aes(attribute, all_cv_sharma)) +        # ggplot2 plot with confidence intervals
  geom_point() +
  geom_errorbar(aes(ymin = all_low_sharma, ymax = all_high_sharma), width = .2) +
  ylab("CV") +
  xlab("Variables") +
  geom_text(aes(label = round(all_cv_sharma,1)), col="blue", hjust = -0.3, size = 3) +
  theme_bw()

ggsave(here::here("analysis/figures/003-cv-sharma.png"),
       width = 4.8,
       height = 4.5,
       units = "in")

#####################
#Make a table of CVs for all variable grouped by site
cv_by_site_df  <-
  df_full_sitename %>%
  select(-SPstage.Stage, -lat_dd, -long_dd, -sitename, -SPstage.Raw_material) %>%
  group_by(full_sitename) %>%
  add_tally() %>%
  filter( n > 1) %>%
  select(-n) %>%
  pivot_longer(-full_sitename,
               names_to = "variable",
               values_to = "value") %>%
  group_by(full_sitename, variable) %>%
  summarise(cv_by_site =   sharma_cv(value),
            cv_low_sharma =  sharma_int_low(value),
            cv_high_sharma = sharma_int_high(value))


## CV for per each site (n=X) with the full site name : label

cv_by_full_site_df_label <-
  df_full_sitename %>%
  select(-SPstage.Stage) %>%
  group_by(full_sitename) %>%
  add_tally() %>%
  filter( n > 1) %>%
  mutate(label = as.factor(paste0(full_sitename,' (N = ', n, ")"))) %>%
  select(-n) %>%
  nest(-full_sitename, -label) %>%
  mutate(cv_by_site = map(data, ~map_df(.x, cv))) %>%
  select (-cv_by_site)


## add label to main dataframe

cv_plot_site_label <- cv_by_site_df %>%
  left_join (cv_by_full_site_df_label) %>%
  select (- data,-full_sitename)


## facet plot

cv_plot_site_label %>%
  ggplot(aes(x = variable, y = cv_by_site, colour = factor(variable))) +
  geom_point() +
  geom_errorbar(aes(ymin = cv_low_sharma,
                    ymax = cv_high_sharma), width = .2) +
  xlab("Site") +
  ylab("Coefficient of Variation on Attributes") +
  facet_wrap( ~ label, scales = "free_y") +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0,130)) +
  theme_bw() +
  labs(color = "Attributes")


# checking SC site values
# Sac <- df_full_sitename %>%
#   filter(full_sitename=="Sachang")
#
# cv(Sac$MW)
# sharma_int_low(Sac$MW)
# sharma_int_high(Sac$MW)
#
# cv(Sac$SL)
# sharma_int_low(Sac$SL)
# sharma_int_high(Sac$SL)


ggsave(here::here("analysis/figures/004-cv-sites.png"),
       width = 8.5,
       height = 4.5,
       units = "in")



############ temporal patterns ###############
# Make a table of CVs for all variable grouped by stage, according to stemmed point chronology
# Exclude the first phase that has only two artefacts
cv_by_stage_df  <-
  df_sitename %>%
  select(-sitename, -SPstage.Raw_material) %>%
  group_by(SPstage.Stage) %>%
  add_tally() %>%
  filter(SPstage.Stage > 1) %>%
  select(-n) %>%
  pivot_longer(-SPstage.Stage,
               names_to = "variable",
               values_to = "value") %>%
  group_by(SPstage.Stage, variable) %>%
  summarise(cv_by_site =   sharma_cv(value),
            cv_low_sharma =  sharma_int_low(value),
            cv_high_sharma = sharma_int_high(value))

####### significant test between the two phases ######
library(cvequality)

# BM testing doing all the CV testing at once

# created nested tibble, nest by variable
df_sitename_nested <-
df_sitename %>%
  select(-sitename, -SPstage.Raw_material) %>%
  group_by(SPstage.Stage) %>%
  add_tally() %>%
  filter(SPstage.Stage > 1) %>%
  select(-n) %>%
  pivot_longer(-SPstage.Stage,
               names_to = "variable",
               values_to = "value") %>%
  ungroup() %>%
  group_by(variable) %>%
  nest()

# apply mslr test to each variable (row with a tibble)
#Krishnamoorthy and Lee’s (2014) modified signed-likelihood ratio test
df_sitename_nested %>%
  mutate(mslr_result = map_df(data, ~mslr_test(nr = 1e4,
                                            .$value,
                                            .$SPstage.Stage))) %>%
  unnest(mslr_result)


# facet plot for two phases
cv_by_stage_df %>%
  ggplot(aes(x = factor(SPstage.Stage), y = cv_by_site, colour = factor(SPstage.Stage))) +
  geom_point() +
  geom_errorbar(aes(ymin = cv_low_sharma,
                    ymax = cv_high_sharma), width = .05) +
  xlab("Attributes") +
  ylab("Coefficient of Variation on Attributes") +
  facet_wrap( ~ variable, scales = "free_y") +
  coord_cartesian(ylim = c(15,45)) +
  theme_bw() +
  labs(color = "Phase")

# bar plot for two phases
# cv_plot_stage %>%
#   select(-SPstage.Raw_material) %>%
#   pivot_longer(cols = -SPstage.Stage, names_to = "group") %>%
#   ggplot(aes(x= SPstage.Stage, y = value, fill = as.factor(SPstage.Stage))) +
#   geom_bar(stat="identity", position=position_dodge()) + labs(x="Attributes", fill= "Phase")+
#   facet_wrap( ~ group) +
#   theme_bw()

ggsave(here::here("analysis/figures/005-cv-phases.png"),
       width = 4.5,
       height = 4.5,
       units = "in")



#### plot CVs by site with number of artefacts showing >5

cv_by_site_df_label_more_than_5 <-
  df_sitename %>%
  unite(col = "site_phase",
        c(sitename,
          SPstage.Stage),
        sep = "_") %>%
  select(-SPstage.Raw_material) %>%
  group_by(site_phase) %>%
  add_tally() %>%
  filter( n > 5) %>%
  mutate(label = as.factor(paste0(site_phase,' (N = ', n, ")"))) %>%
  ungroup() %>%
  select(-n, -site_phase) %>%
  pivot_longer(-label,
               names_to = "variable",
               values_to = "value") %>%
  group_by(label, variable) %>%
  summarise(cv_by_site =   sharma_cv(value),
            cv_low_sharma =  sharma_int_low(value),
            cv_high_sharma = sharma_int_high(value)) %>%
  filter(variable %in% c("ML", "BL","MW", "TW"))

cv_by_site_df_label_more_than_5 %>%
  ggplot(aes(x = factor(variable), y = cv_by_site, colour = factor(variable))) +
  geom_point() +
  geom_errorbar(aes(ymin = cv_low_sharma,
                    ymax = cv_high_sharma), width = .2) +
  xlab("Asseblages") +
  ylab("Coefficient of Variation on Attributes") +
  facet_wrap( ~ label, scales = "free_y") +
  coord_cartesian(ylim = c(15,60)) +
  theme_bw() +
  labs(color = "Attributes")+
  ggtitle("CVs of point attributes by site (with sites >=7 artifacts)")

# #### plot CVs by site with number of artefacts showing >=7
#
# cv_by_site_df_label_more_than_5_ml_mw_bl <-
#   cv_by_site_df_label_more_than_5 %>%
#   pivot_longer(cols = -c(site_phase, label),
#                names_to = "lithic_attribute",
#                values_to = "CV") %>%
#   filter(lithic_attribute %in% c("ML", "BL","MW", "TW"))  %>%
#   mutate(phase = parse_number(str_extract(site_phase, "_."))) %>%
#   group_by(site_phase) %>%
#   mutate(mean_cv = mean(CV))

# bar plot for four assemblages, four attributes
# ggplot(cv_by_site_df_label_more_than_5_ml_mw_bl,
#        aes(x = lithic_attribute,
#            y = CV)) +
#   geom_col() +
#   geom_hline(aes(yintercept = mean_cv),
#              colour = "red") +
#   facet_wrap(~  label) +
#   ggtitle("CVs of point attributes by site (with sites >=7 artifacts)") +
#   theme_minimal()

ggsave(here::here("analysis/figures/006-cv-four-assemblage.png"),
       width = 5,
       height = 4,
       units = "in")


# CV by raw materials
### Make a table of CVs for all variable grouped by raw materials
cv_by_raw_materials <-
  df_full_sitename %>%
  select(-SPstage.Stage, -lat_dd, -long_dd, -sitename, -full_sitename) %>%
  group_by(SPstage.Raw_material) %>%
  add_tally() %>%
  filter( n > 1) %>%
  select(-n) %>%
  pivot_longer(-SPstage.Raw_material,
               names_to = "variable",
               values_to = "value") %>%
  group_by(SPstage.Raw_material, variable) %>%
  summarise(cv_by_materials =   sharma_cv(value),
            cv_low_sharma =  sharma_int_low(value),
            cv_high_sharma = sharma_int_high(value))

######CV for per raw materials (n=X) with the full site name : label
cv_by_raw_materials_df_label <-
  df_sitename %>%
  group_by(SPstage.Raw_material) %>%
  add_tally() %>%
  filter( n > 1) %>%
  mutate(label = as.factor(paste0(SPstage.Raw_material,' (N = ', n, ")"))) %>%
  select(-n) %>%
  nest(-SPstage.Raw_material, -label) %>%
  mutate(cv_by_site = map(data, ~map_df(.x, cv))) %>%
  select (-cv_by_site)

######add label to main dataframe
cv_plot_raw_materials_label <- cv_by_raw_materials %>%
  left_join (cv_by_raw_materials_df_label) %>%
  select (- data,-SPstage.Raw_material)

###### facet plot
cv_plot_raw_materials_label %>%
  ggplot(aes(x = variable, y = cv_by_materials, colour = factor(variable))) +
  geom_point() +
  geom_errorbar(aes(ymin = cv_low_sharma,
                    ymax = cv_high_sharma), width = .2) +
  xlab("Law Materials") +
  ylab("Coefficient of Variation on Attributes") +
  facet_wrap( ~ label, scales = "free_y") +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0,130)) +
  theme_bw() +
  labs(color = "Attributes")


ggsave(here::here("analysis/figures/00000-cv-raw-materials.png"),
       width = 8.5,
       height = 4.5,
       units = "in")

