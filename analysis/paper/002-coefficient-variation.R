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
map(DF, mccv)

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

test_output <- cvCI(c(2,3,4), conf.level = 0.95, method = "sharma", na.rm = FALSE)

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
  select(-SPstage.Stage, -lat_dd, -long_dd, -sitename) %>%
  group_by(full_sitename) %>%
  add_tally() %>%
  filter( n > 1) %>%
  select(-n) %>%
  group_by(full_sitename) %>%
  summarise(cv_by_site =     sharma_cv(c(ML, BL , TL , SL , MW, TW, SW)),
            cv_low_sharma =  sharma_int_low(c(ML, BL , TL , SL , MW, TW, SW)),
            cv_high_sharma = sharma_int_high(c(ML, BL , TL , SL , MW, TW, SW)))





## CV for per each site (n=X) with the full site name : label

cv_by_full_site_df_label <-
  df_full_sitename %>%
  select(-SPstage.Stage) %>%
  group_by(full_sitename) %>%
  add_tally() %>%
  filter( n > 1) %>%
  mutate(label = as.factor(paste0(full_sitename,' (N = ', n, ")"))) %>%   select(-n) %>%
  nest(-full_sitename, -label) %>%
  mutate(cv_by_site = map(data, ~map_df(.x, cv))) %>%
  select (-cv_by_site)

## add label to main dataframe
cv_plot_site_label <- cv_by_site_df %>%
 left_join (cv_by_full_site_df_label) %>%
  select (- data,-full_sitename)

# facetted bar plot per each site, drop legend
site_bar_plot <- cv_plot_site_label %>%
  pivot_longer(cols = -label,
               names_to = "group") %>%
  ggplot(aes(x = group, y = value, fill = label)) +
  geom_col() +
  xlab("Site") +
  ylab("Coefficient of Variation on Attributes") +
  facet_wrap( ~ label) +
  theme(legend.position = "none")

site_bar_plot <- cv_plot_site_label %>%
  pivot_longer(cols = -label,
               names_to = "group") %>%
  ggplot(aes(x = group, y = value, fill = label)) +
  geom_point() +
  geom_errorbar(aes(ymin = min(cv_plot_site_label$group), ymax = max(cv_plot_site_label$group)), width = .2) +
  xlab("Site") +
  ylab("Coefficient of Variation on Attributes") +
  facet_wrap( ~ label) +
  theme(legend.position = "none")

# ggplot(cv_data_sharma_named, aes(attribute, all_cv_sharma)) +        # ggplot2 plot with confidence intervals
#   geom_point() +
#   geom_errorbar(aes(ymin = all_low_sharma, ymax = all_high_sharma), width = .2) +
#   ylab("CV") +
#   xlab("Variables") +
#   geom_text(aes(label = round(all_cv_sharma,1)), col="blue", hjust = -0.3, size = 3) +
#   theme_bw()
# facetted bar plot with (n=X)
#site_bar_plot + scale_fill_discrete(name = "Site name", labels = cv_by_full_site_df_label$label)

ggsave(here::here("analysis/figures/004-cv-sites.png"),
       width = 7,
       height = 4.5,
       units = "in")



############ temporal patterns ###############
# Make a table of CVs for all variable grouped by stage, according to stemmed point chronology
# Exclude the first phase that has only two artefacts
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

# bar plot for two phases
cv_plot_stage %>%
  pivot_longer(cols = -SPstage.Stage, names_to = "group") %>%
  ggplot(aes(x= SPstage.Stage, y = value, fill = as.factor(SPstage.Stage))) +
  geom_bar(stat="identity", position=position_dodge()) + labs(x="Attributes", fill= "Phase")+
  facet_wrap( ~ group) +
  theme_bw()

ggsave(here::here("analysis/figures/004-cv-phases.png"),
       width = 5,
       height = 4.5,
       units = "in")

# bar plot for each attribute
cv_plot_stage %>%
  ggplot(aes(x=SPstage.Stage, y=ML))+
  geom_line()



#### plot CVs by site with number of artefacts showing >10

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

# bar plot for four assemblages, all attributes
cv_ten_att <-
cv_by_site_df_label_more_than_10 %>%
  pivot_longer(cols = -c(site_phase, label),
               names_to = "lithic_attribute",
               values_to = "CV") %>%
  ggplot(aes(x = lithic_attribute,
             y = CV)) +
  geom_col() +
  facet_wrap(~  label) +
  ggtitle("CVs of point attributes by site (with sites >=7 artifacts)")

cv_by_site_df_label_more_than_10_ml_mw_tw <-
  cv_by_site_df_label_more_than_10 %>%
  pivot_longer(cols = -c(site_phase, label),
               names_to = "lithic_attribute",
               values_to = "CV") %>%
  filter(lithic_attribute %in% c("BL", "ML", "MW", "TW"))  %>%
  mutate(phase = parse_number(str_extract(site_phase, "_."))) %>%
  group_by(site_phase) %>%
  mutate(mean_cv = mean(CV))

# bar plot for four assemblages, four attributes
cv_ten_four_att <-
  ggplot(cv_by_site_df_label_more_than_10_ml_mw_tw,
       aes(x = lithic_attribute,
           y = CV)) +
  geom_col() +
  geom_hline(aes(yintercept = mean_cv),
             colour = "red") +
  facet_wrap(~  label) +
  ggtitle("CVs of point attributes by site (with sites >=7 artifacts)") +
  theme_minimal()

#### plot CVs by site with number of artefacts showing >=7

cv_by_site_df_label_more_than_10_ml_mw_bl <-
  cv_by_site_df_label_more_than_10 %>%
  pivot_longer(cols = -c(site_phase, label),
               names_to = "lithic_attribute",
               values_to = "CV") %>%
  filter(lithic_attribute %in% c("ML", "BL","MW", "TW"))  %>%
  mutate(phase = parse_number(str_extract(site_phase, "_."))) %>%
  group_by(site_phase) %>%
  mutate(mean_cv = mean(CV))

# bar plot for four assemblages, four attributes
ggplot(cv_by_site_df_label_more_than_10_ml_mw_bl,
       aes(x = lithic_attribute,
           y = CV)) +
  geom_col() +
  geom_hline(aes(yintercept = mean_cv),
             colour = "red") +
  facet_wrap(~  label) +
  ggtitle("CVs of point attributes by site (with sites >=7 artifacts)") +
  theme_minimal()

library(cowplot)
plot_grid(cv_ten_att, cv_ten_four_att)

