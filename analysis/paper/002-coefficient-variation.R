library(here)
library(tidyverse)
library(tidyr)

##################data treatment###################
#reading landmarks in tps, code from https://gist.github.com/mrdwab/2062329

read.tps = function(data) {
  # Reads the .tps file format produced by TPSDIG
  # (http://life.bio.sunysb.edu/morph/ into a single data frame
  # USAGE: R> read.tps("filename.tps")
  a = readLines(data) # so we can do some searching and indexing
  LM = grep("LM", a) # find the line numbers for LM
  ID.ind = grep("ID", a) # find the line numbers for ID
  # and the ID values, SCALE values, and image names
  ID = gsub("(ID=)(.*)", "\\2", grep("ID", a, value=T))
  SCALE = gsub("(SCALE=)(.*)", "\\2", grep("SCALE", a, value=T))
  images = basename(gsub("(IMAGE=)(.*)", "\\2", a[ID.ind - 1]))
  # FOR EACH LOOP
  skip = LM # set how many lines to skip
  # and how many rows to read
  nrows = as.numeric(gsub("(LM=)(.*)", "\\2", grep("LM", a, value=T)))
  l = length(LM) # number of loops we want

  landmarks = vector("list", l) # create an empty list

  for (i in 1:l) {
    landmarks[i] = list(data.frame(
      read.table(file = data, header = F, skip = LM[i],
                 nrows = nrows[i],
                 col.names = c("X", "Y"),
                 stringsAsFactors = TRUE),
      IMAGE = as.character(images[i]),
      ID = ID[i],
      SCALE = SCALE[i],
      stringsAsFactors = TRUE))
  }
  do.call(rbind, landmarks) # rbind the list items into a data.frame
}

#need to figure out how to read the file inside the project
landmarks <- read.tps(here("analysis/data/raw_data/BG_0_2.tps"))

#getting info for phase of each site
SPstage <- read_csv(here("analysis/data/SPstage.csv"))

#Compare multiple artifacts
#Open multiple tps files. This will need to be changed to read files inside the project
datadir <- here("/analysis/data/raw_data/")
files <- dir(datadir, pattern = "*.tps$|*.TPS$")


#Computing size attributes  of the points by calculating distance between landmarks

#Understanding each landmark
#ML: maximum length, tip to bottom, perpendicular to the length axis 1- 7
#BL: body length: tip to the closest wing, perpendicular to the length axis 1- 3 OR 1-10 depending on the artifact
#TL: Tang length, the closest wing from the tip to bottom, perpendicular to the length axis 9 - 7 depending on the artifact
#SL: max stem length, perpendicular to the length axis. Closer tang curve’s middle point from the tip to the most distant point of the basal end  9-7 depending on the artifact
#MW: mid width: dimension from margin to margin at the mid-point of the length 2 - 11 axis, perpendicular to the width axis.
#TW: tang width: dimension between each wing, perpendicular to the width  axis  3 - 10
#SW: stem width: width of the basal end of the point, 5mm above the end  5 - 8

#distance between landmarks(euclidean)

landmark_dist = function(landmarks, pt1, pt2) {
  sqrt(
    (landmarks$X[pt1] - landmarks$X[pt2])^2 +
      (landmarks$Y[pt1] - landmarks$Y[pt2])^2
  ) * as.numeric(levels(landmarks$SCALE))
}

# distance between two y axis
vertical_dist = function(landmarks, pt1, pt2) {
  d = abs(landmarks$Y[pt1] - landmarks$Y[pt2])
  d * as.numeric(levels(landmarks$SCALE))
}


#Maximum length of sp: ML, distance between landmark 1-7
ML <- landmark_dist(landmarks, 1, 7)

#Body length of sp is 1-3 or 1-10, depending on each artifact.
bl_dist = function(landmarks) {
  d1 = landmark_dist(landmarks, 1, 3)
  d2 = landmark_dist(landmarks, 1, 10)
  max(d1, d2)
}

BL <- bl_dist(landmarks)

# Tang length of sp: TL, distance between landmark 10-7 or 3-7, (y axis)
tl_dist <- function(landmarks) {
  t1 = vertical_dist(landmarks, 10, 7)
  t2 = vertical_dist(landmarks, 3, 7)
  max(t1, t2)
}

TL <- tl_dist

# Maximum length of stem: SL, distance between landmark 9-7 or 4-7 (y axis)
sl_dist <- function(landmarks) {
  s1 = vertical_dist(landmarks, 9, 7)
  s2 = vertical_dist(landmarks, 4, 7)
  max(s1, s2)
}

SL <- sl_dist

# Mid width: MW, distance between landmark 2-11
MW <- landmark_dist(landmarks, 2, 11)

# Maximum width of Tang: TW, distance between landmark 3-10
TW <- landmark_dist(landmarks, 3, 10)

# Stem width: SW, distance between landmark 5-8
SW <- landmark_dist(landmarks, 5, 8)



#Gather attributes from  all artifacts

ML = c()
BL = c()
TL = c()
SL = c()
MW = c()
TW = c()
SW = c()
for (i in 1:length(files)) {
  filename <- paste(datadir, files[i], sep = "")
  # print(filename)
  landmarks <- read.tps(filename)
  ML <- c(ML, landmark_dist(landmarks, 1, 7))
  BL <- c(BL, bl_dist(landmarks))
  TL <- c(TL, tl_dist(landmarks))
  SL <- c(SL, sl_dist(landmarks))
  MW <- c(MW, landmark_dist(landmarks, 2, 11))
  TW <- c(TW, landmark_dist(landmarks, 3, 10))
  SW <- c(SW, landmark_dist(landmarks, 5, 8))
}

#list of all attributes
df <- data.frame(ML, BL, TL, SL, MW, TW, SW, files, SPstage$Stage)

#trim the name of each artifact.
df_sitename <- df %>%
  mutate (files = str_remove_all(files, "_.tps")) %>%
  separate (files, into = c("sitename", "files"), sep = "_") %>%
  select (-files)

#list of all attributes without ID
DF <- df[,1:7]

########################actual analysis##########################

# Coefficient of variation (CV): the ratio of the standard deviation divided by mean

cv <- function(x, ... ) sd(x, ...) / mean(x, ...) *100

# we've got DF a data frame with cols as variables
# rows and specimens, and values as dimensions
# need to add a column of site ID
cvmap <- map(DF, cv)

#Same as map but better view for summary
map_df(DF, cv)

#Make a table of CVs for all variable grouped by site
cv_by_site_df  <-
  df_sitename %>%
  select(-SPstage.Stage) %>%
  group_by(sitename) %>%
  add_tally() %>%
  filter( n > 1) %>%
  select(-n) %>%
  nest(-sitename) %>%
  mutate(cv_by_site = map(data, ~map_df(.x, cv)))

cv_by_site_df_unnest <-
  cv_by_site_df %>%
  unnest(cv_by_site)


cv_plot_site <- cv_by_site_df_unnest %>%
  select(-data)

cv_plot_site %>%
  pivot_longer(cols = -sitename, names_to = "group") %>%
  ggplot(aes(x = group, y = value, fill = sitename)) +
  geom_bar(stat="identity", position=position_dodge())

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


