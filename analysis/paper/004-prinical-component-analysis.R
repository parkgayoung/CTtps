library(here)
library(tidyverse)
library(tidyr)
library(corrr)
library(ggbiplot)
library(ggfortify)

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

### PCA analysis using ggbiplot, ggfortify pkgs


PCA <- prcomp(DF, center = TRUE, scale = TRUE)

# remove stage 1 artefact
xx_pca <-
  df_sitename %>%
  filter(SPstage.Stage != 1) %>%
  select(-sitename, -SPstage.Stage) %>%
  prcomp(., center = TRUE, scale = TRUE)

xx_pca1 <-
  data.frame(xx_pca$x,
             stage = as.factor( df_sitename$SPstage.Stage[df_sitename$SPstage.Stage !=1]))


## CV with site (n=X)
label_phase <-
  df_sitename %>%
  group_by(SPstage.Stage) %>%
  add_tally() %>%
  mutate(label = as.factor(paste0(SPstage.Stage,'\nN = ',n)))


count_phase <- length(which(label_phase[9]>0))
count_phase2 <- length(which(label_phase[9] == 2))
count_phase3 <- length(which(label_phase[9] == 3))


###

#plot for both stage
gg_stage <-
  ggplot(xx_pca1 ,
         aes(PC1,
             PC2,
             colour = stage)) +
  geom_point() +
  stat_ellipse()

# phase 2

xx_pca2 <-
  df_sitename %>%
  filter(SPstage.Stage == 2) %>%
  select(-sitename, -SPstage.Stage) %>%
  prcomp(., center = TRUE, scale = TRUE)



# phase 3

xx_pca3 <-
  df_sitename %>%
  filter(SPstage.Stage == 3) %>%
  select(-sitename, -SPstage.Stage) %>%
  prcomp(., center = TRUE, scale = TRUE)

##### using factoextra pkg
library(factoextra) #http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/

# Color variables by the continuous variable
cv_all <- fviz_pca_var(xx_pca, col.var = "contrib",
                       gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")) + ggtitle(paste("All", "(n =", count_phase,")"))

# phase 2
cv_two <- fviz_pca_var(xx_pca2, col.var = "contrib",
                       gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")) + ggtitle(paste("Phase 2", "(n =", count_phase2,")"))
# phase 3
cv_three <- fviz_pca_var(xx_pca3, col.var = "contrib",
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")) + ggtitle(paste("Phase 3", "(n =", count_phase3,")"))



library(patchwork)
cv_all + cv_two + cv_three + plot_layout(ncol=1)


######### explore SYG6 site which contains more than 1 stage
library(stringr)
SYG6_df <-
  df %>%
  filter(str_detect(files, "SYG6"))


SYG6_df_sitename <- SYG6_df %>%
  mutate (files = str_remove_all(files, "_.tps")) %>%
  separate (files, into = c("sitename", "files"), sep = "_") %>%
  select (-files)

#list of all attributes without ID
SYG6_DF <- SYG6_df[,1:7]

# make n=element number
count_SYG6 <- length(which(SYG6_df$SPstage.Stage == 2))


# using corrr

corrr_SYG6_DF <-correlate(SYG6_DF) %>%
  shave()

cSYG6_DF <- rplot(corrr_SYG6_DF) +
  ggtitle(paste("Correlation SYG6: 2nd Phase", "(n =", count_SYG6,")"))


SYG6_cor <- cor(SYG6_DF)

# CV
map_df(SYG6_DF, cv)

SYG6_PCA <- prcomp(SYG6_cor, center = TRUE, scale = TRUE)

# remove stage 1 artefact
SYG6_xx_pca <-
  SYG6_df_sitename %>%
  select(-sitename, -SPstage.Stage) %>%
  prcomp(., center = TRUE, scale = TRUE)

SYG6_xx_pca1 <-
  data.frame(SYG6_xx_pca$x,
             stage = as.factor( SYG6_df_sitename$SPstage.Stage))

#plot for both stage
SYG6_gg_stage <-
  ggplot(SYG6_xx_pca1 ,
         aes(PC1,
             PC2,
             colour = stage)) +
  geom_point() +
  stat_ellipse()

cv_SYG6 <-fviz_pca_var(SYG6_xx_pca, col.var = "contrib",
                       gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")) +
  ggtitle(paste("SYG6: 2nd Phase", "(n =", count_SYG6,")"))


SYG6_PC1_v <- fviz_contrib(SYG6_xx_pca, choice = "var", axes = 1, top = 7)
SYG6_PC2_v <- fviz_contrib(SYG6_xx_pca, choice = "var", axes = 2, top = 7)

library(gridExtra)
grid.arrange(SYG6_PC1_v, SYG6_PC2_v, nrow =1)


########SYG1_2

SYG1_2df <-
  df %>%
  filter(str_detect(files, "SYG1_2"))


SYG1_2df_sitename <- SYG1_2df %>%
  mutate (files = str_remove_all(files, "_\\d+.tps")) %>%
  dplyr::rename(sitename = files)


#list of all attributes without ID
SYG1_2DF <- SYG1_2df[,1:7]

# make n=element number
count_SYG1_2 <- length(which(SYG1_2df$SPstage.Stage == 3))



# using corrr

corrr_SYG1_2DF <-correlate(SYG1_2DF) %>%
  shave()

cSYG1_2DF <- rplot(corrr_SYG1_2DF) +
  ggtitle(paste("Correlation SYG1_2: 3rd Phase", "(n =", count_SYG1_2,")"))




SYG1_2cor <- cor(SYG1_2DF)

# CV
map_df(SYG1_2DF, cv)

SYG1_2PCA <- prcomp(SYG1_2DF, center = TRUE, scale = TRUE)

# remove stage 1 artefact
SYG1_2xx_pca <-
  SYG1_2df_sitename %>%
  select(-sitename, -SPstage.Stage) %>%
  prcomp(., center = TRUE, scale = TRUE)

SYG1_2xx_pca1 <-
  data.frame(SYG1_2xx_pca$x,
             stage = as.factor(SYG1_2df_sitename$SPstage.Stage))

#plot for both stage
SYG1_2gg_stage <-
  ggplot(SYG1_2xx_pca1 ,
         aes(PC1,
             PC2,
             colour = stage)) +
  geom_point() +
  stat_ellipse()

cv_SYG1_2 <- fviz_pca_var(SYG1_2xx_pca, col.var = "contrib",
                          gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")) +
  ggtitle(paste("SYG1_2: 3rd Phase", "(n =", count_SYG1_2,")"))


SYG1_2PC1_v <- fviz_contrib(SYG1_2xx_pca, choice = "var", axes = 1, top = 7)
SYG1_2PC2_v <- fviz_contrib(SYG1_2xx_pca, choice = "var", axes = 2, top = 7)

grid.arrange(SYG1_2PC1_v, SYG1_2PC2_v, nrow =1)

#########SYG1_1
SYG1_1df <-
  df %>%
  filter(str_detect(files, "SYG1_1"))


SYG1_1df_sitename <- SYG1_1df %>%
  mutate (files = str_remove_all(files, "_\\d+.tps")) %>%
  dplyr::rename(sitename = files)


#list of all attributes without ID
SYG1_1DF <- SYG1_1df[,1:7]

# make n=element number
count_SYG1_1 <- length(which(SYG1_1df$SPstage.Stage == 2))



# using corrr

corrr_SYG1_1DF <-correlate(SYG1_1DF) %>%
  shave()

cSYG1_1DF <- rplot(corrr_SYG1_1DF) +
  ggtitle(paste("Correlation SYG1_1: 2nd Phase", "(n =", count_SYG1_1,")"))



SYG1_1cor <- cor(SYG1_1DF)

# CV
map_df(SYG1_1DF, cv)

SYG1_1PCA <- prcomp(SYG1_1DF, center = TRUE, scale = TRUE)

# remove stage 1 artefact
SYG1_1xx_pca <-
  SYG1_1df_sitename %>%
  select(-sitename, -SPstage.Stage) %>%
  prcomp(., center = TRUE, scale = TRUE)


SYG1_1xx_pca1 <-
  data.frame(SYG1_1xx_pca$x,
             stage = as.factor(SYG1_1df_sitename$SPstage.Stage))

#plot for both stage
SYG1_1gg_stage <-
  ggplot(SYG1_1xx_pca1 ,
         aes(PC1,
             PC2,
             colour = stage)) +
  geom_point() +
  stat_ellipse()

cv_SYG1_1 <- fviz_pca_var(SYG1_1xx_pca, col.var = "contrib",
                          gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")) +
  ggtitle(paste("SYG1_1: 2nd Phase", "(n =", count_SYG1_1,")"))


SYG1_1PC1_v <- fviz_contrib(SYG1_1xx_pca, choice = "var", axes = 1, top = 7)
SYG1_1PC2_v <- fviz_contrib(SYG1_1xx_pca, choice = "var", axes = 2, top = 7)

grid.arrange(SYG1_1PC1_v, SYG1_1PC2_v, nrow =1)

#### three sites at once, using corrr
cSYG6_DF + cSYG1_2DF + cSYG1_1DF + plot_layout(ncol=1)

######## three sites at once, PCA
cv_SYG6 + cv_SYG1_1 + cv_SYG1_2 +  plot_layout(ncol=1)
