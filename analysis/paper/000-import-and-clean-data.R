library(rgl)
library(here)
library(tidyverse)

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
landmarks <- read.tps(here::here("analysis/data/raw_data/original-landmark-data/BG_0_2.tps"))

#getting info for phase of each site
# SPstage <- read_csv(here("analysis/data/SPstage.csv"))
SPstage <- read_excel(here("analysis/data/SPstage_raw_materials.xls"))

#Compare multiple artifacts
#Open multiple tps files. This will need to be changed to read files inside the project
datadir <- here("analysis/data/raw_data/original-landmark-data/")
files <- dir(datadir, pattern = "*.tps$|*.TPS$")

#Getting site info
korean_archaeological_site_locations <- read_csv(here("analysis/data/sitelocation.csv"))

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

hist(MW)
# Maximum width of Tang: TW, distance between landmark 3-10
TW <- landmark_dist(landmarks, 3, 10)

# Stem width: SW, distance between landmark 5-8
SW <- landmark_dist(landmarks, 5, 8)

# get the landmark coordinates
landmarks_list <- vector("list", length = length(files))

files_original <- list.files(datadir, full.names = TRUE)

for (i in 1:length(files_original)) {
  filename <- files_original[i]

  # read in landmark files
  landmarks_list[[i]] <- readLines(here(filename))

  # drop curves it is unnecessary now cuz GP fixed raw data
  landmarks_list[[i]] <- landmarks_list[[i]][c(1:12, 76:77)]
  write_lines(landmarks_list[[i]],
              paste0("analysis/data/raw_data/",
                     "no-curves-landmark-data/",
                     paste0("no_curves_", basename(filename))))

}

# write code to exclude those tps files that lack filename and scale data
# because they wont work for the next lines that import

# this is just for testing and can be deleted.
landmarks_list <- vector("list", length = length(files))
no_curves_files <- list.files(here::here("analysis/data/raw_data/no-curves-landmark-data"),
                              full.names = TRUE)

# no_curves_files_no_bad_ones <- # remove bad files or fix them
#
# for (i in 1:length(no_curves_files)) {
#   filename <- no_curves_files[i]
#   print(filename)
#
#   # read in landmark files
#   landmarks_list[[i]] <- geomorph::readland.tps(filename)
# }

# this is the one we want, it produces the output in a format we can easily work with
x2 <- geomorph::readmulti.tps(no_curves_files)


# attach artefact ID's to landmarks
names(landmarks_list) <- files

# convert list to dataframe
landmarks_list_df <- bind_rows(landmarks_list, .id = "artefact")

library(Momocs)

Out(x2, fac = no_curves_files) %>%
  coo_align() %>%
  coo_center() %>%
  coo_scale() %>%
  pile(transp = 0.9)


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
  landmarks <- read.tps(here(filename))
  ML <- c(ML, landmark_dist(landmarks, 1, 7))
  BL <- c(BL, bl_dist(landmarks))
  TL <- c(TL, tl_dist(landmarks))
  SL <- c(SL, sl_dist(landmarks))
  MW <- c(MW, landmark_dist(landmarks, 2, 11))
  TW <- c(TW, landmark_dist(landmarks, 3, 10))
  SW <- c(SW, landmark_dist(landmarks, 5, 8))
}

#list of all attributes
df <- data.frame(ML, BL, TL, SL, MW, TW, SW, files, SPstage$Stage, SPstage$Raw_material)

#trim the name of each artifact.
df_sitename <- df %>%
  mutate (files = str_remove_all(files, "_.tps")) %>%
  separate (files, into = c("sitename", "files"), sep = "_") %>%
  select (-files)


library(dplyr)
df_full_sitename <-
  left_join(df_sitename, korean_archaeological_site_locations, by = c("sitename"="sitename"))



#list of all attributes without ID
DF <- df[,1:7]


save(df, df_sitename, df_full_sitename, DF, file = "data_main.RData")
# To load the data again
load("data_main.RData")
