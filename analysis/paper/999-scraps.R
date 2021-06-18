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
cv <- function(x, ... ) sd(x, ...) / mean(x, ...) *100

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
