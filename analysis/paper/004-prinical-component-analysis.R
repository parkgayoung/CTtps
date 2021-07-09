library(tidyverse)
library(tidyr)
library(corrr)
library(ggbiplot)
library(ggfortify)

load("data_main.RData")

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
cv_all + cv_two + cv_three + plot_layout(nrow=1)



