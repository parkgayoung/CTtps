
library(tidyr)
library(corrr)
load("data_main.RData")

#using corrr package
corrr_all <- correlate(DF) %>%
  shave()

call<-rplot(corrr_all) +
  ggtitle("Correlation for all phase")

direction_all <- correlate(DF) %>%
  network_plot()

#phase 2

DF2<- df_sitename %>%
  filter(SPstage.Stage == 2) %>%
  select(-sitename, -SPstage.Stage)

corrr_2 <-correlate(DF2) %>%
  shave()

ctwo <- rplot(corrr_2) +
  ggtitle("Correlation for 2nd phase")

direction_2 <- correlate(DF2) %>%
  network_plot()

#phase 3

DF3<-df_sitename %>%
  filter(SPstage.Stage == 3) %>%
  select(-sitename, -SPstage.Stage)

corrr_3 <-correlate(DF3) %>%
  shave()

cthree <- rplot(corrr_3) +
  ggtitle("Correlation for 3rd phase")

direction_3 <-correlate(DF3) %>%
  network_plot()

### three phases at once
library(patchwork)
call + ctwo + cthree + plot_layout(ncol=1)


