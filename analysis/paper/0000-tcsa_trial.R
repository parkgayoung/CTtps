library(readxl)
library(here)

mydata <- read_excel(here("analysis/data/raw_data/TCSA_raw_data.xlsx"))

tcsa <- mydata %>%
  mutate(Tcsa = 0.5 *Width * Thickness)
