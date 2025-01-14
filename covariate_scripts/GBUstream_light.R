## light model 
#
library(lubridate)
library(dplyr)
library(tidyr)
library(tidyverse)

dat <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/saturation/24_GBU_DO_flag_sat.rds")

light_dat <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/climate_dat/24stream_NLDAS_light.rds")%>%
  filter(site=="GBL")

light_dat <- light_dat%>%
  dplyr::select(Year , DOY, Hour, datetime, light) 


GBU_dat <- dat %>%
  left_join(light_dat, by=c("datetime"))


GBU_dat <- GBU_dat  %>%
  fill(light,.direction = "up")

GBU_dat <- GBU_dat  %>%
  fill(light,.direction = "down")

summary(GBU_dat)

saveRDS(GBU_dat, "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/light/24_GBU_DO_flag_sat_light.rds")

