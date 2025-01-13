## light model 
#

dat <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/saturation/24_BWU_DO_flag_sat.rds")

light_dat <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/climate_dat/24stream_NLDAS_light.rds")%>%
  filter(site=="BWL")

light_dat <- light_dat%>%
  dplyr::select(Year , DOY, Hour, datetime, light) 


BWU_dat <- dat %>%
  left_join(light_dat, by=c("datetime"))



BWU_dat <- BWU_dat  %>%
  fill(light,.direction = "up")

BWU_dat <- BWU_dat  %>%
  fill(light,.direction = "down")

summary(BWU_dat)

saveRDS(BWU_dat, "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/light/24_BWU_DO_flag_sat_light.rds")

