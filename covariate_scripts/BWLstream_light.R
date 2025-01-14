## light model 
#

dat <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/saturation/24_BWL_DO_flag_sat.rds")

light_dat <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/climate_dat/24stream_NLDAS_light.rds")

BWL_dat <- dat %>%
  left_join(light_dat, by=c("datetime", "site"))


BWL_dat <- BWL_dat  %>%
  fill(light,.direction = "up")

BWL_dat <- BWL_dat  %>%
  fill(light,.direction = "down")

saveRDS(BWL_dat, "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/light/24_BWL_DO_flag_sat_light.rds")

