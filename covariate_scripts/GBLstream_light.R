## light model 
#

dat <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/saturation/24_GBL_DO_flag_sat.rds")

light_dat <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/climate_dat/24stream_NLDAS_light.rds")

GBL_dat <- dat %>%
  left_join(light_dat, by=c("datetime", "site"))

GBL_dat <- GBL_dat  %>%
  fill(light,.direction = "up")

GBL_dat <- GBL_dat  %>%
  fill(light,.direction = "down")

saveRDS(GBL_dat, "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/light/24_GBL_DO_flag_sat_light.rds")

