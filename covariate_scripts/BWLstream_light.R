## light model 
#25_GBU_DO_flag_sat.rds

## ===============================
### Read in DO data
dat <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/saturation/25_BWL_DO_flag_sat.rds")

names(dat)

## ===============================
## Check initial timestamps: 
# Filter data for a specific date range
data_plot <- dat %>%
  mutate(date = as.Date(rounded_datetime_15min),
         week = week(rounded_datetime_15min)) %>%
  filter(date > as.Date("2021-08-13") & date < as.Date("2021-08-27"))

# Create a dataframe with noon times
noon_lines <- data_plot %>%
  group_by(date) %>%
  summarize(noon_time = as.POSIXct(paste(date, "12:00:00"), tz = "America/Los_Angeles")) %>%
  ungroup()


# Plot with noon vertical lines
data_plot %>%
  ggplot(aes(y = wtr_USGS, x = rounded_datetime_15min)) +
  geom_line(aes(y = wtr_USGS, x = rounded_datetime_15min), col = "red") + 
  geom_line(aes(y = (dischargeCMS*10), x = rounded_datetime_15min), col = "#46b8b8",alpha = .7, size = 0.75) +
  geom_line(aes(y = wtr, x = rounded_datetime_15min), col = "goldenrod",alpha = .7, size = 0.75) +
  geom_vline(data = noon_lines, aes(xintercept = as.numeric(noon_time)), 
             linetype = "dashed", color = "grey50") +
  theme_bw() +
  facet_wrap(~week, scales = "free", ncol = 1) +
  scale_x_datetime(date_breaks = "12 hours", date_labels = "%H:%M") 


## ===============================
## Read in light data
light_dat <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/climate_dat/24stream_NLDAS_light.rds")%>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC")) %>%
  with_tz(tz = "America/Los_Angeles") %>%
  filter(site=="BWL") %>%
  mutate(par= c(light* 2.114))

BWL_dat <- dat %>%
  left_join(light_dat, by=c("rounded_datetime_15min"="datetime"))

BWL_dat2 <- BWL_dat %>%
  mutate(solar.time = calc_solar_time(rounded_datetime_15min, -120.164335),
         calclight = calc_light(solar.time, 39.109784, -120.164335),
         light_in = ifelse(is.na(par), calclight, par))
         
summary(BWL_dat2)

## ===============================
### double check USGS timestamp:
# Filter data for a specific date range
data_plot <- BWL_dat2 %>%
  mutate(date = as.Date(rounded_datetime_15min),
         week = week(rounded_datetime_15min)) %>%
  filter(date > as.Date("2021-08-13") & date < as.Date("2021-08-27"))

# Create a dataframe with noon times
noon_lines <- data_plot %>%
  group_by(date) %>%
  summarize(noon_time = as.POSIXct(paste(date, "12:00:00"), tz = "America/Los_Angeles")) %>%
  ungroup()

# Plot with noon vertical lines
data_plot %>%
  ggplot(aes(y = wtr_USGS, x = solar.time)) +
  geom_line(aes(y = wtr_USGS, x = solar.time), col = "red") + 
  geom_line(aes(y = (dischargeCMS*10), x = solar.time), col = "#46b8b8",alpha = .7, size = 0.75) +
  geom_line(aes(y = wtr, x = solar.time), col = "goldenrod",alpha = .7, size = 0.75) +
  geom_line(aes(y = light_in/100, x = solar.time), col = "yellow",alpha = .7, size = 0.75) +
  # geom_vline(data = noon_lines, aes(xintercept = as.numeric(noon_time)), 
  #            linetype = "dashed", color = "grey50") +
  theme_bw() +
  facet_wrap(~week, scales = "free", ncol = 1) +
  scale_x_datetime(date_breaks = "12 hours", date_labels = "%H:%M") 


# saveRDS(BWL_dat2, "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/light/25_BWL_DO_flag_sat_light.rds")

