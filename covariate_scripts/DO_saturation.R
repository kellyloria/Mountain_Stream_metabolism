
# Load packages 
library(StreamMetabolism)
library(streamMetabolizer)
# double check issues with fxns here: https://github.com/DOI-USGS/streamMetabolizer/tree/main/R
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(rstan)
library(unitted)
library(zoo)
library(lubridate)

##==========================
## Read in baro data 
##==========================
stream_baro_dat <- read.csv("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/climate_dat/stream_NLDAS_baro.csv") %>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC")) %>%
  with_tz(tz = "America/Los_Angeles") %>%
  filter(site=="GBL") %>%
  select(datetime, baro_Pa)

str(stream_baro_dat)

stream_baro_dat %>%
  ggplot(aes(x = datetime, y = baro_Pa, color=as.factor(site))) +
  geom_line() + geom_point(alpha=0.75) + theme_bw() + 
  theme(legend.position = "right") 
  #facet_wrap(~ site)

stream_baro_dat$baro_mmHg <- c(stream_baro_dat$baro_Pa* 0.00750062)



##==========================
## Double check the timestamps of all temp data
##==========================

temp_data_plot <- stream_baro_dat %>%
  mutate(date = as.Date(datetime),
         week = week(datetime)) %>%
  filter(date > as.Date("2023-09-13") & date < as.Date("2023-09-27"))


# Create a dataframe with noon times
noon_lines <- temp_data_plot %>%
  group_by(date) %>%
  summarize(noon_time = as.POSIXct(paste(date, "12:00:00"), tz = "America/Los_Angeles")) %>%
  ungroup()


# Plot with noon vertical lines
#  the pressure is at its lowest around 4 a.m./p.m., and at its highest around 10 a.m./p.m.
temp_data_plot %>%
  ggplot(aes(y = baro_Pa, x = datetime)) +
  geom_line(aes(y = baro_Pa, x = datetime), col = "red") +
  # geom_line(aes(y = wtr, x = datetime), col = "goldenrod", size=1.5, alpha = .7) +
  # geom_point(aes(y = do.obs*1.5, x = datetime), col = "purple", size=1.5) + 
  geom_vline(data = noon_lines, aes(xintercept = as.numeric(noon_time)), 
             linetype = "dashed", color = "grey50") +
  theme_bw() +
  facet_wrap(~week, scales = "free", ncol = 1) +
  scale_x_datetime(date_breaks = "12 hours", date_labels = "%H:%M") +
  labs(color = "Legend")
##==========================

stream_wq_dat <- read.csv("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/WaterQual_dat/YSIWaterQuality.csv") %>%
  mutate(date=as.Date(date, format("%m/%d/%y")))
str(stream_wq_dat)

stream_wq_dat$datetime <- as.POSIXct(paste(stream_wq_dat$date, stream_wq_dat$time), 
                                     format = "%Y-%m-%d %H:%M:%S", 
                                     tz = "America/Los_Angeles") #

# Round the 'datetime' column to the nearest hour and create 'datetime_hr'
stream_wq_dat$datetime <- round_date(stream_wq_dat$datetime, unit = "hour")
stream_wq_dat <- stream_wq_dat%>%
  filter(Site=="BWL" | Site=="GBL" | Site == "GBU" | Site=="BWU") %>%
  dplyr::rename(site="Site") %>%
  select(site, mmHg, datetime)

baro_dat2 <- stream_baro_dat %>%
  left_join(stream_wq_dat, by = c("datetime"))

baro_dat_BW <- baro_dat2 %>%
  filter(site=="BWL") 

bp_BW <- lm(mmHg~baro_mmHg,data=baro_dat_BW)
summary(bp_BW)

# Extract coefficients from the linear model
intercept <- coef(bp_BW)[1]
slope <- coef(bp_BW)[2]

# Adjust baro_mmHg using the linear relationship
baro_dat_BW$corrected_baro_mmHg <- (intercept + 18.3675)+ slope * baro_dat_BW$baro_mmHg

baro_dat_BW %>%
  ggplot(aes(x = corrected_baro_mmHg, y = mmHg, color=as.factor(site))) +
  geom_point(alpha=0.7) + theme_bw() + 
  theme(legend.position = "right")  +
  facet_wrap(~ site)



baro_dat_BW$corrected_baro_mb <- c(baro_dat_BW$corrected_baro_mmHg *1.33322)
range(baro_dat_BW$corrected_baro_mb)

baro_datG <- stream_baro_dat %>%filter(site == "GBL")

baro_G <- baro_datG %>%
  mutate(site = ifelse(site == "GBL", "BWU", site))

baro_dat_BWU <- baro_G %>%
  left_join(stream_wq_dat, by = c("datetime", "site")) %>%
  filter(site=="BWU")

bp_BW <- lm(mmHg~baro_mmHg,data=baro_dat_BWU)
summary(bp_BW)

# Extract coefficients from the linear model
intercept <- (coef(bp_BW)[1])
slope <- coef(bp_BW)[2]

# Adjust baro_mmHg using the linear relationship
baro_dat_BWU$corrected_baro_mmHg <- (intercept-3.5) + slope * baro_dat_BWU$baro_mmHg

baro_dat_BWU %>%
  ggplot(aes(x = corrected_baro_mmHg, y = mmHg, color=as.factor(site))) +
  geom_point(alpha=0.7) + theme_bw() + 
  theme(legend.position = "right")  +
  facet_wrap(~ site)

baro_dat_BWU$corrected_baro_mb <- c(baro_dat_BWU$corrected_baro_mmHg *1.33322)
range(baro_dat_BWU$corrected_baro_mb)


baro_dat_GBL <- baro_dat2 %>% filter(site=="GBL")
bp_GB <- lm(mmHg~baro_mmHg,data=baro_dat_GBL)
summary(bp_GB)

# Extract coefficients from the linear model
intercept <- coef(bp_GB)[1]
slope <- coef(bp_GB)[2]

# Adjust baro_mmHg using the linear relationship
baro_dat_GBL$corrected_baro_mmHg <- intercept + slope * baro_dat_GBL$baro_mmHg

baro_dat_GBL %>%
  ggplot(aes(x = corrected_baro_mmHg, y = mmHg, color=as.factor(site))) +
  geom_point(alpha=0.7) + theme_bw() + 
  theme(legend.position = "right")  +
  facet_wrap(~ site)


baro_dat_GBL$corrected_baro_mb <- c(baro_dat_GBL$corrected_baro_mmHg *1.33322)

baro_dat_GBU <- baro_dat2 %>% filter(site=="GBU")
bp_GB <- lm(mmHg~baro_mmHg,data=baro_dat_GBU)
summary(bp_GB)

# Extract coefficients from the linear model
intercept <- coef(bp_GB)[1]
slope <- coef(bp_GB)[2]

# Adjust baro_mmHg using the linear relationship
baro_dat_GBU$corrected_baro_mmHg <- intercept + slope * baro_dat_GBU$baro_mmHg

baro_dat_GBU %>%
  ggplot(aes(x = corrected_baro_mmHg, y = mmHg, color=as.factor(site))) +
  geom_point(alpha=0.7) + theme_bw() + 
  theme(legend.position = "right")  +
  facet_wrap(~ site)

baro_dat_GBU$corrected_baro_mb <- c(baro_dat_GBU$corrected_baro_mmHg *1.33322)
range(baro_dat_GBU$corrected_baro_mb )

###===========================
# look at what an elevation based range would mean: 
cal_bp <- calc_air_pressure( 
  temp.air = u(33, "degC"),
  elevation = u(1923.288, "m")
)

# roughly 25 range:
# Range BWL (1907.438 m): 793.2117 (-7 C) to 818.9951 (33 C)
# Range BWU (1940.966 m): 789.7994 to 815.9314 

# Range GBL (1902.866 m): 793.6714 (-7 C) to 819.4077 (33 C) 
##  Looks good [1] range(baro_dat_GBL$corrected_baro_mb) 796.8397 819.2899
# Range GBU (1923.288 m): 791.5936 to 817.5426 
## Looks good 799.7196 816.1404

###===========================
baro_dat_GBU
baro_dat_GBU
baro_dat_BWU
baro_dat_BW

baro_dat_C <-rbind(baro_dat_GBU, baro_dat_GBU, baro_dat_BWU, baro_dat_BW)

baro_dat_C %>%
  ggplot(aes(x = datetime, y = corrected_baro_mb, color=as.factor(site))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right")  # + facet_wrap(~ site)


# saveRDS(baro_dat_C, "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/climate_dat/24_corrected_baro_record.rds")

###===========================
BWL_spc <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/WaterQual_dat/24_BWL_SPCv2.rds")

est_salinity <- function(spc, temp) {
  # Correct SPC to 25°C
  spc_25 <- spc / (1 + 0.0191 * (temp - 25))
  
  # Calculate salinity in PSU for freshwater to brackish conditions
  salinity <- 0.5 * 10^-3 * spc_25
  return(salinity)
}

BWL_spc$sal_PSU <- est_salinity(BWL_spc$SPC, BWL_spc$wtr)

hist(BWL_spc$sal_PSU)

BWL_spc<- BWL_spc%>%
  select(datetime,SPC, sal_PSU)


## DO 
# saveRDS(BWL_DO_Q, "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/24_BWL_DO_flag_recordv2.rds")


###===========================
## create an empty vector of times 
# Define the start and end dates
start_datetime <- as.POSIXct("2020-09-27 07:05:00", tz = "America/Los_Angeles")
end_datetime <- as.POSIXct("2024-08-16 06:55:00", tz = "America/Los_Angeles")

# Create a sequence of datetime values in 15-minute intervals
datetime_seq <- seq(from = start_datetime, to = end_datetime, by = "5 mins")

# Convert to a dataframe
datetime_df <- data.frame(datetime = datetime_seq)



start_datetime_GB <- as.POSIXct("2021-03-11 07:05:00", tz = "America/Los_Angeles")
end_datetime_GB <- as.POSIXct("2024-10-16 06:55:00", tz = "America/Los_Angeles")

# Create a sequence of datetime values in 15-minute intervals
datetime_seq_GB <- seq(from = start_datetime_GB, to = end_datetime_GB, by = "5 mins")

# Convert to a dataframe
datetime_df_GB <- data.frame(datetime = datetime_seq_GB)




###===========================
## Read in DO data 
###===========================

BW_DO_dat <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/25_BWL_DO_flag_flow.rds")%>%
  filter(wtr_outlier<1 & DO_outlier <1 ) %>%
  mutate(
    maintenance_flag = ifelse(is.na(maintenance_flag), 0, maintenance_flag),
    fouling_flag = ifelse(is.na(fouling_flag), 0, fouling_flag)
  )
  

BW_DO_dat_15 <- BW_DO_dat %>%
filter(datetime > as.POSIXct("2021-04-29 12:15:00")) %>%  
# Step 1: Create a column for rounding to the nearest 15 minutes
mutate(
  rounded_datetime_15min = as.POSIXct(
    round(as.numeric(datetime) / (15 * 60)) * (15 * 60),
    origin = "1970-01-01",
    tz = "America/Los_Angeles")) %>%
  group_by(rounded_datetime_15min, maintenance_flag, fouling_flag) %>%
  summarize(
    do.obs = mean(do.obs, na.rm = TRUE),
    wtr = mean(wtr, na.rm = TRUE),
    dischargeCMS = mean(adjusted_dischargeCMS, na.rm = TRUE),
    depth = mean(adjusted_depth, na.rm = TRUE),
    wtr_USGS = mean(wtr_USGS, na.rm = TRUE),
    w = mean(w, na.rm = TRUE),
    v = mean(v_estimation, na.rm = TRUE))
    

BW_dat_15q <-BW_DO_dat_15 %>%
  left_join(BWL_spc, by = c("rounded_datetime_15min" = "datetime")) %>%
  left_join(stream_baro_dat, by = c("rounded_datetime_15min" = "datetime"))

names(BW_dat_15q)

str(BW_dat_15q)

BW_dat_15q %>%
  ggplot(aes(x = rounded_datetime_15min, y = sal_PSU)) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right")

BW_dat_15q %>%
  ggplot(aes(x = rounded_datetime_15min, y = do.obs)) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right")

BW_dat_15q$sal_PSUi <- na.approx(BW_dat_15q$sal_PSU, x = BW_dat_15q$rounded_datetime_15min, na.rm = FALSE)

BW_dat_15q1 %>%
  ggplot(aes(x = rounded_datetime_15min, y = sal_PSUi)) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right")

BW_datq<- BW_dat_15q1 %>%
  fill(sal_PSUi,.direction = "up")%>%
  fill(baro_mmHg,.direction = "up")%>%
  dplyr::ungroup()

summary(BW_datq)

BW_datq1<- BW_datq %>%
  fill(sal_PSUi,.direction = "down")%>%
  fill(baro_mmHg,.direction = "down")%>%
  dplyr::ungroup()

summary(BW_datq1)

BW_datqb<-  BW_datq1%>% 
  mutate(
  maintenance_flag = ifelse(is.na(maintenance_flag), 0, maintenance_flag),
  fouling_flag = ifelse(is.na(fouling_flag), 0, fouling_flag))



library(dplyr)
library(zoo)
library(lubridate)


BW_datq_a <- BW_datqb %>%
  filter(rounded_datetime_15min > as.POSIXct("2021-04-29 12:30:00"))%>%
  group_by(rounded_datetime_15min, maintenance_flag, fouling_flag) %>%
  summarize(
    do.obs = mean(do.obs, na.rm = TRUE),
    wtr = mean(wtr, na.rm = TRUE),
    dischargeCMS = mean(dischargeCMS, na.rm = TRUE),
    depth = mean(depth, na.rm = TRUE),
    wtr_USGS = mean(wtr_USGS, na.rm = TRUE),
    w = mean(w, na.rm = TRUE),
    v = mean(v, na.rm = TRUE),
    SPC= mean(SPC, na.rm = TRUE),
    sal_PSUi= mean(sal_PSUi, na.rm = TRUE),
    baro_mmHg = mean(baro_mmHg, na.rm = TRUE))


# Specify columns to infill
columns_to_infill <- c("do.obs", "wtr", "dischargeCMS", "depth", 
                       "wtr_USGS", "w", "v", "SPC", "baro_mmHg", "sal_PSUi")

# Define the rolling 3-hour window infill function
infill_rolling_avg <- function(data, columns, window_size = 3 * 60) {
  data %>%
    arrange(rounded_datetime_15min) %>%  # Ensure data is sorted by time
    group_by(rounded_datetime_15min) %>%
    mutate(across(all_of(columns), 
                  ~ ifelse(is.na(.), 
                           rollapplyr(.x, width = window_size / 15, 
                                      FUN = mean, fill = NA, na.rm = TRUE, partial = TRUE), 
                           .x))) %>%
    ungroup()
}

# Apply the infill function to the dataset
BW_datq2 <- infill_rolling_avg(BW_datq_a, columns_to_infill)

summary(BW_datq2)
summary(BW_datq_a)


BW_datq2$corrected_baro_mb <-c(BW_datq2$baro_mmHg*1.33322)
BW_datq3 <- BW_datq2 %>% 
  mutate(
    sal_PSUi = ifelse(sal_PSUi<0.0001, 0.0001, sal_PSUi))

summary(BW_datq3)

###########
###########

BW_datq3$DO.sat <- calc_DO_sat(BW_datq3$wtr, 
                               BW_datq3$corrected_baro_mb,
                               BW_datq3$sal_PSUi, 
                          model = "garcia-benson") 
hist(BW_datq3$DO.sat)
?calc_DO_sat()

BW_datq3$year <- year(BW_datq3$rounded_datetime_15min)
BW_datq3$yday <- yday(BW_datq3$rounded_datetime_15min)


BWU_datq4<- BW_datq3 %>%
  fill(dischargeCMS,.direction = "down")%>%
  fill(depth,.direction = "down")%>%
  dplyr::ungroup()

summary(BWU_datq4)



BWU_datq4 <- BWU_datq4 %>%
  mutate(DO_supersaturated = ifelse(do.obs > DO.sat, TRUE, FALSE))


BWU_datq4 %>%
  ggplot(aes(x = yday, y = DO.sat, color=as.factor(DO_supersaturated))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right") + 
  facet_grid(year~.)

BW_datq3 %>%
  filter(dischargeCMS<9)%>%
  ggplot(aes(x = dischargeCMS, y = DO.sat, color=as.factor(DO_supersaturated), 
             shape=as.factor(discharge_T50_flag))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right") +  facet_grid(year~.)


summary(BW_datq3)

# saveRDS(BWU_datq4, "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/saturation/25_BWL_DO_flag_sat.rds")


###===========================
### BWU
###===========================

BWU_spc <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/WaterQual_dat/24_BWU_SPCv2.rds")
BWU_spc$sal_PSU <- est_salinity(BWU_spc$SPC, BWU_spc$wtr)

hist(BWU_spc$sal_PSU)

BWU_spc<- BWU_spc%>%
  select(datetime, SPC, sal_PSU)

###===========================
## Read in DO data 
###===========================

# "R:\Users\kloria\Documents\Stream_Metab_24\Core_sites\offset_DO_dat\24_BWU_DO_flag_recordv3.rds"
# "R:\Users\kloria\Documents\Stream_Metab_24\Core_sites\offset_DO_dat\25_BWU_DO_flag_flow.rds"
BWU_DO_dat <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/25_BWU_DO_flag_flow.rds")%>%
  filter(wtr_outlier<1 & DO_outlier <1 ) %>%
  mutate(
    maintenance_flag = ifelse(is.na(maintenance_flag), 0, maintenance_flag),
    fouling_flag = ifelse(is.na(fouling_flag), 0, fouling_flag)
  )

BWU_DO_dat_15 <- BWU_DO_dat %>%
  filter(datetime > as.POSIXct("2021-04-29 12:15:00")) %>%  
  # Step 1: Create a column for rounding to the nearest 15 minutes
  mutate(
    rounded_datetime_15min = as.POSIXct(
      round(as.numeric(datetime) / (15 * 60)) * (15 * 60),
      origin = "1970-01-01",
      tz = "America/Los_Angeles")) %>%
  group_by(rounded_datetime_15min, maintenance_flag, fouling_flag) %>%
  summarize(
    do.obs = mean(do.obs, na.rm = TRUE),
    wtr = mean(wtr, na.rm = TRUE),
    dischargeCMS = mean(adjusted_dischargeCMS, na.rm = TRUE),
    depth = mean(adjusted_depth, na.rm = TRUE),
    wtr_USGS = mean(wtr_USGS, na.rm = TRUE),
    w = mean(w, na.rm = TRUE),
    v = mean(v_estimation, na.rm = TRUE))


BWU_dat_15q <-BWU_DO_dat_15 %>%
  left_join(BWU_spc, by = c("rounded_datetime_15min" = "datetime")) %>%
  left_join(stream_baro_dat, by = c("rounded_datetime_15min" = "datetime"))

names(BWU_dat_15q)

str(BWU_dat_15q)

BWU_dat_15q %>%
  ggplot(aes(x = rounded_datetime_15min, y = sal_PSU)) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right")

###

BWU_dat_15q$sal_PSUi <- na.approx(BWU_dat_15q$sal_PSU, x = BWU_dat_15q$rounded_datetime_15min, na.rm = FALSE)

BWU_dat_15q %>%
  ggplot(aes(x = rounded_datetime_15min, y = sal_PSUi)) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right")

BWU_datq<- BWU_dat_15q %>%
  fill(sal_PSUi,.direction = "up")%>%
  fill(baro_mmHg,.direction = "up")%>%
  dplyr::ungroup()

summary(BWU_datq)

BWU_datq1<- BWU_datq %>%
  fill(sal_PSUi,.direction = "down")%>%
  fill(baro_mmHg,.direction = "down")%>%
  dplyr::ungroup()

summary(BWU_datq1)

BWU_datqb<-  BWU_datq1%>% 
  mutate(
    maintenance_flag = ifelse(is.na(maintenance_flag), 0, maintenance_flag),
    fouling_flag = ifelse(is.na(fouling_flag), 0, fouling_flag))

library(dplyr)
library(zoo)
library(lubridate)

BWU_datq_a <- BWU_datqb %>%
  filter(rounded_datetime_15min > as.POSIXct("2021-04-29 12:30:00"))%>%
  group_by(rounded_datetime_15min, maintenance_flag, fouling_flag) %>%
  summarize(
    do.obs = mean(do.obs, na.rm = TRUE),
    wtr = mean(wtr, na.rm = TRUE),
    dischargeCMS = mean(dischargeCMS, na.rm = TRUE),
    depth = mean(depth, na.rm = TRUE),
    wtr_USGS = mean(wtr_USGS, na.rm = TRUE),
    w = mean(w, na.rm = TRUE),
    v = mean(v, na.rm = TRUE),
    SPC= mean(SPC, na.rm = TRUE),
    sal_PSUi= mean(sal_PSUi, na.rm = TRUE),
    baro_mmHg = mean(baro_mmHg, na.rm = TRUE))


# Specify columns to infill
columns_to_infill <- c("do.obs", "wtr", "dischargeCMS", "depth", 
                       "wtr_USGS", "w", "v", "SPC", "baro_mmHg", "sal_PSUi")

# Define the rolling 3-hour window infill function
infill_rolling_avg <- function(data, columns, window_size = 3 * 60) {
  data %>%
    arrange(rounded_datetime_15min) %>%  # Ensure data is sorted by time
    group_by(rounded_datetime_15min) %>%
    mutate(across(all_of(columns), 
                  ~ ifelse(is.na(.), 
                           rollapplyr(.x, width = window_size / 15, 
                                      FUN = mean, fill = NA, na.rm = TRUE, partial = TRUE), 
                           .x))) %>%
    ungroup()
}

# Apply the infill function to the dataset
BWU_datq2 <- infill_rolling_avg(BWU_datq_a, columns_to_infill)


summary(BWU_datq2)
summary(BWU_datq_a)


BWU_datq3<- BWU_datq2 %>%
  fill(dischargeCMS,.direction = "down")%>%
  fill(depth,.direction = "down")%>%
  dplyr::ungroup()

summary(BWU_datq3)


BWU_datq3$corrected_baro_mb <-c(BWU_datq3$baro_mmHg*1.33322)
BWU_datq4 <- BWU_datq3 %>% 
  mutate(
    sal_PSUi = ifelse(sal_PSUi<0.0001, 0.0001, sal_PSUi),
    sal_PSUi = ifelse(is.na(sal_PSUi), 0.0001, sal_PSUi)
    )

BWU_datq4$corrected_baro_mb <-c(BWU_datq4$baro_mmHg*1.33322)


summary(BWU_datq4)


###########

BWU_datq4$DO.sat <- calc_DO_sat(BWU_datq4$wtr, 
                               BWU_datq4$corrected_baro_mb,
                               BWU_datq4$sal_PSUi, 
                               model = "garcia-benson") 
hist(BWU_datq4$DO.sat)


BWU_datq4$year <- year(BWU_datq4$rounded_datetime_15min)
BWU_datq4$yday <- yday(BWU_datq4$rounded_datetime_15min)


BWU_datq4 <- BWU_datq4 %>%
  mutate(DO_supersaturated = ifelse(do.obs > DO.sat, TRUE, FALSE))

BWU_datq4 %>%
  ggplot(aes(x = yday, y = DO.sat, color=as.factor(DO_supersaturated))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right") + 
  facet_grid(year~.)


summary(BWU_datq4)

# saveRDS(BWU_datq4, "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/saturation/25_BWU_DO_flag_sat.rds")


###===========================
### GBL
###===========================

GBL_spc <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/WaterQual_dat/24_GBL_SPCv2.rds")

GBL_spc$sal_PSU <- est_salinity(GBL_spc$SPC, GBL_spc$wtr)

hist(GBL_spc$sal_PSU)

GBL_spc<- GBL_spc%>%
  select(datetime, SPC, sal_PSU)

###===========================
## Read in DO data 
###===========================

GBL_DO_dat <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/25_GBL_DO_flag_flow.rds")%>%
  filter(wtr_outlier<1 & DO_outlier <1 ) %>%
  mutate(
    maintenance_flag = ifelse(is.na(maintenance_flag), 0, maintenance_flag),
    fouling_flag = ifelse(is.na(fouling_flag), 0, fouling_flag)
  )

GBL_DO_dat_15 <- GBL_DO_dat %>%
 # filter(datetime > as.POSIXct("2021-02-29 12:15:00")) %>%  
  # Step 1: Create a column for rounding to the nearest 15 minutes
  mutate(
    rounded_datetime_15min = as.POSIXct(
      round(as.numeric(datetime) / (15 * 60)) * (15 * 60),
      origin = "1970-01-01",
      tz = "America/Los_Angeles")) %>%
  group_by(rounded_datetime_15min, maintenance_flag, fouling_flag) %>%
  summarize(
    do.obs = mean(do.obs, na.rm = TRUE),
    wtr = mean(wtr, na.rm = TRUE),
    dischargeCMS = mean(adjusted_dischargeCMS, na.rm = TRUE),
    depth = mean(adjusted_depth, na.rm = TRUE),
    wtr_USGS = mean(wtr_USGS, na.rm = TRUE),
    w = mean(w, na.rm = TRUE),
    v = mean(v_estimation, na.rm = TRUE))


GBL_dat_15q <-GBL_DO_dat_15 %>%
  left_join(GBL_spc, by = c("rounded_datetime_15min" = "datetime")) %>%
  left_join(stream_baro_dat, by = c("rounded_datetime_15min" = "datetime"))

names(GBL_dat_15q)

str(GBL_dat_15q)

GBL_dat_15q %>%
  ggplot(aes(x = rounded_datetime_15min, y = sal_PSU)) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right")

###

GBL_dat_15q$sal_PSUi <- na.approx(GBL_dat_15q$sal_PSU, x = GBL_dat_15q$rounded_datetime_15min, na.rm = FALSE)

GBL_dat_15q %>%
  ggplot(aes(x = rounded_datetime_15min, y = sal_PSUi)) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right")

GBL_datq<- GBL_dat_15q %>%
  fill(sal_PSUi,.direction = "up")%>%
  fill(baro_mmHg,.direction = "up")%>%
  dplyr::ungroup()

summary(GBL_datq)

GBL_datq1<- GBL_datq %>%
  fill(sal_PSUi,.direction = "down")%>%
  fill(baro_mmHg,.direction = "down")%>%
  dplyr::ungroup()

summary(GBL_datq1)

GBL_datqb<-  GBL_datq1%>% 
  mutate(
    maintenance_flag = ifelse(is.na(maintenance_flag), 0, maintenance_flag),
    fouling_flag = ifelse(is.na(fouling_flag), 0, fouling_flag))

library(dplyr)
library(zoo)
library(lubridate)

GBL_datq_a <- GBL_datqb %>%
  filter(rounded_datetime_15min > as.POSIXct("2021-04-29 12:30:00"))%>%
  group_by(rounded_datetime_15min, maintenance_flag, fouling_flag) %>%
  summarize(
    do.obs = mean(do.obs, na.rm = TRUE),
    wtr = mean(wtr, na.rm = TRUE),
    dischargeCMS = mean(dischargeCMS, na.rm = TRUE),
    depth = mean(depth, na.rm = TRUE),
    wtr_USGS = mean(wtr_USGS, na.rm = TRUE),
    w = mean(w, na.rm = TRUE),
    v = mean(v, na.rm = TRUE),
    SPC= mean(SPC, na.rm = TRUE),
    sal_PSUi= mean(sal_PSUi, na.rm = TRUE),
    baro_mmHg = mean(baro_mmHg, na.rm = TRUE))


# Specify columns to infill
columns_to_infill <- c("do.obs", "wtr", "dischargeCMS", "depth", 
                       "wtr_USGS", "w", "v", "SPC", "baro_mmHg", "sal_PSUi")

# Define the rolling 3-hour window infill function
infill_rolling_avg <- function(data, columns, window_size = 3 * 60) {
  data %>%
    arrange(rounded_datetime_15min) %>%  # Ensure data is sorted by time
    group_by(rounded_datetime_15min) %>%
    mutate(across(all_of(columns), 
                  ~ ifelse(is.na(.), 
                           rollapplyr(.x, width = window_size / 15, 
                                      FUN = mean, fill = NA, na.rm = TRUE, partial = TRUE), 
                           .x))) %>%
    ungroup()
}

# Apply the infill function to the dataset
GBL_datq2 <- infill_rolling_avg(GBL_datq_a, columns_to_infill)


summary(GBL_datq2)
summary(GBL_datq_a)


GBL_datq3<- GBL_datq2 %>%
  fill(dischargeCMS,.direction = "down")%>%
  fill(depth,.direction = "down")%>%
  dplyr::ungroup()

summary(GBL_datq3)


GBL_datq4 <- GBL_datq3 %>% 
  mutate(
    sal_PSUi = ifelse(sal_PSUi<0.0001, 0.0001, sal_PSUi),
    sal_PSUi = ifelse(is.na(sal_PSUi), 0.0001, sal_PSUi)
  )

GBL_datq4$corrected_baro_mb <-c(GBL_datq4$baro_mmHg*1.33322)


summary(GBL_datq4)


###########

GBL_datq4$DO.sat <- calc_DO_sat(GBL_datq4$wtr, 
                                GBL_datq4$corrected_baro_mb,
                                GBL_datq4$sal_PSUi, 
                                model = "garcia-benson") 
hist(GBL_datq4$DO.sat)


GBL_datq4$year <- year(GBL_datq4$rounded_datetime_15min)
GBL_datq4$yday <- yday(GBL_datq4$rounded_datetime_15min)


GBL_datq4 <- GBL_datq4 %>%
  mutate(DO_supersaturated = ifelse(do.obs > DO.sat, TRUE, FALSE))

GBL_datq4 %>%
  ggplot(aes(x = yday, y = DO.sat, color=as.factor(DO_supersaturated))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right") + 
  facet_grid(year~.)


summary(GBL_datq4)

# saveRDS(GBL_datq4, "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/saturation/25_GBL_DO_flag_sat.rds")



###===========================
### GBU
###===========================

GBU_spc <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/WaterQual_dat/24_GBL_SPCv2.rds")

GBU_spc$sal_PSU <- est_salinity(GBU_spc$SPC, GBU_spc$wtr)

hist(GBU_spc$sal_PSU)

GBU_spc<- GBU_spc%>%
  select(datetime, SPC, sal_PSU)

###===========================
## Read in DO data 
###===========================
GBU_DO_dat <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/25_GBU_DO_flag_flow.rds")%>%
  filter(wtr_outlier<1 & DO_outlier <1 ) %>%
  mutate(
    maintenance_flag = ifelse(is.na(maintenance_flag), 0, maintenance_flag),
    fouling_flag = ifelse(is.na(fouling_flag), 0, fouling_flag)
  )

GBU_DO_dat_15 <- GBU_DO_dat %>%
  # filter(datetime > as.POSIXct("2021-02-29 12:15:00")) %>%  
  # Step 1: Create a column for rounding to the nearest 15 minutes
  mutate(
    rounded_datetime_15min = as.POSIXct(
      round(as.numeric(datetime) / (15 * 60)) * (15 * 60),
      origin = "1970-01-01",
      tz = "America/Los_Angeles")) %>%
  group_by(rounded_datetime_15min, maintenance_flag, fouling_flag) %>%
  summarize(
    do.obs = mean(do.obs, na.rm = TRUE),
    wtr = mean(wtr, na.rm = TRUE),
    dischargeCMS = mean(adjusted_dischargeCMS, na.rm = TRUE),
    depth = mean(adjusted_depth, na.rm = TRUE),
    wtr_USGS = mean(wtr_USGS, na.rm = TRUE),
    w = mean(w, na.rm = TRUE),
    v = mean(v_estimation, na.rm = TRUE))


GBU_dat_15q <-GBU_DO_dat_15 %>%
  left_join(GBU_spc, by = c("rounded_datetime_15min" = "datetime")) %>%
  left_join(stream_baro_dat, by = c("rounded_datetime_15min" = "datetime"))

names(GBU_dat_15q)

str(GBU_dat_15q)

GBU_dat_15q %>%
  ggplot(aes(x = rounded_datetime_15min, y = sal_PSU)) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right")

###

GBU_dat_15q$sal_PSUi <- na.approx(GBU_dat_15q$sal_PSU, x = GBU_dat_15q$rounded_datetime_15min, na.rm = FALSE)

GBU_dat_15q %>%
  ggplot(aes(x = rounded_datetime_15min, y = sal_PSUi)) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right")

GBU_datq<- GBU_dat_15q %>%
  fill(sal_PSUi,.direction = "up")%>%
  fill(baro_mmHg,.direction = "up")%>%
  dplyr::ungroup()

summary(GBU_datq)

GBU_datq1<- GBU_datq %>%
  fill(sal_PSUi,.direction = "down")%>%
  fill(baro_mmHg,.direction = "down")%>%
  dplyr::ungroup()

summary(GBU_datq1)

GBU_datqb<-  GBU_datq1%>% 
  mutate(
    maintenance_flag = ifelse(is.na(maintenance_flag), 0, maintenance_flag),
    fouling_flag = ifelse(is.na(fouling_flag), 0, fouling_flag))

library(dplyr)
library(zoo)
library(lubridate)

GBU_datq_a <- GBU_datqb %>%
  filter(rounded_datetime_15min > as.POSIXct("2021-04-29 12:30:00"))%>%
  group_by(rounded_datetime_15min, maintenance_flag, fouling_flag) %>%
  summarize(
    do.obs = mean(do.obs, na.rm = TRUE),
    wtr = mean(wtr, na.rm = TRUE),
    dischargeCMS = mean(dischargeCMS, na.rm = TRUE),
    depth = mean(depth, na.rm = TRUE),
    wtr_USGS = mean(wtr_USGS, na.rm = TRUE),
    w = mean(w, na.rm = TRUE),
    v = mean(v, na.rm = TRUE),
    SPC= mean(SPC, na.rm = TRUE),
    sal_PSUi= mean(sal_PSUi, na.rm = TRUE),
    baro_mmHg = mean(baro_mmHg, na.rm = TRUE))


# Specify columns to infill
columns_to_infill <- c("do.obs", "wtr", "dischargeCMS", "depth", 
                       "wtr_USGS", "w", "v", "SPC", "baro_mmHg", "sal_PSUi")

# Define the rolling 3-hour window infill function
infill_rolling_avg <- function(data, columns, window_size = 3 * 60) {
  data %>%
    arrange(rounded_datetime_15min) %>%  # Ensure data is sorted by time
    group_by(rounded_datetime_15min) %>%
    mutate(across(all_of(columns), 
                  ~ ifelse(is.na(.), 
                           rollapplyr(.x, width = window_size / 15, 
                                      FUN = mean, fill = NA, na.rm = TRUE, partial = TRUE), 
                           .x))) %>%
    ungroup()
}

# Apply the infill function to the dataset
GBU_datq2 <- infill_rolling_avg(GBU_datq_a, columns_to_infill)


summary(GBU_datq2)
summary(GBU_datq_a)


GBU_datq3<- GBU_datq2 %>%
  fill(dischargeCMS,.direction = "down")%>%
  fill(depth,.direction = "down")%>%
  dplyr::ungroup()

summary(GBU_datq3)


GBU_datq4 <- GBU_datq3 %>% 
  mutate(
    sal_PSUi = ifelse(sal_PSUi<0.0001, 0.0001, sal_PSUi),
    sal_PSUi = ifelse(is.na(sal_PSUi), 0.0001, sal_PSUi)
  )

GBU_datq4$corrected_baro_mb <-c(GBU_datq4$baro_mmHg*1.33322)


summary(GBU_datq4)


###########

GBU_datq4$DO.sat <- calc_DO_sat(GBU_datq4$wtr, 
                                GBU_datq4$corrected_baro_mb,
                                GBU_datq4$sal_PSUi, 
                                model = "garcia-benson") 
hist(GBU_datq4$DO.sat)


GBU_datq4$year <- year(GBU_datq4$rounded_datetime_15min)
GBU_datq4$yday <- yday(GBU_datq4$rounded_datetime_15min)


GBU_datq4 <- GBU_datq4 %>%
  mutate(DO_supersaturated = ifelse(do.obs > DO.sat, TRUE, FALSE))

GBU_datq4 %>%
  ggplot(aes(x = yday, y = DO.sat, color=as.factor(DO_supersaturated))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right") + 
  facet_grid(year~.)


summary(GBU_datq4)

# saveRDS(GBU_datq4, "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/saturation/25_GBU_DO_flag_sat.rds")

