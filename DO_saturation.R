
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
  select(site, datetime, baro_Pa)

str(stream_baro_dat)

stream_baro_dat %>%
  ggplot(aes(x = datetime, y = baro_Pa, color=as.factor(site))) +
  geom_line() + geom_point(alpha=0.75) + theme_bw() + 
  theme(legend.position = "right") 
  #facet_wrap(~ site)

stream_baro_dat$baro_mmHg <- c(stream_baro_dat$baro_Pa* 0.00750062)

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
  left_join(stream_wq_dat, by = c("datetime", "site"))

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
start_datetime <- as.POSIXct("2020-09-27 07:15:00", tz = "America/Los_Angeles")
end_datetime <- as.POSIXct("2024-08-16 07:00:00", tz = "America/Los_Angeles")

# Create a sequence of datetime values in 15-minute intervals
datetime_seq <- seq(from = start_datetime, to = end_datetime, by = "15 mins")

# Convert to a dataframe
datetime_df <- data.frame(datetime = datetime_seq)


###===========================
## Read in DO data 
###===========================

BW_DO_dat <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/flagged/24_BWL_DO_flag_recordv2.rds")

BW_df <- datetime_df %>%
  left_join(BW_DO_dat,  by = c("datetime"))


# Round datetime to the nearest 15 minutes
BW_DO_dat1 <- BW_df %>%
  mutate(
    datetime = round_date(datetime, unit = "15 minutes")
  ) %>%
  # Group by rounded datetime and aggregate selected columns to their nearest value within each 15-minute interval
  group_by(datetime, site, maintenance_flag, DO_outlier, wtr_outlier, 
           fouling_flag, discharge_D50_flag, 
           discharge_D100_flag, discharge_T25_flag, discharge_T10_flag, discharge_T50_flag, 
           discharge_T100_flag, discharge_T150_flag, wtr_D50_flag, wtr_D100_flag) %>%
  summarise(
    do.obs = mean(do.obs, na.rm = TRUE),
    wtr = mean(wtr, na.rm = TRUE),
    dischargeCMS = mean(dischargeCMS, na.rm = TRUE),
    stage_m = mean(stage_m, na.rm = TRUE),
    depth = mean(depth, na.rm = TRUE),
    v = mean(v, na.rm = TRUE),
    .groups = "drop"
  )

names(BW_DO_dat1)

BW_DO_dat1 <- BW_DO_dat1%>%
  select(site, datetime, do.obs, wtr, dischargeCMS, stage_m, depth, v,
         maintenance_flag, DO_outlier, wtr_outlier, 
         fouling_flag, discharge_D50_flag, 
         discharge_D100_flag, discharge_T25_flag, discharge_T10_flag, discharge_T50_flag, 
         discharge_T100_flag, discharge_T150_flag)


BW_dat <-BW_DO_dat1 %>%
  left_join(BWL_spc, by = c("datetime")) %>%
  left_join(baro_dat_BW, by = c("datetime", "site"))

names(BW_dat)

str(BW_dat)

BW_dat %>%
  ggplot(aes(x = datetime, y = sal_PSU, color=as.factor(site))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right")

BW_dat %>%
  ggplot(aes(x = datetime, y = do.obs, color=as.factor(site))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right")

BW_dat$sal_PSUi <- na.approx(BW_dat$sal_PSU, x = BW_dat$datetime, na.rm = FALSE)


BW_datq <- BW_dat %>%
  arrange(datetime) %>%
  mutate(
    # Calculate 6-month rolling mean of SPC (in days, approx. 182 days)
    SPC_6mo_mean = rollapply(sal_PSU, width = 24, FUN = mean, 
                             fill = NA, align = "right", na.rm = TRUE),
    # Fill NA values in SPC with the 6-month rolling mean where SPC is NA
    sal_PSU2 = ifelse(is.na(sal_PSU), SPC_6mo_mean, sal_PSU)
  ) %>%
  select(-SPC_6mo_mean) # Optionally, remove the helper column

BW_dat %>%
  ggplot(aes(x = datetime, y = sal_PSUi, color=as.factor(site))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right")


BW_datq %>%
  ggplot(aes(x = datetime, y = sal_PSU2, color=as.factor(site))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right")

BW_datq<- BW_datq %>%
  fill(sal_PSUi,.direction = "up")%>%
  fill(dischargeCMS,.direction = "up")%>% 
  fill(stage_m,.direction = "up")%>% 
  fill(depth,.direction = "up")%>%
  fill(sal_PSU2,.direction = "up")%>%
  fill(corrected_baro_mb,.direction = "up")%>%
  dplyr::ungroup()

summary(BW_datq)

BW_datq<- BW_datq %>%
  fill(sal_PSUi,.direction = "down")%>%
  fill(dischargeCMS,.direction = "down")%>% 
  fill(stage_m,.direction = "down")%>% 
  fill(depth,.direction = "down")%>%
  fill(sal_PSU2,.direction = "down")%>%
  fill(corrected_baro_mb,.direction = "down")%>%
  dplyr::ungroup()

BW_datq$DO.sat <- calc_DO_sat(BW_datq$wtr, 
                             BW_datq$corrected_baro_mb,
                          BW_datq$sal_PSU2, 
                          model = "garcia-benson") 
hist(BW_datq$DO.sat)
?calc_DO_sat()

BW_datq$year <- year(BW_datq$datetime)
BW_datq$yday <- yday(BW_datq$datetime)


BW_datq <- BW_datq %>%
  mutate(DO_supersaturated = ifelse(do.obs > DO.sat, TRUE, FALSE))


BW_datq %>%
  ggplot(aes(x = yday, y = DO.sat, color=as.factor(DO_supersaturated))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right") + 
  facet_grid(year~.)

BW_datq %>%
  filter(dischargeCMS<9)%>%
  ggplot(aes(x = dischargeCMS, y = DO.sat, color=as.factor(DO_supersaturated), 
             shape=as.factor(discharge_T50_flag))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right") +  facet_grid(year~.)


share_sat <- BW_datq %>%
  #filter(maintenance_flag<1 & fouling_flag <1)%>%
  ggplot(aes(x = yday, y = DO.sat, color=as.factor(DO_supersaturated), 
             shape=as.factor(discharge_T50_flag))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right") +  facet_grid(year~.)

# saveRDS(BW_datq, "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/saturation/24_BWL_DO_flag_sat.rds")


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

BWU_DO_dat <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/24_BWU_DO_flag_recordv3.rds")

BWU_df <- datetime_df %>%
  left_join(BWU_DO_dat,  by = c("datetime"))

BWU_df <- BWU_df%>%
  filter(datetime > as.POSIXct("2021-06-29 08:05:00 PDT"))

###
BWU_DO_dat1 %>%
  ggplot(aes(x = datetime , y = do.obs)) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right")


# Round datetime to the nearest 15 minutes
BWU_DO_dat1 <- BWU_df %>%
  mutate(
    datetime = round_date(datetime, unit = "15 minutes")
  ) %>%
  # Group by rounded datetime and aggregate selected columns to their nearest value within each 15-minute interval
  dplyr::group_by(datetime, site, maintenance_flag, DO_outlier, wtr_outlier, 
           fouling_flag, discharge_D50_flag, 
           discharge_D100_flag, discharge_T25_flag, discharge_T10_flag, discharge_T50_flag, 
           discharge_T100_flag, discharge_T150_flag) %>%
  summarise(
    do.obs = mean(do.obs, na.rm = TRUE),
    wtr = mean(wtr, na.rm = TRUE),
    dischargeCMS = mean(dischargeCMS, na.rm = TRUE),
    stage_m = mean(stage_m, na.rm = TRUE),
    depth = mean(depth, na.rm = TRUE),
    v = mean(v, na.rm = TRUE),
    .groups = "drop"
  )

names(BWU_DO_dat1)

BWU_DO_dat1 <- BWU_DO_dat1%>%
  select(site, datetime, do.obs, wtr, dischargeCMS, stage_m, depth, v,
         maintenance_flag, DO_outlier, wtr_outlier, 
         fouling_flag, discharge_D50_flag, 
         discharge_D100_flag, discharge_T25_flag, discharge_T10_flag, discharge_T50_flag, 
         discharge_T100_flag, discharge_T150_flag)


BWU_dat <-BWU_DO_dat1 %>%
  left_join(BWU_spc, by = c("datetime")) %>%
  left_join(baro_dat_BWU, by = c("datetime", "site"))

names(BWU_dat)

str(BWU_dat)

BWU_dat %>%
  ggplot(aes(x = datetime, y = sal_PSU, color=as.factor(site))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right")

BWU_dat %>%
  ggplot(aes(x = datetime, y = do.obs, color=as.factor(site))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right")

BWU_dat$sal_PSUi <- na.approx(BWU_dat$sal_PSU, x = BWU_dat$datetime, na.rm = FALSE)


BWU_datq <- BWU_dat %>%
  arrange(datetime) %>%
  mutate(
    # Calculate 6-month rolling mean of SPC (in days, approx. 182 days)
    SPC_6mo_mean = rollapply(sal_PSU, width = 24, FUN = mean, 
                             fill = NA, align = "right", na.rm = TRUE),
    # Fill NA values in SPC with the 6-month rolling mean where SPC is NA
    sal_PSU2 = ifelse(is.na(sal_PSU), SPC_6mo_mean, sal_PSU)
  ) %>%
  select(-SPC_6mo_mean) # Optionally, remove the helper column

BWU_dat %>%
  ggplot(aes(x = datetime, y = sal_PSUi, color=as.factor(site))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right")


BWU_datq %>%
  ggplot(aes(x = datetime, y = sal_PSU2, color=as.factor(site))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right")

BWU_datq1 <- BWU_datq  %>%
  fill(sal_PSUi,.direction = "up")%>%
  fill(dischargeCMS,.direction = "up")%>% 
  fill(stage_m,.direction = "up")%>% 
  fill(depth,.direction = "up")%>%
  fill(sal_PSU2,.direction = "up")%>%
  fill(corrected_baro_mb,.direction = "up")%>%
  dplyr::ungroup()

summary(BWU_datq)

BWU_datq<- BWU_datq1 %>%
  fill(sal_PSUi,.direction = "down")%>%
  fill(dischargeCMS,.direction = "down")%>% 
  fill(stage_m,.direction = "down")%>% 
  fill(depth,.direction = "down")%>%
  fill(sal_PSU2,.direction = "down")%>%
  fill(corrected_baro_mb,.direction = "down")%>%
  dplyr::ungroup()

BWU_datq$DO.sat <- calc_DO_sat(BWU_datq$wtr, 
                              BWU_datq$corrected_baro_mb,
                              BWU_datq$sal_PSU2, 
                              model = "garcia-benson") 
hist(BWU_datq$DO.sat)
?calc_DO_sat()

BWU_datq$year <- year(BWU_datq$datetime)
BWU_datq$yday <- yday(BWU_datq$datetime)


BWU_datq <- BWU_datq %>%
  mutate(DO_supersaturated = ifelse(do.obs > DO.sat, TRUE, FALSE))


BWU_datq %>%
  ggplot(aes(x = yday, y = DO.sat, color=as.factor(DO_supersaturated))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right") + 
  facet_grid(year~.)

BWU_datq %>%
  filter(dischargeCMS<9)%>%
  ggplot(aes(x = dischargeCMS, y = DO.sat, color=as.factor(DO_supersaturated), 
             shape=as.factor(discharge_T50_flag))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right") +  facet_grid(year~.)


share_sat <- BWU_datq %>%
  #filter(maintenance_flag<1 & fouling_flag <1)%>%
  ggplot(aes(x = yday, y = DO.sat, color=as.factor(DO_supersaturated), 
             shape=as.factor(discharge_T50_flag))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right") +  facet_grid(year~.)

# saveRDS(BWU_datq, "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/saturation/24_BWU_DO_flag_sat.rds")




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

GBL_DO_dat <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/flagged/24_GBL_DO_flag_recordv2.rds")

datetime_df <- datetime_df %>%
  filter(datetime> as.POSIXct("2021-03-12 06:00:00"))


GBL_df <- datetime_df %>%
  left_join(GBL_DO_dat,  by = c("datetime"))

# Round datetime to the nearest 15 minutes
GBL_DO_dat1 <- GBL_df %>%
  mutate(
    datetime = round_date(datetime, unit = "15 minutes")
  ) %>%
  # Group by rounded datetime and aggregate selected columns to their nearest value within each 15-minute interval
  group_by(datetime, site, maintenance_flag, DO_outlier, wtr_outlier, 
           fouling_flag, 
           discharge_D50_flag, 
           discharge_D100_flag, 
           discharge_T25_flag, 
           discharge_T50_flag,
           discharge_T100_flag) %>%
  summarise(
    do.obs = mean(do.obs, na.rm = TRUE),
    wtr = mean(wtr, na.rm = TRUE),
    dischargeCMS = mean(dischargeCMS, na.rm = TRUE),
    stage_m = mean(stage_m, na.rm = TRUE),
    depth = mean(depth, na.rm = TRUE),
    v = mean(v, na.rm = TRUE),
    .groups = "drop"
  )


GBL_DO_dat1 <- GBL_DO_dat1%>%
  select(site, datetime, do.obs, wtr, dischargeCMS, stage_m, depth, v,
         maintenance_flag, DO_outlier, wtr_outlier, 
         fouling_flag, discharge_D50_flag, 
         discharge_D100_flag, discharge_T25_flag, discharge_T50_flag, 
         discharge_T100_flag)


GBL_dat <-GBL_DO_dat1 %>%
  left_join(GBL_spc, by = c("datetime")) %>%
  left_join(baro_dat_GBL, by = c("datetime", "site"))

names(GBL_dat)
str(GBL_dat)

GBL_dat %>%
  ggplot(aes(x = datetime, y = sal_PSU, color=as.factor(site))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right")

GBL_dat %>%
  ggplot(aes(x = datetime, y = do.obs, color=as.factor(site))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right")

GBL_dat$sal_PSUi <- na.approx(GBL_dat$sal_PSU, x = GBL_dat$datetime, na.rm = FALSE)


GBL_datq <- GBL_dat %>%
  arrange(datetime) %>%
  mutate(
    # Calculate 6-month rolling mean of SPC (in days, approx. 182 days)
    SPC_6mo_mean = rollapply(sal_PSU, width = 24, FUN = mean, 
                             fill = NA, align = "right", na.rm = TRUE),
    # Fill NA values in SPC with the 6-month rolling mean where SPC is NA
    sal_PSU2 = ifelse(is.na(sal_PSU), SPC_6mo_mean, sal_PSU)
  ) %>%
  select(-SPC_6mo_mean) # Optionally, remove the helper column

GBL_datq %>%
  ggplot(aes(x = datetime, y = sal_PSUi, color=as.factor(site))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right")


GBL_datq %>%
  ggplot(aes(x = datetime, y = sal_PSU2, color=as.factor(site))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right")

GBL_datq <- GBL_datq  %>%
  fill(sal_PSUi,.direction = "up")%>%
  fill(dischargeCMS,.direction = "up")%>% 
  fill(stage_m,.direction = "up")%>% 
  fill(depth,.direction = "up")%>%
  fill(sal_PSU2,.direction = "up")%>%
  fill(corrected_baro_mb,.direction = "up")%>%
  dplyr::ungroup()

summary(GBL_datq)

GBL_datq<- GBL_datq %>%
  fill(sal_PSUi,.direction = "down")%>%
  fill(dischargeCMS,.direction = "down")%>% 
  fill(stage_m,.direction = "down")%>% 
  fill(depth,.direction = "down")%>%
  fill(sal_PSU2,.direction = "down")%>%
  fill(corrected_baro_mb,.direction = "down")%>%
  dplyr::ungroup()

GBL_datq$DO.sat <- calc_DO_sat(GBL_datq$wtr, 
                               GBL_datq$corrected_baro_mb,
                               GBL_datq$sal_PSU2, 
                               model = "garcia-benson") 
hist(GBL_datq$DO.sat)


GBL_datq$year <- year(GBL_datq$datetime)
GBL_datq$yday <- yday(GBL_datq$datetime)

GBL_datq <- GBL_datq %>%
  mutate(DO_supersaturated = ifelse(do.obs > DO.sat, TRUE, FALSE))


GBL_datq %>%
  ggplot(aes(x = yday, y = DO.sat, color=as.factor(DO_supersaturated))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right") + 
  facet_grid(year~.)

GBL_datq %>%
  filter(dischargeCMS<9)%>%
  ggplot(aes(x = dischargeCMS, y = DO.sat, color=as.factor(DO_supersaturated), 
             shape=as.factor(discharge_T50_flag))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right") +  facet_grid(year~.)


share_sat <- GBL_datq %>%
  #filter(maintenance_flag<1 & fouling_flag <1)%>%
  ggplot(aes(x = yday, y = DO.sat, color=as.factor(DO_supersaturated), 
             shape=as.factor(discharge_T50_flag))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right") +  facet_grid(year~.)

# saveRDS(GBL_datq, "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/saturation/24_GBL_DO_flag_sat.rds")




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

GBU_DO_dat <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/flagged/24_GBU_DO_flag_recordv2.rds")

GBU_df <- datetime_df %>%
  left_join(GBU_DO_dat,  by = c("datetime"))

# Round datetime to the nearest 15 minutes
GBU_DO_dat1 <- GBU_df %>%
  mutate(
    datetime = round_date(datetime, unit = "15 minutes")
  ) %>%
  # Group by rounded datetime and aggregate selected columns to their nearest value within each 15-minute interval
  group_by(datetime, site, maintenance_flag, DO_outlier, wtr_outlier, 
           fouling_flag, 
           discharge_D50_flag, 
           discharge_D100_flag, 
           discharge_T25_flag, 
           discharge_T50_flag,
           discharge_T100_flag) %>%
  summarise(
    do.obs = mean(do.obs, na.rm = TRUE),
    wtr = mean(wtr, na.rm = TRUE),
    dischargeCMS = mean(dischargeCMS, na.rm = TRUE),
    stage_m = mean(stage_m, na.rm = TRUE),
    depth = mean(depth, na.rm = TRUE),
    v = mean(v, na.rm = TRUE),
    .groups = "drop"
  )


GBU_DO_dat1 <- GBU_DO_dat1%>%
  select(site, datetime, do.obs, wtr, dischargeCMS, stage_m, depth, v,
         maintenance_flag, DO_outlier, wtr_outlier, 
         fouling_flag, discharge_D50_flag, 
         discharge_D100_flag, discharge_T25_flag, discharge_T50_flag, 
         discharge_T100_flag)


GBU_dat <-GBU_DO_dat1 %>%
  left_join(GBU_spc, by = c("datetime")) %>%
  left_join(baro_dat_GBU, by = c("datetime", "site"))

names(GBU_dat)
str(GBU_dat)

GBU_dat %>%
  ggplot(aes(x = datetime, y = sal_PSU, color=as.factor(site))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right")

GBU_dat %>%
  ggplot(aes(x = datetime, y = do.obs, color=as.factor(site))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right")

GBU_dat$sal_PSUi <- na.approx(GBU_dat$sal_PSU, x = GBU_dat$datetime, na.rm = FALSE)


GBU_datq <- GBU_dat %>%
  arrange(datetime) %>%
  mutate(
    # Calculate 6-month rolling mean of SPC (in days, approx. 182 days)
    SPC_6mo_mean = rollapply(sal_PSU, width = 24, FUN = mean, 
                             fill = NA, align = "right", na.rm = TRUE),
    # Fill NA values in SPC with the 6-month rolling mean where SPC is NA
    sal_PSU2 = ifelse(is.na(sal_PSU), SPC_6mo_mean, sal_PSU)
  ) %>%
  select(-SPC_6mo_mean) # Optionally, remove the helper column

GBU_datq %>%
  ggplot(aes(x = datetime, y = sal_PSUi, color=as.factor(site))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right")


GBU_datq %>%
  ggplot(aes(x = datetime, y = sal_PSU2, color=as.factor(site))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right")

GBU_datq <- GBU_datq  %>%
  fill(sal_PSUi,.direction = "up")%>%
  fill(dischargeCMS,.direction = "up")%>% 
  fill(stage_m,.direction = "up")%>% 
  fill(depth,.direction = "up")%>%
  fill(sal_PSU2,.direction = "up")%>%
  fill(corrected_baro_mb,.direction = "up")%>%
  dplyr::ungroup()

summary(GBU_datq)

GBU_datq<- GBU_datq %>%
  fill(sal_PSUi,.direction = "down")%>%
  fill(dischargeCMS,.direction = "down")%>% 
  fill(stage_m,.direction = "down")%>% 
  fill(depth,.direction = "down")%>%
  fill(sal_PSU2,.direction = "down")%>%
  fill(corrected_baro_mb,.direction = "down")%>%
  dplyr::ungroup()

GBU_datq$DO.sat <- calc_DO_sat(GBU_datq$wtr, 
                               GBU_datq$corrected_baro_mb,
                               GBU_datq$sal_PSU2, 
                               model = "garcia-benson") 
hist(GBU_datq$DO.sat)


GBU_datq$year <- year(GBU_datq$datetime)
GBU_datq$yday <- yday(GBU_datq$datetime)

GBU_datq <- GBU_datq %>%
  mutate(DO_supersaturated = ifelse(do.obs > DO.sat, TRUE, FALSE))


GBU_datq %>%
  ggplot(aes(x = yday, y = DO.sat, color=as.factor(DO_supersaturated))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right") + 
  facet_grid(year~.)

GBL_datq %>%
  filter(dischargeCMS<9)%>%
  ggplot(aes(x = dischargeCMS, y = DO.sat, color=as.factor(DO_supersaturated), 
             shape=as.factor(discharge_T50_flag))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right") +  facet_grid(year~.)


share_sat <- GBL_datq %>%
  #filter(maintenance_flag<1 & fouling_flag <1)%>%
  ggplot(aes(x = yday, y = DO.sat, color=as.factor(DO_supersaturated), 
             shape=as.factor(discharge_T50_flag))) +
  geom_point(alpha=0.1) + theme_bw() + 
  theme(legend.position = "right") +  facet_grid(year~.)

# saveRDS(GBU_datq, "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/saturation/24_GBU_DO_flag_sat_v2.rds")


