library(dataRetrieval)
library(lubridate)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(plotly)

##==========================
## Read in DO data from GBL
##===========================
## BWL
BWL_DO <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/25_BWL_DO_flag_record.rds")

NS_plot_DO <- plot_ly(data = BWL_DO, x = ~datetime, y = ~Dissolved_Oxygen_offset, type = 'scatter',  mode = 'markers',  color = ~Concat_Date) %>%
  layout(xaxis = list(title = "Date and Time"),
         yaxis = list(title = "DO (mgL)"))

##==========================
## Flag sediment burial 
BWL_DO <- BWL_DO %>%
  mutate(fouling_flag = if_else(
    (datetime > as.POSIXct("2021-11-11 00:00:00") & datetime < as.POSIXct("2021-11-13 00:00:00")) | 
      (datetime > as.POSIXct("2022-05-08 00:00:00") & datetime < as.POSIXct("2022-06-09 00:00:00")) |
      (datetime > as.POSIXct("2023-03-07 00:00:00") & datetime < as.POSIXct("2023-03-08 00:00:00")) |
      (datetime > as.POSIXct("2023-04-05 00:00:00") & datetime < as.POSIXct("2023-04-07 00:00:00")) |
      (datetime > as.POSIXct("2023-05-06 00:00:00") & datetime < as.POSIXct("2023-06-04 00:00:00")) |
      (datetime > as.POSIXct("2024-06-05 00:00:00") & datetime < as.POSIXct("2024-06-26 00:00:00")), 
    1, 0))

temp_plot <- BWL_DO %>% filter(fouling_flag<1) %>%
  ggplot(aes(x = datetime, y = Dissolved_Oxygen_offset, colour = as.factor(maintenance_flag))) +
  geom_point(size=1, alpha=0.75) + theme_bw()

##==========================
## Bring in Streamflow
siteNo_BW <- "10336660"

# Define the parameter codes for flow and stage
pCode_flow <- "00060"
pCode_stage <- "00065"
pCode_temp <-"00010"

# Set the start date to "1980-01-01" and end date to today (current date)
start.date <- "2021-04-20"
end.date <- "2024-08-10"

HRflow_data_BW <- readNWISuv(
  siteNumbers = siteNo_BW,
  parameterCd = c(pCode_temp, pCode_flow, pCode_stage),
  startDate = start.date,
  endDate = end.date
) %>%
  select(
    datetime = "dateTime",
    dischargeCFS = "X_00060_00000",
    wtr_USGS = "X_.AquaTroll._00010_00000",
    stageF = "X_00065_00000"
  ) %>%
  mutate(datetime = as.POSIXct(datetime, tz = "UTC") %>% 
           with_tz("America/Los_Angeles")) # Convert to Pacific Time

### mutate for SI units: 
HRflow_BW <- HRflow_data_BW %>%
  mutate(
    dischargeCMS= c(dischargeCFS*0.0283168),
    scale_Q= c((dischargeCFS*0.0283168)/29.00787),
    stage_m= c(stageF*0.3048)
  ) 

## ===============================
### double check USGS timestamp:
# Filter data for a specific date range
HRflow_data_plot <- HRflow_BW %>%
  mutate(date = as.Date(datetime),
         week = week(datetime)) %>%
  filter(date > as.Date("2021-07-24") & date < as.Date("2021-10-01"))

# Create a dataframe with noon times
noon_lines <- HRflow_data_plot %>%
  group_by(date) %>%
  summarize(noon_time = as.POSIXct(paste(date, "12:00:00"), tz = "America/Los_Angeles")) %>%
  ungroup()

# Plot with noon vertical lines
HRflow_data_plot %>%
  ggplot(aes(y = wtr_USGS, x = datetime)) +
  geom_line(aes(y = wtr_USGS, x = datetime), col = "red") + 
  geom_line(aes(y = (discharge_USGS_cfs*10), x = datetime), col = "#46b8b8",alpha = .7, size = 0.75) +
  geom_vline(data = noon_lines, aes(xintercept = as.numeric(noon_time)), 
             linetype = "dashed", color = "grey50") +
  theme_bw() +
  facet_wrap(~week, scales = "free", ncol = 1) +
  scale_x_datetime(date_breaks = "12 hours", date_labels = "%H:%M") 



## ===============================
### read in morphology data:
morph_dat <- read.csv("R:/Users/kloria/Documents/Stream_Metab_24/DO_calibration_dat/stream_morphology.csv") %>%
  mutate(
    datetime = as.POSIXct(Date, format="%Y-%m-%dT%H:%M:%SZ"),
    date= as.Date(datetime)) %>%
  filter(Site=="BWL")

hist(morph_dat$Discharge_cms)


morphQ <- morph_dat %>%
  left_join(HRflow_BW, by = c("datetime"))

plot(morphQ$dischargeCMS, morphQ$Discharge_cms)

plot(morphQ$w, morphQ$stage_m)

plot(morphQ$w, morphQ$z)

plot(morphQ$stage_m, morphQ$z)


stage_lm <- lm(dischargeCMS ~ Discharge_cms, data=morphQ)
summary(stage_lm)


stage_lm <- lm(stage_m ~ z, data=morphQ)
summary(stage_lm)

## effective mean depth and stage
HRflow_na <- HRflow_BW %>%
  filter(stage_m < 2.9)

plot(HRflow_na$z, HRflow_na$stage_m)

stage_lm <- lm(stage_m ~ z, data=HRflow_na)
summary(stage_lm)


HRflow_na$depth <- c(HRflow_na$stage_m * c(0.49898/2.5))
plot(HRflow_na$z, HRflow_na$depth)

HRflow_BW$depth <- c(HRflow_BW$stage_m * c(0.49898/2.5))



HRflow_BW1 <- HRflow_BW %>%
  dplyr::select(datetime, dischargeCFS, stageF, dischargeCMS, scale_Q, stage_m,
                depth,
                #Site, 
                Discharge, Temp, z, w, v)

### merge DO and flow data
BWL_DO1 <- BWL_DO%>%
  left_join(HRflow_BW1, by = c("datetime"))
str(BWL_DO1)
range((BWL_DO1$date))

BWL_DOQ <- as.data.frame(BWL_DO1)

##==========================
## Infill observations btwn the 15 minutes
BWL_DOQ1<- BWL_DOQ %>%
  fill(dischargeCFS,.direction = "down")%>%
  fill(dischargeCMS,.direction = "down")%>% 
  fill(scale_Q,.direction = "down")%>% 
  fill(depth,.direction = "down")%>%
  dplyr::ungroup()

str(BWL_DOQ)

BWL_DOQ2 <- BWL_DOQ1 %>%
  fill(dischargeCFS,.direction = "up")%>%
  fill(dischargeCMS,.direction = "up")%>% 
  fill(scale_Q,.direction = "up")%>% 
  fill(depth,.direction = "up")%>%
  dplyr::ungroup()


##==========================
##  Summarize daily discharge statistics
BWL_sum <- BWL_DOQ2 %>%
  filter(fouling_flag<1 & maintenance_flag<1) %>%
  dplyr::group_by(date) %>%
  summarise(
    day_mean_flow = mean(dischargeCMS, na.rm = TRUE),
    day_mean_temp = mean(wtr, na.rm = TRUE))

Flow_sum <- BWL_DOQ2 %>%
  summarise(
    day_mean_flow = mean(dischargeCMS, na.rm = TRUE))

BWL_DO_Q <- BWL_DOQ2 %>%
  left_join(BWL_sum, by = "date") %>%
  mutate(
    discharge_D50_flag = if_else(dischargeCMS > c(1.5 * day_mean_flow), 1, 0),
    discharge_D100_flag = if_else(dischargeCMS > c(2 * day_mean_flow), 1, 0),
    discharge_T25_flag = if_else(dischargeCMS > c(1.25 * 1.065283), 1, 0),
    discharge_T10_flag = if_else(dischargeCMS > c(1.1 * 1.065283), 1, 0),
    discharge_T50_flag = if_else(dischargeCMS > c(1.5 * 1.065283), 1, 0),
    discharge_T100_flag = if_else(dischargeCMS > c(2 * 1.065283), 1, 0),
    discharge_T150_flag = if_else(dischargeCMS > c(2.5 * 1.065283), 1, 0),
    wtr_D50_flag = if_else(wtr < c(0.5 * day_mean_temp), 1, 0),
    wtr_D100_flag = if_else(wtr < c(2 * day_mean_temp), 1, 0))

temp_plot <- BWL_DO_Q %>% filter(datetime > as.POSIXct("2023-09-02 00:00:00")) %>%
  ggplot(aes(x = datetime, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_D100_flag))) +
  geom_point(size=1, alpha=0.75) + theme_bw()

temp_plot <- BWL_DO_Q %>% filter(datetime > as.POSIXct("2023-09-02 00:00:00")) %>%
  ggplot(aes(x = dischargeCMS, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_D100_flag))) +
  geom_point(size=1, alpha=0.75) + theme_bw()

temp_plot <- BWL_DO_Q %>% filter(datetime > as.POSIXct("2023-09-02 00:00:00")) %>%
  ggplot(aes(x = dischargeCMS, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_T25_flag))) +
  geom_point(size=1, alpha=0.75) + theme_bw()

temp_plot <- BWL_DO_Q %>% filter(datetime > as.POSIXct("2023-09-02 00:00:00")) %>%
  ggplot(aes(x = dischargeCMS, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_T100_flag))) +
  geom_point(size=1, alpha=0.75) + theme_bw()

temp_plot <- BWL_DO_Q %>% filter(datetime > as.POSIXct("2023-09-02 00:00:00")) %>%
  ggplot(aes(x = dischargeCMS, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_T150_flag))) +
  geom_point(size=1, alpha=0.75) + theme_bw()

### save data:
# saveRDS(BWL_DO_Q, "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/24_BWL_DO_flag_recordv2.rds")


## BWU 
BWU_DO <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/24_BWU_DO_flag_record.rds")

NS_plot_DO <- plot_ly(data = BWU_DO, x = ~datetime, y = ~Dissolved_Oxygen_offset, type = 'scatter',  mode = 'markers',  color = ~Concat_Date) %>%
  layout(xaxis = list(title = "Date and Time"),
         yaxis = list(title = "DO (mgL)"))

##==========================
## Flag sediment burial 
BWU_DO <- BWU_DO %>%
  mutate(fouling_flag = if_else(
    (datetime > as.POSIXct("2021-08-01 00:00:00") & datetime < as.POSIXct("2021-08-12 00:00:00")) | 
      (datetime > as.POSIXct("2022-07-08 00:00:00") & datetime < as.POSIXct("2022-07-09 00:00:00")) | 
      (datetime > as.POSIXct("2022-09-11 00:00:00") & datetime < as.POSIXct("2022-10-13 00:00:00")) |
      (datetime > as.POSIXct("2022-11-21 00:00:00") & datetime < as.POSIXct("2023-01-29 00:00:00"))|
      (datetime > as.POSIXct("2023-04-25 00:00:00") & datetime < as.POSIXct("2023-04-27 00:00:00")), 
    1, 0))

temp_plot <- BWU_DO %>% filter(fouling_flag<1) %>%
  ggplot(aes(x = datetime, y = Dissolved_Oxygen_offset, colour = as.factor(maintenance_flag))) +
  geom_point(size=1, alpha=0.75) + theme_bw()

morph_dat <- read.csv("R:/Users/kloria/Documents/Stream_Metab_24/DO_calibration_dat/stream_morphology.csv") %>%
  mutate(
    datetime = as.POSIXct(Date, format="%Y-%m-%dT%H:%M:%SZ"),
    date= as.Date(datetime)) %>%
  filter(Site=="BWU")

HRflow_BW <- HRflow_BW %>%
  left_join(morph_dat, by = c("datetime"))

## effective mean depth and stage
HRflow_na <- HRflow_BW %>%
  filter(stage_m < 2.9)

plot(HRflow_na$z, HRflow_na$stage_m)

stage_lm <- lm(stage_m ~ z, data=HRflow_na)
summary(stage_lm)


HRflow_na$depth <- c(HRflow_na$stage_m * c(0.1111111))
plot(HRflow_na$z, HRflow_na$depth)

HRflow_BW$depth <- c(HRflow_BW$stage_m * c(0.1111111))
HRflow_BW$dischargeCMS <- c(HRflow_BW$dischargeCMS/(2))


HRflow_BW1 <- HRflow_BW %>%
  dplyr::select(datetime, dischargeCFS, stageF, dischargeCMS, scale_Q, stage_m,
                depth)

### merge DO and flow data
BWU_DO1 <- BWU_DO%>%
  left_join(HRflow_BW1, by = c("datetime"))
str(BWU_DO1)


## Infill observations btwn the 15 minutes
BWU_DOQ1<- BWU_DO1 %>%
  fill(dischargeCFS,.direction = "down")%>%
  fill(dischargeCMS,.direction = "down")%>% 
  fill(scale_Q,.direction = "down")%>% 
  fill(depth,.direction = "down")%>%
  dplyr::ungroup()

str(BWU_DOQ1)

BWU_DOQ2<- BWU_DOQ1 %>%
  fill(dischargeCFS,.direction = "up")%>%
  fill(dischargeCMS,.direction = "up")%>% 
  fill(scale_Q,.direction = "up")%>% 
  fill(depth,.direction = "up")%>%
  dplyr::ungroup()

summary(BWU_DOQ2)

##==========================
##  Summarize daily discharge statistics
BWU_sum <- BWU_DOQ2 %>%
  filter(fouling_flag<1 & maintenance_flag<1) %>%
  dplyr::group_by(date) %>%
  summarise(
    day_mean_flow = mean(dischargeCMS, na.rm = T),
    day_mean_temp = mean(wtr, na.rm = TRUE))

BWU_sumT <- BWU_DOQ2 %>%
  filter(fouling_flag<1 & maintenance_flag<1) %>%
  summarise(
    mean_flow = mean(dischargeCMS, na.rm = TRUE),
    mean_temp = mean(wtr, na.rm = TRUE))

BWU_DO_Q <- BWU_DOQ2 %>%
  left_join(BWU_sum, by = "date") %>%
  mutate(
    discharge_D50_flag = if_else(dischargeCMS > c(1.5 * day_mean_flow), 1, 0),
    discharge_D100_flag = if_else(dischargeCMS > c(2 * day_mean_flow), 1, 0),
    discharge_T25_flag = if_else(dischargeCMS > c(1.25 * 0.779), 1, 0),
    discharge_T10_flag = if_else(dischargeCMS > c(1.1 * 0.779), 1, 0),
    discharge_T50_flag = if_else(dischargeCMS > c(1.5 * 0.779), 1, 0),
    discharge_T100_flag = if_else(dischargeCMS > c(2 * 0.779), 1, 0),
    discharge_T150_flag = if_else(dischargeCMS > c(2.5 * 0.779), 1, 0))

temp_plot <- BWU_DO_Q %>% filter(fouling_flag<1 & maintenance_flag<1) %>%
  ggplot(aes(x = datetime, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_T50_flag))) +
  geom_point(size=1, alpha=0.75) + theme_bw()


temp_plot <- BWU_DO_Q %>% filter(fouling_flag<1 & maintenance_flag<1) %>%
  ggplot(aes(x = dischargeCMS, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_50_flag))) +
  geom_point(size=1, alpha=0.75) + theme_bw()

# saveRDS(BWU_DO_Q, "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/24_BWU_DO_flag_recordv3.rds")


