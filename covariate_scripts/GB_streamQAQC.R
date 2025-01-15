library(dataRetrieval)
library(lubridate)
library(ggplot2)
library(dplyr)
library(tidyverse)

##==========================
## create date sequence: 
##==========================
## create an empty vector of times 
#   Define the start and end dates
start_datetime <- as.POSIXct("2021-03-11 00:00:00", tz = "America/Los_Angeles")
end_datetime <- as.POSIXct("2024-10-16 00:00:00", tz = "America/Los_Angeles")

# Create a sequence of datetime values in 15-minute intervals
datetime_seq <- seq(from = start_datetime, to = end_datetime, by = "5 mins")

# Convert to a dataframe
datetime_df <- data.frame(datetime = datetime_seq)

##==========================
## Read in DO data from GBL
##===========================
## GBL
GBL_DO <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/25_GBL_DO_flag_record.rds")

NS_plot_DO <- plot_ly(data = GBL_DO, x = ~datetime, y = ~Dissolved_Oxygen_offset, type = 'scatter',  mode = 'markers',  color = ~Concat_Date) %>%
  layout(xaxis = list(title = "Date and Time"),
         yaxis = list(title = "DO (mgL)"))

##==========================
## Flag sediment burial 
GBL_DOQ <- GBL_DO %>%
  mutate(fouling_flag = if_else(
      (datetime > as.POSIXct("2021-07-09 00:00:00") & datetime < as.POSIXct("2021-07-22 00:00:00")) | 
      (datetime > as.POSIXct("2021-08-02 00:00:00") & datetime < as.POSIXct("2021-08-18 00:00:00")) |
      (datetime > as.POSIXct("2021-10-17 00:00:00") & datetime < as.POSIXct("2021-10-22 00:00:00")) |
      (datetime > as.POSIXct("2022-06-01 00:00:00") & datetime < as.POSIXct("2022-06-03 00:00:00")) | 
      (datetime > as.POSIXct("2022-07-04 00:00:00") & datetime < as.POSIXct("2022-07-12 00:00:00")) | 
      (datetime > as.POSIXct("2022-08-21 00:00:00") & datetime < as.POSIXct("2022-09-20 00:00:00")) | 
      (datetime > as.POSIXct("2022-10-18 00:00:00") & datetime < as.POSIXct("2022-11-04 00:00:00")) |
      (datetime > as.POSIXct("2022-11-28 00:00:00") & datetime < as.POSIXct("2022-12-12 00:00:00")) |
      (datetime > as.POSIXct("2022-12-24 00:00:00") & datetime < as.POSIXct("2023-01-24 00:00:00")) |
      (datetime > as.POSIXct("2023-03-06 00:00:00") & datetime < as.POSIXct("2023-03-08 00:00:00")) |
      (datetime > as.POSIXct("2023-04-09 00:00:00") & datetime < as.POSIXct("2023-04-27 00:00:00")) |
      (datetime > as.POSIXct("2023-06-23 00:00:00") & datetime < as.POSIXct("2023-07-11 00:00:00")) |
      (datetime > as.POSIXct("2023-09-08 00:00:00") & datetime < as.POSIXct("2023-09-11 00:00:00")) |
      (datetime > as.POSIXct("2024-09-09 00:00:00") & datetime < as.POSIXct("2024-09-20 00:00:00")), 
    1, 0))

temp_plot <- GBL_DOQ %>% filter(fouling_flag<1) %>%
  ggplot(aes(x = datetime, y = Dissolved_Oxygen_offset, colour = as.factor(maintenance_flag))) +
  geom_point(size=1, alpha=0.75) + theme_bw()

##==========================
## Bring in Streamflow
siteNo_GB <- "10336730" #siteNo_BW <- "10336660"

# Define the parameter codes for flow and stage
pCode_flow <- "00060"
pCode_stage <- "00065"
pCode_temp <-"00010"


# Set the start date to "1980-01-01" and end date to today (current date)
start.date <- "2021-03-11"
end.date <- "2024-10-16"

HRflow_data_GB <- readNWISuv(
  siteNumbers = siteNo_GB,
  parameterCd = c(pCode_temp, pCode_flow, pCode_stage),
  startDate = start.date,
  endDate = end.date
) %>%
  select(
    datetime = "dateTime",
    dischargeCFS = "X_00060_00000",
    wtr_USGS = "X_00010_00000",
    stageF = "X_00065_00000"
  ) %>%
  mutate(datetime = as.POSIXct(datetime, tz = "UTC") %>% 
           with_tz("America/Los_Angeles")) # Convert to Pacific Time

### mutate for SI units: 
HRflow_GB <- HRflow_data_GB %>%
  mutate(
    dischargeCMS= c(dischargeCFS*0.0283168),
    scale_Q= c((dischargeCFS*0.0283168)/10.64485),
    stage_m= c(stageF*0.3048)) 

## ===============================
### double check USGS timestamp:
# Filter data for a specific date range
HRflow_data_plot <- HRflow_GB %>%
  mutate(date = as.Date(datetime),
         week = week(datetime)) %>%
  filter(date > as.Date("2023-08-13") & date < as.Date("2023-08-27"))

# Create a dataframe with noon times
noon_lines <- HRflow_data_plot %>%
  group_by(date) %>%
  summarize(noon_time = as.POSIXct(paste(date, "12:00:00"), tz = "America/Los_Angeles")) %>%
  ungroup()

# Plot with noon vertical lines
HRflow_data_plot %>%
  ggplot(aes(y = wtr_USGS, x = datetime)) +
  geom_line(aes(y = wtr_USGS, x = datetime), col = "red") + 
  geom_line(aes(y = (dischargeCFS*5), x = datetime), col = "#46b8b8",alpha = .7, size = 0.75) +
  geom_vline(data = noon_lines, aes(xintercept = as.numeric(noon_time)), 
             linetype = "dashed", color = "grey50") +
  theme_bw() +
  facet_wrap(~week, scales = "free", ncol = 1) +
  scale_x_datetime(date_breaks = "12 hours", date_labels = "%H:%M") 


##==========================
## Read in morphology observations
##==========================
morph_dat_GBL <- read.csv("R:/Users/kloria/Documents/Stream_Metab_24/stream_morphology_core.csv") %>%
  mutate(datetime = as.POSIXct(datetime, format="%Y-%m-%dT%H:%M:%SZ"),
         date= as.Date(datetime)) %>% filter(Site=="GBL")

## Join with flow 
HRflow_GBL <- HRflow_GB %>%
  left_join(morph_dat_GBL, by = c("datetime"))

HRflow_GBL_plot <- morph_dat_GBL %>%
  left_join(HRflow_GB, by = c("datetime"))

plot(HRflow_GBL$dischargeCMS, HRflow_GBL$Discharge_cms)
plot(HRflow_GBL$stage_m, HRflow_GBL$z)

## quick linear models 
flow_lm <- lm(Discharge_cms ~ dischargeCMS, data=HRflow_GBL)
summary(flow_lm)

## flow plot
flow_summary <- summary(flow_lm)
slope <- round(flow_summary$coefficients[2, 1], 3) 
r_squared <- round(flow_summary$r.squared, 2) 
p_value <- round(flow_summary$coefficients[2, 4], 4)  

## Dynamic Adjustment Based on Slope
HRflow_GBL$adjusted_dischargeCMS <- c(HRflow_GBL$dischargeCMS * slope)
hist(HRflow_GBL$adjusted_dischargeCMS)
hist(HRflow_GBL$dischargeCMS)

GBLflow_plot <- ggplot(HRflow_GBL_plot, aes(x = dischargeCMS, y = Discharge_cms, color=method)) +
  geom_point(alpha=0.75, size =2) + geom_smooth(method="lm", se=F, color="grey50") +
  labs(title = "GBL", x = "USGS flow (cms)", y = "Observed flow (cms)") + theme_bw() +
  annotate(
    "text",
    x = Inf, y = Inf, 
    label = paste("Slope:", slope, "\nR²:", r_squared, "\np-value:", p_value),
    hjust = 1.1, vjust = 4.5,
    size = 3,
    color = "black")

GBLflow_plot

### depth 
stage_lm <- lm(z ~stage_m, data=HRflow_GBL)
summary(stage_lm)

## flow plot
depth_summary <- summary(stage_lm)
slope <- round(depth_summary$coefficients[2, 1], 3) 
r_squared <- round(depth_summary$r.squared, 2) 
p_value <- round(depth_summary$coefficients[2, 4], 4)  

## Dynamic Adjustment Based on Slope
HRflow_GBL$adjusted_depth <- c(HRflow_GBL$stage_m * slope)
hist(HRflow_GBL$stage_m)
hist(HRflow_GBL$adjusted_depth)

GBLdepth_plot <- ggplot(HRflow_GBL_plot, aes(x = stage_m, y = z, color=method)) +
  geom_point(alpha=0.75, size =2) + geom_smooth(method="lm", se=F, color="grey50") +
  labs(title = "GBL", x = "USGS stage (m)", y = "Obsersed depth (m)") + theme_bw() +
  annotate(
    "text",
    x = Inf, y = Inf, 
    label = paste("Slope:", slope, "\nR²:", r_squared, "\np-value:", p_value),
    hjust = 1.1, vjust = 4.5,
    size = 3,
    color = "black")

GBLdepth_plot


#### 
names(HRflow_GBL)

HRflow_GBL <- HRflow_GBL %>%
  dplyr::select(datetime, adjusted_dischargeCMS, adjusted_depth, dischargeCFS, 
                dischargeCMS,
                stageF, stage_m, z, wtr_USGS, stage_m, w, v, v_estimation)

### merge DO and flow data
GBL_dat <- datetime_df %>%
  left_join(GBL_DOQ, by = c("datetime")) %>%
  left_join(HRflow_GBL, by = c("datetime"))

range((GBL_dat$date))
str(GBL_dat)

GBL_datQ <- as.data.frame(GBL_dat)

##==========================
## Double check the timestamps of all temp data
##==========================

temp_data_plot <- GBL_dat %>%
  mutate(date = as.Date(datetime),
         week = week(datetime)) %>%
  filter(date > as.Date("2023-09-13") & date < as.Date("2023-09-27"))


# Create a dataframe with noon times
noon_lines <- temp_data_plot %>%
  group_by(date) %>%
  summarize(noon_time = as.POSIXct(paste(date, "12:00:00"), tz = "America/Los_Angeles")) %>%
  ungroup()


# Plot with noon vertical lines
temp_data_plot %>%
  ggplot(aes(y = wtr_USGS, x = datetime)) +
  geom_line(aes(y = wtr_USGS, x = datetime), col = "red") + 
  geom_line(aes(y = wtr, x = datetime), col = "goldenrod", size=1.5, alpha = .7) +
  geom_point(aes(y = do.obs*1.5, x = datetime), col = "purple", size=1.5) + 
  geom_line(aes(y = (dischargeCFS*10), x = datetime), col = "#46b8b8",alpha = .7, size = 0.75) +
  # geom_line(aes(y = (dischargeCMS*250), x = datetime), col = "darkblue",alpha = .7, size = 0.75) +
  geom_vline(data = noon_lines, aes(xintercept = as.numeric(noon_time)), 
             linetype = "dashed", color = "grey50") +
  theme_bw() +
  facet_wrap(~week, scales = "free", ncol = 1) +
  scale_x_datetime(date_breaks = "12 hours", date_labels = "%H:%M") +
  labs(color = "Legend")


##===========================
## Things look good time stamp wise ## 
## pausing here to create a new infil script 

# saveRDS(GBL_datQ, "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/25_GBL_DO_flag_flow.rds")



##==========================
## Read in DO data from GBU
##===========================
## GBL
GBU_DO <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/25_GBU_DO_flag_record.rds")

NS_plot_DO <- plot_ly(data = GBU_DO, x = ~datetime, y = ~Dissolved_Oxygen_offset, type = 'scatter',  mode = 'markers',  color = ~Concat_Date) %>%
  layout(xaxis = list(title = "Date and Time"),
         yaxis = list(title = "DO (mgL)"))

##==========================
## Flag sediment burial 
GBU_DO <- GBU_DO %>%
  mutate(fouling_flag = if_else(
    (datetime > as.POSIXct("2021-11-29 00:00:00") & datetime < as.POSIXct("2021-12-02 00:00:00")) | 
      (datetime > as.POSIXct("2022-06-01 00:00:00") & datetime < as.POSIXct("2022-06-03 00:00:00")) | 
      (datetime > as.POSIXct("2022-09-11 00:00:00") & datetime < as.POSIXct("2022-09-15 00:00:00")) | 
      (datetime > as.POSIXct("2022-09-24 00:00:00") & datetime < as.POSIXct("2022-10-02 00:00:00")) | 
      (datetime > as.POSIXct("2022-10-21 00:00:00") & datetime < as.POSIXct("2022-10-25 00:00:00")),
    1, 0))

temp_plot <- GBU_DO %>% filter(fouling_flag<1) %>%
  ggplot(aes(x = datetime, y = Dissolved_Oxygen_offset, colour = as.factor(maintenance_flag))) +
  geom_point(size=1, alpha=0.75) + theme_bw()



##==========================
## Read in morphology observations
##==========================

morph_dat_GBU <- read.csv("R:/Users/kloria/Documents/Stream_Metab_24/stream_morphology_core.csv") %>%
  mutate(datetime = as.POSIXct(datetime, format="%Y-%m-%dT%H:%M:%SZ"),
         date= as.Date(datetime)) %>% filter(Site=="GBU")

## Join with flow 
HRflow_GBU_plot <- morph_dat_GBU %>%
  left_join(HRflow_GB, by = c("datetime"))

HRflow_GBU <- HRflow_GB %>%
  left_join(morph_dat_GBU, by = c("datetime"))

plot(HRflow_GBU$dischargeCMS, HRflow_GBU$Discharge_cms)
plot(HRflow_GBU$stage_m, HRflow_GBU$z)

## quick linear models 
flow_lm <- lm(dischargeCMS ~ Discharge_cms, data=HRflow_GBU)
summary(flow_lm)

## flow plot
flow_summary <- summary(flow_lm)
slope <- round(flow_summary$coefficients[2, 1], 3) 
r_squared <- round(flow_summary$r.squared, 3) 
p_value <- round(flow_summary$coefficients[2, 4], 3)  

## Dynamic Adjustment Based on Slope
HRflow_GBU$adjusted_dischargeCMS <- c(HRflow_GBU$dischargeCMS * slope)
hist(HRflow_GBU$adjusted_dischargeCMS)
hist(HRflow_GBU$dischargeCMS)

GBUflow_plot <- ggplot(HRflow_GBU_plot, aes(x = dischargeCMS, y = Discharge_cms, color=method)) +
  geom_point(alpha=0.75, size =2) + geom_smooth(method="lm", se=F, color="grey50") +
  labs(title = "GBU", x = "USGS flow (cms)", y = "Observed flow (cms)") + theme_bw() +
  annotate(
    "text",
    x = Inf, y = Inf, 
    label = paste("Slope:", slope, "\nR²:", r_squared, "\np-value:", p_value),
    hjust = 1.1, vjust = 4.5,
    size = 3,
    color = "black")
GBUflow_plot

## depth

stage_lm <- lm(z~stage_m, data=HRflow_GBU)
summary(stage_lm)

## flow plot
depth_summary <- summary(stage_lm)
slope <- round(depth_summary$coefficients[2, 1], 3) 
r_squared <- round(depth_summary$r.squared, 3) 
p_value <- round(depth_summary$coefficients[2, 4], 3)  


## Dynamic Adjustment Based on Slope
HRflow_GBU$adjusted_depth <- c(HRflow_GBU$stage_m * slope)
hist(HRflow_GBU$stage_m)
hist(HRflow_GBU$adjusted_depth)

GBUdepth_plot <- ggplot(HRflow_GBU_plot, aes(x = stage_m, y = z, color=method)) +
  geom_point(alpha=0.75, size =2) + geom_smooth(method="lm", se=F, color="grey50") +
  labs(title = "GBU", x = "USGS stage (m)", y = "Obsersed depth (m)") + theme_bw() +
  annotate(
    "text",
    x = Inf, y = Inf, 
    label = paste("Slope:", slope, "\nR²:", r_squared, "\np-value:", p_value),
    hjust = 1.1, vjust = 4.5,
    size = 3,
    color = "black")

GBUdepth_plot


#### 
names(HRflow_GBU)

HRflow_GBU <- HRflow_GBU %>%
  dplyr::select(datetime, adjusted_dischargeCMS, adjusted_depth, dischargeCFS, 
                dischargeCMS,
                stageF, stage_m, z, wtr_USGS, stage_m, w, v, v_estimation)

### merge DO and flow data
GBU_dat <- datetime_df %>%
  left_join(GBU_DO, by = c("datetime")) %>%
  left_join(HRflow_GBU, by = c("datetime"))

range((GBU_dat$date))
str(GBU_dat)

GBU_datQ <- as.data.frame(GBU_dat)

##==========================
## Double check the timestamps of all temp data
##==========================

temp_data_plot <- GBU_dat %>%
  mutate(date = as.Date(datetime),
         week = week(datetime)) %>%
  filter(date > as.Date("2023-09-13") & date < as.Date("2023-09-27"))


# Create a dataframe with noon times
noon_lines <- temp_data_plot %>%
  group_by(date) %>%
  summarize(noon_time = as.POSIXct(paste(date, "12:00:00"), tz = "America/Los_Angeles")) %>%
  ungroup()


# Plot with noon vertical lines
temp_data_plot %>%
  ggplot(aes(y = wtr_USGS, x = datetime)) +
  geom_line(aes(y = wtr_USGS, x = datetime), col = "red") + 
  geom_line(aes(y = wtr, x = datetime), col = "goldenrod", size=1.5, alpha = .7) +
  geom_point(aes(y = do.obs*1.5, x = datetime), col = "purple", size=1.5) + 
  geom_line(aes(y = (dischargeCFS*10), x = datetime), col = "#46b8b8",alpha = .7, size = 0.75) +
  # geom_line(aes(y = (dischargeCMS*250), x = datetime), col = "darkblue",alpha = .7, size = 0.75) +
  geom_vline(data = noon_lines, aes(xintercept = as.numeric(noon_time)), 
             linetype = "dashed", color = "grey50") +
  theme_bw() +
  facet_wrap(~week, scales = "free", ncol = 1) +
  scale_x_datetime(date_breaks = "12 hours", date_labels = "%H:%M") +
  labs(color = "Legend")


##===========================
## Things look good time stamp wise ## 
## pausing here to create a new infil script 

# saveRDS(GBU_datQ, "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/25_GBU_DO_flag_flow.rds")













# # Set the start date to "1980-01-01" and end date to today (current date)
# start.date <- "2021-03-10"
# end.date <- "2024-10-10"  # Use
# ## Download hourly flow data for both sites separately and combine
# HRflow_data_GB <- readNWISuv(siteNumbers = siteNo_GB, parameterCd = c("00060", "00065"), startDate = start.date, endDate = end.date) 
# 
# HRflow_GB <- HRflow_data_GB %>%
#   select(datetime = "dateTime", dischargeCFS = "X_00060_00000", stageF = "X_00065_00000") %>%
#   mutate(
#     dischargeCMS= c(dischargeCFS*0.0283168),
#     scale_Q= c((dischargeCFS*0.0283168)/10.64485),
#     stage_m= c(stageF*0.3048)
#     ) 
# 
# 
# morph_dat <- read.csv("R:/Users/kloria/Documents/Stream_Metab_24/DO_calibration_dat/stream_morphology.csv") %>%
#   mutate(
#     datetime = as.POSIXct(Date, format="%Y-%m-%dT%H:%M:%SZ"),
#     date= as.Date(datetime)) %>%
#   filter(Site=="GBL")
# 
# HRflow_GB <- HRflow_GB %>%
#   left_join(morph_dat, by = c("datetime"))
#  
# ## effective mean depth and stage
# HRflow_na <- HRflow_GB 
# plot(HRflow_na$z, HRflow_na$stage_m)
# 
# stage_lm <- lm(z ~ stage_m, data=HRflow_na)
# summary(stage_lm)
# 
# plot(HRflow_na$z, HRflow_na$stage_m)
# 
# HRflow_na$depth <- c(HRflow_na$stage_m * 0.09446)
# plot(HRflow_na$z, HRflow_na$depth)
# 
# HRflow_GB$depth <- c(HRflow_GB$stage_m * 0.09446)
# 
# 
# #0.43423
# 
# hist(HRflow_GB$z)
# hist(HRflow_GB$stage_m)
# hist(HRflow_GB$depth)
# 
# 
# 
# HRflow_GB <- HRflow_GB%>%
#   dplyr::select(datetime, dischargeCFS, stageF, dischargeCMS, scale_Q, stage_m,
#                 #depth,
#                 #Site, 
#                 Discharge, Temp, z, w, v)
# 
# ### merge DO and flow data
# GBL_DOQ <- GBL_DOQ%>%
#   left_join(HRflow_GB, by = c("datetime"))
# str(GBL_DOQ)
# range((GBL_DOQ$date))
# 
# GBL_DOQ <- as.data.frame(GBL_DOQ)
# str(GBL_DOQ)
# 
# ##==========================
# ## Infill observations btwn the 15 minutes
# GBL_DOQ1<- GBL_DOQ %>%
#   fill(dischargeCFS,.direction = "down")%>%
#   fill(dischargeCMS,.direction = "down")%>% 
#   fill(scale_Q,.direction = "down")%>% 
#   fill(depth,.direction = "down")%>%
#   dplyr::ungroup()
# 
# str(GBL_DOQ1)
# 
# GBL_DOQ2<- GBL_DOQ1 %>%
#   fill(dischargeCFS,.direction = "up")%>%
#   fill(dischargeCMS,.direction = "up")%>% 
#   fill(scale_Q,.direction = "up")%>% 
#   fill(depth,.direction = "up")%>%
#   dplyr::ungroup()
# 
# summary(GBL_DOQ2)
# 
# ##==========================
# ##  Summarize daily discharge statistics
# GBL_sum <- GBL_DOQ2 %>%
#   filter(fouling_flag<1 & maintenance_flag<1) %>%
#   dplyr::group_by(date) %>%
#   summarise(
#     day_mean_flow = mean(dischargeCMS, na.rm = TRUE),
#     day_mean_temp = mean(wtr, na.rm = TRUE))
# 
# Flow_sum <- HRflow_GB %>%
#   summarise(
#     day_mean_flow = mean(dischargeCMS, na.rm = TRUE))
#     
# GBL_DO_Q <- GBL_DOQ2 %>%
#   left_join(GBL_sum, by = "date") %>%
#   mutate(
#     discharge_D50_flag = if_else(dischargeCMS > c(1.5 * day_mean_flow), 1, 0),
#     discharge_D100_flag = if_else(dischargeCMS > c(2 * day_mean_flow), 1, 0),
#     discharge_T25_flag = if_else(dischargeCMS > c(1.25 * 0.07115973), 1, 0),
#     discharge_T50_flag = if_else(dischargeCMS > c(1.5 * 0.07115973), 1, 0),
#     discharge_T100_flag = if_else(dischargeCMS > c(2 * 0.07115973), 1, 0),
#     wtr_D50_flag = if_else(wtr < c(0.5 * day_mean_temp), 1, 0),
#     wtr_D100_flag = if_else(wtr < c(2 * day_mean_temp), 1, 0))
# 
# temp_plot <- GBL_DO_Q %>% filter(datetime > as.POSIXct("2023-09-02 00:00:00")) %>%
#   ggplot(aes(x = datetime, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_D100_flag))) +
#   geom_point(size=1, alpha=0.75) + theme_bw()
# 
# temp_plot <- GBL_DO_Q %>% filter(datetime > as.POSIXct("2023-09-02 00:00:00")) %>%
#   ggplot(aes(x = dischargeCMS, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_D100_flag))) +
#   geom_point(size=1, alpha=0.75) + theme_bw()
# 
# temp_plot <- GBL_DO_Q %>% filter(datetime > as.POSIXct("2023-09-02 00:00:00")) %>%
#   ggplot(aes(x = dischargeCMS, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_T25_flag))) +
#   geom_point(size=1, alpha=0.75) + theme_bw()
# 
# temp_plot <- GBL_DO_Q %>% filter(datetime > as.POSIXct("2023-09-02 00:00:00")) %>%
#   ggplot(aes(x = dischargeCMS, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_T50_flag))) +
#   geom_point(size=1, alpha=0.75) + theme_bw()
# 
# temp_plot <- GBL_DO_Q %>% filter(datetime > as.POSIXct("2023-09-02 00:00:00")) %>%
#   ggplot(aes(x = dischargeCMS, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_T50_flag))) +
#   geom_point(size=1, alpha=0.75) + theme_bw()
# 
# temp_plot <- GBL_DO_Q %>% filter(datetime > as.POSIXct("2023-09-02 00:00:00")) %>%
#   ggplot(aes(x = dischargeCMS, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_T100_flag))) +
#   geom_point(size=1, alpha=0.75) + theme_bw()
# 
# temp_plot <- GBL_DO_Q %>% filter(datetime > as.POSIXct("2023-09-02 00:00:00")) %>%
#   ggplot(aes(x = wtr, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_T100_flag))) +
#   geom_point(size=1, alpha=0.75) + theme_bw()
# 
# temp_plot <- GBL_DO_Q %>% 
#   filter(wtr_D100_flag<2) %>%
#   ggplot(aes(x = wtr, y = Dissolved_Oxygen_offset, colour = as.factor(wtr_D50_flag))) +
#   geom_point(size=1, alpha=0.75) + theme_bw()
# 
# hist(GBL_DO_Q$dischargeCMS)
# 
# 
# ### save data:
# saveRDS(GBL_DO_Q, "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/24_GBL_DO_flag_record.rds")
# 
# 
# ####
# # 1. bring in Baro data
# # 3. Phil's light model 
# 
# ################
# ################
# 
# ## GBU 
# GBU_DO <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/24_GBU_DO_flag_record.rds")
# 
# NS_plot_DO <- plot_ly(data = GBU_DO, x = ~datetime, y = ~Dissolved_Oxygen_offset, type = 'scatter',  mode = 'markers',  color = ~Concat_Date) %>%
#   layout(xaxis = list(title = "Date and Time"),
#          yaxis = list(title = "DO (mgL)"))
# ##==========================
# ## Flag sediment burial 
# GBU_DOQ <- GBU_DO %>%
#   mutate(fouling_flag = if_else(
#     (datetime > as.POSIXct("2022-09-11 00:00:00") & datetime < as.POSIXct("2022-09-16 00:00:00")) | 
#       (datetime > as.POSIXct("2022-09-23 00:00:00") & datetime < as.POSIXct("2022-10-04 00:00:00")), 
#     1, 0))
# 
# temp_plot <- GBL_DOQ %>% filter(fouling_flag<1) %>%
#   ggplot(aes(x = datetime, y = Dissolved_Oxygen_offset, colour = as.factor(maintenance_flag))) +
#   geom_point(size=1, alpha=0.75) + theme_bw()
# 
# 
# morph_datU <- read.csv("R:/Users/kloria/Documents/Stream_Metab_24/DO_calibration_dat/stream_morphology.csv") %>%
#   mutate(
#     datetime = as.POSIXct(Date, format="%Y-%m-%dT%H:%M:%SZ"),
#     date= as.Date(datetime)) %>%
#   filter(Site=="GBU")
# 
# 
# GBLQ <- c(0.098, 0.152, NA, 0.071, NA)
# GBU <- c(0.098, 0.111, 0.060, 0.064, NA)
# 
# 
# # Create a data frame with these values
# GB_flow <- data.frame(GBLQ, GBU)
# 
# # Fit a linear model with GBU as the dependent variable and GBLQ as the independent variable
# lm_Q <- lm(GBU ~ GBLQ, data = GB_flow)
# 
# # View the summary to get the slope (coefficient of GBLQ)
# summary(lm_Q)
# 
# # Extract the slope from the model
# slope <- coef(lm_Q)["GBLQ"]
# 
# # 0.531746 
# 
# # View the corrected data
# print(GB_flow)
# 
# 
# HRflow_GBU <- HRflow_GB %>%
#   left_join(morph_datU, by = c("datetime"))
# 
# ### merge DO and flow data # 0.7963961
# HRflow_GBU1 <- HRflow_GBU %>%
#   mutate(
#     dischargeCFS = c(dischargeCFS * 0.7963961), 
#     dischargeCMS = c(dischargeCMS*0.7963961),
#     scale_Q = c(scale_Q*0.7963961), 
#     stage_m =  c(stage_m*0.531746))
# 
# HRflow_GBU2 <- HRflow_GBU1%>%
#   dplyr::select(datetime, dischargeCFS, stageF, dischargeCMS, scale_Q, stage_m,
#                 Site, Discharge, Temp, z, w, v)
# 
# GBU_DOQ1 <- GBU_DOQ %>%
#   left_join(HRflow_GBU2,by = c("datetime"))
# 
# 
# 
# str(GBU_DOQ1)
# range((GBU_DOQ$date))
# hist(GBU_DOQ$dischargeCMS)
# 
# 
# # dept relationship 
# GBU_DOQ1$depth <- c(GBU_DOQ1$stage_m * 0.1437)
# 
# 
# ##
# 
# 
# ## Infill observations btwn the 15 minutes
# GBU_DOQ1<- GBU_DOQ1 %>%
#   fill(dischargeCFS,.direction = "down")%>%
#   fill(dischargeCMS,.direction = "down")%>% 
#   fill(scale_Q,.direction = "down")%>% 
#   fill(depth,.direction = "down")%>%
#   dplyr::ungroup()
# 
# str(GBL_DOQ1)
# 
# GBU_DOQ2<- GBU_DOQ1 %>%
#   fill(dischargeCFS,.direction = "up")%>%
#   fill(dischargeCMS,.direction = "up")%>% 
#   fill(scale_Q,.direction = "up")%>% 
#   fill(depth,.direction = "up")%>%
#   dplyr::ungroup()
# 
# summary(GBU_DOQ2)
# 
# 
# ##==========================
# ##  Summarize daily discharge statistics
# GBU_sum <- GBU_DOQ2 %>%
#   #filter(fouling_flag<1 & maintenance_flag<1) %>%
#   dplyr::group_by(date) %>%
#   summarise(
#     day_mean_flow = mean(dischargeCMS, na.rm = T),
#     day_mean_temp = mean(wtr, na.rm = TRUE))
# 
# GBU_sumT <- GBU_DOQ2 %>%
#   filter(fouling_flag<1 & maintenance_flag<1) %>%
#   summarise(
#     mean_flow = mean(dischargeCMS, na.rm = TRUE),
#     mean_temp = mean(wtr, na.rm = TRUE))
# 
# GBU_DO_Q <- GBU_DOQ2 %>%
#   left_join(GBU_sum, by = "date") %>%
#   mutate(discharge_D50_flag = if_else(dischargeCMS > c(1.5 * day_mean_flow), 1, 0), 
#          discharge_D100_flag = if_else(dischargeCMS > c(2 * day_mean_flow), 1, 0),
#          discharge_T25_flag = if_else(dischargeCMS > c(1.25 * 0.0482), 1, 0),
#          discharge_T50_flag = if_else(dischargeCMS > c(1.5 * 0.0482), 1, 0),
#          discharge_T100_flag = if_else(dischargeCMS > c(2 * 0.0482), 1, 0),
#          wtr_75_flag = if_else(wtr > c(1.75 * day_mean_temp), 1, 0))
# 
# temp_plot <- GBU_DO_Q %>% filter(fouling_flag<1 & maintenance_flag<1) %>%
#   ggplot(aes(x = datetime, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_T50_flag))) +
#   geom_point(size=1, alpha=0.75) + theme_bw()
# 
# 
# temp_plot <- GBU_DO_Q %>% filter(fouling_flag<1 & maintenance_flag<1) %>%
#   ggplot(aes(x = dischargeCMS, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_50_flag))) +
#   geom_point(size=1, alpha=0.75) + theme_bw()
# 
# # saveRDS(GBU_DO_Q, "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/24_GBU_DO_flag_recordv2.rds")
# 
