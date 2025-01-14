library(dataRetrieval)
library(lubridate)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(plotly)


##==========================
## create date sequence: 
## create an empty vector of times 
# Define the start and end dates
start_datetime <- as.POSIXct("2020-09-27 00:00:00", tz = "America/Los_Angeles")
end_datetime <- as.POSIXct("2024-08-17 00:00:00", tz = "America/Los_Angeles")

# Create a sequence of datetime values in 15-minute intervals
datetime_seq <- seq(from = start_datetime, to = end_datetime, by = "5 mins")

# Convert to a dataframe
datetime_df <- data.frame(datetime = datetime_seq)

##==========================
## Read in DO data from BWL
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
  filter(date > as.Date("2021-08-13") & date < as.Date("2021-08-27"))

# Create a dataframe with noon times
noon_lines <- HRflow_data_plot %>%
  group_by(date) %>%
  summarize(noon_time = as.POSIXct(paste(date, "12:00:00"), tz = "America/Los_Angeles")) %>%
  ungroup()

# Plot with noon vertical lines
HRflow_data_plot %>%
  ggplot(aes(y = wtr_USGS, x = datetime)) +
  geom_line(aes(y = wtr_USGS, x = datetime), col = "red") + 
  geom_line(aes(y = (dischargeCFS*10), x = datetime), col = "#46b8b8",alpha = .7, size = 0.75) +
  geom_vline(data = noon_lines, aes(xintercept = as.numeric(noon_time)), 
             linetype = "dashed", color = "grey50") +
  theme_bw() +
  facet_wrap(~week, scales = "free", ncol = 1) +
  scale_x_datetime(date_breaks = "12 hours", date_labels = "%H:%M") 

## ===============================
### read in morphology data:
##==========================
## Read in morphology observations
morph_dat_BWL <- read.csv("R:/Users/kloria/Documents/Stream_Metab_24/stream_morphology_core.csv") %>%
  mutate(datetime = as.POSIXct(datetime, format="%Y-%m-%dT%H:%M:%SZ"),
         date= as.Date(datetime)) %>% filter(Site=="BWL")

## Join with flow 
HRflow_BWL <- HRflow_BW %>%
  left_join(morph_dat_BWL, by = c("datetime"))

plot(HRflow_BWL$dischargeCMS, HRflow_BWL$Discharge_cms)
plot(HRflow_BWL$stage_m, HRflow_BWL$z)

## quick linear models 
flow_lm <- lm(Discharge_cms~dischargeCMS, data=HRflow_BWL)
summary(flow_lm)

## flow plot
flow_summary <- summary(flow_lm)
slope <- round(flow_summary$coefficients[2, 1], 3) 
r_squared <- round(flow_summary$r.squared, 2) 
p_value <- round(flow_summary$coefficients[2, 4], 4)

## Dynamic Adjustment Based on Slope
HRflow_BWL$adjusted_dischargeCMS <- c(HRflow_BWL$dischargeCMS * slope)
hist(HRflow_BWL$adjusted_dischargeCMS)
hist(HRflow_BWL$dischargeCMS)

BWLflow_plot <- ggplot(HRflow_BWL_plot, aes(x = dischargeCMS, y = Discharge_cms, color=method)) +
  geom_point(alpha=0.75, size =2) +
  geom_smooth(method="lm", se=F, color="grey50") +
  labs(title = "BWL", x = "USGS flow (cms)", y = "Observed flow (cms)") + theme_bw() +
  annotate(
    "text",
    x = Inf, y = Inf, 
    label = paste("Slope:", slope, "\nR²:", r_squared, "\np-value:", p_value),
    hjust = 1.1, vjust = 4.5,
    size = 3,
    color = "black"
  )
BWLflow_plot

### depth :
stage_lm <- lm(z ~ stage_m, data=HRflow_BWL)
summary(stage_lm)

## flow plot
depth_summary <- summary(stage_lm)
slope <- round(depth_summary$coefficients[2, 1], 3) 
r_squared <- round(depth_summary$r.squared, 2) 
p_value <- round(depth_summary$coefficients[2, 4], 4)  

## Dynamic Adjustment Based on Slope
HRflow_BWL$adjusted_depth <- c(HRflow_BWL$stage_m * slope)
hist(HRflow_BWL$stage_m)
hist(HRflow_BWL$adjusted_depth)

BWLdepth_plot <- ggplot(HRflow_BWL_plot, aes(x = stage_m, y = z, color=method)) +
  geom_point(alpha=0.75, size =2) + geom_smooth(method="lm", se=F, color="grey50") +
  labs(title = "BWL", x = "USGS stage (m)", y = "Obsersed depth (m)") + theme_bw() +
  annotate(
    "text",
    x = Inf, y = Inf, 
    label = paste("Slope:", slope, "\nR²:", r_squared, "\np-value:", p_value),
    hjust = 1.1, vjust = 4.5,
    size = 3,
    color = "black")
BWLdepth_plot


#### 
names(HRflow_BWL)

HRflow_BWL <- HRflow_BWL %>%
  dplyr::select(datetime, adjusted_dischargeCMS, adjusted_depth, dischargeCFS, 
                dischargeCMS,
                stageF, stage_m, z, wtr_USGS, stage_m, w, v, v_estimation)

### merge DO and flow data
BWL_dat <- datetime_df %>%
  left_join(BWL_DO, by = c("datetime")) %>%
  left_join(HRflow_BWL, by = c("datetime"))

range((BWL_dat$date))
str(BWL_dat)

BWL_datQ <- as.data.frame(BWL_dat)

##==========================
## Double check the timestamps of all temp data
##==========================

temp_data_plot <- BWL_datQ %>%
  mutate(date = as.Date(datetime),
         week = week(datetime)) %>%
  filter(date > as.Date("2022-09-13") & date < as.Date("2022-09-27"))


# Create a dataframe with noon times
noon_lines <- temp_data_plot %>%
  group_by(date) %>%
  summarize(noon_time = as.POSIXct(paste(date, "12:00:00"), tz = "America/Los_Angeles")) %>%
  ungroup()


# Plot with noon vertical lines
temp_data_plot %>%
  ggplot(aes(y = wtr_USGS, x = datetime)) +
  geom_line(aes(y = wtr_USGS, x = datetime), col = "red") + 
  geom_line(aes(y = wtr, x = datetime), col = "goldenrod", size=1.5) +
  geom_point(aes(y = do.obs*2, x = datetime), col = "purple", size=1.5) + 
  geom_line(aes(y = (dischargeCFS*2), x = datetime), col = "#46b8b8",alpha = .7, size = 0.75) +
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

# saveRDS(BWL_datQ, "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/25_BWL_DO_flag_flow.rds")


##==========================
### BWU
##==========================

##==========================
## Read in DO data from BWL
##===========================
## BWU
BWU_DO <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/25_BWU_DO_flag_record.rds")

NS_plot_DO <- plot_ly(data = BWU_DO, x = ~datetime, y = ~Dissolved_Oxygen_offset, type = 'scatter',  mode = 'markers',  color = ~Concat_Date) %>%
  layout(xaxis = list(title = "Date and Time"),
         yaxis = list(title = "DO (mgL)"))

##==========================
## Flag sediment burial 
##==========================
BWU_DO <- BWU_DO %>%
  mutate(fouling_flag = if_else(
    (datetime > as.POSIXct("2022-09-15 00:00:00") & datetime < as.POSIXct("2021-10-04 00:00:00")) | 
      (datetime > as.POSIXct("2022-11-21 00:00:00") & datetime < as.POSIXct("2023-01-03 00:00:00")),
    1, 0))

temp_plot <- BWU_DO %>% filter(fouling_flag<1) %>%
  ggplot(aes(x = datetime, y = Dissolved_Oxygen_offset, colour = as.factor(maintenance_flag))) +
  geom_point(size=1, alpha=0.75) + theme_bw()

##==========================
## Read in morphology observations
##==========================
morph_dat_BWU <- read.csv("R:/Users/kloria/Documents/Stream_Metab_24/stream_morphology_core.csv") %>%
  mutate(datetime = as.POSIXct(datetime, format="%Y-%m-%dT%H:%M:%SZ"),
         date= as.Date(datetime)) %>% filter(Site=="BWU")

## Join with flow 
HRflow_BWU <- HRflow_BW %>%
  left_join(morph_dat_BWU, by = c("datetime"))

HRflow_BWU_plot <- morph_dat_BWU %>%
  left_join(HRflow_BW, by = c("datetime"))

plot(HRflow_BWU$dischargeCMS, HRflow_BWU$Discharge_cms)
plot(HRflow_BWU$stage_m, HRflow_BWU$z)

## quick linear models 
flow_lm <- lm(Discharge_cms ~dischargeCMS, data=HRflow_BWU)
summary(flow_lm)

## flow plot
flow_summary <- summary(flow_lm)
slope <- round(flow_summary$coefficients[2, 1], 3) 
r_squared <- round(flow_summary$r.squared, 2) 
p_value <- round(flow_summary$coefficients[2, 4], 3)  

## Dynamic Adjustment Based on Slope
HRflow_BWU$adjusted_dischargeCMS <- c(HRflow_BWU$dischargeCMS * slope)
hist(HRflow_BWU$adjusted_dischargeCMS)
hist(HRflow_BWU$dischargeCMS)

BWUflow_plot <- ggplot(HRflow_BWU_plot, aes(x = dischargeCMS, y = Discharge_cms, color=method)) +
  geom_point(alpha=0.75, size =2) + geom_smooth(method="lm", se=F, color="grey50") +
  labs(title = "BWU", x = "USGS flow (cms)", y = "Observed flow (cms)") + theme_bw() +
  annotate(
    "text",
    x = Inf, y = Inf, 
    label = paste("Slope:", slope, "\nR²:", r_squared, "\np-value:", p_value),
    hjust = 1.1, vjust = 4.5,
    size = 3,
    color = "black"
  )
BWUflow_plot

stage_lm <- lm(z ~stage_m, data=HRflow_BWU)
summary(stage_lm)

## depth plot
depth_summary <- summary(stage_lm)
slope <- round(depth_summary$coefficients[2, 1], 3) 
r_squared <- round(depth_summary$r.squared, 2) 
p_value <- round(depth_summary$coefficients[2, 4], 4)  

## Dynamic Adjustment Based on Slope
HRflow_BWU$adjusted_depth <- c(HRflow_BWU$stage_m * slope)
hist(HRflow_BWU$stage_m)
hist(HRflow_BWU$adjusted_depth)

BWUdepth_plot <- ggplot(HRflow_BWU_plot, aes(x = stage_m, y = z, color=method)) +
  geom_point(alpha=0.75, size =2) + geom_smooth(method="lm", se=F, color="grey50") +
  labs(title = "BWU", x = "USGS stage (m)", y = "Obsersed depth (m)") + theme_bw() +
  annotate(
    "text",
    x = Inf, y = Inf, 
    label = paste("Slope:", slope, "\nR²:", r_squared, "\np-value:", p_value),
    hjust = 1.1, vjust = 4.5,
    size = 3,
    color = "black")

BWUdepth_plot

#### 
names(HRflow_BWU)

HRflow_BWU <- HRflow_BWU %>%
  dplyr::select(datetime, adjusted_dischargeCMS, adjusted_depth, dischargeCFS, 
                dischargeCMS,
                stageF, stage_m, z, wtr_USGS, stage_m, w, v, v_estimation)

### merge DO and flow data
BWU_dat <- datetime_df %>%
  filter(datetime > as.POSIXct("2021-06-29 00:0:00") & 
           datetime < as.POSIXct("2023-10-26 00:0:00"))%>%
  left_join(BWU_DO, by = c("datetime")) %>%
  left_join(HRflow_BWU, by = c("datetime"))

range((BWU_dat$date))
str(BWU_dat)

BWU_datQ <- as.data.frame(BWU_dat)

##==========================
## Double check the timestamps of all temp data
##==========================
temp_data_plot <- BWU_datQ %>%
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
  geom_line(aes(y = wtr, x = datetime), col = "goldenrod", size=1.5) +
  geom_point(aes(y = do.obs*2, x = datetime), col = "purple", size=1.5) + 
  geom_line(aes(y = (dischargeCFS*2), x = datetime), col = "#46b8b8",alpha = .7, size = 0.75) +
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

# saveRDS(BWU_datQ, "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/25_BWU_DO_flag_flow.rds")













# OLD CODE BELOW
# # ##===========================
# 
# ##===========================
# ## Average data at 15 minutes 
# str(BWL_datQ)
# 
# BWL_15m <- BWL_datQ %>%
#   # Step 1: Create a column for rounding to the nearest 15 minutes
#   mutate(
#     rounded_datetime_15min = as.POSIXct(
#       round(as.numeric(datetime) / (15 * 60)) * (15 * 60),
#       origin = "1970-01-01",
#       tz = "America/Los_Angeles"
#     )
#   ) %>%
#   # Step 2: Group by 15-minute intervals and compute averages
#   group_by(rounded_datetime_15min) %>%
#   summarize(
#     avg_DO.obs = mean(do.obs, na.rm = TRUE),
#     avg_temp.water = mean(wtr, na.rm = TRUE),
#     avg_depth = mean(depth, na.rm = TRUE),
#     avg_discharge = mean(dischargeCMS, na.rm = TRUE),
#     avg_wtr_USGS = mean(wtr_USGS, na.rm =TRUE),
#     .groups = "drop"
#   ) %>%
#   # Step 3: Merge 15-minute averages back into the original data
#   right_join(
#     BWL_datQ %>% 
#       mutate(
#         rounded_datetime_15min = as.POSIXct(
#           round(as.numeric(datetime) / (15 * 60)) * (15 * 60),
#           origin = "1970-01-01",
#           tz = "America/Los_Angeles"
#         )
#       ),
#     by = "rounded_datetime_15min"
#   ) %>%
#   arrange(datetime) %>%
#   # Step 4: Ensure all 15-minute intervals are present
#   mutate(
#     date_only = as.Date(datetime)
#   ) %>%
#   group_by(date_only) %>%
#   complete(rounded_datetime_15min = seq(min(rounded_datetime_15min), max(rounded_datetime_15min), by = "15 min")) %>%
#   ungroup() %>%
#   # Step 5: Fill missing 15-minute values with computed averages
#   mutate(
#     DO.obs = ifelse(is.na(do.obs), avg_DO.obs, do.obs),
#     temp.water = ifelse(is.na(wtr), avg_temp.water, wtr),
#     wtr_USGS = ifelse(is.na(wtr_USGS), avg_wtr_USGS, wtr_USGS),
#     depth = ifelse(is.na(depth), avg_depth, depth),
#     discharge = ifelse(is.na(dischargeCMS), avg_discharge, dischargeCMS)
#     ) %>%
#   # Step 6: Drop intermediate columns
#   select(-c(avg_DO.obs, avg_temp.water, avg_depth, avg_discharge, date_only))
# 
# # Display column names and summary of the processed data
# names(BWL_15m)
# summary(BWL_15m)
# 


# 
# 
# ##==========================
# ## Infill observations btwn the 15 minutes
# BWL_DOQ1<- BWL_DOQ %>%
#   fill(dischargeCFS,.direction = "down")%>%
#   fill(dischargeCMS,.direction = "down")%>% 
#   fill(scale_Q,.direction = "down")%>% 
#   fill(depth,.direction = "down")%>%
#   dplyr::ungroup()
# 
# 
# 
# str(BWL_DOQ)
# 
# BWL_DOQ2 <- BWL_DOQ1 %>%
#   fill(dischargeCFS,.direction = "up")%>%
#   fill(dischargeCMS,.direction = "up")%>% 
#   fill(scale_Q,.direction = "up")%>% 
#   fill(depth,.direction = "up")%>%
#   dplyr::ungroup()
# 
# 
# ##==========================
# ##  Summarize daily discharge statistics
# BWL_sum <- BWL_DOQ2 %>%
#   filter(fouling_flag<1 & maintenance_flag<1) %>%
#   dplyr::group_by(date) %>%
#   summarise(
#     day_mean_flow = mean(dischargeCMS, na.rm = TRUE),
#     day_mean_temp = mean(wtr, na.rm = TRUE))
# 
# Flow_sum <- BWL_DOQ2 %>%
#   summarise(
#     day_mean_flow = mean(dischargeCMS, na.rm = TRUE))
# 
# BWL_DO_Q <- BWL_DOQ2 %>%
#   left_join(BWL_sum, by = "date") %>%
#   mutate(
#     discharge_D50_flag = if_else(dischargeCMS > c(1.5 * day_mean_flow), 1, 0),
#     discharge_D100_flag = if_else(dischargeCMS > c(2 * day_mean_flow), 1, 0),
#     discharge_T25_flag = if_else(dischargeCMS > c(1.25 * 1.065283), 1, 0),
#     discharge_T10_flag = if_else(dischargeCMS > c(1.1 * 1.065283), 1, 0),
#     discharge_T50_flag = if_else(dischargeCMS > c(1.5 * 1.065283), 1, 0),
#     discharge_T100_flag = if_else(dischargeCMS > c(2 * 1.065283), 1, 0),
#     discharge_T150_flag = if_else(dischargeCMS > c(2.5 * 1.065283), 1, 0),
#     wtr_D50_flag = if_else(wtr < c(0.5 * day_mean_temp), 1, 0),
#     wtr_D100_flag = if_else(wtr < c(2 * day_mean_temp), 1, 0))
# 
# temp_plot <- BWL_DO_Q %>% filter(datetime > as.POSIXct("2023-09-02 00:00:00")) %>%
#   ggplot(aes(x = datetime, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_D100_flag))) +
#   geom_point(size=1, alpha=0.75) + theme_bw()
# 
# temp_plot <- BWL_DO_Q %>% filter(datetime > as.POSIXct("2023-09-02 00:00:00")) %>%
#   ggplot(aes(x = dischargeCMS, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_D100_flag))) +
#   geom_point(size=1, alpha=0.75) + theme_bw()
# 
# temp_plot <- BWL_DO_Q %>% filter(datetime > as.POSIXct("2023-09-02 00:00:00")) %>%
#   ggplot(aes(x = dischargeCMS, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_T25_flag))) +
#   geom_point(size=1, alpha=0.75) + theme_bw()
# 
# temp_plot <- BWL_DO_Q %>% filter(datetime > as.POSIXct("2023-09-02 00:00:00")) %>%
#   ggplot(aes(x = dischargeCMS, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_T100_flag))) +
#   geom_point(size=1, alpha=0.75) + theme_bw()
# 
# temp_plot <- BWL_DO_Q %>% filter(datetime > as.POSIXct("2023-09-02 00:00:00")) %>%
#   ggplot(aes(x = dischargeCMS, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_T150_flag))) +
#   geom_point(size=1, alpha=0.75) + theme_bw()
# 
# ### save data:
# # saveRDS(BWL_DO_Q, "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/24_BWL_DO_flag_recordv2.rds")
# 
# 
# ## BWU 
# BWU_DO <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/24_BWU_DO_flag_record.rds")
# 
# NS_plot_DO <- plot_ly(data = BWU_DO, x = ~datetime, y = ~Dissolved_Oxygen_offset, type = 'scatter',  mode = 'markers',  color = ~Concat_Date) %>%
#   layout(xaxis = list(title = "Date and Time"),
#          yaxis = list(title = "DO (mgL)"))
# 
# ##==========================
# ## Flag sediment burial 
# BWU_DO <- BWU_DO %>%
#   mutate(fouling_flag = if_else(
#     (datetime > as.POSIXct("2021-08-01 00:00:00") & datetime < as.POSIXct("2021-08-12 00:00:00")) | 
#       (datetime > as.POSIXct("2022-07-08 00:00:00") & datetime < as.POSIXct("2022-07-09 00:00:00")) | 
#       (datetime > as.POSIXct("2022-09-11 00:00:00") & datetime < as.POSIXct("2022-10-13 00:00:00")) |
#       (datetime > as.POSIXct("2022-11-21 00:00:00") & datetime < as.POSIXct("2023-01-29 00:00:00"))|
#       (datetime > as.POSIXct("2023-04-25 00:00:00") & datetime < as.POSIXct("2023-04-27 00:00:00")), 
#     1, 0))
# 
# temp_plot <- BWU_DO %>% filter(fouling_flag<1) %>%
#   ggplot(aes(x = datetime, y = Dissolved_Oxygen_offset, colour = as.factor(maintenance_flag))) +
#   geom_point(size=1, alpha=0.75) + theme_bw()
# 
# morph_dat <- read.csv("R:/Users/kloria/Documents/Stream_Metab_24/DO_calibration_dat/stream_morphology.csv") %>%
#   mutate(
#     datetime = as.POSIXct(Date, format="%Y-%m-%dT%H:%M:%SZ"),
#     date= as.Date(datetime)) %>%
#   filter(Site=="BWU")
# 
# HRflow_BW <- HRflow_BW %>%
#   left_join(morph_dat, by = c("datetime"))
# 
# ## effective mean depth and stage
# HRflow_na <- HRflow_BW %>%
#   filter(stage_m < 2.9)
# 
# plot(HRflow_na$z, HRflow_na$stage_m)
# 
# stage_lm <- lm(stage_m ~ z, data=HRflow_na)
# summary(stage_lm)
# 
# 
# HRflow_na$depth <- c(HRflow_na$stage_m * c(0.1111111))
# plot(HRflow_na$z, HRflow_na$depth)
# 
# HRflow_BW$depth <- c(HRflow_BW$stage_m * c(0.1111111))
# HRflow_BW$dischargeCMS <- c(HRflow_BW$dischargeCMS/(2))
# 
# 
# HRflow_BW1 <- HRflow_BW %>%
#   dplyr::select(datetime, dischargeCFS, stageF, dischargeCMS, scale_Q, stage_m,
#                 depth)
# 
# ### merge DO and flow data
# BWU_DO1 <- BWU_DO%>%
#   left_join(HRflow_BW1, by = c("datetime"))
# str(BWU_DO1)
# 
# 
# ## Infill observations btwn the 15 minutes
# BWU_DOQ1<- BWU_DO1 %>%
#   fill(dischargeCFS,.direction = "down")%>%
#   fill(dischargeCMS,.direction = "down")%>% 
#   fill(scale_Q,.direction = "down")%>% 
#   fill(depth,.direction = "down")%>%
#   dplyr::ungroup()
# 
# str(BWU_DOQ1)
# 
# BWU_DOQ2<- BWU_DOQ1 %>%
#   fill(dischargeCFS,.direction = "up")%>%
#   fill(dischargeCMS,.direction = "up")%>% 
#   fill(scale_Q,.direction = "up")%>% 
#   fill(depth,.direction = "up")%>%
#   dplyr::ungroup()
# 
# summary(BWU_DOQ2)
# 
# ##==========================
# ##  Summarize daily discharge statistics
# BWU_sum <- BWU_DOQ2 %>%
#   filter(fouling_flag<1 & maintenance_flag<1) %>%
#   dplyr::group_by(date) %>%
#   summarise(
#     day_mean_flow = mean(dischargeCMS, na.rm = T),
#     day_mean_temp = mean(wtr, na.rm = TRUE))
# 
# BWU_sumT <- BWU_DOQ2 %>%
#   filter(fouling_flag<1 & maintenance_flag<1) %>%
#   summarise(
#     mean_flow = mean(dischargeCMS, na.rm = TRUE),
#     mean_temp = mean(wtr, na.rm = TRUE))
# 
# BWU_DO_Q <- BWU_DOQ2 %>%
#   left_join(BWU_sum, by = "date") %>%
#   mutate(
#     discharge_D50_flag = if_else(dischargeCMS > c(1.5 * day_mean_flow), 1, 0),
#     discharge_D100_flag = if_else(dischargeCMS > c(2 * day_mean_flow), 1, 0),
#     discharge_T25_flag = if_else(dischargeCMS > c(1.25 * 0.779), 1, 0),
#     discharge_T10_flag = if_else(dischargeCMS > c(1.1 * 0.779), 1, 0),
#     discharge_T50_flag = if_else(dischargeCMS > c(1.5 * 0.779), 1, 0),
#     discharge_T100_flag = if_else(dischargeCMS > c(2 * 0.779), 1, 0),
#     discharge_T150_flag = if_else(dischargeCMS > c(2.5 * 0.779), 1, 0))
# 
# temp_plot <- BWU_DO_Q %>% filter(fouling_flag<1 & maintenance_flag<1) %>%
#   ggplot(aes(x = datetime, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_T50_flag))) +
#   geom_point(size=1, alpha=0.75) + theme_bw()
# 
# 
# temp_plot <- BWU_DO_Q %>% filter(fouling_flag<1 & maintenance_flag<1) %>%
#   ggplot(aes(x = dischargeCMS, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_50_flag))) +
#   geom_point(size=1, alpha=0.75) + theme_bw()

# saveRDS(BWU_DO_Q, "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/24_BWU_DO_flag_recordv3.rds")


