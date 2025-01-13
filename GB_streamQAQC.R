library(dataRetrieval)
library(lubridate)
library(ggplot2)
library(dplyr)
library(tidyverse)

##==========================
## Read in DO data from GBL
##===========================
## GBL
GBL_DO <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/24_GBL_DO_flag_record.rds")

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

# Set the start date to "1980-01-01" and end date to today (current date)
start.date <- "2021-03-10"
end.date <- "2024-10-10"  # Use
## Download hourly flow data for both sites separately and combine
HRflow_data_GB <- readNWISuv(siteNumbers = siteNo_GB, parameterCd = c("00060", "00065"), startDate = start.date, endDate = end.date) 

HRflow_GB <- HRflow_data_GB %>%
  select(datetime = "dateTime", dischargeCFS = "X_00060_00000", stageF = "X_00065_00000") %>%
  mutate(
    dischargeCMS= c(dischargeCFS*0.0283168),
    scale_Q= c((dischargeCFS*0.0283168)/10.64485),
    stage_m= c(stageF*0.3048)
    ) 


morph_dat <- read.csv("R:/Users/kloria/Documents/Stream_Metab_24/DO_calibration_dat/stream_morphology.csv") %>%
  mutate(
    datetime = as.POSIXct(Date, format="%Y-%m-%dT%H:%M:%SZ"),
    date= as.Date(datetime)) %>%
  filter(Site=="GBL")

HRflow_GB <- HRflow_GB %>%
  left_join(morph_dat, by = c("datetime"))
 
## effective mean depth and stage
HRflow_na <- HRflow_GB 
plot(HRflow_na$z, HRflow_na$stage_m)

stage_lm <- lm(z ~ stage_m, data=HRflow_na)
summary(stage_lm)

plot(HRflow_na$z, HRflow_na$stage_m)

HRflow_na$depth <- c(HRflow_na$stage_m * 0.09446)
plot(HRflow_na$z, HRflow_na$depth)

HRflow_GB$depth <- c(HRflow_GB$stage_m * 0.09446)


#0.43423

hist(HRflow_GB$z)
hist(HRflow_GB$stage_m)
hist(HRflow_GB$depth)



HRflow_GB <- HRflow_GB%>%
  dplyr::select(datetime, dischargeCFS, stageF, dischargeCMS, scale_Q, stage_m,
                #depth,
                #Site, 
                Discharge, Temp, z, w, v)

### merge DO and flow data
GBL_DOQ <- GBL_DOQ%>%
  left_join(HRflow_GB, by = c("datetime"))
str(GBL_DOQ)
range((GBL_DOQ$date))

GBL_DOQ <- as.data.frame(GBL_DOQ)
str(GBL_DOQ)

##==========================
## Infill observations btwn the 15 minutes
GBL_DOQ1<- GBL_DOQ %>%
  fill(dischargeCFS,.direction = "down")%>%
  fill(dischargeCMS,.direction = "down")%>% 
  fill(scale_Q,.direction = "down")%>% 
  fill(depth,.direction = "down")%>%
  dplyr::ungroup()

str(GBL_DOQ1)

GBL_DOQ2<- GBL_DOQ1 %>%
  fill(dischargeCFS,.direction = "up")%>%
  fill(dischargeCMS,.direction = "up")%>% 
  fill(scale_Q,.direction = "up")%>% 
  fill(depth,.direction = "up")%>%
  dplyr::ungroup()

summary(GBL_DOQ2)

##==========================
##  Summarize daily discharge statistics
GBL_sum <- GBL_DOQ2 %>%
  filter(fouling_flag<1 & maintenance_flag<1) %>%
  dplyr::group_by(date) %>%
  summarise(
    day_mean_flow = mean(dischargeCMS, na.rm = TRUE),
    day_mean_temp = mean(wtr, na.rm = TRUE))

Flow_sum <- HRflow_GB %>%
  summarise(
    day_mean_flow = mean(dischargeCMS, na.rm = TRUE))
    
GBL_DO_Q <- GBL_DOQ2 %>%
  left_join(GBL_sum, by = "date") %>%
  mutate(
    discharge_D50_flag = if_else(dischargeCMS > c(1.5 * day_mean_flow), 1, 0),
    discharge_D100_flag = if_else(dischargeCMS > c(2 * day_mean_flow), 1, 0),
    discharge_T25_flag = if_else(dischargeCMS > c(1.25 * 0.07115973), 1, 0),
    discharge_T50_flag = if_else(dischargeCMS > c(1.5 * 0.07115973), 1, 0),
    discharge_T100_flag = if_else(dischargeCMS > c(2 * 0.07115973), 1, 0),
    wtr_D50_flag = if_else(wtr < c(0.5 * day_mean_temp), 1, 0),
    wtr_D100_flag = if_else(wtr < c(2 * day_mean_temp), 1, 0))

temp_plot <- GBL_DO_Q %>% filter(datetime > as.POSIXct("2023-09-02 00:00:00")) %>%
  ggplot(aes(x = datetime, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_D100_flag))) +
  geom_point(size=1, alpha=0.75) + theme_bw()

temp_plot <- GBL_DO_Q %>% filter(datetime > as.POSIXct("2023-09-02 00:00:00")) %>%
  ggplot(aes(x = dischargeCMS, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_D100_flag))) +
  geom_point(size=1, alpha=0.75) + theme_bw()

temp_plot <- GBL_DO_Q %>% filter(datetime > as.POSIXct("2023-09-02 00:00:00")) %>%
  ggplot(aes(x = dischargeCMS, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_T25_flag))) +
  geom_point(size=1, alpha=0.75) + theme_bw()

temp_plot <- GBL_DO_Q %>% filter(datetime > as.POSIXct("2023-09-02 00:00:00")) %>%
  ggplot(aes(x = dischargeCMS, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_T50_flag))) +
  geom_point(size=1, alpha=0.75) + theme_bw()

temp_plot <- GBL_DO_Q %>% filter(datetime > as.POSIXct("2023-09-02 00:00:00")) %>%
  ggplot(aes(x = dischargeCMS, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_T50_flag))) +
  geom_point(size=1, alpha=0.75) + theme_bw()

temp_plot <- GBL_DO_Q %>% filter(datetime > as.POSIXct("2023-09-02 00:00:00")) %>%
  ggplot(aes(x = dischargeCMS, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_T100_flag))) +
  geom_point(size=1, alpha=0.75) + theme_bw()

temp_plot <- GBL_DO_Q %>% filter(datetime > as.POSIXct("2023-09-02 00:00:00")) %>%
  ggplot(aes(x = wtr, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_T100_flag))) +
  geom_point(size=1, alpha=0.75) + theme_bw()

temp_plot <- GBL_DO_Q %>% 
  filter(wtr_D100_flag<2) %>%
  ggplot(aes(x = wtr, y = Dissolved_Oxygen_offset, colour = as.factor(wtr_D50_flag))) +
  geom_point(size=1, alpha=0.75) + theme_bw()

hist(GBL_DO_Q$dischargeCMS)


### save data:
saveRDS(GBL_DO_Q, "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/24_GBL_DO_flag_record.rds")


####
# 1. bring in Baro data
# 3. Phil's light model 

################
################

## GBU 
GBU_DO <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/24_GBU_DO_flag_record.rds")

NS_plot_DO <- plot_ly(data = GBU_DO, x = ~datetime, y = ~Dissolved_Oxygen_offset, type = 'scatter',  mode = 'markers',  color = ~Concat_Date) %>%
  layout(xaxis = list(title = "Date and Time"),
         yaxis = list(title = "DO (mgL)"))
##==========================
## Flag sediment burial 
GBU_DOQ <- GBU_DO %>%
  mutate(fouling_flag = if_else(
    (datetime > as.POSIXct("2022-09-11 00:00:00") & datetime < as.POSIXct("2022-09-16 00:00:00")) | 
      (datetime > as.POSIXct("2022-09-23 00:00:00") & datetime < as.POSIXct("2022-10-04 00:00:00")), 
    1, 0))

temp_plot <- GBL_DOQ %>% filter(fouling_flag<1) %>%
  ggplot(aes(x = datetime, y = Dissolved_Oxygen_offset, colour = as.factor(maintenance_flag))) +
  geom_point(size=1, alpha=0.75) + theme_bw()


morph_datU <- read.csv("R:/Users/kloria/Documents/Stream_Metab_24/DO_calibration_dat/stream_morphology.csv") %>%
  mutate(
    datetime = as.POSIXct(Date, format="%Y-%m-%dT%H:%M:%SZ"),
    date= as.Date(datetime)) %>%
  filter(Site=="GBU")


GBLQ <- c(0.098, 0.152, NA, 0.071, NA)
GBU <- c(0.098, 0.111, 0.060, 0.064, NA)


# Create a data frame with these values
GB_flow <- data.frame(GBLQ, GBU)

# Fit a linear model with GBU as the dependent variable and GBLQ as the independent variable
lm_Q <- lm(GBU ~ GBLQ, data = GB_flow)

# View the summary to get the slope (coefficient of GBLQ)
summary(lm_Q)

# Extract the slope from the model
slope <- coef(lm_Q)["GBLQ"]

# 0.531746 

# View the corrected data
print(GB_flow)


HRflow_GBU <- HRflow_GB %>%
  left_join(morph_datU, by = c("datetime"))

### merge DO and flow data # 0.7963961
HRflow_GBU1 <- HRflow_GBU %>%
  mutate(
    dischargeCFS = c(dischargeCFS * 0.7963961), 
    dischargeCMS = c(dischargeCMS*0.7963961),
    scale_Q = c(scale_Q*0.7963961), 
    stage_m =  c(stage_m*0.531746))

HRflow_GBU2 <- HRflow_GBU1%>%
  dplyr::select(datetime, dischargeCFS, stageF, dischargeCMS, scale_Q, stage_m,
                Site, Discharge, Temp, z, w, v)

GBU_DOQ1 <- GBU_DOQ %>%
  left_join(HRflow_GBU2,by = c("datetime"))



str(GBU_DOQ1)
range((GBU_DOQ$date))
hist(GBU_DOQ$dischargeCMS)


# dept relationship 
GBU_DOQ1$depth <- c(GBU_DOQ1$stage_m * 0.1437)


##


## Infill observations btwn the 15 minutes
GBU_DOQ1<- GBU_DOQ1 %>%
  fill(dischargeCFS,.direction = "down")%>%
  fill(dischargeCMS,.direction = "down")%>% 
  fill(scale_Q,.direction = "down")%>% 
  fill(depth,.direction = "down")%>%
  dplyr::ungroup()

str(GBL_DOQ1)

GBU_DOQ2<- GBU_DOQ1 %>%
  fill(dischargeCFS,.direction = "up")%>%
  fill(dischargeCMS,.direction = "up")%>% 
  fill(scale_Q,.direction = "up")%>% 
  fill(depth,.direction = "up")%>%
  dplyr::ungroup()

summary(GBU_DOQ2)


##==========================
##  Summarize daily discharge statistics
GBU_sum <- GBU_DOQ2 %>%
  #filter(fouling_flag<1 & maintenance_flag<1) %>%
  dplyr::group_by(date) %>%
  summarise(
    day_mean_flow = mean(dischargeCMS, na.rm = T),
    day_mean_temp = mean(wtr, na.rm = TRUE))

GBU_sumT <- GBU_DOQ2 %>%
  filter(fouling_flag<1 & maintenance_flag<1) %>%
  summarise(
    mean_flow = mean(dischargeCMS, na.rm = TRUE),
    mean_temp = mean(wtr, na.rm = TRUE))

GBU_DO_Q <- GBU_DOQ2 %>%
  left_join(GBU_sum, by = "date") %>%
  mutate(discharge_D50_flag = if_else(dischargeCMS > c(1.5 * day_mean_flow), 1, 0), 
         discharge_D100_flag = if_else(dischargeCMS > c(2 * day_mean_flow), 1, 0),
         discharge_T25_flag = if_else(dischargeCMS > c(1.25 * 0.0482), 1, 0),
         discharge_T50_flag = if_else(dischargeCMS > c(1.5 * 0.0482), 1, 0),
         discharge_T100_flag = if_else(dischargeCMS > c(2 * 0.0482), 1, 0),
         wtr_75_flag = if_else(wtr > c(1.75 * day_mean_temp), 1, 0))

temp_plot <- GBU_DO_Q %>% filter(fouling_flag<1 & maintenance_flag<1) %>%
  ggplot(aes(x = datetime, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_T50_flag))) +
  geom_point(size=1, alpha=0.75) + theme_bw()


temp_plot <- GBU_DO_Q %>% filter(fouling_flag<1 & maintenance_flag<1) %>%
  ggplot(aes(x = dischargeCMS, y = Dissolved_Oxygen_offset, colour = as.factor(discharge_50_flag))) +
  geom_point(size=1, alpha=0.75) + theme_bw()

# saveRDS(GBU_DO_Q, "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/24_GBU_DO_flag_recordv2.rds")

