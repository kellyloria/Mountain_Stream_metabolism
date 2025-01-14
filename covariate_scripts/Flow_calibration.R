
library(dataRetrieval)
library(lubridate)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)

## flow corrections

##==========================
## Read USGS flow data BWL
##===========================
## Bring in Streamflow
siteNo_BW <- "10336660"
# Define the parameter codes for flow and stage
pCode_flow <- "00060"
pCode_stage <- "00065"
# Set the start date to "1980-01-01" and end date to today (current date)
start.date <- "2020-09-26"
end.date <- "2024-08-18"  # Use
## Download hourly flow data for both sites separately and combine
HRflow_data_BW <- readNWISuv(siteNumbers = siteNo_BW, parameterCd = c("00060", "00065"), startDate = start.date, endDate = end.date) 

HRflow_BW <- HRflow_data_BW %>%
  select(datetime = "dateTime", dischargeCFS = "X_00060_00000", stageF = "X_00065_00000") %>%
  mutate(dischargeCMS= c(dischargeCFS*0.0283168), # meter conversion
    scale_Q= c((dischargeCFS*0.0283168)/29.00787),  # scale for km2
    stage_m= c(stageF*0.3048))  # meter conversion 

##==========================
## Read in morphology observations
morph_dat_BWL <- read.csv("R:/Users/kloria/Documents/Stream_Metab_24/stream_morphology_core.csv") %>%
  mutate(datetime = as.POSIXct(datetime, format="%Y-%m-%dT%H:%M:%SZ"),
    date= as.Date(datetime)) %>% filter(Site=="BWL")

## Join with flow 
HRflow_BWL <- HRflow_BW %>%
  left_join(morph_dat_BWL, by = c("datetime"))

HRflow_BWL_plot <- morph_dat_BWL %>%
  left_join(HRflow_BW, by = c("datetime"))

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


##==========================
### BWU
##==========================
## Read in morphology observations
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


##==========================
## Read USGS flow data BWL
##===========================
## Bring in Streamflow
siteNo_GB <- "10336730"
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
  mutate(dischargeCMS= c(dischargeCFS*0.0283168),
    scale_Q= c((dischargeCFS*0.0283168)/10.64485),
    stage_m= c(stageF*0.3048))

##==========================
## Read in morphology observations
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


##==========================
## Read in morphology observations
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


flow_grid <- ggarrange(BWLflow_plot,
                      BWUflow_plot, 
                      GBLflow_plot,
                      GBUflow_plot,
                      ncol = 2, nrow = 2,
                      common.legend = TRUE, 
                      legend = "bottom")

# ggsave("R:/Users/kloria/Documents/Stream_Metab_24/figures/Core_sites_flowcorrections.png", plot = flow_grid, width = 10, height = 6, units = "in")


depth_grid <- ggarrange(BWLdepth_plot,
                        BWUdepth_plot, 
                       GBLdepth_plot,
                       GBUdepth_plot,
                       ncol = 2, nrow = 2,
                       common.legend = TRUE, 
                       legend = "bottom")
# ggsave("R:/Users/kloria/Documents/Stream_Metab_24/figures/Core_sites_depthorrections.png", plot = depth_grid, width = 10, height = 6, units = "in")


## temp end.