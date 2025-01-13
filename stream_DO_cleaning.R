##===============================================================================
## Created  10/31/2024 by KAL
##===============================================================================

## KAL's temporary path reminders: 
## setwd("/Users/kellyloria/Documents/LittoralMetabModeling")
## PC: setwd("R:/Users/kloria/Documents/Stream_Metab_24")

# install.packages("remotes")
# remotes::install_github("nrlottig/nrlmetab")

lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate", "zoo",
         "tidyverse","data.table","xts","dygraphs",
         "nrlmetab"), require, character.only=T)
##==========================
## Read in DO data
##===========================
# 1. Read in raw aggregated DO data
files <- list.files("./Core_sites/raw_merged_dat", pattern = "\\.rds$", full.names = TRUE)
rawdatDO <- do.call(rbind, lapply(files, readRDS))
str(rawdatDO)

unique(rawdatDO$Sensor)

####
####

##===========================
## 1. OFF SET by CALIBRATION 
##===========================
### Stream miniDOT calibration (ref senosr Cat_7450-666671)
ms_DOcal <-read.csv("./DO_calibration_dat/DO_cal_20241010.csv")
unique(ms_DOcal$serial)

###
rawdatDO_offset <- rawdatDO %>%
  mutate(offset = case_when( 
    Sensor == "7450-617000" ~ -0.0484,
    Sensor == "7450-686243" ~ -0.1173, #
    Sensor == "7450-162475" ~ -0.0484,
    Sensor == "7450-714094" ~ -0.0484,
    Sensor == "7450-666671" ~ -0.0001, #
    Sensor == "7450-099447" ~ -0.0484,
    Sensor == "7450-559438" ~ 0.0649, #
    Sensor == "7450-547404" ~ 0.1197,
    Sensor == "7450-559438" ~ 0.0484,
    Sensor == "7450-411027" ~ -0.2416, #
    TRUE ~ (0.0484))) # average offset value if the sensor was not found for calibration.

DO_offset <- rawdatDO_offset %>%
  mutate(Dissolved_Oxygen_offset = (do.obs + offset),
         date = as.Date(datetime)) 

hist(DO_offset$Dissolved_Oxygen_offset)

# saveRDS(DO_offset, file = "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/24_offset_DO_record.rds")

##===========================
## 2. Add in cleaning date obs
clean_dat <- read.csv("./DO_calibration_dat/Sensor_cleaning_record.csv")%>%
  mutate(date= as.Date(date, format("%Y-%m-%d")),
  datetime = as.POSIXct(paste(date, "11:00:00"), format = "%Y-%m-%d %H:%M:%S", tz = "America/Los_Angeles"))%>%
  select("site", "datetime", "cleaned", "plate_add")
str(clean_dat)

DO_offset <- DO_offset%>%
  left_join(clean_dat, by = c("site","datetime"))
str(DO_offset)

# saveRDS(DO_offset1, file = "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/24_offset_DO_record.rds")


##===========================
# 2. Create the plot using plot_ly


## BWL
NS_plot_BWL <- plot_ly(data = DO_offset_out%>%filter(site=="BWL"), x = ~datetime, y = ~Dissolved_Oxygen_offset, type = 'scatter',  mode = 'markers',  color = ~Sensor) %>%
  layout(xaxis = list(title = "Date and Time"),
         yaxis = list(title = "DO (mgL)"))


filtered_BWL <- DO_offset %>%
  filter(site == "BWL")  %>%
  filter(datetime > as.POSIXct("2020-09-27 9:00:00"))%>%
  filter(datetime < as.POSIXct("2021-04-29 00:00:00") | datetime > as.POSIXct("2021-04-29 12:00:00"))%>%
  filter(datetime < as.POSIXct("2021-05-12 04:30:00") | datetime > as.POSIXct("2021-05-12 08:00:00")) %>%
  filter(datetime < as.POSIXct("2021-05-26 06:00:00") | datetime > as.POSIXct("2021-05-26 07:30:00")) %>%
  filter(datetime < as.POSIXct("2021-06-02 06:30:00") | datetime > as.POSIXct("2021-06-02 07:50:00")) %>%
  filter(datetime < as.POSIXct("2021-11-12 04:30:00") | datetime > as.POSIXct("2021-11-12 10:30:00")) %>%
  filter(datetime < as.POSIXct("2023-03-07 00:30:00") | datetime > as.POSIXct("2023-03-09 00:00:00")) %>%
  filter(datetime < as.POSIXct("2023-03-07 00:30:00") | datetime > as.POSIXct("2023-03-09 00:00:00")) %>%
  filter(datetime < as.POSIXct("2023-04-05 00:30:00") | datetime > as.POSIXct("2023-04-06 00:00:00")) %>%
  filter(datetime < as.POSIXct("2023-05-07 00:00:00") | datetime > as.POSIXct("2023-05-07 02:00:00")) %>%
  filter(datetime < as.POSIXct("2023-05-19 00:00:00") | datetime > as.POSIXct("2023-05-19 02:00:00")) %>%
  filter(datetime < as.POSIXct("2023-06-12 03:00:00") | datetime > as.POSIXct("2023-06-12 05:00:00")) %>%
  filter(datetime < as.POSIXct("2023-07-12 03:00:00") | datetime > as.POSIXct("2023-07-12 05:00:00")) %>%
  filter(datetime < as.POSIXct("2023-08-23 03:00:00") | datetime > as.POSIXct("2023-08-23 09:00:00")) %>%
  filter(datetime < as.POSIXct("2024-06-26 03:00:00") | datetime > as.POSIXct("2024-06-26 12:00:00")) %>%
  filter(datetime < as.POSIXct("2024-06-20 00:00:00") | datetime > as.POSIXct("2024-06-26 12:00:00")) %>%
  filter(datetime < as.POSIXct("2024-08-15 00:00:00"))

NS_plot_BWLF <- plot_ly(data = filtered_BWL, x = ~datetime, y = ~Dissolved_Oxygen_offset, type = 'scatter',  mode = 'markers',  color = ~Sensor) %>%
  layout(xaxis = list(title = "Date and Time"),
         yaxis = list(title = "DO (mgL)"))


## BWU
filtered_BWU <- DO_offset %>%
  filter(site == "BWU")  %>%
  filter(datetime > as.POSIXct("2021-06-30 10:00:00")) %>%
  filter(datetime < as.POSIXct("2021-09-07 00:00:00") | datetime > as.POSIXct("2021-09-22 12:00:00"))%>%
  #filter(datetime < as.POSIXct("2022-07-08 00:00:00") | datetime > as.POSIXct("2022-07-10 00:00:00"))%>%
  #filter(datetime < as.POSIXct("2022-09-01 08:00:00") | datetime > as.POSIXct("2022-09-01 11:00:00")) %>%
  filter(datetime < as.POSIXct("2023-10-25 00:00:00"))
  
NS_plot_BWU <- plot_ly(data = filtered_BWU, x = ~datetime, y = ~Dissolved_Oxygen_offset, type = 'scatter',  mode = 'markers',  color = ~Sensor) %>%
  layout(xaxis = list(title = "Date and Time"),
         yaxis = list(title = "temp"))

## GBL
NS_plot_GBL <- plot_ly(data = DO_offset%>%filter(site=="GBL"), x = ~datetime, y = ~Dissolved_Oxygen_offset, type = 'scatter',  mode = 'markers',  color = ~Sensor) %>%
  layout(xaxis = list(title = "Date and Time"),
         yaxis = list(title = "DO (mgL)"))


filtered_GBL1 <- DO_offset %>%
  filter(site == "GBL")  %>%
  filter(wtr < 17)  %>%
  filter(datetime > as.POSIXct("2021-03-13 00:00:00")) %>%
  filter(datetime < as.POSIXct("2021-05-28 03:00:00") | datetime > as.POSIXct("2021-05-28 06:00:00"))%>%
  filter(datetime < as.POSIXct("2021-06-03 02:00:00") | datetime > as.POSIXct("2021-06-03 06:00:00"))%>%
  filter(datetime < as.POSIXct("2021-07-21 18:00:00") | datetime > as.POSIXct("2021-07-22 20:00:00"))%>%
  filter(datetime < as.POSIXct("2021-08-05 03:00:00") | datetime > as.POSIXct("2021-08-05 06:00:00")) %>%
  filter(datetime < as.POSIXct("2021-10-20 08:00:00") | datetime > as.POSIXct("2021-10-20 11:00:00")) %>%
  filter(datetime < as.POSIXct("2021-12-01 02:00:00") | datetime > as.POSIXct("2021-12-01 11:00:00")) %>%
  filter(datetime < as.POSIXct("2022-01-27 01:00:00") | datetime > as.POSIXct("2022-01-27 03:00:00")) %>%
  filter(datetime < as.POSIXct("2022-04-07 02:00:00") | datetime > as.POSIXct("2022-04-07 04:00:00")) %>%
  filter(datetime < as.POSIXct("2022-07-12 09:00:00") | datetime > as.POSIXct("2022-07-12 12:00:00")) %>%
  filter(datetime < as.POSIXct("2022-11-04 02:00:00") | datetime > as.POSIXct("2022-11-04 12:00:00")) %>%
  filter(datetime < as.POSIXct("2023-01-04 02:00:00") | datetime > as.POSIXct("2023-01-04 04:00:00")) %>%
  filter(datetime < as.POSIXct("2023-01-23 00:00:00") | datetime > as.POSIXct("2023-01-23 09:00:00")) %>%
  filter(datetime < as.POSIXct("2023-03-27 02:00:00") | datetime > as.POSIXct("2023-03-27 08:00:00")) %>%
  filter(datetime < as.POSIXct("2023-06-15 02:00:00") | datetime > as.POSIXct("2023-06-15 08:00:00")) %>%
  filter(datetime < as.POSIXct("2023-07-10 02:00:00") | datetime > as.POSIXct("2023-07-10 08:00:00")) %>%
  filter(datetime < as.POSIXct("2023-08-03 02:00:00") | datetime > as.POSIXct("2023-08-03 08:00:00")) %>%
  filter(datetime < as.POSIXct("2023-09-05 02:00:00") | datetime > as.POSIXct("2023-09-05 08:00:00")) %>%
  filter(datetime < as.POSIXct("2023-09-18 02:00:00") | datetime > as.POSIXct("2023-09-18 08:00:00")) %>%
  filter(datetime < as.POSIXct("2023-10-19 03:00:00") | datetime > as.POSIXct("2023-10-19 06:00:00")) %>%
  filter(datetime < as.POSIXct("2024-01-17 03:00:00") | datetime > as.POSIXct("2024-01-17 06:00:00")) %>%
  filter(datetime < as.POSIXct("2024-04-17 04:00:00") | datetime > as.POSIXct("2024-04-17 06:00:00")) %>%
  filter(datetime < as.POSIXct("2024-05-29 05:00:00") | datetime > as.POSIXct("2024-05-29 07:00:00")) %>%
  filter(datetime < as.POSIXct("2024-07-10 04:00:00") | datetime > as.POSIXct("2024-07-10 07:00:00")) %>%
  filter(datetime < as.POSIXct("2024-08-07 04:00:00") | datetime > as.POSIXct("2024-08-07 07:00:00")) %>%
  filter(datetime < as.POSIXct("2024-10-09 00:00:00"))

NS_plot_GBL1 <- plot_ly(data = filtered_GBL1, x = ~datetime, y = ~wtr, type = 'scatter',  mode = 'markers',  color = ~Sensor) %>%
  layout(xaxis = list(title = "Date and Time"),
         yaxis = list(title = "Temp ()"))

## GBU
NS_plot_GBU <- plot_ly(data = DO_offset%>%filter(site=="GBU"), x = ~datetime, y = ~wtr, type = 'scatter',  mode = 'markers',  color = ~Sensor) %>%
  layout(xaxis = list(title = "Date and Time"),
         yaxis = list(title = "Temp ()"))


filtered_GBU <- DO_offset %>%
  filter(site == "GBU")  %>%
  filter(wtr < 17.9)  %>%
  filter(datetime > as.POSIXct("2021-03-13 00:00:00")) %>%
  filter(datetime < as.POSIXct("2021-08-05 03:00:00") | datetime > as.POSIXct("2021-08-05 06:00:00")) %>%
  filter(datetime < as.POSIXct("2021-04-21 03:00:00") | datetime > as.POSIXct("2021-04-21 06:00:00")) %>%
  filter(datetime < as.POSIXct("2021-05-28 03:00:00") | datetime > as.POSIXct("2021-05-28 06:00:00")) %>%
  filter(datetime < as.POSIXct("2021-07-22 06:00:00") | datetime > as.POSIXct("2021-07-22 08:00:00")) %>%
  filter(datetime < as.POSIXct("2021-08-05 03:00:00") | datetime > as.POSIXct("2021-08-05 06:00:00")) %>%
  filter(datetime < as.POSIXct("2021-10-20 07:00:00") | datetime > as.POSIXct("2021-10-20 11:00:00")) %>%
  filter(datetime < as.POSIXct("2021-12-01 00:00:00") | datetime > as.POSIXct("2021-12-02 00:00:00")) %>%
  filter(datetime < as.POSIXct("2023-03-06 00:00:00") | datetime > as.POSIXct("2023-03-08 00:00:00")) %>%
  filter(datetime < as.POSIXct("2023-05-08 03:00:00") | datetime > as.POSIXct("2023-05-08 05:00:00")) %>%
  filter(datetime < as.POSIXct("2023-07-10 05:30:00") | datetime > as.POSIXct("2023-07-10 70:00:00")) %>%
  filter(datetime < as.POSIXct("2023-09-18 03:30:00") | datetime > as.POSIXct("2023-09-18 05:00:00")) %>%
  filter(datetime < as.POSIXct("2024-10-09 00:00:00"))


NS_plot_GBU <- plot_ly(data = filtered_GBU, x = ~datetime, y = ~wtr, type = 'scatter',  mode = 'markers',  color = ~Sensor) %>%
  layout(xaxis = list(title = "Date and Time"),
         yaxis = list(title = "Temp ()"))



DO_offset_BWL <- DO_offset %>%
  filter(site == "BWL")%>%
  mutate(maintenance_flag = if_else(
      (datetime < as.POSIXct("2020-09-27 09:00:00")) |
        (datetime > as.POSIXct("2021-04-29 00:00:00") & datetime < as.POSIXct("2021-04-29 12:00:00")) |
        (datetime > as.POSIXct("2021-05-12 04:30:00") & datetime < as.POSIXct("2021-05-12 08:00:00")) |
        (datetime > as.POSIXct("2021-05-26 06:00:00") & datetime < as.POSIXct("2021-05-26 07:30:00")) |
        (datetime > as.POSIXct("2021-06-02 06:30:00") & datetime < as.POSIXct("2021-06-02 07:50:00")) |
        (datetime > as.POSIXct("2021-11-12 04:30:00") & datetime < as.POSIXct("2021-11-12 10:30:00")) |
        (datetime > as.POSIXct("2023-03-07 00:30:00") & datetime < as.POSIXct("2023-03-09 00:00:00")) |
        (datetime > as.POSIXct("2023-04-05 00:30:00") & datetime < as.POSIXct("2023-04-06 00:00:00")) |
        (datetime > as.POSIXct("2023-05-07 00:00:00") & datetime < as.POSIXct("2023-05-07 02:00:00")) |
        (datetime > as.POSIXct("2023-05-19 00:00:00") & datetime < as.POSIXct("2023-05-19 02:00:00")) |
        (datetime > as.POSIXct("2023-06-12 03:00:00") & datetime < as.POSIXct("2023-06-12 05:00:00")) |
        (datetime > as.POSIXct("2023-07-12 03:00:00") & datetime < as.POSIXct("2023-07-12 05:00:00")) |
        (datetime > as.POSIXct("2023-08-23 03:00:00") & datetime < as.POSIXct("2023-08-23 09:00:00")) |
        (datetime > as.POSIXct("2024-06-26 03:00:00") & datetime < as.POSIXct("2024-06-26 12:00:00")) |
        (datetime > as.POSIXct("2024-08-15 00:00:00")), 
      1, 0))

temp_plot <- DO_offset_BWL %>%
  ggplot(aes(x = datetime, y = wtr, colour = as.factor(maintenance_flag))) +
  geom_point(size=1, alpha=0.75) + theme_bw() + facet_grid(site~.)


DO_offset_BWU <- DO_offset %>%
  filter(site == "BWU") %>%
  mutate(maintenance_flag = if_else(
    (datetime < as.POSIXct("2021-06-30 10:00:00")) |
        (datetime > as.POSIXct("2021-09-07 00:00:00") & datetime < as.POSIXct("2021-09-22 12:00:00")) |
        (datetime > as.POSIXct("2022-07-08 00:00:00") & datetime < as.POSIXct("2022-07-10 00:00:00")) |
        (datetime > as.POSIXct("2022-09-01 08:00:00") & datetime < as.POSIXct("2022-09-01 11:00:00")) |
        (datetime > as.POSIXct("2023-10-25 00:00:00")), 
      1, 0))

temp_plot <- DO_offset_BWU %>%
  ggplot(aes(x = datetime, y = wtr, colour = as.factor(maintenance_flag))) +
  geom_point(size=1, alpha=0.75) + theme_bw() + facet_grid(site~.)


DO_offset_GBL <- DO_offset %>%
  filter(site == "GBL" & wtr<17 ) %>%
  mutate(maintenance_flag = if_else(
      (datetime < as.POSIXct("2021-03-13 00:00:00")) |
        (datetime > as.POSIXct("2021-05-28 03:00:00") & datetime < as.POSIXct("2021-05-28 06:00:00")) | 
        (datetime > as.POSIXct("2021-06-03 02:00:00") & datetime < as.POSIXct("2021-06-03 06:00:00")) |
        (datetime > as.POSIXct("2021-07-21 18:00:00") & datetime < as.POSIXct("2021-07-22 20:00:00")) |
        (datetime > as.POSIXct("2021-08-05 03:00:00") & datetime < as.POSIXct("2021-08-05 06:00:00")) | 
        (datetime > as.POSIXct("2021-10-20 08:00:00") & datetime < as.POSIXct("2021-10-20 11:00:00")) |
        (datetime > as.POSIXct("2021-12-01 02:00:00") & datetime < as.POSIXct("2021-12-01 11:00:00")) |
        (datetime > as.POSIXct("2022-01-27 01:00:00") & datetime < as.POSIXct("2022-01-27 03:00:00")) |
        (datetime > as.POSIXct("2022-04-07 02:00:00") & datetime < as.POSIXct("2022-04-07 04:00:00")) |
        (datetime > as.POSIXct("2022-07-12 09:00:00") & datetime < as.POSIXct("2022-07-12 12:00:00")) |
        (datetime > as.POSIXct("2022-11-04 02:00:00") & datetime < as.POSIXct("2022-11-04 12:00:00")) |
        (datetime > as.POSIXct("2023-01-04 02:00:00") & datetime < as.POSIXct("2023-01-04 04:00:00")) |
        (datetime > as.POSIXct("2023-01-23 00:00:00") & datetime < as.POSIXct("2023-01-23 09:00:00")) |
        (datetime > as.POSIXct("2023-03-27 02:00:00") & datetime < as.POSIXct("2023-03-27 08:00:00")) |
        (datetime > as.POSIXct("2023-06-15 02:00:00") & datetime < as.POSIXct("2023-06-15 08:00:00")) |
        (datetime > as.POSIXct("2023-07-10 02:00:00") & datetime < as.POSIXct("2023-07-10 08:00:00")) |
        (datetime > as.POSIXct("2023-08-03 02:00:00") & datetime < as.POSIXct("2023-08-03 08:00:00")) |
        (datetime > as.POSIXct("2023-09-05 02:00:00") & datetime < as.POSIXct("2023-09-05 08:00:00")) |
        (datetime > as.POSIXct("2023-09-18 02:00:00") & datetime < as.POSIXct("2023-09-18 08:00:00")) |
        (datetime > as.POSIXct("2023-10-19 03:00:00") & datetime < as.POSIXct("2023-10-19 06:00:00")) |
        (datetime > as.POSIXct("2024-01-17 03:00:00") & datetime < as.POSIXct("2024-01-17 06:00:00")) |
        (datetime > as.POSIXct("2024-04-17 04:00:00") & datetime < as.POSIXct("2024-04-17 06:00:00")) |
        (datetime > as.POSIXct("2024-05-29 05:00:00") & datetime < as.POSIXct("2024-05-29 07:00:00")) |
        (datetime > as.POSIXct("2024-07-10 04:00:00") & datetime < as.POSIXct("2024-07-10 07:00:00")) |
        (datetime > as.POSIXct("2024-08-07 04:00:00") & datetime < as.POSIXct("2024-08-07 07:00:00")) |
        (datetime > as.POSIXct("2024-10-09 00:00:00")), 
      1, 0))

temp_plot <- DO_offset_GBL %>%
  ggplot(aes(x = datetime, y = wtr, colour = as.factor(maintenance_flag))) +
  geom_point(size=1, alpha=0.75) + theme_bw() + facet_grid(site~.)



DO_offset_GBU <- DO_offset %>%
  filter(site == "GBU" &  wtr < 17.9) %>%
  mutate(maintenance_flag = if_else( 
      (datetime < as.POSIXct("2021-03-13 00:00:00")) |
      (datetime > as.POSIXct("2021-08-05 03:00:00") & datetime < as.POSIXct("2021-08-05 06:00:00")) |
      (datetime > as.POSIXct("2021-04-21 03:00:00") & datetime < as.POSIXct("2021-04-21 06:00:00")) |
      (datetime > as.POSIXct("2021-05-28 03:00:00") & datetime < as.POSIXct("2021-05-28 06:00:00")) |
      (datetime > as.POSIXct("2021-07-22 06:00:00") & datetime < as.POSIXct("2021-07-22 08:00:00")) |
      (datetime > as.POSIXct("2021-08-05 03:00:00") & datetime < as.POSIXct("2021-08-05 06:00:00")) |
      (datetime > as.POSIXct("2021-10-20 07:00:00") & datetime < as.POSIXct("2021-10-20 11:00:00")) |
      (datetime > as.POSIXct("2021-12-01 00:00:00") & datetime < as.POSIXct("2021-12-02 00:00:00")) |
      (datetime > as.POSIXct("2023-03-06 00:00:00") & datetime < as.POSIXct("2023-03-08 00:00:00")) |
      (datetime > as.POSIXct("2023-05-08 03:00:00") & datetime < as.POSIXct("2023-05-08 05:00:00")) |
      (datetime > as.POSIXct("2023-07-10 05:30:00") & datetime < as.POSIXct("2023-07-10 70:00:00")) |
      (datetime > as.POSIXct("2023-09-18 03:30:00") & datetime < as.POSIXct("2023-09-18 05:00:00")) |
      (datetime > as.POSIXct("2024-10-09 00:00:00")), 
      1, 0)) 

DO_offset_GBU <- DO_offset %>%
  filter(site == "GBU" & wtr < 17.9) %>%
  mutate(maintenance_flag = if_else(
    (datetime > as.POSIXct("2021-03-13 00:00:00") & datetime < as.POSIXct("2021-08-05 03:00:00")) |
      (datetime > as.POSIXct("2021-08-05 03:00:00") & datetime < as.POSIXct("2021-08-05 06:00:00")) |
      (datetime > as.POSIXct("2021-04-21 03:00:00") & datetime < as.POSIXct("2021-04-21 06:00:00")) |
      (datetime > as.POSIXct("2021-05-28 03:00:00") & datetime < as.POSIXct("2021-05-28 06:00:00")) |
      (datetime > as.POSIXct("2021-07-22 06:00:00") & datetime < as.POSIXct("2021-07-22 08:00:00")) |
      (datetime > as.POSIXct("2021-10-20 07:00:00") & datetime < as.POSIXct("2021-10-20 11:00:00")) |
      (datetime > as.POSIXct("2021-12-01 00:00:00") & datetime < as.POSIXct("2021-12-02 00:00:00")) |
      (datetime > as.POSIXct("2023-03-06 00:00:00") & datetime < as.POSIXct("2023-03-08 00:00:00")) |
      (datetime > as.POSIXct("2023-05-08 03:00:00") & datetime < as.POSIXct("2023-05-08 05:00:00")) |
      (datetime > as.POSIXct("2023-07-10 05:30:00") & datetime < as.POSIXct("2023-07-10 07:00:00")) |
      (datetime > as.POSIXct("2023-09-18 03:30:00") & datetime < as.POSIXct("2023-09-18 05:00:00")) |
      (datetime < as.POSIXct("2024-10-09 00:00:00")),
    1, 0
  ))



temp_plot <- DO_offset_GBU %>%
  ggplot(aes(x = datetime, y = wtr, colour = as.factor(maintenance_flag))) +
  geom_point(size=1, alpha=0.75) + theme_bw() + facet_grid(site~.)




# Define a function to identify if the current value is outside bounds from the rolling window
identify_outliers_any <- function(x) {
  q1 <- quantile(x, 0.15, na.rm = TRUE)
  q3 <- quantile(x, 0.85, na.rm = TRUE)
  iqr <- q3 - q1
  lower_bound <- q1 - c(3 * iqr)
  upper_bound <- q3 + c(3 * iqr)
    current_value <- x[length(x)]  
  is_outlier <- current_value < lower_bound | current_value > upper_bound
  
  return(is_outlier)
}

# Apply the outlier detection within each day
DO_flag_BWU <- DO_offset_BWU %>%
  #filter(datetime > as.POSIXct("2022-05-01 00:00:00") & datetime < as.POSIXct("2022-06-01 00:00:00")) %>%
  arrange(datetime) %>%
  #group_by(day = as.Date(datetime)) %>%
  mutate(
    DO_outlier = rollapply(Dissolved_Oxygen_offset, width = 36, 
                           FUN = identify_outliers_any, fill = NA, align = "right"),
    wtr_outlier = rollapply(wtr, width = 24, 
                            FUN = identify_outliers_any, fill = NA, align = "right")
  ) %>%
  ungroup()

# Convert logical outlier columns to character labels if desired
DO_flag_BWU <- DO_flag_BWU %>%
  mutate(
    DO_outlier = ifelse(DO_outlier, 1, 0),
    wtr_outlier = ifelse(wtr_outlier, 1, 0)
  )

# Quick visualization
DO_plot <- DO_flag_BWU %>%
  ggplot(aes(x = datetime, y = Dissolved_Oxygen_offset, colour = as.factor(DO_outlier))) +
  geom_point(size=1, alpha=0.75) + theme_bw() +
  labs(colour = "Outlier Status")
DO_plot


wtemp_plot <- DO_flag_BWL %>%
  ggplot(aes(x = datetime, y = wtr, colour = wtr_outlier)) +
  geom_point(size=1, alpha=0.75) + theme_bw() +
  labs(colour = "Outlier Status")
wtemp_plot



# saveRDS(DO_flag_BWU, file = "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/24_BWU_DO_flag_record.rds")



# 
# # Define a function to detect outliers using a 1.5*IQR rule for a single centered value in the rolling window
# identify_outlier_center <- function(x) {
#   q1 <- quantile(x, 0.25, na.rm = TRUE)
#   q3 <- quantile(x, 0.75, na.rm = TRUE)
#   iqr <- q3 - q1
#   lower_bound <- q1 - c(1.5 * iqr)
#   upper_bound <- q3 + c(1.5 * iqr)
#   # Return TRUE if the center value of the window is an outlier
#   x[ceiling(length(x) / 2)] < lower_bound | x[ceiling(length(x) / 2)] > upper_bound
# }
# 
# 
# ####
# ####
# ####
# 
# # Define a function to identify any values in the window that exceed the bounds
# identify_outliers_in_window <- function(x) {
#   q1 <- quantile(x, 0.25, na.rm = TRUE)
#   q3 <- quantile(x, 0.75, na.rm = TRUE)
#   iqr <- q3 - q1
#   lower_bound <- q1 - 1.5 * iqr
#   upper_bound <- q3 + 1.5 * iqr
#   
#   # Identify any values in the window that are outside the bounds
#   outliers <- x[x < lower_bound | x > upper_bound]
#   
#   # Return the outliers if there are any; otherwise, return NA
#   if (length(outliers) > 0) {
#     return(outliers)
#   } else {
#     return(NA)
#   }
# }
# 
# 
# 
# # Apply the outlier detection within each day
# DO_flag_BWL <- DO_offset_BWL %>%
#   filter(datetime > as.POSIXct("2022-04-01 00:00:00") & datetime < as.POSIXct("2022-06-01 00:00:00")) %>%
#   arrange(datetime) %>%
#   group_by(day = as.Date(datetime)) %>%
#   mutate(
#     DO_outlier = rollapply(Dissolved_Oxygen_offset, width = 8, 
#                            FUN = identify_outliers_in_window, fill = NA, align = "right"),
#     wtr_outlier = rollapply(wtr, width = 8, 
#                             FUN = identify_outliers_in_window, fill = NA, align = "right")
#   ) %>%
#   ungroup()
# 
# # Convert logical outlier columns to character labels if desired
# DO_flag_BWL <- DO_flag_BWL %>%
#   mutate(
#     DO_outlier = ifelse(DO_outlier, "DO Outlier", NA),
#     wtr_outlier = ifelse(wtr_outlier, "WTR Outlier", NA)
#   )
# 
# # Quick visualization
# DO_plot <- DO_flag_BWL %>%
#   ggplot(aes(x = datetime, y = Dissolved_Oxygen_offset, colour = DO_outlier)) +
#   geom_point(size=1, alpha=0.75) + theme_bw() +
#   labs(colour = "Outlier Status")
# DO_plot
# 
# 
# temp_plot <- DO_flag_BWL %>%
#   ggplot(aes(x = datetime, y = wtr, colour = wtr_outlier)) +
#   geom_point(size=1, alpha=0.75) + theme_bw() +
#   labs(colour = "Outlier Status")
# temp_plot
# 
# 
# # 
# # 
# # # calculate % saturation based on bucket: 
# # ###function to estimate oxygen saturation. 
# # # From Garcia and Gordon 1992 L&O.  Takes temp in deg C and bp in mm Hg. Yes, I know taht this calculationis a slight approximation because we do not account for water density, but the difference is much smaller (0.01 mg/L at 20 deg) than anyone's ability to calibrate an oxygen sonde.
# # # Define the function to calculate oxygen saturation
# # osat <- function(temp, bp) {
# #   sato <- (exp(2.00907 + 3.22014 * (log((298.15 - temp) / (273.15 + temp))) + 
# #                  4.0501 * (log((298.15 - temp) / (273.15 + temp))) ^ 2 + 
# #                  4.94457 * (log((298.15 - temp) / (273.15 + temp))) ^ 3 - 
# #                  0.256847 * (log((298.15 - temp) / (273.15 + temp))) ^ 4 + 
# #                  3.88767 * (log((298.15 - temp) / (273.15 + temp))) ^ 5)) * 
# #     1.4276 * bp / 760
# #   return(sato)
# # }
# # 
# # # Apply the function to calculate oxygen saturation for each row
# # df <- ns_DOcal2 %>%
# #   mutate(oxygen_saturation = round(osat(Temperature_mean, bp),4))
# # 
