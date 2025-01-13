# GBU stream DO
# only updated lower 2025-01-12

#=========== Preliminaries
rm(list=ls())
# load packages
library(tidyverse)
library(lubridate)
library(plotly)
library(devtools)
library(nrlmetab)
 # remotes::install_github("nrlottig/nrlmetab")
library(zoo)
library(suncalc)
library(readxl)
library(patchwork)
library(gridExtra)
library(dplyr)
library(reshape)
library(lubridate)
library(ggplot2)
library(scales)

library(dplyr)
library(readr)
library(purrr)
library(stringr)


#=========================================== 
# Get and process high frequency sensor data

# Directory with files
data_dir <- "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/GBU_DO/RawDO_GBU"

# Column names to assign after skipping the first 9 rows
column_names <- c("Unix", "UTC_Date_Time","PCT_datetime", "Battery", "Temperature", "DO", "DO_sat", "Q", "Sensor", "Concat_Date")

# List and process files
rawdat <- list.files(data_dir, full.names = TRUE) %>%
  # Map over each file to read, extract sensor and concatenation date, and bind data
  map_df(~{
    # Read in the full file (including header lines for sensor and concat date extraction)
    file_data <- read_csv(.x, col_names = FALSE)
    
    # Extract the sensor number from the first line
    sensor_line <- file_data$X1[2]
    sensor_number <- str_extract(sensor_line, "\\d{4}-\\d{6}")
    
    # Extract the concatenation date from line 3
    concat_date_line <- file_data$X1[3]
    concat_date <- str_extract(concat_date_line, "\\d{4}\\w{3}\\d{2} \\d{2}:\\d{2}:\\d{2}")
    
    # Read data starting from line 10, skipping headers
    data_start <- read_csv(.x, skip = 9, col_names = column_names)
    
    # Add sensor and concatenation date columns
    data_start %>%
      mutate(Sensor = sensor_number, Concat_Date = concat_date)
  })


range(rawdat$PCT_datetime)
unique(rawdat$Sensor)
summary(rawdat)


# Convert PCT_datetime to correct for daylight savings
rawdat1 <- rawdat %>%
  mutate(
    # Specify that the current datetime is in UTC
    UTC_Date_Time = force_tz(UTC_Date_Time, tzone = "UTC"),
    
    # Convert from UTC to Pacific Time (automatically adjusts for daylight savings)
    PCT_datetime_adjusted = with_tz(UTC_Date_Time, tzone = "America/Los_Angeles")
  )


range(na.omit(rawdat1$PCT_datetime_adjusted))
range((rawdat1$PCT_datetime_adjusted))

# time sieries for data cleaning 
do.ts <- rawdat1 %>%
  dplyr::rename(datetime = "PCT_datetime_adjusted", do.obs = "DO", wtr="Temperature") %>%
  #mutate(datetime = as.POSIXct((datetime), format ="%Y-%m-%d %H:%M:%S"))%>%
  dplyr::select("datetime", "do.obs", "wtr", "Battery", "Q", "Sensor", "Concat_Date")

do.ts$datetime <-  (as.POSIXct(round_date(
  as.POSIXct(do.ts$datetime, format="%Y-%m-%d %H:%M:%S"), unit="5 minutes")))

### save concated sensor records:
do.ts$site <- "GBU"
# write.csv(x = do.ts, file = "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/raw_merged_dat/25_GBU_rawDO_record.csv", row.names = TRUE)
# saveRDS(do.ts, file = "R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/raw_merged_dat/25_GBU_rawDO_record.rds")


tempcheck <- subset(do.ts, 
                    datetime > '2021-09-16 00:00:00' & datetime < '2021-09-17 :00:00') # sensor quality minimum threshold

qplot(datetime, do.obs, data = tempcheck, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(breaks = scales::date_breaks("1 hours")) +
  geom_vline(xintercept = as.POSIXct(c("2021-09-16 00:00:00", 
                                       "2021-09-16 12:00:00", 
                                       "2021-09-16 24:00:00")), color = "yellow")


qplot(datetime, wtr, data = tempcheck, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(breaks = date_breaks("1 hours")) +
  geom_vline(xintercept = as.POSIXct(c("2021-09-16 00:00:00", 
                                       "2021-09-16 12:00:00", 
                                       "2021-09-16 24:00:00")), color = "yellow")


do.ts1 <- do.ts %>%
  dplyr::select("datetime", "do.obs")

###clean DO data
ggplot(data = do.ts1,aes(x=datetime,y=do.obs)) + geom_line()
describe.ts(do.ts)
do.ts.clean <- trim.outliers(do.ts, width = 7,sd.dev = 3) # usually use 3, but the spikes might be poor.  
ggplot(data = do.ts,aes(x=datetime,y=do.obs)) + geom_line() + 
  geom_line(data=do.ts.clean,aes(x=datetime,y=do.obs),col="red")
do.ts.avg <- aggregate.data(data = do.ts.clean,time.step = 15)

######
wtr.ts <- do.ts %>%
  dplyr::select("datetime", "wtr")


###clean wtr data
ggplot(data = wtr.ts,aes(x=datetime,y=wtr)) + geom_line()
describe.ts(wtr.ts)
wtr.ts.clean <- trim.outliers(wtr.ts, width = 7, sd.dev = 3) 
ggplot(data = wtr.ts,aes(x=datetime,y=wtr)) + geom_line() + 
  geom_line(data=wtr.ts.clean,aes(x=datetime,y=wtr),col="red")
wtr.ts.avg <- aggregate.data(data = wtr.ts.clean,time.step = 15)

### meta data 
do.ts_met <-do.ts %>%
  select(datetime, Battery, Q, Sensor, Concat_Date)

dat_out <- do.ts.avg  %>% 
  full_join(wtr.ts.avg, by = "datetime") %>%
  full_join(do.ts_met, by = "datetime")

dat_out<- na.omit(dat_out)
ggplot(data = dat_out,aes(x=datetime,y=wtr)) + geom_line() 
ggplot(data = dat_out,aes(x=datetime,y=do.obs)) + geom_line() 
  
####
# Create flags 
# days with large variation in water temperature
# days with high flow
# days where sensor was downloaded 

dat_out$date <- as.Date(dat_out$datetime)
str(dat_out)


daily_temp_range <- dat_out %>%
  group_by(date) %>%
  summarize(temp_range = max(wtr, na.rm = TRUE) - min(wtr, na.rm = TRUE))

hist(daily_temp_range$temp_range)
ggplot(data = daily_temp_range,aes(x=date,y=temp_range)) + geom_line() 


daily_DO_range <- dat_out %>%
  group_by(date) %>%
  summarize(DO_range = max(do.obs, na.rm = TRUE) - min(do.obs, na.rm = TRUE))

hist(daily_DO_range$DO_range)
ggplot(data = daily_DO_range,aes(x=date,y=DO_range)) + geom_point() 

plot(x=daily_temp_range$temp_range,y=daily_DO_range$DO_range)
































### climate data ###
# 2. merge in weather station data:
baro_dat <- read.csv("/Users/kellyloria/Documents/UNR/MSMmetab/climate/stream_NLDAS_baro.csv")

baro_datQ <- baro_dat %>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC")) %>%
  with_tz(tz = "America/Los_Angeles") %>%
  select(site, datetime, baro_Pa)

light_dat <- read.csv("/Users/kellyloria/Documents/UNR/MSMmetab/climate/stream_NLDAS_light.csv")
light_datQ <- light_dat %>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC")) %>%
  with_tz(tz = "America/Los_Angeles") %>%
  dplyr::select(site, datetime, light)

# 
dat2 <- do.ts.avg  %>% full_join(wtr.ts.avg, by = "datetime")
dat2$site <- "BWL"

summary(dat2)
# okay lets figure out these NA's 
dat2t<- na.omit(dat2)
summary(dat2t)

qplot(datetime, do.obs, data = dat2t, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(breaks = date_breaks("1000 hours"))

qplot(datetime, wtr, data = dat2t, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(breaks = date_breaks("1000 hours"))

dat_Q <- dat2t %>%
  left_join(light_datQ, by = c("datetime", "site")) %>%
  left_join(baro_datQ, by = c("datetime", "site")) 

# infil missing values
dat_Q1 <- dat_Q %>%
  fill(light, .direction = "down")%>%
  fill(baro_Pa,.direction = "down")%>% 
  mutate(baro= (baro_Pa*0.01)) # %>% dplyr::ungroup()

summary(dat_Q1)

dat_Q1 <- dat_Q1 %>%
  fill(light, .direction = "up")%>%
  fill(baro,.direction = "up")


## Add in flow
library(dataRetrieval)

siteNo <- "10336660"
pCode <- c("00060", "00065")
start.date <- "2021-04-29"
end.date <- "2024-01-01"

BWflow <- readNWISuv(siteNumbers = siteNo,
                     parameterCd = pCode,
                     startDate = start.date,
                     endDate = end.date)

flow.ts <- BWflow %>% 
  dplyr::rename(datetime = "dateTime", dischargeCFS = "X_00060_00000", gageHF= "X_00065_00000") %>%
  mutate(datetime = as.POSIXct((datetime), format ="%Y-%m-%d %H:%M:%S"))%>%
  dplyr::select("datetime", "dischargeCFS", "gageHF")

dat_Q2 <-  dat_Q1 %>%
  left_join(flow.ts, by = c("datetime")) 

summary(dat_Q2)

# infil missing values
dat_Q3 <- dat_Q2 %>%
  fill(dischargeCFS, .direction = "down")%>%
  fill(gageHF,.direction = "down")

qplot(datetime, light, data = dat_Q3, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(breaks = date_breaks("1000 hours"))

qplot(datetime, dischargeCFS, data = dat_Q3, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(breaks = date_breaks("1000 hours"))


dat <- dat_Q3 %>%
  dplyr::select(site,
                datetime,
                do.obs,
                wtr,
                light, 
                baro,
                dischargeCFS,
                gageHF) #

summary(dat)
## write.csv(x = dat, file = "/Users/kellyloria/Documents/UNR/MSMmetab/23_CleanDat/24_BWLInputs.csv", row.names = TRUE)

# saveRDS(dat, file = "/Users/kellyloria/Documents/UNR/MSMmetab/23_CleanDat/24_BWL_Inputs.rds")



##
##
##
#
#######################
##
###
#=========================
# UPPER site?
#=========================
###
##

rawdat = list.files(paste("/Users/kellyloria/Documents/UNR/MSMmetab/Blackwood/DodatL/",sep=""), full.names = T) %>%
  lapply(read_csv) %>%
  bind_rows() 


do.ts <- rawdat %>%
  dplyr::rename(datetime = "UTC_Date_&_Time", do.obs = "Dissolved Oxygen") %>%
  mutate(datetime = as_datetime(datetime, "America/Los_Angeles")) %>% 
  #mutate(datetime = as.POSIXct((datetime), format ="%Y-%m-%d %H:%M:%S"))%>%
  dplyr::select("datetime", "do.obs", "Battery", "Q")

do.ts$datetime <-  (as.POSIXct(round_date(
  as.POSIXct(do.ts$datetime, format="%Y-%m-%d %H:%M:%S"), unit="5 minutes")))

do.ts <- subset(do.ts, 
                datetime > "2021-04-30 00:00:00",
                Q>0.7) # sensor quality minimum threshold
#do.ts$datetime <- do.ts$datetime - hours(4)
range(do.ts$datetime)


tempcheck <- subset(do.ts, 
                    datetime > '2022-09-16 00:00:00' & datetime < '2022-09-17 :00:00') # sensor quality minimum threshold

qplot(datetime, do.obs, data = tempcheck, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(breaks = date_breaks("1 hours")) +
  geom_vline(xintercept = as.POSIXct(c("2022-09-16 00:00:00", 
                                       "2022-09-16 12:00:00", 
                                       "2022-09-16 24:00:00")), color = "yellow")

###
wtr.ts <- rawdat %>% 
  dplyr::rename(datetime = "UTC_Date_&_Time", wtr = "Temperature") %>%
  mutate(datetime = as_datetime(datetime, "America/Los_Angeles")) %>% 
  #mutate(datetime = as.POSIXct((datetime), format ="%Y-%m-%d %H:%M:%S"))%>%
  dplyr::select("datetime", "wtr")

wtr.ts$datetime <-  (as.POSIXct(round_date(
  as.POSIXct(wtr.ts$datetime, format="%Y-%m-%d %H:%M:%S"), unit="5 minutes")))

tempcheck <- subset(wtr.ts, 
                    datetime > '2022-09-16 00:00:00' & datetime < '2022-09-17 :00:00') # sensor quality minimum threshold

qplot(datetime, wtr, data = tempcheck, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(breaks = date_breaks("1 hours")) +
  geom_vline(xintercept = as.POSIXct(c("2022-09-16 00:00:00", 
                                       "2022-09-16 12:00:00", 
                                       "2022-09-16 24:00:00")), color = "yellow")


###clean DO data
ggplot(data = do.ts,aes(x=datetime,y=do.obs)) + geom_line()
describe.ts(do.ts)
do.ts.clean <- trim.outliers(do.ts,width = 8,sd.dev = 6) # usually use 3, but the spikes might be poor.  
ggplot(data = do.ts,aes(x=datetime,y=do.obs)) + geom_line() + 
  geom_line(data=do.ts.clean,aes(x=datetime,y=do.obs),col="red")
do.ts.avg21 <- aggregate.data(data = do.ts.clean21,time.step = 15)



######
###clean wtr data
ggplot(data = wtr.ts,aes(x=datetime,y=wtr)) + geom_line()
describe.ts(wtr.ts)
wtr.ts.clean <- trim.outliers(wtr.ts,width = 12, sd.dev = 3) 
ggplot(data = wtr.ts,aes(x=datetime,y=wtr)) + geom_line() + 
  geom_line(data=wtr.ts.clean,aes(x=datetime,y=wtr),col="red")
wtr.ts.avg21 <- aggregate.data(data = wtr.ts.clean21,time.step = 15)

##
##
##
# add in climate:
climate <- read_csv("/Users/kellyloria/Documents/UNR/MSMmetab/D9413.2023-09-01.csv", skip=) %>% 
  dplyr::rename(datetime='Date_Time',
                Atemp='air_temp_set_1',
                wspeed='wind_speed_set_1',
                baro='pressure_set_1d')%>%
  mutate(datetime = as_datetime(datetime, "America/Los_Angeles")) %>% 
  mutate(year=year(datetime),
         yday=yday(datetime),
         hour=hour(datetime)) %>%
  dplyr::select(datetime,
                Atemp,
                baro,
                wspeed,
                hour,
                yday) # datetime in UTC

climate$datetime <- (as.POSIXct(round_date(
  as.POSIXct(climate$datetime, format="%Y-%m-%d %H:%M:%S"), unit="5 minutes")))

tempcheck <- subset(climate, 
                    datetime > '2021-09-16 00:00:00' & datetime < '2021-09-17 :00:00') # sensor quality minimum threshold

qplot(datetime, Atemp, data = tempcheck, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(breaks = date_breaks("1 hours")) +
  geom_vline(xintercept = as.POSIXct(c("2021-09-16 00:00:00", 
                                       "2021-09-16 12:00:00", 
                                       "2021-09-16 24:00:00")), color = "yellow")

baro.ts <- climate %>%
  dplyr::select(datetime,
                baro) #

ggplot(data = baro.ts, aes(x=datetime,y=baro)) + geom_line()
describe.ts(baro.ts)
baro.ts.clean <- trim.outliers(baro.ts,width = 10,sd.dev = 5) 
ggplot(data = baro.ts,aes(x=datetime,y=baro)) + geom_line() + 
  geom_line(data=baro.ts.clean,aes(x=datetime,y=baro),col="red")
baro.ts.avg <- aggregate.data(data = baro.ts.clean,time.step = 15)

###
clim2 <- read_csv("/Users/kellyloria/Documents/UNR/MSMmetab/HMDC1.2023-09-01.csv") %>% 
  dplyr::rename(datetime='Date_Time',
                Atemp='air_temp_set_1',
                par='solar_radiation_set_1',
                wspeed='wind_speed_set_1')%>%
  # mutate(datetime = as.POSIXct((datetime), format ="%Y-%m-%d %H:%M:%S")) %>%
  mutate(datetime = as_datetime(datetime, "America/Los_Angeles")) %>% 
  mutate(year=year(datetime),
         yday=yday(datetime),
         hour=hour(datetime)) %>%
  dplyr::select(datetime,
                Atemp,
                par,
                wspeed,
                hour,
                yday) # datetime in UTC

#clim2$datetime <- clim2$datetime - c(60*60*8)
## Check the dates:
tempcheck <- subset(clim2, 
                    datetime > '2021-09-16 00:00:00' & datetime < '2021-09-17 :00:00') # sensor quality minimum threshold

qplot(datetime, par, data = tempcheck, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(breaks = date_breaks("1 hours")) +
  geom_vline(xintercept = as.POSIXct(c("2021-09-16 00:00:00", 
                                       "2021-09-16 12:00:00", 
                                       "2021-09-16 24:00:00")), color = "yellow")


clim2 <- clim2 %>%
  dplyr::select(datetime,
                par) #

clim2$datetime <- (as.POSIXct(round_date(
  as.POSIXct(clim2$datetime, format="%Y-%m-%d %H:%M:%S"), unit="5 minutes")))



ggplot(data = clim2,aes(x=datetime,y=par)) + geom_line()
describe.ts(clim2)
par.ts.clean <- trim.outliers(clim2,width = 10,sd.dev = 5) 
ggplot(data = clim2,aes(x=datetime,y=par)) + geom_line() + 
  geom_line(data=par.ts.clean,aes(x=datetime,y=par),col="red")
par.ts.avg <- aggregate.data(data = par.ts.clean,time.step = 15)


### infill ##
###
clim3 <- read_csv("/Users/kellyloria/Documents/UNR/MSMmetab/F9917.2023-09-01.csv") %>% 
  dplyr::rename(datetime='Date_Time',
                Atemp='air_temp_set_1',
                par='solar_radiation_set_1',
                baro='pressure_set_1d')%>%
  # mutate(datetime = as.POSIXct((datetime), format ="%Y-%m-%d %H:%M:%S")) %>%
  mutate(datetime = as_datetime(datetime, "America/Los_Angeles")) %>% 
  mutate(year=year(datetime),
         yday=yday(datetime),
         hour=hour(datetime)) %>%
  dplyr::select(datetime,
                Atemp,
                par,
                baro,
                hour,
                yday) # datetime in UTC

#clim2$datetime <- clim2$datetime - c(60*60*8)
## Check the dates:
tempcheck <- subset(clim3, 
                    datetime > '2023-06-16 00:00:00' & datetime < '2023-06-17 :00:00') # sensor quality minimum threshold

qplot(datetime, par, data = tempcheck, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(breaks = date_breaks("1 hours")) +
  geom_vline(xintercept = as.POSIXct(c("2023-06-16 00:00:00", 
                                       "2023-06-16 12:00:00", 
                                       "2023-06-16 24:00:00")), color = "yellow")

clim3$datetime <- (as.POSIXct(round_date(
  as.POSIXct(clim3$datetime, format="%Y-%m-%d %H:%M:%S"), unit="5 minutes")))

clim3 <- clim3 %>%
  dplyr::select(datetime,
                par, 
                baro) #



ggplot(data = clim3,aes(x=datetime,y=par)) + geom_line()
describe.ts(clim2)
par.ts.clean <- trim.outliers(clim2,width = 10,sd.dev = 5) 
ggplot(data = clim2,aes(x=datetime,y=par)) + geom_line() + 
  geom_line(data=par.ts.clean,aes(x=datetime,y=par),col="red")
par.ts.avg <- aggregate.data(data = par.ts.clean,time.step = 15)




# 
dat <- do.ts.clean%>% #note in this case I am not using the drift correction -KL
  full_join(wtr.ts.clean) %>% 
  full_join(clim2) %>%
  full_join(baro.ts.clean) 

summary(dat)

## Add in flow
library(dataRetrieval)

siteNo <- "10336660"
pCode <- c("00060", "00065")
start.date <- "2021-04-29"
end.date <- "2023-09-01"

BWflow <- readNWISuv(siteNumbers = siteNo,
                     parameterCd = pCode,
                     startDate = start.date,
                     endDate = end.date)

flow.ts <- BWflow %>% 
  dplyr::rename(datetime = "dateTime", dischargeCFS = "X_00060_00000", gageHF= "X_00065_00000") %>%
  mutate(datetime = as.POSIXct((datetime), format ="%Y-%m-%d %H:%M:%S"))%>%
  dplyr::select("datetime", "dischargeCFS", "gageHF")

flow.ts$datetime <- (as.POSIXct(round_date(
  as.POSIXct(flow.ts$datetime, format="%Y-%m-%d %H:%M:%S"), unit="5 minutes")))


## need to build a flow based rating curve here ##
datv1 <- dat %>%
  full_join(flow.ts) 

dat_unique <- datv1 %>%
  group_by(datetime) %>%
  summarize(do.obs = mean(do.obs, na.rm=T),
            wtr = mean(wtr, na.rm=T),
            par = mean(par, na.rm=T),
            baro = mean(baro, na.rm=T),
            dischargeCFS = mean(dischargeCFS, na.rm=T),
            gageHF = mean(gageHF, na.rm=T))

#dat_unique <- dat %>% distinct(datetime, .keep_all = TRUE)

dat1<- dat_unique[!is.na(dat_unique$do.obs), ]
dat2<- dat1[!is.na(dat1$wtr), ]

summary(dat2)

# Calculate the rolling average and store it in 'gageM'
dat2$gageM <- rollapply(dat2$gageHF, width = 10,
                           FUN = function(x) mean(x, na.rm = TRUE), by = 1,
                           by.column = TRUE, partial = TRUE, fill = NA, align = "center")
# Loop through each row and fill missing 'gageHF' values with 'gageM'
for (i in 1:nrow(dat2)) {
  if (is.na(dat2$gageHF[i])) {
    dat2$gageHF[i] <- dat2$gageM[i]
  }
}
# Summary of the updated 'dat21_2' dataset
summary(dat2)

qplot(datetime, gageHF, data = dat2, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(breaks = date_breaks("1000 hours"))


# Calculate the rolling average and store it in 'gageM'
dat2$dischargeM <- rollapply(dat2$dischargeCFS, width = 10,
                        FUN = function(x) mean(x, na.rm = TRUE), by = 1,
                        by.column = TRUE, partial = TRUE, fill = NA, align = "center")
# Loop through each row and fill missing 'gageHF' values with 'gageM'
for (i in 1:nrow(dat2)) {
  if (is.na(dat2$dischargeCFS[i])) {
    dat2$dischargeCFS[i] <- dat2$dischargeM[i]
  }
}
# Summary of the updated 'dat21_2' dataset
summary(dat2)

qplot(datetime, dischargeCFS, data = dat2, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(breaks = date_breaks("1000 hours"))

# Calculate the rolling average and store it in 'gageM'
dat2$baroM <- rollapply(dat2$baro, width = 10,
                             FUN = function(x) mean(x, na.rm = TRUE), by = 1,
                             by.column = TRUE, partial = TRUE, fill = NA, align = "center")
# Loop through each row and fill missing 'gageHF' values with 'gageM'
for (i in 1:nrow(dat2)) {
  if (is.na(dat2$baro[i])) {
    dat2$baro[i] <- dat2$baroM[i]
  }
}

qplot(datetime, baro, data = dat2, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(breaks = date_breaks("1000 hours"))

# Summary of the updated 'dat21_2' dataset
summary(dat2)

datsave<- dat2

baro_infil <-left_join(clim3, climate, 
                      by=c('datetime'))

bar_reg <- lm(baro_infil$baro.x ~ baro_infil$baro.y, data=baro_infil)
summary(bar_reg)

baro_infil$baroM<- ((baro_infil$baro.x) * (summary(bar_reg)$coef[2,1]))
hist(baro_infil$baroM) 

for (i in 1:nrow(dat2)) {
  if (is.na(dat2$baro[i])) {
    dat2$baro[i] <- baro_infil$baro.x[i]
  }
}
summary(dat2)

qplot(datetime, baro, data = dat2, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(breaks = date_breaks("1000 hours"))

# Calculate the rolling average and store it in 'gageM'
dat2$parM <- rollapply(dat2$par, width = 100,
                        FUN = function(x) mean(x, na.rm = TRUE), by = 1,
                        by.column = TRUE, partial = TRUE, fill = NA, align = "center")
# Loop through each row and fill missing 'gageHF' values with 'gageM'
for (i in 1:nrow(dat2)) {
  if (is.na(dat2$par[i])) {
    dat2$par[i] <- dat2$parM[i]
  }
}

par_infil <- clim3 %>%
  full_join(clim2)

par_infil <-left_join(clim3, clim2, 
                      by=c('datetime'))

par_reg <- lm(par_infil$par.x ~ par_infil$par.y, data=par_infil)
summary(par_reg)

par_infil$parM<- ((par_infil$par.x) * summary(par_reg)$coef[2,1])
hist(par_infil$parM) 
hist(par_infil$par.x)
hist(par_infil$par.y)

for (i in 1:nrow(dat2)) {
  if (is.na(dat2$par[i])) {
    dat2$par[i] <- par_infil$parM[i]
  }
}
# Summary of the updated 'dat21_2' dataset
summary(dat2)

qplot(datetime, parM, data = dat2, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(breaks = date_breaks("1000 hours"))

dat3 <- dat2 %>%
  dplyr::select(datetime,
                do.obs,
                wtr,
                par, 
                baro,
                dischargeCFS,
                gageHF) #


# write.csv(x = dat3, file = "/Users/kellyloria/Documents/UNR/MSMmetab/23_CleanDat/23_BWUInputs.csv", row.names = TRUE)
