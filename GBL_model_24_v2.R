# single stage
# 2024-12-09

## ---------------------------
rm(list=ls())

# Load packages
library(StreamMetabolism)
library(streamMetabolizer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(rstan)
library(unitted)
library(zoo)
library(lubridate)
library(dataRetrieval)
library(scales)
library(StanHeaders)

set.seed(2021)

#  install.packages("remotes")
#  remotes::install_github("appling/unitted")

#remotes::install_github("USGS-R/streamMetabolizer", force = TRUE)

## ---------------------------
## Read in DO data from miniDOT deployments 03/25-10/01
##
# PC path  setwd("R:/Users/kloria/Documents/2023_StreamMetab")

#"R:\Users\kloria\Documents\Stream_Metab_24\Core_sites\offset_DO_dat\saturation\24_BWL_DO_flag_sat.rds"

site <- "GBL"
dat <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/light/24_GBL_DO_flag_sat_light.rds") %>%
  filter(maintenance_flag<1 & maintenance_flag<1 & datetime > as.POSIXct("2021-03-12 00:00:00") & datetime < as.POSIXct("2024-08-01 00:00:00") & 
           do.obs>3.85 & wtr>0.01)
summary(dat)


dat <- dat %>%
  # Step 1: Compute hourly averages
  mutate(
    rounded_datetime = as.POSIXct(round(as.numeric(datetime) / 3600) * 3600, origin = "1970-01-01", tz = "America/Los_Angeles")
  ) %>%
  group_by(rounded_datetime) %>%
  summarize(
    hourly_DO.obs = mean(do.obs, na.rm = TRUE),
    hourly_DO.sat = mean(DO.sat, na.rm = TRUE),
    hourly_temp.water = mean(wtr, na.rm = TRUE),
    hourly_depth = mean(depth, na.rm = TRUE),
    hourly_discharge = mean(dischargeCMS, na.rm = TRUE),
    hourly_light = mean(light, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Step 2: Merge hourly averages back into the original 15-minute data
  right_join(
    dat %>% mutate(rounded_datetime = as.POSIXct(round(as.numeric(datetime) / 3600) * 3600, origin = "1970-01-01", tz = "America/Los_Angeles")),
    by = "rounded_datetime"
  ) %>%
  arrange(datetime) %>%
  # Step 3: Ensure all 15-minute intervals are present
  mutate(
    date_only = as.Date(datetime),
    datetime_15min = floor_date(datetime, unit = "15 minutes")
  ) %>%
  filter(!is.na(datetime_15min)) %>% # Ensure no NAs in datetime_15min
  group_by(date_only) %>%
  complete(datetime_15min = seq(min(datetime_15min), max(datetime_15min), by = "15 min")) %>% # Ensure valid "by" argument
  ungroup() %>%
  # Step 4: Fill missing 15-minute values with hourly averages
  mutate(
    DO.obs = ifelse(is.na(do.obs), hourly_DO.obs, do.obs),
    DO.sat = ifelse(is.na(DO.sat), hourly_DO.sat, DO.sat),
    temp.water = ifelse(is.na(wtr), hourly_temp.water, wtr),
    depth = ifelse(is.na(depth), hourly_depth, depth),
    discharge = ifelse(is.na(dischargeCMS), hourly_discharge, dischargeCMS),
    light = ifelse(is.na(light), hourly_light, light)
  ) %>%
  select(-c(hourly_DO.obs, hourly_DO.sat, hourly_temp.water, hourly_depth, hourly_discharge, hourly_light))

names(dat)

summary(dat)



dat$solar.time <- calc_solar_time(dat$datetime_15min, -119.9399057)

plot(dat$solar.time, dat$datetime_15min)

names(dat)

## Compile Data

dat <- dat%>%
  dplyr::select(solar.time, DO.obs, DO.sat, depth, temp.water, light, discharge) 

## Check for NAs in time series
which(is.na(dat$solar.time)) 
dat<- na.omit((dat)) 

## mle model needs: solar.time, DO.obs, DO.sat, depth, temp.water, light
names(dat)
dat <- dat[,c("solar.time","DO.obs","DO.sat","depth","temp.water","light","discharge")]


# Check for duplicated timestamps
if (any(duplicated(dat$solar.time))) {
  cat("Duplicates detected. Removing duplicates...\n")
  # Remove duplicates
  dat <- dat[!duplicated(dat$solar.time), ]
} else {
  cat("No duplicates detected.\n")
}

# Verify that duplicates are removed
any(duplicated(dat$solar.time))
# 
# 
# # 
# HighQ_test <- subset(dat, solar.time > "2023-05-10 00:00:00" & solar.time < "2023-06-18 20:00:00")
# 
# lowQ_test <- subset(dat, solar.time > "2023-09-15 00:00:00" & solar.time < "2023-10-05 20:00:00")
# 
# 
# lowQ2_test <- subset(dat, solar.time > "2022-08-05 00:00:00" & solar.time < "2022-09-05 20:00:00")

Summer_test <- subset(dat, solar.time > "2022-05-05 00:00:00" & solar.time < "2022-10-20 20:00:00")


# Check timestamp format with light or temp
qplot(solar.time,  DO.obs, data = Summer_test, geom="point") +
  theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))+
  scale_x_datetime(labels = date_format("%m/%d %H:%M"),
                   breaks = date_breaks("180 hours"))


# dat<- Summer_test


# library(dplyr)
# library(tidyverse)
# # Remove flows that violate model assumptions
# dat_cleaned <- dat %>%
#   mutate(date = as.Date(solar.time)) %>%
#   dplyr::group_by(date) %>%
#   # Calculate the daily average discharge
#   mutate(daily_avg_discharge = mean(discharge, na.rm = TRUE)) %>%
#   # Identify days where any discharge is double the daily average
#   mutate(flag_double_discharge = any(discharge > 2.05 * daily_avg_discharge)) %>%
#   # Ungroup to remove grouping
#   dplyr::ungroup() %>%
#   # Filter out flagged days
#   filter(!flag_double_discharge) %>%
#   # Remove intermediate columns
#   dplyr::select(-daily_avg_discharge, -flag_double_discharge)

#####################################################
#Visualize the data    
#####################################################
# dat <- dat_cleaned

dat %>% unitted::v() %>%
  mutate(DO.pctsat = 100 * (DO.obs / DO.sat)) %>%
  dplyr::select(solar.time, starts_with('DO')) %>%
  gather(type, DO.value, starts_with('DO')) %>%
  mutate(units=ifelse(type == 'DO.pctsat', 'DO\n(% sat)', 'DO\n(mg/L)')) %>%
  ggplot(aes(x=solar.time, y=DO.value, color=type)) + geom_line() + 
  facet_grid(units ~ ., scale='free_y') + theme_bw() +
  scale_color_discrete('variable')

labels <- c(depth='depth\n(m)', temp.water='water temp\n(deg C)', light='PAR\n(umol m^-2 s^-1)', discharge='Q\n(cms)')
dat %>% unitted::v() %>%
  dplyr::select(solar.time, depth, temp.water, light, discharge) %>%
  gather(type, value, depth, temp.water, light, discharge) %>%
  mutate(
    type=ordered(type, levels=c('depth','temp.water','light','discharge')),
    units=ordered(labels[type], unname(labels))) %>%
  ggplot(aes(x=solar.time, y=value, color=type)) + geom_line() + 
  facet_grid(units ~ ., scale='free_y') + theme_bw() +
  scale_color_discrete('variable')


dat2<- as.data.frame(dat)

dat3 <- dat2 %>%
  dplyr::select(solar.time, DO.obs, DO.sat, depth, temp.water, light, discharge)

  

##====================================================
## Model the data 
##====================================================

### ** Pause **
### code for binned K600 model:
# ## Set bayes specs
bayes_name_new <- mm_name(type='bayes', pool_K600="binned", err_obs_iid=TRUE, err_proc_iid = TRUE, ode_method = "trapezoid", deficit_src='DO_mod', engine='stan')
bayes_specs_new <- specs(bayes_name_new)
## Based on range of log daily Q (readjust based on range of your discharge, but then remember to reset number of nodes_meanlog and sdlog)
mean(log(dat3$discharge))
median(log(dat3$discharge))
range(log(dat3$discharge))
hist(log(dat3$discharge))
sd(log(dat3$discharge))
##  log-transformed discharge for setting K600 priors
logQ <- log(dat3$discharge)
# look at log(Q) range to define node centers
logQ_range <- range(logQ, na.rm = TRUE)
bayes_specs_new$K600_lnQ_nodes_centers <- seq(logQ_range[1]+1.40, logQ_range[2], length.out = 6)
## Based on Pete Raymond's data for small headwater streams
## (might leave at default values but make sure to adjust number of nodes)
bayes_specs_new$K600_lnQ_nodes_meanlog <- c(rep(-3.83, 6))
bayes_specs_new$K600_lnQ_nodes_sdlog <- c(rep(1.56, 6))
## Change sigma if need to constrain var
bayes_specs_new$K600_daily_sigma_sigma <- 0.05
# bayes_specs_new$GPP_daily_lower <- c(0)
# bayes_specs_new$ER_daily_upper <- c(0)

# added stuff
bayes_specs_new$n_chains <- c(4)
bayes_specs_new$n_cores <- c(4)
bayes_specs_new$burnin_steps <- c(5000)
bayes_specs_new$saved_steps <- c(2500)

bayes_specs_new$stan_control <- list(adapt_delta = 0.99, max_treedepth = 18)


### ** Pause **
### code for implementing measured priors on K600 
## Compute log-transformed discharge for setting K600 priors
logQ <- log(dat3$discharge)

# Inspect log(Q) range to define node centers
logQ_range <- range(logQ, na.rm = TRUE)
hist(logQ, breaks = 20, main = "Histogram of log(Q)")

# Set bayes specs
bayes_name_new <- mm_name(type = 'bayes', pool_K600 = "binned",
                          err_obs_iid = TRUE, err_proc_iid = TRUE,
                          ode_method = "trapezoid", deficit_src = 'DO_mod', engine = 'stan')
bayes_specs_new <- specs(bayes_name_new)

## # Customize Stan settings
bayes_specs_new$stan_control <- list(adapt_delta = 0.99, max_treedepth = 18)
# defaults at adapt_delta = 0.90, max_treedepth = 10
#
bayes_specs_new$n_chains <- 4
bayes_specs_new$n_cores <- 4
bayes_specs_new$burnin_steps <- 5000
bayes_specs_new$saved_steps <- 2500

## Could get less certain ranges for the upper and lower limits
#  147.03x + 0.8307
logK600 = c((147.03*(dat3$discharge)) + 0.8307) #y = mx+b from measured flow and K
logK600_range <- range(logK600, na.rm = TRUE)
median(logK600)
hist(logK600, breaks = 20, main = "Histogram of log(K600)")
summary(logK600)
sd(logK600)
# plot(dat$discharge, logK)

# Adjust priors for K600
bayes_specs_new$K600_lnQ_nodes_centers <- seq(logQ_range[1]+1.85, logQ_range[2], length.out = 6)
bayes_specs_new$K600_lnQ_nodes_meanlog <- c(4.86, # median
                                            4.35, # measured
                                            3.65, # measured
                                            3.08, # measured
                                            2.99, # 1 quartile
                                            0.90) # 1 min
bayes_specs_new$K600_lnQ_nodes_sdlog <- rep(0.5, 6)
bayes_specs_new$K600_daily_sigma_sigma <- 0.05
# 
# # Constrain GPP and ER values
# # bayes_specs_new$GPP_daily_lower <- 0
# # bayes_specs_new$ER_daily_upper <- 0

##==============
# Run the model
##==============
dat_metab_GBL <- metab(bayes_specs_new, data = dat3)

# extract fit:
dat_fit_GBL <- get_fit(dat_metab_GBL)

## Visualize
DOplot <-plot_DO_preds((dat_metab_GBL))
DOplot
# 
ggsave(plot = DOplot, filename = paste("R:/Users/kloria/Documents/Stream_Metab_24/24_metabolism_output/Full_figure/DO_preds_GBL_bayesmodel_mK600_filter_250106.png",sep=""),width=10,height=4,dpi=300)

metabplot<- plot_metab_preds(predict_metab(dat_metab_GBL))
# 
ggsave(plot = metabplot, filename = paste("R:/Users/kloria/Documents/Stream_Metab_24/24_metabolism_output/Full_figure/metabplot_GBL_bayesmodel_mK600_filter_250106.png",sep=""),width=10,height=4,dpi=300)


## Check binning
Binning <- function(fit_Site, Site){
  SM_output <- fit_Site$daily
  SM_day <- get_data_daily(Site)
  SM_KQbin <- fit_Site$KQ_binned
  SM_specs <- get_specs(Site)
  
  day <- data.frame(SM_day$discharge.daily, SM_output$K600_daily_50pct, rep('daily', dim(SM_output)[1]))
  colnames(day)<-c('Q', 'K600', 'Group')
  
  nodes<-data.frame(exp(as.numeric(as.character(SM_specs$K600_lnQ_nodes_centers))), exp(SM_KQbin$lnK600_lnQ_nodes_50pct), rep('node', dim(SM_KQbin)[1]))
  colnames(nodes)<-c('Q', 'K600', 'Group')
  KQ_plot<-rbind(day,nodes)
  
  ggplot(data=KQ_plot, aes(x=log(Q), y=K600, group=Group, colour=Group)) + 
    geom_point(size=3) +
    #geom_line() + 
    scale_color_manual(name="K-Q",
                       breaks = c("daily", "node"),
                       values=c("grey", "purple"),
                       labels=c("Daily","Bin")) +
    ylab("K600") +
    xlab("logQ") +
    theme_bw() +
    theme(legend.position = "top")
}

binplot<- Binning(dat_fit_GBL, dat_metab_GBL)

ggsave(plot = binplot, filename = paste("R:/Users/kloria/Documents/Stream_Metab_24/24_metabolism_output/Full_figure/binplot_GBL_bayesmodel_mK600_filter_250106.png",sep=""),width=5,height=4,dpi=300)

# ggsave(plot = binplot, filename = paste("R:/Users/kloria/Documents/2023_StreamMetab/figures/binplot_GBL_K600_measured_prior_filter_241225.png",sep=""),width=5,height=4,dpi=300)


get_fit(dat_metab_GBL)$overall %>%
  dplyr::select(ends_with('Rhat')) # might be best rhat 

## optional save details
site <- "GBL"
rundate <- format(Sys.Date(), "%y%m%d")
runmeta_k600 <- "K600_m_prior_filter"
runmeta_light <- "NLDAS"


## Save info:
writefiles <- function(data, data2, path = "R:/Users/kloria/Documents/Stream_Metab_24/24_metabolism_output/") {
  for (i in seq_along(data)) {
    filename = paste(path,site,"_",runmeta_k600,"_",runmeta_light,"_",names(data)[i],rundate,".csv", sep = "")
    write.csv(data[[i]], filename)
  }
  
  write.csv(unlist(get_specs(data2)), paste(path,site,"_",runmeta_k600,"_",runmeta_light,"_",rundate,"specs.csv", sep = ""))
  write.csv(get_data_daily(data2), paste(path,site,"_",runmeta_k600,"_",runmeta_light,"_",rundate,"datadaily.csv", sep = ""))
  write.csv(get_data(data2), paste(path,site,"_",runmeta_k600,"_",runmeta_light,"_","mod_and_obs.csv", sep = ""))
}

## Create new folder for site and write csv info
writefiles(dat_fit_GBL, dat_metab_GBL)


getwd()


# saveRDS(dat_metab_GBL, file = "R:/Users/kloria/Documents/Stream_Metab_24/24_metabolism_output/Full_out/GBL_Summer_K600_measured_prior_noGPPprior_241225.rds")



###
###

# plot for k600
# R:/Users/kloria/Documents/Stream_Metab_24/24_metabolism_output/BWL_HighFlow_test_pool_K600_binned_NLDAS_241204datadaily.csv
GBL <- read.csv("R:/Users/kloria/Documents/Stream_Metab_24/24_metabolism_output/GBL_K600_est_prior_filter_NLDAS_daily250102.csv")

GBL$date <- as.Date(GBL$date, origin="2021-01-01")
GBL$site <- "GBL"
GBL$shore <- "west"

Odd_plot <- ggplot(GBL, aes(x = K600_daily_mean, color = shore, fill = shore)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 15) +
  scale_fill_manual(values = alpha(c("#3283a8"), 0.2)) +
  labs(subtitle = "K600 measured no GPP/ER prior") +
  scale_color_manual(values = alpha(c("#3283a8"), 0.9)) + theme_bw() + 
  geom_vline(data = GBL, aes(xintercept = mean(na.omit(K600_daily_mean))), linetype = "dashed") 

Gplot_sp <- ggplot(GBL, aes(x = K600_daily_mean, y = (ER_mean*-1))) +  
  geom_point(shape= 17, col = alpha(c("#3283a8"),0.75)) +
  geom_smooth(method ="lm", se=F)+  facet_grid(.~site)+
  theme_bw()


GB_lm <- lm(K600_daily_mean~(ER_mean*-1), data=GBL)
summary(GB_lm)

# Extract R-squared value
r_gsquared <- summary(GB_lm)$r.squared
GB_plot <- Gplot_sp + 
  geom_text(aes(x = mean((GBL$K600_daily_mean), na.rm=T), y = max((GBL$ER_mean*-1), na.rm=T), 
                label = paste("R2=", round(r_gsquared, 3))),
            hjust = 1, vjust = 0, size = 3, col = "blue",
            parse = T, check_overlap = T, na.rm = T)
library(ggpubr)

k_grid <- ggarrange(GB_plot,
                    Odd_plot,
                    ncol = 2, nrow = 1,
                    common.legend = TRUE, 
                    legend = "bottom",
                    widths = c(0.55, 0.45))

#
ggsave(plot = k_grid, filename = paste("R:/Users/kloria/Documents/Stream_Metab_24/24_metabolism_output/Full_out/ER_K600plot_GBL_K600_est_prior_filter_NLDAS_daily250102.png",sep=""),width=6,height=3,dpi=300)

dat3_sum <- dat3%>%
  mutate(date = as.Date(solar.time))%>%
  dplyr::group_by(date)%>%
  dplyr::summarise(
    DO.obs=mean(DO.obs, na.rm=T),
    DO.sat=mean(DO.sat, na.rm=T), 
    depth=mean(depth , na.rm=T), 
    temp.water=mean(temp.water, na.rm=T),
    discharge=mean(discharge, na.rm=T))


str(GBL)


GBL2 <- GBL%>%
  left_join(dat3_sum, by= c("date"))

Gplot_DOsat <- ggplot(GBL2, aes(y = K600_daily_mean, x = DO.sat)) +  
  labs(subtitle = "2021 summer test K600 measured") +
  geom_point(shape= 17, col = alpha(c("#3283a8"),0.75)) +
  #geom_smooth(method ="lm", se=F)+  
  facet_grid(.~site)+
  theme_bw()


Gplot_temp <- ggplot(GBL2, aes(y = K600_daily_mean, x = temp.water)) +  
  geom_point(shape= 17, col = alpha(c("#3283a8"),0.75)) +
  #geom_smooth(method ="lm", se=F)+ 
  facet_grid(.~site)+
  theme_bw()

Gplot_Q <- ggplot(GBL2, aes(y = K600_daily_mean, x = discharge)) +  
  geom_point(shape= 17, col = alpha(c("#3283a8"),0.75)) +
  #geom_smooth(method ="lm", se=F)+  
  facet_grid(.~site)+
  theme_bw()


GB_lmq <- lm(K600_daily_mean~discharge, data=GBL2)
summary(GB_lmq)

# Extract R-squared value
r_gsquaredq <- summary(GB_lmq)$r.squared
Gplot_Q2 <- Gplot_Q + 
  geom_text(aes(y = mean((GBL2$K600_daily_mean), na.rm=T), x = max((GBL2$discharge), na.rm=T), 
                label = paste("R2=", round(r_gsquaredq, 3))),
            hjust = 1, vjust = 0, size = 3, col = "blue",
            parse = T, check_overlap = T, na.rm = T)

Gplot_z <- ggplot(GBL2, aes(y = K600_daily_mean, x = depth)) +  
  geom_point(shape= 17, col = alpha(c("#3283a8"),0.75)) +
  #geom_smooth(method ="lm", se=F)+  
  facet_grid(.~site)+
  theme_bw()

GB_lmz <- lm(K600_daily_mean~depth, data=GBL2)
summary(GB_lmz)
r_gsquaredz <- summary(GB_lmz)$r.squared

Gplot_z2 <- Gplot_z + 
  geom_text(aes(y = mean((GBL2$K600_daily_mean), na.rm=T), x = max((GBL2$depth), na.rm=T), 
                label = paste("R2=", round(r_gsquaredz, 3))),
            hjust = 1, vjust = 0, size = 3, col = "blue",
            parse = T, check_overlap = T, na.rm = T)

k_e_grid <- ggarrange(Gplot_DOsat,
                    Gplot_temp,
                    Gplot_Q2,
                    Gplot_z2,
                    ncol = 1, nrow = 4,
                    common.legend = TRUE, 
                    legend = "bottom")


ggsave(plot = k_e_grid, filename = paste("R:/Users/kloria/Documents/Stream_Metab_24/24_metabolism_output/Full_out/GBL_K600_conditions_prior_filter_NLDAS_daily250102.png",sep=""),width=4,height=12,dpi=300)



