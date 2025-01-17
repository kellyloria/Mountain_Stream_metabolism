
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




BWL_dat <- readRDS("R:/Users/kloria/Documents/Stream_Metab_24/Core_sites/offset_DO_dat/25_BWL_DO_flag_flow.rds")


