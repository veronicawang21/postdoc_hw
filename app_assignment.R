
##################################################################
### Post-doc application assignment: wildfire PM2.5 and Covid ####
##################################################################

# clear environment
rm(list=ls())

# import libraries needed
library(tidyverse)
library(lubridate)
library(data.table)
library(remotes)
install_github("nyiuab/NBZIMM", force=T, build_vignettes=F)
library(NBZIMM)
library(splines)

# read in data
firecovid <- read.csv("moddat_Feb2021.csv")

# Without a given data dictionary, I will assume that the following:
# 'FIPS': US FIPS county code
# 'date': Date
# 'cases': Number of COVID cases
# 'cumu_cases': uncertain when the start date is -- so will not use
# 'deaths': Number of COVID deaths
# 'cumu_deaths': uncertain when the start date is -- so will not use
# 'population': County population
# 'pm25': Daily average PM2.5
# 'Long': Longitude of County centroid
# 'Lat' : Longitude of County centroid
# 'month': Month of year extracted from 'date'
# 'day': Day of Month extracted from 'date'
# 'md': month-date
# 'County': County name
# 'State': State abbreviation
# 'StateFIPS': US FIPS state code
# 'dayofweek': day of the week
# 'date_num': number of days since 2020-03-15
# 'date_str': date as string
# 'pm25_history': median of 2016-2019 pm
# 'wildfire': whether it was a wildfire day
# 'pm_wildfire': pm25-pm_ambient
# 'pm_ambient': either pm25 or pm25_history, whichever has a higher concentration
# 
# 
# it is unclear what the following variables represented, and therefore was not used:
# 'cumu_cases', 'cumu_deaths', 'tmmx', 'rmax', 'sph', 'hazardmap', 
# 'relative_change_feb', 'ratio_travelers', 'pm25_history_raw', 'pm25_raw'


##################################################################
################ Data prep #######################################
##################################################################

# aggregate data to state level and create exposure histories
agg_dat <- firecovid %>%
  dplyr::select(FIPS:cases, deaths, population, pm25, 
         County:date_str, wildfire, pm_wildfire, pm_ambient) %>% # keep columns needed
  na.omit %>%                                                    # death counts are unavailable for some counties (remove these counties) 
  mutate(dat = as.Date(date)) %>%
  as.data.table()

## create exposure history for pm2.5
cols = c("pm25")
agg_firecovid <- agg_dat
anscols = paste(cols, "lag", 1:28, sep="_")
agg_pm <- agg_firecovid[, (anscols) := shift(.SD, 1:28, "lag"), .SDcols=cols, by=c("State", "County")]

## create exposure history for wildfire pm2.5
cols = c("pm_wildfire")
agg_firecovid <- agg_dat
anscols = paste(cols, "lag", 1:28, sep="_")
agg_pm_fire <- agg_firecovid[, (anscols) := shift(.SD, 1:28, "lag"), .SDcols=cols, by=c("State", "County")]

## create exposure history for ambient pm
cols = c("pm_ambient")
agg_firecovid <- agg_dat
anscols = paste(cols, "lag", 1:28, sep="_")
agg_pm_amb <- agg_firecovid[, (anscols) := shift(.SD, 1:28, "lag"), .SDcols=cols, by=c("State", "County")]

## merge data
agg_hist <- agg_dat %>%
  left_join(agg_pm) %>%
  left_join(agg_pm_fire) %>%
  left_join(agg_pm_amb) %>%
  na.omit() %>%
  as.data.frame() %>%
  rowwise() %>%
  mutate(pm_week1 = mean(pm25_lag_1:pm25_lag_7),
         pm_week2 = mean(pm25_lag_8:pm25_lag_14),
         pm_week3 = mean(pm25_lag_15:pm25_lag_21),
         pm_week4 = mean(pm25_lag_22:pm25_lag_28),
         pm_fire_week1 = mean(pm_wildfire_lag_1:pm_wildfire_lag_7),
         pm_fire_week2 = mean(pm_wildfire_lag_8:pm_wildfire_lag_14),
         pm_fire_week3 = mean(pm_wildfire_lag_15:pm_wildfire_lag_21),
         pm_fire_week4 = mean(pm_wildfire_lag_22:pm_wildfire_lag_28),
         pm_amb_week1 = mean(pm_ambient_lag_1:pm_ambient_lag_7),
         pm_amb_week2 = mean(pm_ambient_lag_8:pm_ambient_lag_14),
         pm_amb_week3 = mean(pm_ambient_lag_15:pm_ambient_lag_21),
         pm_amb_week4 = mean(pm_ambient_lag_22:pm_ambient_lag_28))

# create dataframe by state
## California
agg_ca <- agg_hist %>%
  filter(State=="CA")

## Oregon
agg_or <- agg_hist %>%
  filter(State=="OR")

## Washington
agg_wa <- agg_hist %>%
  filter(State=="WA")

##################################################################
################ Analysis ########################################
##################################################################

# Zero-inflated negative binomial models

## for PM2.5
form_pm_cases <- "cases ~ pm_week1 + pm_week2 + pm_week3 + pm_week4 + ns(date_num, 4) + offset(log(population))"
form_pm_deaths <- "deaths ~ pm_week1 + pm_week2 + pm_week3 + pm_week4 + ns(date_num, 4) + offset(log(population))"

## for wildfire PM2.5
form_pmfire_cases <- "cases ~ pm_week1 + pm_week2 + pm_week3 + pm_week4 + ns(date_num, 4) + offset(log(population))"
form_pmfire_deaths <- "deaths ~ pm_week1 + pm_week2 + pm_week3 + pm_week4 + ns(date_num, 4) + offset(log(population))"

## for ambient PM2.5
form_pmamb_cases <- "cases ~ pm_week1 + pm_week2 + pm_week3 + pm_week4 + ns(date_num, 4) + offset(log(population))"
form_pmamb_deaths <- "deaths ~ pm_week1 + pm_week2 + pm_week3 + pm_week4 + ns(date_num, 4) + offset(log(population))"

## California
######################################

### Cases
#### PM2.5
mod_ca_cases_pm <- glmm.zinb(as.formula(form_pmfire_cases), data = agg_ca, random = ~ 1|County, zi_fixed = ~1, zi_random =NULL)
summary(mod_ca_cases_pm)

#### wildfire PM2.5
mod_ca_cases_pmfire <- glmm.zinb(as.formula(form_pmfire_cases), data = agg_ca, random = ~ 1|County, zi_fixed = ~1, zi_random =NULL)
summary(mod_ca_cases_pmfire)

#### ambient PM2.5
mod_ca_cases_pmamb <- glmm.zinb(as.formula(form_pmamb_cases), data = agg_ca, random = ~ 1|County, zi_fixed = ~1, zi_random =NULL)
summary(mod_ca_cases_pmamb)

### Deaths
#### PM2.5
mod_ca_deaths_pm <- glmm.zinb(as.formula(form_pmfire_deaths), data = agg_ca, random = ~ 1|County, zi_fixed = ~1, zi_random =NULL)
summary(mod_ca_deaths_pm)

#### wildfire PM2.5
mod_ca_deaths_pmfire <- glmm.zinb(as.formula(form_pmfire_deaths), data = agg_ca, random = ~ 1|County, zi_fixed = ~1, zi_random =NULL)
summary(mod_ca_deaths_pmfire)

#### ambient PM2.5
mod_ca_deaths_pmamb <- glmm.zinb(as.formula(form_pmamb_deaths), data = agg_ca, random = ~ 1|County, zi_fixed = ~1, zi_random =NULL)
summary(mod_ca_deaths_pmamb)


## Oregon
######################################

### Cases
#### PM2.5
mod_or_cases_pm <- glmm.zinb(as.formula(form_pmfire_cases), data = agg_or, random = ~ 1|County, zi_fixed = ~1, zi_random =NULL)
summary(mod_or_cases_pm)

#### wildfire PM2.5
mod_or_cases_pmfire <- glmm.zinb(as.formula(form_pmfire_cases), data = agg_or, random = ~ 1|County, zi_fixed = ~1, zi_random =NULL)
summary(mod_or_cases_pmfire)

#### ambient PM2.5
mod_or_cases_pmamb <- glmm.zinb(as.formula(form_pmamb_cases), data = agg_or, random = ~ 1|County, zi_fixed = ~1, zi_random =NULL)
summary(mod_or_cases_pmamb)

### Deaths
#### PM2.5
mod_or_deaths_pm <- glmm.zinb(as.formula(form_pmfire_deaths), data = agg_or, random = ~ 1|County, zi_fixed = ~1, zi_random =NULL)
summary(mod_or_deaths_pm)

#### wildfire PM2.5
mod_or_deaths_pmfire <- glmm.zinb(as.formula(form_pmfire_deaths), data = agg_or, random = ~ 1|County, zi_fixed = ~1, zi_random =NULL)
summary(mod_or_deaths_pmfire)

#### ambient PM2.5
mod_or_deaths_pmamb <- glmm.zinb(as.formula(form_pmamb_deaths), data = agg_or, random = ~ 1|County, zi_fixed = ~1, zi_random =NULL)
summary(mod_or_deaths_pmamb)


## Washington
######################################

### Cases
#### PM2.5
mod_wa_cases_pm <- glmm.zinb(as.formula(form_pmfire_cases), data = agg_wa, random = ~ 1|County, zi_fixed = ~1, zi_random =NULL)
summary(mod_wa_cases_pm)

#### wildfire PM2.5
mod_wa_cases_pmfire <- glmm.zinb(as.formula(form_pmfire_cases), data = agg_wa, random = ~ 1|County, zi_fixed = ~1, zi_random =NULL)
summary(mod_wa_cases_pmfire)

#### ambient PM2.5
mod_wa_cases_pmamb <- glmm.zinb(as.formula(form_pmamb_cases), data = agg_wa, random = ~ 1|County, zi_fixed = ~1, zi_random =NULL)
summary(mod_wa_cases_pmamb)

### Deaths
#### PM2.5
mod_wa_deaths_pm <- glmm.zinb(as.formula(form_pmfire_deaths), data = agg_wa, random = ~ 1|County, zi_fixed = ~1, zi_random =NULL)
summary(mod_wa_deaths_pm)

#### wildfire PM2.5
mod_wa_deaths_pmfire <- glmm.zinb(as.formula(form_pmfire_deaths), data = agg_wa, random = ~ 1|County, zi_fixed = ~1, zi_random =NULL)
summary(mod_wa_deaths_pmfire)

#### ambient PM2.5
mod_wa_deaths_pmamb <- glmm.zinb(as.formula(form_pmamb_deaths), data = agg_wa, random = ~ 1|County, zi_fixed = ~1, zi_random =NULL)
summary(mod_wa_deaths_pmamb)
