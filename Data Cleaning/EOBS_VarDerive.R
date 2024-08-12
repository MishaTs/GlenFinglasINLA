#########################################################################################################################
#### Climate data - Glen Finglas ####
#### Original Author: Stuart Smith ####
#### Editor: Misha Tseitlin ####
#### Date: 21 May 2024 ####
#########################################################################################################################
#clear system & package libraries
#rm(list=ls())

library(tidyverse)
library(sf)
library(terra)
# pre-empt namespace conflicts
select <- dplyr::select
#########################################################################################################################

#########################################################################################################################
#### Glen Finglas spatial area ####
#########################################################################################################################
plotBounds <- st_read("./Spatial/boundaries/newplotsCleaned.shp")
plotBounds <- st_transform(plotBounds, "EPSG:4326") # reproject to WGS84

#write to shapefile
st_write(obj = plotBounds,
         dsn="./Spatial/DEM/newplots_wgs84.shp")

#Bounding box
bb <- c(xmin = -4.46, xmax = -4.37, ymin = 56.26, ymax = 56.31)
# create bounding box manually to export for visualisations
#bbShp <- st_as_sfc(st_bbox(c(xmin = -4.46, xmax = -4.37, ymin = 56.26, ymax = 56.31), crs = st_crs(4326)))
#st_write(obj = bbShp,
#         dsn="./Spatial/Bounding Box/StudyArea_WGS84.shp")
gfCrop<- st_crop(plotBounds,bb)
ggplot(gfCrop) + geom_sf()

#########################################################################################################################
# E-OBS - European climate data
# https://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php#datafiles
# E-OBS -  0.1 degrees regular grid ensemble mean
# TG daily mean temperature - checked
# TN daily minimum temperature - checked
# TX daily maximum temperature - checked
# RR daily precipitation sum - checked
# PP daily averaged sea level pressure
# HU daily averaged relative humidity
# FG daily mean wind speed - checked
# QQ daily mean global radiation - checked
#########################################################################################################################

###########################################################################
##### E-OBS temporal selection ####
###########################################################################

# In-situ climatic variables
# Total rainfall (mm) since previous NVC monitoring period
# Average wind speed (knots) since previous NVC monitoring period
# Mean DAILY temperature since previous NVC sampling
# Mean MIN temperature since previous NVC sampling
# Mean MAX temperature since previous NVC sampling
# Mean VARiance in temperature since previous NVC sampling - using daily temperature

# Seasonal measures
# Mean autumn-winter temperature since previous NVC sampling- using daily temperature
# Mean spring-summer temperature since previous NVC to sampling (survey within summer period)- using daily temperature

# Relative climate variables based on long-term trends/thresholds
# Growing degree days - seasonal/annual
# Growing degree days above long-term mean (1950-2022) sum of degrees above mean since previous monitoring period
# Calculated for the following seasons and temperature variables
# Spring degree days minimum temperature
# Summer degree days minimum temperature
# Autumn degree days minimum temperature
# Winter degree days minimum temperature
# Spring degree days max temperature
# Summer degree days max temperature
# Autumn degree days max temperature
# Winter degree days max temperature

# Growing degree days above 5.5C sum of degrees above threshold since previous monitoring period
# Growing degree days daily temperature
# Growing degree days minimum temperature
# Frost days - number of days temperatures <0C - use min temperature 

# Other calculations not used
# Heatwaves - max temperature >25C for 3 consecutive days (Met office heatwave) - too few occurrences
# Very Warm Days (VWD) - number of occurrences 3 consecutive days above 95% quantile for spring/summer/winter - use max temp long-term quantile 1950-2022
# Too few occurrences

# Possible environmental variables?
# Evapo-transpiration - how to calculate? 
# Snow line- cannot find data source? 

# Seasonal dates - UK met office https://www.metoffice.gov.uk/weather/learn-about/weather/seasons/
# meteorological (not astronomical) seasons
# Spring March - May
# Summer June - August
# Autumn September - November
# Winter December - February

# Import processed climate data
TG<-read.csv("./Climate/EOBS/tg.csv", header=TRUE, sep=',',strip.white=T, dec='.',na.strings=c("","NA"))
TN<-read.csv("./Climate/EOBS/tn.csv", header=TRUE, sep=',',strip.white=T, dec='.',na.strings=c("","NA"))
TX<-read.csv("./Climate//EOBS/tx.csv", header=TRUE, sep=',',strip.white=T, dec='.',na.strings=c("","NA"))
TGfull<-read.csv("./Climate//EOBS/tg_full.csv", header=TRUE, sep=',',strip.white=T, dec='.',na.strings=c("","NA"))
TXfull<-read.csv("./Climate//EOBS/tx_full.csv", header=TRUE, sep=',',strip.white=T, dec='.',na.strings=c("","NA"))
TNfull<-read.csv("./Climate//EOBS/tn_full.csv", header=TRUE, sep=',',strip.white=T, dec='.',na.strings=c("","NA"))
RR<-read.csv("./Climate//EOBS/rr.csv", header=TRUE, sep=',',strip.white=T, dec='.',na.strings=c("","NA"))
RRfull<-read.csv("./Climate//EOBS/rr_full.csv", header=TRUE, sep=',',strip.white=T, dec='.',na.strings=c("","NA"))
FG<-read.csv("./Climate//EOBS/fg.csv", header=TRUE, sep=',',strip.white=T, dec='.',na.strings=c("","NA"))
FGfull<-read.csv("./Climate//EOBS/fg_full.csv", header=TRUE, sep=',',strip.white=T, dec='.',na.strings=c("","NA"))
QQ<-read.csv("./Climate//EOBS/qq.csv", header=TRUE, sep=',',strip.white=T, dec='.',na.strings=c("","NA"))
QQfull<-read.csv("./Climate//EOBS/qq_full.csv", header=TRUE, sep=',',strip.white=T, dec='.',na.strings=c("","NA"))
HU<-read.csv("./Climate//EOBS/hu.csv", header=TRUE, sep=',',strip.white=T, dec='.',na.strings=c("","NA"))
HUfull<-read.csv("./Climate//EOBS/hu_full.csv", header=TRUE, sep=',',strip.white=T, dec='.',na.strings=c("","NA"))
PP<-read.csv("./Climate//EOBS/pp.csv", header=TRUE, sep=',',strip.white=T, dec='.',na.strings=c("","NA"))
PPfull<-read.csv("./Climate//EOBS/pp_full.csv", header=TRUE, sep=',',strip.white=T, dec='.',na.strings=c("","NA"))

# harmonise spatially
# originally up to 6 decimal points (usually), but they're all equal at 4
# so, round everything to 4 to standardise
TGmerge <- TGfull %>% mutate(xS = round(x,3),
                             yS = round(y,4))
TXmerge <- TXfull %>% mutate(xS = round(x,3),
                             yS = round(y,4))
TNmerge <- TNfull %>% mutate(xS = round(x,3),
                             yS = round(y,4))
RRmerge <- RRfull %>% mutate(xS = round(x,3),
                             yS = round(y,4))
FGmerge <- FGfull %>% mutate(xS = round(x,3),
                             yS = round(y,4))
QQmerge <- QQfull %>% mutate(xS = round(x,3),
                             yS = round(y,4))
HUmerge <- HUfull %>% mutate(xS = round(x,3),
                             yS = round(y,4))
PPmerge <- PPfull %>% mutate(xS = round(x,3),
                             yS = round(y,4))
fullClim <- full_join(QQmerge, TGmerge, by = c("xS", "yS", "Date")) %>% 
  full_join(., TXmerge, by = c("xS", "yS", "Date")) %>% 
  full_join(., TNmerge, by = c("xS", "yS", "Date")) %>% 
  full_join(., RRmerge, by = c("xS", "yS", "Date")) %>% 
  full_join(., FGmerge, by = c("xS", "yS", "Date")) %>% 
  full_join(., HUmerge, by = c("xS", "yS", "Date")) %>% 
  full_join(., PPmerge, by = c("xS", "yS", "Date")) %>% 
  rename(charDate = "Date") %>% 
  mutate(Date = as.Date(charDate),
         Yr = year(Date)) %>% 
  select(c("xS", "yS", "Date", "Yr",
           "mean_temp_C", "max_temp_C", "min_temp_C", 
           "precip", "wind_speed", "solar_radiation", "humid", "pressure")) %>% 
  rename(x = "xS",
         y = "yS")

# confirm that everything passed through OK -- it did
#sum(!is.na(fullClim$mean_temp_C)) - nrow(TGmerge)
#sum(!is.na(fullClim$max_temp_C)) - nrow(TXmerge)
#sum(!is.na(fullClim$min_temp_C)) - nrow(TNmerge)
#sum(!is.na(fullClim$precip)) - nrow(RRmerge)
#sum(!is.na(fullClim$humid)) - nrow(HUmerge)
#sum(!is.na(fullClim$pressure)) - nrow(PPmerge)
#sum(!is.na(fullClim$wind_speed)) - nrow(FGmerge)
#sum(!is.na(fullClim$solar_radiation)) - nrow(QQmerge)
# slight row mismatch
nrow(fullClim) - nrow(QQmerge)
# this is because solar radiation is missing some observations in 2020 (over 5 days in April, May, and June)
#table(fullClim %>% filter(is.na(solar_radiation)) %>% select(Date))

# make this wide with data for each block group
# climate data agrees on the cell distribution of E&F blocks except for wind speed which uses slightly different cells
# rather than stick too closely to the actual values, we'll assume here that E&F are 13% in a new cell and 87% in the A&B cell
# this will accentuate contrast and improve separability in calculations
# see lines 370-400 in "/JHIGit/Climate/GF_EOBS_DataImport.R" for the full workings
fullClim <- fullClim %>% mutate(spatID = paste(x, y, sep = "_"),
                                # A&B: -4.383_56.2749 * 100%
                                # C&D: -4.45_56.2749 * 100%
                                # E&F: -4.383_56.2749 * 87% + -4.383_56.3249 * 13%
                                block = factor(spatID, 
                                              #levels = ,
                                              labels = c("AB", "EFRaw", "CD")))

# easier to handle this with a wide dataset
climWide <- fullClim %>% select(-c("x","y","spatID")) %>% 
  pivot_wider(
  names_from = block,
  names_glue = "{.value}Plot{block}",
  values_from = colnames(fullClim)[5:12])

# not efficient but decently readable
# it takes a moment to run
climWide <- climWide %>% rowwise() %>% mutate(mean_temp_CPlotEF = ifelse(is.na(mean_temp_CPlotEFRaw), # if data missing from upper cell
                                                                         mean_temp_CPlotAB, # then, use the lower cell only
                                                                         # otherwise, take the weighted mean
                                                                         weighted.mean(c_across(c(mean_temp_CPlotAB, 
                                                                                                  mean_temp_CPlotEFRaw)),
                                                                                       # heavily weight AB
                                                                                       # slightly weight raw EF
                                                                                       c(0.87, 0.13))),
                                              # repeat for all other variables
                                              max_temp_CPlotEF = ifelse(is.na(max_temp_CPlotEFRaw), 
                                                                        max_temp_CPlotAB,
                                                                        weighted.mean(c_across(c(max_temp_CPlotAB, 
                                                                                                 max_temp_CPlotEFRaw)),
                                                                                      c(0.87, 0.13))),
                                              min_temp_CPlotEF = ifelse(is.na(min_temp_CPlotEFRaw), 
                                                                        min_temp_CPlotAB,
                                                                        weighted.mean(c_across(c(min_temp_CPlotAB, 
                                                                                                 min_temp_CPlotEFRaw)),
                                                                                      c(0.87, 0.13))),
                                              precipPlotEF = ifelse(is.na(precipPlotEFRaw), 
                                                                    precipPlotAB,
                                                                    weighted.mean(c_across(c(precipPlotAB, 
                                                                                             precipPlotEFRaw)),
                                                                                  c(0.87, 0.13))),
                                              wind_speedPlotEF = ifelse(is.na(wind_speedPlotEFRaw), 
                                                                        wind_speedPlotAB,
                                                                        weighted.mean(c_across(c(wind_speedPlotAB, 
                                                                                                 wind_speedPlotEFRaw)),
                                                                                      c(0.87, 0.13))),
                                              solar_radiationPlotEF = ifelse(is.na(solar_radiationPlotEFRaw), 
                                                                             solar_radiationPlotAB,
                                                                             weighted.mean(c_across(c(solar_radiationPlotAB, 
                                                                                                      solar_radiationPlotEFRaw)),
                                                                                           c(0.87, 0.13))),
                                              humidPlotEF = ifelse(is.na(humidPlotEFRaw), 
                                                                   humidPlotAB,
                                                                   weighted.mean(c_across(c(humidPlotAB, 
                                                                                            humidPlotEFRaw)),
                                                                                 c(0.87, 0.13))),
                                              pressurePlotEF = ifelse(is.na(pressurePlotEFRaw), 
                                                                      pressurePlotAB,
                                                                      weighted.mean(c_across(c(pressurePlotAB, 
                                                                                               pressurePlotEFRaw)),
                                                                                    c(0.87, 0.13))))
# remove raw helper column
climWide <- climWide %>% select(-ends_with("EFRaw"))
# convert back to long format
climCalc <- climWide %>% 
  # first step, get every variable and block into its own row
  pivot_longer(
    cols = !c("Date","Yr"),
    names_to = c("var", "block"),
    names_sep = "Plot",
    values_to = "val") %>% 
  # second step, widen back out to give each variable its own column
  pivot_wider(
    names_from = var,
    values_from = val
  ) %>% 
  # rudimentary season attribution by meteorological seaason (not astronomical)
  mutate(Mth = month(Date),
         Season = ifelse(Mth < 3 | Mth == 12, "Winter", 
                         ifelse(Mth < 6, "Spring",
                                ifelse(Mth < 9, "Summer", 
                                       "Autumn"))))

# so what's happening with missing wind speed?
# no clear seasonality or bias from visual inspection
ggplot(data = climWide %>% filter(Yr >= 1980), aes(x = Date, y = wind_speedPlotAB)) + 
  geom_line() + facet_wrap(~Yr, scales = "free_x")
#library(forecast)
#windSpeedAB <- auto.arima(climWide$wind_speedPlotAB, seasonal = TRUE,
#                          stepwise = FALSE, approximation = FALSE, 
                          #parallel = TRUE, 
#                          trace = TRUE)
climTest <- climWide %>% mutate(Mth = month(Date),
                                Season = ifelse(Mth < 3 | Mth == 12, "Winter", 
                                                ifelse(Mth < 6, "Spring",
                                                       ifelse(Mth < 9, "Summer", 
                                                              "Autumn"))))
# no significant difference between June, July and August within a season: this is fine!
# but...there is difference between seasons so make sure to weight them equally
windFactors <- lm(wind_speedPlotAB ~ relevel(as.factor(Mth), ref = 6), data = climTest)
summary(windFactors)

# do some summarising by season 
# just tendency, no built-in uncertainty for now
climGF <- climCalc %>% 
  # rudimentary for now, but assume all surveying concludes by August
  # only untrue for NVC samples in 2006, 2007, 2008, 2011, 2012, and 2023: never past August though
  # Pinframe samples are weird -- they start as early as Jan and end as late as August
  mutate(sampleYr = ifelse(Mth < 9, Yr, Yr + 1)) %>% 
  # autumn "2002" is actually Sep - Nov 2001, so this is fine
  filter(sampleYr >= 2002) %>% group_by(sampleYr, block, Season) %>% 
  summarise(
    # average mean temperature
    meanT = mean(mean_temp_C, na.rm = TRUE),
    # absolute max temperature
    maxT = max(max_temp_C, na.rm = TRUE),
    # absolute min temperature
    minT = min(min_temp_C, na.rm = TRUE),
    # get the spread of temperature values
    rangeT = maxT - minT,
    # average rainfall
    precip = mean(precip, na.rm = TRUE),
    # average wind speed
    windSpeed = mean(wind_speed, na.rm = TRUE),
    # average solar radiation
    solarRad = mean(solar_radiation, na.rm = TRUE),
    # average humidity
    humid = mean(humid, na.rm = TRUE),
    # average pressure
    pressure = mean(pressure, na.rm = TRUE),
    # Growing degree days (>5.5C) - sum of degrees above threshold - daily temp
    growDaysMin = sum(mean_temp_C > 0, na.rm = TRUE),
    # Growing degree days (>5.5C) - sum of degrees above threshold - min temp
    growDays = sum(min_temp_C > 5.5, na.rm = TRUE),
    #Frost days - number of days temperatures <0C - use min temperature 
    frostDays = sum(min_temp_C < 0, na.rm = TRUE))
  
# prepare for exporting
climExport <- climGF %>% 
  pivot_wider(
    names_from = Season,
    names_glue = "{.value}_{Season}",
    values_from = colnames(climGF)[4:15]
  ) %>% 
  # derive annual measures too
  # months have slightly different number of days, but this weights all seasons equally
  mutate(meanT = mean(c_across(starts_with("meanT")), na.rm = TRUE),
         maxT = max(c_across(starts_with("maxT")), na.rm = TRUE),
         minT = min(c_across(starts_with("minT")), na.rm = TRUE),
         rangeT = maxT - minT,
         precip = mean(c_across(starts_with("precip")), na.rm = TRUE),
         windSpeed = mean(c_across(starts_with("windSpeed")), na.rm = TRUE),
         solarRad = mean(c_across(starts_with("solarRad")), na.rm = TRUE),
         humid = mean(c_across(starts_with("humid")), na.rm = TRUE),
         pressure = mean(c_across(starts_with("pressure")), na.rm = TRUE),
         growDaysMin = sum(c_across(starts_with("growDaysMin")), na.rm = TRUE),
         growDays = sum(c_across(starts_with("growDays")), na.rm = TRUE),
         frostDays = sum(c_across(starts_with("frostDays")), na.rm = TRUE)) %>% 
  filter(sampleYr < 2024)

write_rds(climExport, "./Climate/climBasic.rds")



# if we want to do sampling time-specific stuff, then NVC data is needed
# Import NVC sampling dates
NVC <- read_rds("./Vegetation/vegCleaned.rds")

# Climate variable date columns
NVCsampling <- NVC %>% group_by(Year, Block) %>% 
  summarise(firstNVCDay = min(nvcDatFull, na.rm = TRUE),
            # last day is roughly either the same as or less than a month after the first sampling date
            lastNVCDay = max(nvcDatFull, na.rm = TRUE),
            firstPinDay = min(pinDat, na.rm = TRUE),
            lastPinDay = max(pinDat, na.rm = TRUE)) %>% 
  # remove Inf and -Inf that occur with null values
  mutate(across(everything(),
                ~ifelse(. == Inf | . == -Inf, NA, .))) %>% 
  # reconvert to Date format
  mutate(firstNVCDay = as.Date(firstNVCDay),
         lastNVCDay = as.Date(lastNVCDay),
         firstPinDay = as.Date(firstPinDay),
         lastPinDay = as.Date(lastPinDay),
         daysNVC = as.numeric(lastNVCDay - firstNVCDay) + 1,
         daysPin = as.numeric(lastPinDay - firstPinDay) + 1)

# get climate variables for the survey periods
# these are averages to minimise any bias from uneven survey periods (between 215 and 1 day)
climSample <- NVCsampling %>% rowwise() %>% 
  mutate(nvcMeanT = mean(as.vector(climCalc %>% 
                                     filter(Date >= (firstNVCDay - 1) & 
                                              Date <= lastNVCDay & 
                                              str_detect(block, Block)) %>% 
                                     select(mean_temp_C))$mean_temp_C, 
                         na.rm = TRUE),
         nvcMaxT = mean(as.vector(climCalc %>% 
                                    filter(Date >= (firstNVCDay - 1) & 
                                             Date <= lastNVCDay & 
                                             str_detect(block, Block)) %>% 
                                    select(max_temp_C))$max_temp_C, 
                        na.rm = TRUE),
         nvcMinT = mean(as.vector(climCalc %>% 
                                    filter(Date >= (firstNVCDay - 1) & 
                                             Date <= lastNVCDay & 
                                             str_detect(block, Block)) %>% 
                                    select(min_temp_C))$min_temp_C, 
                        na.rm = TRUE),
         nvcPrecip = mean(as.vector(climCalc %>% 
                                      filter(Date >= (firstNVCDay - 1) & 
                                               Date <= lastNVCDay & 
                                               str_detect(block, Block)) %>% 
                                      select(precip))$precip, 
                          na.rm = TRUE),
         # NA for 2023 -- can be imputed with ARIMA if needed
         nvcWindSpeed = mean(as.vector(climCalc %>% 
                                         filter(Date >= (firstNVCDay - 1) & 
                                                  Date <= lastNVCDay & 
                                                  str_detect(block, Block)) %>% 
                                         select(wind_speed))$wind_speed, 
                             na.rm = TRUE),
         nvcSolarRad = mean(as.vector(climCalc %>% 
                                        filter(Date >= (firstNVCDay - 1) & 
                                                 Date <= lastNVCDay & 
                                                 str_detect(block, Block)) %>% 
                                        select(solar_radiation))$solar_radiation, 
                            na.rm = TRUE),
         nvcHumid = mean(as.vector(climCalc %>% 
                                     filter(Date >= (firstNVCDay - 1) & 
                                              Date <= lastNVCDay & 
                                              str_detect(block, Block)) %>% 
                                     select(humid))$humid, 
                         na.rm = TRUE),
         nvcPressure = mean(as.vector(climCalc %>% 
                                        filter(Date >= (firstNVCDay - 1) & 
                                                 Date <= lastNVCDay & 
                                                 str_detect(block, Block)) %>% 
                                        select(pressure))$pressure, 
                            na.rm = TRUE),
         pinMeanT = mean(as.vector(climCalc %>% 
                                     filter(Date >= (firstPinDay - 1) & 
                                              Date <= lastPinDay & 
                                              str_detect(block, Block)) %>% 
                                     select(mean_temp_C))$mean_temp_C, 
                         na.rm = TRUE),
         pinMaxT = mean(as.vector(climCalc %>% 
                                    filter(Date >= (firstPinDay - 1) & 
                                             Date <= lastPinDay & 
                                             str_detect(block, Block)) %>% 
                                    select(max_temp_C))$max_temp_C, 
                        na.rm = TRUE),
         pinMinT = mean(as.vector(climCalc %>% 
                                    filter(Date >= (firstPinDay - 1) & 
                                             Date <= lastPinDay & 
                                             str_detect(block, Block)) %>% 
                                    select(min_temp_C))$min_temp_C, 
                        na.rm = TRUE),
         pinPrecip = mean(as.vector(climCalc %>% 
                                      filter(Date >= (firstPinDay - 1) & 
                                               Date <= lastPinDay & 
                                               str_detect(block, Block)) %>% 
                                      select(precip))$precip, 
                          na.rm = TRUE),
         # NA for 2023 -- can be imputed with ARIMA if needed
         pinWindSpeed = mean(as.vector(climCalc %>% 
                                         filter(Date >= (firstPinDay - 1) & 
                                                  Date <= lastPinDay & 
                                                  str_detect(block, Block)) %>% 
                                         select(wind_speed))$wind_speed, 
                             na.rm = TRUE),
         pinSolarRad = mean(as.vector(climCalc %>% 
                                        filter(Date >= (firstPinDay - 1) & 
                                                 Date <= lastPinDay & 
                                                 str_detect(block, Block)) %>% 
                                        select(solar_radiation))$solar_radiation, 
                            na.rm = TRUE),
         pinHumid = mean(as.vector(climCalc %>% 
                                     filter(Date >= (firstPinDay - 1) & 
                                              Date <= lastPinDay & 
                                              str_detect(block, Block)) %>% 
                                     select(humid))$humid, 
                         na.rm = TRUE),
         pinPressure = mean(as.vector(climCalc %>% 
                                        filter(Date >= (firstPinDay - 1) & 
                                                 Date <= lastPinDay & 
                                                 str_detect(block, Block)) %>% 
                                        select(pressure))$pressure, 
                            na.rm = TRUE)) %>% 
  # remove NAs but leave out date columns since ifelse strips as.Date() formatting
  mutate(across(daysNVC:pinPressure,
                ~ifelse(is.nan(.), NA, .)))
# export in case it's needed
write_rds(climSample, "./Climate/climSample.rds")

###########################################################################
##### END ####
###########################################################################