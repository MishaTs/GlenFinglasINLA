#########################################################################################################################
#### Climate data - Glen Finglas ####
#### Original Author: Stuart Smith ####
#### Editor: Misha Tseitlin ####
#### Date: 3 July 2024 ####
#########################################################################################################################
#clear system & package libraries
#rm(list=ls())

library(tidyverse)
library(devtools)
library(raster)
library(sf)
library(terra)
library(exactextractr)

#########################################################################################################################

#########################################################################################################################
#### Glen Finglas spatial area ####
#########################################################################################################################
plotBounds <- st_read("./Spatial/boundaries/newplotsCleaned.shp")
plotBounds <- st_transform(plotBounds, "EPSG:4326") # reproject to WGS84

#write to shapefile
st_write(obj = plotBounds,
         dsn="./Spatial/DTM/newplots_wgs84.shp")

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
#### Process E-OBS climate data ####
#devtools::install_github("RetoSchmucki/climateExtract")
library(climateExtract)
ecad_version # 29.0 newest version 2024-06-19
#########################################################################################################################

#getwd()
# assume that we're working from a local version of the git
# "/[Personal_Path_to_Folder]/GlenFinglas/"
setwd("./Climate/Raw/")

# downloads from online don't work without increasing timeout
#options(timeout = max(1000, getOption("timeout")))

# TG daily mean temperature
TG <- climateExtract::extract_nc_value(first_year = 2000, 
                                       last_year = 2024,
                                       local_file = TRUE, # source data from online (if running for the first time)
                                       file_path = NULL,
                                       spatial_extent = gfCrop,
                                       clim_variable = "mean temp",
                                       statistic = "mean", # getting "spread" for uncertainty might be a good idea too
                                       grid_size = 0.1,
                                       ecad_v = NULL,
                                       write_raster = TRUE,
                                       out = "raster_mean_temp.tiff",
                                       return_data = TRUE)

# original non-terra version
#rbkTG <- raster::brick("raster_mean_temp.grd")

# new terra version
rbkTG <- terra::rast("raster_mean_temp.tiff")
format(object.size(TG), "MB")
format(object.size(rbkTG), "MB") # file would be too large, so this is a pointer to the .nc

#daily_TG<- extract(rbk,GF,sp=T, method='simple')
#write.csv(daily_TG, "GF_TG.csv",row.names=F) #,sep = ",",dec = ".",col.names = TRUE, row.names=F)

# TN daily minimum temperature
TN <- climateExtract::extract_nc_value(first_year = 2000, 
                                       last_year = 2024,
                                       local_file = TRUE, # Source the data from file online = eaiser
                                       file_path = NULL,
                                       spatial_extent = gfCrop,
                                       clim_variable = "min temp",
                                       statistic = "mean",
                                       grid_size = 0.1,
                                       ecad_v = NULL,
                                       write_raster = TRUE,
                                       out = "raster_min_temp.tiff",
                                       return_data = TRUE)

#rbkTN <- raster::brick("raster_min_temp.grd")
# new terra version
rbkTN <- terra::rast("raster_min_temp.tiff")
format(object.size(TN), "MB") # "0.3 Mb"
format(object.size(rbkTN), "MB") # file would be too large, so this is a pointer to the .nc

# TX daily max temp
TX <- climateExtract::extract_nc_value(first_year = 2000, 
                                       last_year = 2024,
                                       local_file = TRUE, # Source the data from file online = easier
                                       file_path = NULL,
                                       spatial_extent = gfCrop,
                                       clim_variable = "max temp",
                                       statistic = "mean",
                                       grid_size = 0.1,
                                       ecad_v = NULL,
                                       write_raster = TRUE,
                                       out = "raster_max_temp.tiff",
                                       return_data = TRUE)

rbkTX <- terra::rast("raster_max_temp.tiff")
format(object.size(TX), "MB") # "0.3 Mb"
format(object.size(rbkTX), "MB")

# RR daily precipitation sum
RR <- climateExtract::extract_nc_value(first_year = 2000, 
                                       last_year = 2024,
                                       local_file = FALSE, # Source the data from file online = easier
                                       file_path = NULL,
                                       spatial_extent = gfCrop,
                                       clim_variable = "precipitation",
                                       statistic = "mean",
                                       grid_size = 0.1,
                                       ecad_v = NULL,
                                       write_raster = TRUE,
                                       out = "raster_precip.tiff",
                                       return_data = TRUE)

rbkRR <- terra::rast("raster_precip.tiff")
format(object.size(RR), "MB") # "0.3 Mb"
format(object.size(rbkRR), "MB")

# QQ daily mean global radiation
QQ <- climateExtract::extract_nc_value(first_year = 2000, 
                                       last_year = 2024,
                                       local_file = FALSE, # Source the data from file online = easier
                                       file_path = NULL,
                                       spatial_extent = gfCrop,
                                       clim_variable = "global radiation",
                                       statistic = "mean",
                                       grid_size = 0.1,
                                       ecad_v = NULL,
                                       write_raster = TRUE,
                                       out = "raster_solar_radiation.tiff",
                                       return_data = TRUE)

rbkQQ <- terra::rast("raster_solar_radiation.tiff")
format(object.size(QQ), "MB") # "0.3 Mb"
format(object.size(rbkQQ), "MB")

# FG wind speed
# this is only updated to 2023-06-30 as of EOBS v29.0..........so some thinking will need to be done
# downloading from online does not work through this package
FG <- climateExtract::extract_nc_value(first_year = 2000, 
                                       last_year = 2024,
                                       local_file = TRUE, # Source the data from file online doesn't work here
                                       file_path = "fg_ens_mean_0.1deg_reg_v28.0e.nc", # cannot manually click this time either
                                       spatial_extent = gfCrop,
                                       clim_variable = "wind speed",
                                       statistic = "mean",
                                       grid_size = 0.1,
                                       ecad_v = 28.0, # :'( data availability
                                       write_raster = TRUE,
                                       out = "raster_wind_speed.tiff",
                                       return_data = TRUE)

rbkFG <- terra::rast("raster_wind_speed.tiff")
format(object.size(FG), "MB") # "0.3 Mb"
format(object.size(rbkFG), "MB")

# HU daily averaged relative humidity
HU <- climateExtract::extract_nc_value(first_year = 2000, 
                                       last_year = 2024,
                                       local_file = FALSE, # Source the data from file online = easier
                                       file_path = NULL,
                                       spatial_extent = gfCrop,
                                       clim_variable = "relative humidity",
                                       statistic = "mean",
                                       grid_size = 0.1,
                                       ecad_v = NULL,
                                       write_raster = TRUE,
                                       out = "raster_relative_humidity.tiff",
                                       return_data = TRUE)

rbkHU <- terra::rast("raster_relative_humidity.tiff")
format(object.size(HU), "MB") # "0.3 Mb"
format(object.size(rbkHU), "MB")

# PP daily averaged sea level pressure
PP <- climateExtract::extract_nc_value(first_year = 2000, 
                                       last_year = 2024,
                                       local_file = FALSE, # Source the data from file online = easier
                                       file_path = NULL,
                                       spatial_extent = gfCrop,
                                       clim_variable = "sea level pressure",
                                       statistic = "mean",
                                       grid_size = 0.1,
                                       ecad_v = NULL,
                                       write_raster = TRUE,
                                       out = "raster_sea_pressure.tiff",
                                       return_data = TRUE)

rbkPP <- terra::rast("raster_sea_pressure.tiff")
format(object.size(PP), "MB") # "0.3 Mb"
format(object.size(rbkPP), "MB")

#### Getting full data -- not needed for this go-around ####

# TN daily minimum temperature - FULL - 1950 to 2022
#TNfull<-climateExtract::extract_nc_value(first_year = 1950, 
#                                         last_year = 2022,
#                                         local_file = T, # Source the data from file online = eaiser
#                                         file_path = NULL,
#                                         spatial_extent = gfCrop,
#                                         clim_variable = "min temp",
#                                         statistic = "mean",
#                                         grid_size = 0.1,
#                                         ecad_v = NULL,
#                                         write_raster = TRUE,
#                                         out = "raster_min_temp_full.grd",
#                                         return_data = TRUE)

#rbkTNfull<-raster::brick("raster_min_temp_full.grd")
#format(object.size(TNfull), "MB") # "0.8 Mb"
#format(object.size(rbkTNfull), "MB") # "2.3 Mb"

# TX daily max temp - FULL - 1950 to 2022
#TXfull<-climateExtract::extract_nc_value(first_year = 1950, 
#                                         last_year = 2022,
#                                         local_file = T, # Source the data from file online = easier
#                                         file_path = NULL,
#                                         spatial_extent = gfCrop,
#                                         clim_variable = "max temp",
#                                         statistic = "mean",
#                                         grid_size = 0.1,
#                                         ecad_v = NULL,
#                                         write_raster = TRUE,
#                                         out = "raster_max_temp_full.grd",
#                                         return_data = TRUE)

#rbkTXfull<-raster::brick("raster_max_temp_full.grd")
#format(object.size(TXfull), "MB") # "0.8 Mb"
#format(object.size(rbkTXfull), "MB") # "2.3 Mb"

# RR daily precipitation sum - FULL - 1950 to 2022
#RRfull<-climateExtract::extract_nc_value(first_year = 1950, 
#                                         last_year = 2022,
#                                         local_file = T, # Source the data from file online = easier
#                                         file_path = NULL,
#                                         spatial_extent = gfCrop,
#                                         clim_variable = "precip",
#                                         statistic = "sum",
#                                         grid_size = 0.1,
#                                         ecad_v = NULL,
#                                         write_raster = TRUE,
#                                         out = "raster_precip_full.grd",
#                                         return_data = TRUE)

#rbkRRfull<-raster::brick("raster_precip_full.grd")
#format(object.size(RRfull), "MB") # "0.8 Mb"
#format(object.size(rbkRRfull), "MB") # "2.3 Mb"

# QQ daily mean global radiation  - FULL - 1950 to 2022
#QQfull<-climateExtract::extract_nc_value(first_year = 1950, 
#                                         last_year = 2022,
#                                         local_file = T, # Source the data from file online = easier
#                                         file_path = NULL,
#                                         spatial_extent = gfCrop,
#                                         clim_variable = "solar_radiation",
#                                         statistic = "sum",
#                                         grid_size = 0.1,
#                                         ecad_v = NULL,
#                                         write_raster = TRUE,
#                                         out = "solar_radiation_full.grd",
#                                         return_data = TRUE)

#rbkQQfull<-raster::brick("solar_radiation_full.grd")
#format(object.size(QQfull), "MB") # "0.8 Mb"
#format(object.size(rbkQQfull), "MB") # "2.3 Mb"

# FG wind speed  - FULL - 1950 to 2022
#FGfull<-climateExtract::extract_nc_value(first_year = 1950, 
#                                         last_year = 2022,
#                                         local_file = T, # Source the data from file online = easier
#                                         file_path = NULL,
#                                         spatial_extent = gfCrop,
#                                         clim_variable = "wind_speed",
#                                         statistic = "mean",
#                                         grid_size = 0.1,
#                                         ecad_v = NULL,
#                                         write_raster = TRUE,
#                                         out = "wind_speed_full.grd",
#                                         return_data = TRUE)

#rbkFGfull<-raster::brick("wind_speed_full.grd")
#format(object.size(FGfull), "MB") # "0.5 Mb"
#format(object.size(rbkFGfull), "MB") # "1.3 Mb"





#### Checking selected pixels -- skipped ####
# Daily temperature averages 
#monthly_avg_temp_R<-temporal_aggregate(x = rbkTG,
#                                       agg_function = "mean",
#                                       variable_name = "average temp",
#                                       time_step = "monthly")
#par(mfrow=c(1,1))
#raster::plot(monthly_avg_temp_R[["2002.09"]])
#plot(GF,add=T)
# x = -4.4501395912103, y = 56.274860511038
# x = -4.38347292481033, y = 56.274860511038

#monthly_avg_precip_R<-temporal_aggregate(x = rbkRR,
#                                      agg_function = "sum",
#                                      variable_name = "precip",
#                                       time_step = "monthly")
#par(mfrow=c(1,1))
#raster::plot(monthly_avg_precip_R[["2002.09"]])
#plot(GF,add=T)

monthly_avg_solar_R<-temporal_aggregate(x = rbkQQ,
                                        agg_function = "mean",
                                        variable_name = "solar_radiation",
                                        time_step = "monthly")
par(mfrow=c(1,1))
raster::plot(monthly_avg_solar_R[["2002.09"]])
plot(GF,add=T)

monthly_avg_windspeed_R<-temporal_aggregate(x = rbkFG,
                                            agg_function = "mean",
                                            variable_name = "wind_speed",
                                            time_step = "monthly")
par(mfrow=c(1,1))
raster::plot(monthly_avg_windspeed_R[["2002.09"]])
plot(GF,add=T)

# Check identify cell correcly
#library(terra)
#r<-monthly_avg_temp_R[["2002.09"]]
#y <- terra::ifel(r < 11.42, 0, r)
#y <- terra::ifel(y > 11.7, 0, y)
#as.data.frame(y, xy = TRUE)
#5 -4.450140 56.27486 11.44533 # desired cells
#6 -4.383473 56.27486 11.50500

# Visualize data overlay Glen Finglas plots
#filename <- paste0("/Users/stuartsmith/Documents/JHI/Projects/WPS/WP36/362_Trophic_Interactions/Outputs/GlenFinglas/Vegetation/Climate/", "GF_TG",".jpeg" )
#jpeg(filename,width= 12, height = 12,units ="cm",bg ="transparent", res = 800)
#raster::plot(monthly_avg_temp_R[["2002.09"]])
#plot(GF,add=T)
#dev.off() # 6 grid cells, but we only need 2 

#########################################################################################################
#### Glen Finglas E-OBS climate ####
#########################################################################################################

# fix working directory
setwd('..')
#setwd("./JHIGit/")

#### quick aside to get overlap percentages ####
# E&F is a bit more complicated, so calculate the area weighting here
# the area of interest is
plots <- st_transform(st_read("./Spatial/boundaries/newplotsCleaned.shp"), "EPSG:4326") %>% 
  mutate(plot = str_sub(Plot_ID, 2),
         block = str_sub(Plot_ID, 0, 1),
         blockGroup = factor(block, levels = LETTERS[1:6], labels = c("A&B", "A&B", "C&D","C&D", "E&F", "E&F")))
# confirm that these are ok by a plot
#plot(rbkTG[[1]])
#plot(plots, add = T)
#extract the area of each cell that is contained within each polygon
overlap <- exact_extract(rbkTG[[1]], plots, coverage_area = TRUE)
#add polygon names that the results will be grouped by
names(overlap) <- plots$blockGroup
#bind the list output into a df and calculate the proportion cover for each category
overlapCalcs <- bind_rows(overlap, .id = "region") %>%
  group_by(region, value) %>%
  summarize(total_area = sum(coverage_area)) %>%
  group_by(region) %>%
  mutate(proportion = total_area/sum(total_area))
# kind of a "hack" solution, but cells are always the same, so it's fine
# precision checks for other rasters show the same thing
# so, the maths come out to...
print(paste0("E&F Blocks are ", round(overlapCalcs$proportion[3]*100,2), 
             "% in Cell 2 and ", round(overlapCalcs$proportion[4]*100,2), "% in Cell 3"))

# TG - daily temperature
daily_temp_df <- as.data.frame(rbkTG, xy = TRUE)
daily_temp <- reshape2::melt(daily_temp_df, id.vars=c("x","y"),
                             value.name = "mean_temp_C", variable.name='date_yy_mm_dd', na.rm = T)
dim(daily_temp) # 52596     4

# Retain only three grid cells for Glen Finglas
levels(as.factor(daily_temp$x)) # "-4.51680625761026" "-4.4501395912103"  "-4.38347292481033"
levels(as.factor(daily_temp$y)) # "56.274860511038"  "56.3248605108792"
# three cells
# Cell 1: x = -4.4501395912103, y = 56.2748605110379, C&D Blocks
# Cell 2: x = -4.38347292481033, y = 56.2748605110379, A&B Blocks, E&F Blocks
# Cell 3: x = -4.38347292481033, y = 56.3248605108792, E&F Blocks

#### return to raster processing ####
# only weather within the study period is used: see original script for full data
daily_temp$xy<-as.factor(with(daily_temp, paste(x,y, sep="_")))
TGgf <- daily_temp %>% filter(xy %in% c("-4.4501395912103_56.2748605110379", 
                                        "-4.38347292481033_56.2748605110379",
                                        "-4.38347292481033_56.3248605108792")) %>% 
  select(c(x, y, date_yy_mm_dd, mean_temp_C))%>% 
  # date processing
  mutate(Date = as.Date(date_yy_mm_dd ,"%Y-%m-%d"),
         Yr = format(as.Date(Date), "%Y"))
# not sure what the code is below but doesn't look necessary -- no matches for "X" in date
#sum(grepl("X",TGgf$date_yy_mm_dd))
#TGgf$date_yy_mm_dd <- gsub("X","",TGgf$date_yy_mm_dd)
# Export tg - daily temperature
write.csv(TGgf,file="./Climate/EOBS/tg.csv", row.names = F)


# TN - min temperature
min_temp_df <- as.data.frame(rbkTN, xy = TRUE)
min_temp <- reshape2::melt(min_temp_df, id.vars=c("x","y"),
                           value.name = "min_temp_C", variable.name='date_yy_mm_dd',na.rm = T)
dim(min_temp) # 52596     4
# Retain only three grid cells for Glen Finglas
min_temp$xy <- as.factor(with(min_temp, paste(x,y, sep="_")))
# x coordinates change slightly here (move just slightly more south)
levels(min_temp$xy)
TNgf <- min_temp %>% filter(xy %in% c("-4.4501395912103_56.2748605109631", 
                                      "-4.38347292481033_56.2748605109631",
                                      "-4.38347292481033_56.3248605108047")) %>% 
  select(c(x, y, date_yy_mm_dd, min_temp_C)) %>% 
  # date processing
  mutate(Date = as.Date(date_yy_mm_dd ,"%Y-%m-%d"),
         Yr = format(as.Date(Date), "%Y"))
# plot to see really high-level trends
TNgf %>% group_by(Yr) %>% summarise(min_temp=mean(min_temp_C))%>%
  ggplot(aes(x=Yr, y=min_temp)) + geom_point()
# more granular look at long-term trends
TNgf %>% group_by(Date) %>% summarise(min_temp=mean(min_temp_C))%>%
  ggplot(aes(x=Date, y=min_temp)) + geom_point()
# Export tn - min temperature
write.csv(TNgf, file="./Climate/EOBS/tn.csv", row.names = F)


# TX - max temperature
max_temp_df <- as.data.frame(rbkTX, xy = TRUE)
max_temp<-reshape2::melt(max_temp_df, id.vars=c("x","y"),
                         value.name = "max_temp_C", variable.name='date_yy_mm_dd',na.rm = T)
dim(max_temp) # 52596     4
# Retain only three grid cells for Glen Finglas
max_temp$xy<-as.factor(with(max_temp, paste(x,y, sep="_")))
# same coordinates as min_temp
# always worth a manual check if you're not looking at this script in 07/2024
#levels(max_temp$xy)
TXgf <- max_temp %>% filter(xy %in% c("-4.4501395912103_56.2748605109631", 
                                      "-4.38347292481033_56.2748605109631",
                                      "-4.38347292481033_56.3248605108047")) %>% 
  select(c(x, y, date_yy_mm_dd, max_temp_C)) %>% 
  # date processing
  mutate(Date = as.Date(date_yy_mm_dd ,"%Y-%m-%d"),
         Yr = format(as.Date(Date), "%Y"))
# plotting
TXgf %>% group_by(Yr) %>% summarise(max_temp=mean(max_temp_C))%>%
  ggplot(aes(x=Yr, y=max_temp)) + geom_point()
TXgf %>% group_by(Date) %>% summarise(max_temp=mean(max_temp_C))%>%
  ggplot(aes(x=Date, y=max_temp)) + geom_point()
# Export tx - max temperature
write.csv(TXgf, file="./Climate/EOBS/tx.csv", row.names = F)


# RR - precipitation
precip_df <- as.data.frame(rbkRR, xy = TRUE)
precip<-reshape2::melt(precip_df, id.vars=c("x","y"),
                       value.name = "precip", variable.name='date_yy_mm_dd',na.rm = T)
dim(precip) # 52596     4
# Retain only three grid cells for Glen Finglas
precip$xy<-as.factor(with(precip, paste(x,y, sep="_")))
# levels change again here for the latitude
levels(precip$xy)
RRgf <- precip %>% filter(xy %in% c("-4.45013959120935_56.274859277992", 
                                    "-4.38347292480938_56.274859277992",
                                    "-4.38347292480938_56.3248592757786")) %>% 
  select(c(x, y, date_yy_mm_dd, precip)) %>% 
  # date processing
  mutate(Date = as.Date(date_yy_mm_dd ,"%Y-%m-%d"),
         Yr = format(as.Date(Date), "%Y"))
# Export RR - precipitation
write.csv(RRgf, file="./Climate/EOBS/rr.csv", row.names = F)


# QQ - solar radiation
solar_df <- as.data.frame(rbkQQ, xy = TRUE)
solar<-reshape2::melt(solar_df, id.vars=c("x","y"),
                      value.name = "solar_radiation", variable.name='date_yy_mm_dd',na.rm = T)
# row mismatch is bc QQ contains data for 01/2024 as well
dim(solar) # 52752     4
tail(solar)
# Retain only three grid cells for Glen Finglas
solar$xy<-as.factor(with(solar, paste(x,y, sep="_")))
# both x and y change this time
levels(solar$xy)
QQgf <- solar %>% filter(xy %in% c("-4.4501395912_56.2748605259", 
                                   "-4.3834729248_56.2748605259",
                                   "-4.3834729248_56.3248605257")) %>% 
  select(c(x, y, date_yy_mm_dd, solar_radiation)) %>% 
  # date processing
  mutate(Date = as.Date(date_yy_mm_dd ,"%Y-%m-%d"),
         Yr = format(as.Date(Date), "%Y"))
# see what's up here with the dimension mismatch
#QQgf2 <- QQgf %>% mutate(xM = round(x, 2),
#                         yM = round(y, 2))
#TXgf2 <- TXgf %>% mutate(xM = round(x, 2),
#                         yM = round(y, 2))
#
#testMerge <- left_join(QQgf2, TXgf2, by = c("xM", "yM", "Date"))
# Export QQ -solar radiation
write.csv(QQgf, file="./Climate/EOBS/qq.csv", row.names = F)


# FG - wind speed
windspeed_df <- as.data.frame(rbkFG, xy = TRUE)
windspeed<-reshape2::melt(windspeed_df, id.vars=c("x","y"),
                          value.name = "wind_speed", variable.name='date_yy_mm_dd',na.rm = T)
dim(windspeed) # 51483     4
tail(windspeed)
# Retain only three grid cells for Glen Finglas
windspeed$xy<-as.factor(with(windspeed, paste(x,y, sep="_")))
levels(windspeed$xy)
FGgf <- windspeed %>% filter(xy %in% c("-4.4501_56.2749", 
                                       "-4.38343333333333_56.2749",
                                       "-4.38343333333333_56.3249")) %>% 
  select(c(x, y, date_yy_mm_dd, wind_speed)) %>% 
  # date processing
  mutate(Date = as.Date(date_yy_mm_dd ,"%Y-%m-%d"),
         Yr = format(as.Date(Date), "%Y"))
# Export FG - windspeed
write.csv(FGgf, file="./Climate/EOBS/fg.csv", row.names = F)


# HU - relative humidity
humid_df <- as.data.frame(rbkHU, xy = TRUE)
humid <- reshape2::melt(precip_df, id.vars=c("x","y"),
                        value.name = "humid", variable.name='date_yy_mm_dd',na.rm = T)
dim(humid) # 52596     4
# Retain only three grid cells for Glen Finglas
humid$xy <- as.factor(with(humid, paste(x,y, sep="_")))
# levels change again here for the latitude
levels(humid$xy)
HUgf <- humid %>% filter(xy %in% c("-4.45013959120935_56.274859277992", 
                                   "-4.38347292480938_56.274859277992",
                                   "-4.38347292480938_56.3248592757786")) %>% 
  select(c(x, y, date_yy_mm_dd, humid)) %>% 
  # date processing
  mutate(Date = as.Date(date_yy_mm_dd ,"%Y-%m-%d"),
         Yr = format(as.Date(Date), "%Y"))
# Export HU - relative humidity
write.csv(HUgf, file="./Climate/EOBS/hu.csv", row.names = F)

# PP - sea pressure
sea_pressure_df <- as.data.frame(rbkPP, xy = TRUE)
sea_pressure <- reshape2::melt(sea_pressure_df, id.vars=c("x","y"),
                               value.name = "pressure", variable.name='date_yy_mm_dd',na.rm = T)
dim(sea_pressure) # 52596     4
# Retain only three grid cells for Glen Finglas
sea_pressure$xy <- as.factor(with(sea_pressure, paste(x,y, sep="_")))
# levels change again here for the latitude
levels(sea_pressure$xy)
PPgf <- sea_pressure %>% filter(xy %in% c("-4.45013959120917_56.2748605122776", 
                                          "-4.3834729248092_56.2748605122776",
                                          "-4.3834729248092_56.3248605121161")) %>% 
  select(c(x, y, date_yy_mm_dd, pressure)) %>% 
  # date processing
  mutate(Date = as.Date(date_yy_mm_dd ,"%Y-%m-%d"),
         Yr = format(as.Date(Date), "%Y"))
# Export PP - sea pressure
write.csv(PPgf, file="./Climate/EOBS/pp.csv", row.names = F)

#########################################################################################################
#### Glen Finglas Full Time Period (1950 - 2024) E-OBS climate -- repeat previous steps ####
#########################################################################################################
# assume that we're working from a local version of the git
# "/[Personal_Path_to_Folder]/GlenFinglas/"
setwd("./Climate/Raw/")

# downloads from online don't work without increasing timeout
#options(timeout = max(1000, getOption("timeout")))


# TG_full daily mean temperature
TG_full <- climateExtract::extract_nc_value(#first_year = 2000, 
  last_year = 2024,
  local_file = TRUE, # source data from online (if running for the first time)
  file_path = "tg_ens_mean_0.1deg_reg_v29.0e.nc",
  spatial_extent = gfCrop,
  clim_variable = "mean temp",
  statistic = "mean", # getting "spread" for uncertainty might be a good idea too
  grid_size = 0.1,
  ecad_v = NULL,
  write_raster = TRUE,
  out = "raster_full_mean_temp.tiff",
  return_data = TRUE)

# new terra version
rbkTG_full <- terra::rast("raster_full_mean_temp.tiff")
format(object.size(TG_full), "MB")
format(object.size(rbkTG_full), "MB") # file would be too large, so this is a pointer to the .nc

#daily_TG_full<- extract(rbk,GF,sp=T, method='simple')
#write.csv(daily_TG_full, "GF_TG_full.csv",row.names=F) #,sep = ",",dec = ".",col.names = TRUE, row.names=F)

# TN_full daily minimum temperature
TN_full <- climateExtract::extract_nc_value(#first_year = 2000, 
  last_year = 2024,
  local_file = TRUE, # Source the data from file online = eaiser
  file_path = "tn_ens_mean_0.1deg_reg_v29.0e.nc",
  spatial_extent = gfCrop,
  clim_variable = "min temp",
  statistic = "mean",
  grid_size = 0.1,
  ecad_v = NULL,
  write_raster = TRUE,
  out = "raster_full_min_temp.tiff",
  return_data = TRUE)

#rbkTN_full <- raster::brick("raster_full_min_temp.grd")
# new terra version
rbkTN_full <- terra::rast("raster_full_min_temp.tiff")
format(object.size(TN_full), "MB") # "0.3 Mb"
format(object.size(rbkTN_full), "MB") # file would be too large, so this is a pointer to the .nc

# TX_full daily max temp
TX_full <- climateExtract::extract_nc_value(#first_year = 2000, 
  last_year = 2024,
  local_file = TRUE, # Source the data from file online = easier
  file_path = "tx_ens_mean_0.1deg_reg_v29.0e.nc",
  spatial_extent = gfCrop,
  clim_variable = "max temp",
  statistic = "mean",
  grid_size = 0.1,
  ecad_v = NULL,
  write_raster = TRUE,
  out = "raster_full_max_temp.tiff",
  return_data = TRUE)

rbkTX_full <- terra::rast("raster_full_max_temp.tiff")
format(object.size(TX_full), "MB") # "0.3 Mb"
format(object.size(rbkTX_full), "MB")

# RR_full daily precipitation sum
RR_full <- climateExtract::extract_nc_value(#first_year = 2000, 
  last_year = 2024,
  local_file = TRUE, # Source the data from file online = easier
  file_path = "rr_ens_mean_0.1deg_reg_v29.0e.nc",
  spatial_extent = gfCrop,
  clim_variable = "precipitation",
  statistic = "mean",
  grid_size = 0.1,
  ecad_v = NULL,
  write_raster = TRUE,
  out = "raster_full_precip.tiff",
  return_data = TRUE)

rbkRR_full <- terra::rast("raster_full_precip.tiff")
format(object.size(RR_full), "MB") # "0.3 Mb"
format(object.size(rbkRR_full), "MB")

# QQ_full daily mean global radiation
QQ_full <- climateExtract::extract_nc_value(#first_year = 2000, 
  last_year = 2024,
  local_file = TRUE, # Source the data from file online = easier
  file_path = "qq_ens_mean_0.1deg_reg_v29.0e.nc",
  spatial_extent = gfCrop,
  clim_variable = "global radiation",
  statistic = "mean",
  grid_size = 0.1,
  ecad_v = NULL,
  write_raster = TRUE,
  out = "raster_full_solar_radiation.tiff",
  return_data = TRUE)

rbkQQ_full <- terra::rast("raster_full_solar_radiation.tiff")
format(object.size(QQ_full), "MB") # "0.3 Mb"
format(object.size(rbkQQ_full), "MB")

# FG_full wind speed
# this is only updated to 2023-06-30 as of EOBS v29.0..........so some thinking will need to be done
# downloading from online does not work through this package
FG_full <- climateExtract::extract_nc_value(#first_year = 2000, 
  last_year = 2024,
  local_file = TRUE, # Source the data from file online doesn't work here
  file_path = "fg_ens_mean_0.1deg_reg_v28.0e.nc", # cannot manually click this time either
  spatial_extent = gfCrop,
  clim_variable = "wind speed",
  statistic = "mean",
  grid_size = 0.1,
  ecad_v = 28.0, # :'( data availability
  write_raster = TRUE,
  out = "raster_full_wind_speed.tiff",
  return_data = TRUE)

rbkFG_full <- terra::rast("raster_full_wind_speed.tiff")
format(object.size(FG_full), "MB") # "0.3 Mb"
format(object.size(rbkFG_full), "MB")

# HU_full daily averaged relative humidity
HU_full <- climateExtract::extract_nc_value(#first_year = 2000, 
  last_year = 2024,
  local_file = TRUE, # Source the data from file online = easier
  file_path = "hu_ens_mean_0.1deg_reg_v29.0e.nc",
  spatial_extent = gfCrop,
  clim_variable = "relative humidity",
  statistic = "mean",
  grid_size = 0.1,
  ecad_v = NULL,
  write_raster = TRUE,
  out = "raster_full_relative_humidity.tiff",
  return_data = TRUE)

rbkHU_full <- terra::rast("raster_full_relative_humidity.tiff")
format(object.size(HU_full), "MB") # "0.3 Mb"
format(object.size(rbkHU_full), "MB")

# PP daily averaged sea level pressure
PP_full <- climateExtract::extract_nc_value(#first_year = 2000, 
  last_year = 2024,
  local_file = TRUE, # Source the data from file online = easier
  file_path = "pp_ens_mean_0.1deg_reg_v29.0e.nc",
  spatial_extent = gfCrop,
  clim_variable = "sea level pressure",
  statistic = "mean",
  grid_size = 0.1,
  ecad_v = NULL,
  write_raster = TRUE,
  out = "raster_full_sea_pressure.tiff",
  return_data = TRUE)

rbkPP_full <- terra::rast("raster_full_sea_pressure.tiff")
format(object.size(PP_full), "MB") # "0.3 Mb"
format(object.size(rbkPP_full), "MB")

#########################################################################################################
#### Glen Finglas E-OBS climate ####
#########################################################################################################

# fix working directory
setwd('..')

#### quick aside to get overlap percentages ####
# E&F is a bit more complicated, so calculate the area weighting here
# the area of interest is
plots <- st_transform(st_read("./Spatial/boundaries/newplotsCleaned.shp"), "EPSG:4326") %>% 
  mutate(plot = str_sub(Plot_ID, 2),
         block = str_sub(Plot_ID, 0, 1),
         blockGroup = factor(block, levels = LETTERS[1:6], labels = c("A&B", "A&B", "C&D","C&D", "E&F", "E&F")))
# confirm that these are ok by a plot
#plot(rbkTG_full[[1]])
#plot(plots, add = T)
#extract the area of each cell that is contained within each polygon
overlap <- exact_extract(rbkTG_full[[1]], plots, coverage_area = TRUE)
#add polygon names that the results will be grouped by
names(overlap) <- plots$blockGroup
#bind the list output into a df and calculate the proportion cover for each category
overlapCalcs <- bind_rows(overlap, .id = "region") %>%
  group_by(region, value) %>%
  summarize(total_area = sum(coverage_area)) %>%
  group_by(region) %>%
  mutate(proportion = total_area/sum(total_area))
# kind of a "hack" solution, but cells are always the same, so it's fine
# precision checks for other rasters show the same thing
# so, the maths come out to...
print(paste0("E&F Blocks are ", round(overlapCalcs$proportion[3]*100,2), 
             "% in Cell 2 and ", round(overlapCalcs$proportion[4]*100,2), "% in Cell 3"))

# TG_full - daily temperature
daily_temp_full_df <- as.data.frame(rbkTG_full, xy = TRUE)
daily_temp_full <- reshape2::melt(daily_temp_full_df, id.vars=c("x","y"),
                                  value.name = "mean_temp_C", variable.name='date_yy_mm_dd', na.rm = T)
dim(daily_temp_full) # 52596     4

# Retain only three grid cells for Glen Finglas
levels(as.factor(daily_temp_full$x)) # "-4.51680625761026" "-4.4501395912103"  "-4.38347292481033"
levels(as.factor(daily_temp_full$y)) # "56.274860511038"  "56.3248605108792"
# three cells
# Cell 1: x = -4.4501395912103, y = 56.2748605110379, C&D Blocks
# Cell 2: x = -4.38347292481033, y = 56.2748605110379, A&B Blocks, E&F Blocks
# Cell 3: x = -4.38347292481033, y = 56.3248605108792, E&F Blocks

#### return to raster processing ####
# only weather within the study period is used: see original script for full data
daily_temp_full$xy<-as.factor(with(daily_temp_full, paste(x,y, sep="_")))
TG_fullgf <- daily_temp_full %>% filter(xy %in% c("-4.4501395912103_56.2748605110379", 
                                                  "-4.38347292481033_56.2748605110379",
                                                  "-4.38347292481033_56.3248605108792")) %>% 
  select(c(x, y, date_yy_mm_dd, mean_temp_C))%>% 
  # date processing
  mutate(Date = as.Date(date_yy_mm_dd ,"%Y-%m-%d"),
         Yr = format(as.Date(Date), "%Y"))
# not sure what the code is below but doesn't look necessary -- no matches for "X" in date
#sum(grepl("X",TG_fullgf$date_yy_mm_dd))
#TG_fullgf$date_yy_mm_dd <- gsub("X","",TG_fullgf$date_yy_mm_dd)
# Export tg - daily temperature
write.csv(TG_fullgf,file="./Climate/EOBS/tg_full.csv", row.names = F)


# TN_full - min temperature
min_temp_full_df <- as.data.frame(rbkTN_full, xy = TRUE)
min_temp_full <- reshape2::melt(min_temp_full_df, id.vars=c("x","y"),
                                value.name = "min_temp_C", variable.name='date_yy_mm_dd',na.rm = T)
dim(min_temp) # 52596     4
# Retain only three grid cells for Glen Finglas
min_temp_full$xy <- as.factor(with(min_temp_full, paste(x,y, sep="_")))
# x coordinates change slightly here (move just slightly more south)
levels(min_temp_full$xy)
TN_fullgf <- min_temp_full %>% filter(xy %in% c("-4.4501395912103_56.2748605109631", 
                                                "-4.38347292481033_56.2748605109631",
                                                "-4.38347292481033_56.3248605108047")) %>% 
  select(c(x, y, date_yy_mm_dd, min_temp_C)) %>% 
  # date processing
  mutate(Date = as.Date(date_yy_mm_dd ,"%Y-%m-%d"),
         Yr = format(as.Date(Date), "%Y"))
# plot to see really high-level trends
TN_fullgf %>% group_by(Yr) %>% summarise(min_temp=mean(min_temp_C))%>%
  ggplot(aes(x=Yr, y=min_temp)) + geom_point()
# more granular look at long-term trends
TN_fullgf %>% group_by(Date) %>% summarise(min_temp=mean(min_temp_C))%>%
  ggplot(aes(x=Date, y=min_temp)) + geom_point()
# Export tn - min temperature
write.csv(TN_fullgf, file="./Climate/EOBS/tn_full.csv", row.names = F)


# TX_full - max temperature
max_temp_full_df <- as.data.frame(rbkTX_full, xy = TRUE)
max_temp_full <-reshape2::melt(max_temp_full_df, id.vars=c("x","y"),
                               value.name = "max_temp_C", variable.name='date_yy_mm_dd',na.rm = T)
dim(max_temp_full) # 52596     4
# Retain only three grid cells for Glen Finglas
max_temp_full$xy<-as.factor(with(max_temp_full, paste(x,y, sep="_")))
# same coordinates as min_temp
# always worth a manual check if you're not looking at this script in 07/2024
#levels(max_temp_full$xy)
TX_fullgf <- max_temp_full %>% filter(xy %in% c("-4.4501395912103_56.2748605109631", 
                                                "-4.38347292481033_56.2748605109631",
                                                "-4.38347292481033_56.3248605108047")) %>% 
  select(c(x, y, date_yy_mm_dd, max_temp_C)) %>% 
  # date processing
  mutate(Date = as.Date(date_yy_mm_dd ,"%Y-%m-%d"),
         Yr = format(as.Date(Date), "%Y"))
# plotting
TX_fullgf %>% group_by(Yr) %>% summarise(max_temp=mean(max_temp_C))%>%
  ggplot(aes(x=Yr, y=max_temp)) + geom_point()
TX_fullgf %>% group_by(Date) %>% summarise(max_temp=mean(max_temp_C))%>%
  ggplot(aes(x=Date, y=max_temp)) + geom_point()
# Export tx - max temperature
write.csv(TX_fullgf, file="./Climate/EOBS/tx_full.csv", row.names = F)


# RR_full - precipitation
precip_full_df <- as.data.frame(rbkRR_full, xy = TRUE)
precip_full <-reshape2::melt(precip_full_df, id.vars=c("x","y"),
                             value.name = "precip", variable.name='date_yy_mm_dd',na.rm = T)
dim(precip_full) # 52596     4
# Retain only three grid cells for Glen Finglas
precip_full$xy<-as.factor(with(precip_full, paste(x,y, sep="_")))
# levels change again here for the latitude
levels(precip_full$xy)
RR_fullgf <- precip_full %>% filter(xy %in% c("-4.45013959120935_56.274859277992", 
                                              "-4.38347292480938_56.274859277992",
                                              "-4.38347292480938_56.3248592757786")) %>% 
  select(c(x, y, date_yy_mm_dd, precip)) %>% 
  # date processing
  mutate(Date = as.Date(date_yy_mm_dd ,"%Y-%m-%d"),
         Yr = format(as.Date(Date), "%Y"))
# Export RR_full - precipitation
write.csv(RR_fullgf, file="./Climate/EOBS/rr_full.csv", row.names = F)


# QQ_full - solar radiation
solar_full_df <- as.data.frame(rbkQQ_full, xy = TRUE)
solar_full <- reshape2::melt(solar_full_df, id.vars=c("x","y"),
                             value.name = "solar_radiation", variable.name='date_yy_mm_dd',na.rm = T)
# row mismatch is bc QQ_full contains data for 01/2024 as well
dim(solar_full) # 52752     4
tail(solar_full)
# Retain only three grid cells for Glen Finglas
solar_full$xy<-as.factor(with(solar_full, paste(x,y, sep="_")))
# both x and y change this time
levels(solar_full$xy)
QQ_fullgf <- solar_full %>% filter(xy %in% c("-4.4501395912_56.2748605259", 
                                             "-4.3834729248_56.2748605259",
                                             "-4.3834729248_56.3248605257")) %>% 
  select(c(x, y, date_yy_mm_dd, solar_radiation)) %>% 
  # date processing
  mutate(Date = as.Date(date_yy_mm_dd ,"%Y-%m-%d"),
         Yr = format(as.Date(Date), "%Y"))
# see what's up here with the dimension mismatch
#QQ_fullgf2 <- QQ_fullgf %>% mutate(xM = round(x, 2),
#                         yM = round(y, 2))
#TX_fullgf2 <- TX_fullgf %>% mutate(xM = round(x, 2),
#                         yM = round(y, 2))
#
#testMerge <- left_join(QQ_fullgf2, TX_fullgf2, by = c("xM", "yM", "Date"))
# Export QQ_full -solar radiation
write.csv(QQ_fullgf, file="./Climate/EOBS/qq_full.csv", row.names = F)


# FG_full - wind speed
windspeed_full_df <- as.data.frame(rbkFG_full, xy = TRUE)
windspeed_full <- reshape2::melt(windspeed_full_df, id.vars=c("x","y"),
                                 value.name = "wind_speed", variable.name='date_yy_mm_dd',na.rm = T)
dim(windspeed_full) # 51483     4
tail(windspeed_full)
# Retain only three grid cells for Glen Finglas
windspeed_full$xy<-as.factor(with(windspeed_full, paste(x,y, sep="_")))
levels(windspeed_full$xy)
FG_fullgf <- windspeed_full %>% filter(xy %in% c("-4.4501_56.2749", 
                                                 "-4.38343333333333_56.2749",
                                                 "-4.38343333333333_56.3249")) %>% 
  select(c(x, y, date_yy_mm_dd, wind_speed)) %>% 
  # date processing
  mutate(Date = as.Date(date_yy_mm_dd ,"%Y-%m-%d"),
         Yr = format(as.Date(Date), "%Y"))
# Export FG_full - windspeed
write.csv(FG_fullgf, file="./Climate/EOBS/fg_full.csv", row.names = F)


# HU_full - relative humidity
humid_full_df <- as.data.frame(rbkHU_full, xy = TRUE)
humid_full <- reshape2::melt(precip_full_df, id.vars=c("x","y"),
                             value.name = "humid", variable.name='date_yy_mm_dd',na.rm = T)
dim(humid_full) # 52596     4
# Retain only three grid cells for Glen Finglas
humid_full$xy <- as.factor(with(humid_full, paste(x,y, sep="_")))
# levels change again here for the latitude
levels(humid_full$xy)
HU_fullgf <- humid_full %>% filter(xy %in% c("-4.45013959120935_56.274859277992", 
                                             "-4.38347292480938_56.274859277992",
                                             "-4.38347292480938_56.3248592757786")) %>% 
  select(c(x, y, date_yy_mm_dd, humid)) %>% 
  # date processing
  mutate(Date = as.Date(date_yy_mm_dd ,"%Y-%m-%d"),
         Yr = format(as.Date(Date), "%Y"))
# Export HU_full - relative humidity
write.csv(HU_fullgf, file="./Climate/EOBS/hu_full.csv", row.names = F)

# PP - sea pressure
sea_pressure_full_df <- as.data.frame(rbkPP_full, xy = TRUE)
sea_pressure_full <- reshape2::melt(sea_pressure_full_df, id.vars=c("x","y"),
                                    value.name = "pressure", variable.name='date_yy_mm_dd',na.rm = T)
dim(sea_pressure_full) # 52596     4
# Retain only three grid cells for Glen Finglas
sea_pressure_full$xy <- as.factor(with(sea_pressure_full, paste(x,y, sep="_")))
# levels change again here for the latitude
levels(sea_pressure_full$xy)
PPgf_full <- sea_pressure_full %>% filter(xy %in% c("-4.45013959120917_56.2748605122776", 
                                                    "-4.3834729248092_56.2748605122776",
                                                    "-4.3834729248092_56.3248605121161")) %>% 
  select(c(x, y, date_yy_mm_dd, pressure)) %>% 
  # date processing
  mutate(Date = as.Date(date_yy_mm_dd ,"%Y-%m-%d"),
         Yr = format(as.Date(Date), "%Y"))
# Export PP - sea pressure
write.csv(PPgf_full, file="./Climate/EOBS/pp_full.csv", row.names = F)