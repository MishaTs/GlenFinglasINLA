# libraries
library(tidyverse)
library(terra)

# read in raster
demFull <- terra::rast("./Spatial/DEM/EEA10_INSP_2011_03_31.tif")
plot(demFull)

##### Derive DEM Variables #####
# get everything possible for now
slope <- demFull %>% terrain(v = "slope",
                             neighbors = 8)
aspect <- demFull %>% terrain(v = "aspect",
                              neighbors = 8)
# two alternative measures of roughness
# Roughness is the difference between the maximum and the minimum value of a cell and its 8 surrounding cells.
roughness <- demFull %>% terrain(v = "roughness")
# TRI (Terrain Ruggedness Index) is the mean of the absolute differences between the value of a cell and its 8 surrounding cells.
tri <- demFull %>% terrain(v = "TRI")
# TRIriley (TRI according to Riley et al., 2007) returns the square root of summed squared differences between the value of a cell and its 8 surrounding cells
triRiley <- demFull %>% terrain(v = "TRIriley")
# TRIrmsd computes the square root of the mean of the squared differences between these cells (TRIriley)
triRmsd <- demFull %>% terrain(v = "TRIrmsd")
# TPI (Topographic Position Index) is the difference between the value of a cell and the mean value of its 8 surrounding cells
# less granular TRI rm sd
tpi <- demFull %>% terrain(v = "TPI")
plot(tpi)
# flow direction
# flowdir returns the "flow direction" (of water), that is the direction of the greatest drop in elevation (or the smallest rise if all neighbors are higher). They are encoded as powers of 2 (0 to 7). 
flowDir <- demFull %>% terrain(v = "flowdir")
plot(flowDir)

##### Process Vector Geodata #####
# handle point data
points81 <- st_read("./Spatial/81points/81PointsCorrected.shp")
pts81 <- st_transform(points81, "EPSG:4326") # reproject to WGS84
# maximum ID length
#max(str_length(unique(pts81$ID)))
# manually created points have _ attached; doesn't include manually changed points
pts81 <- pts81 %>% mutate(id = str_extract(ID, "[[:alnum:]]{0,7}")) %>% 
  select(-c("X","Y","ID")) %>% 
  rowid_to_column("ID")
pts81$lat <- st_coordinates(pts81)[,1]
pts81$long <- st_coordinates(pts81)[,2]
ggplot(pts81, aes(x = lat, y = long)) + geom_point()

# handle plot boundaries
plotBounds <- st_read("./Spatial/boundaries/newplotsCleaned.shp")
plotBounds <- st_transform(plotBounds, "EPSG:4326") # reproject to WGS84
plotBounds <- plotBounds %>% mutate(plot = str_sub(Plot_ID, 2),
                                    block = str_sub(Plot_ID, 0, 1),
                                    blockGroup = factor(block, 
                                                        levels = LETTERS[1:6],
                                                        labels = c("A&B", "A&B", "C&D","C&D", "E&F", "E&F"))) %>%
  select(-ID) %>% rowid_to_column("ID")


##### Vectorise Raster Variables For Points #####
# first make compatible with terra
pts81Vect <- vect(pts81)
# then, extract
pts81Elev <- terra::extract(demFull,
                            pts81Vect,
                            na.rm = T,
                            # "simple": values for the cell a point falls in are returned
                            # "bilinear" returned values are interpolated from the values of the four nearest raster cells
                            method="simple",
                            ID = TRUE) %>% rename(elev = "EEA10_INSP_2011_03_31")
# create a final covariate dataset
pts81Cov <- left_join(pts81, pts81Elev, by = "ID")

# repeat for slope
pts81Slope <- terra::extract(slope,
                             pts81Vect,
                             na.rm = T,
                             # "simple": values for the cell a point falls in are returned
                             # "bilinear" returned values are interpolated from the values of the four nearest raster cells
                             method="simple",
                             ID = TRUE)
pts81Cov <- left_join(pts81Cov, pts81Slope, by = "ID")

# repeat for aspect
pts81Aspect <- terra::extract(aspect,
                              pts81Vect,
                              na.rm = T,
                              # "simple": values for the cell a point falls in are returned
                              # "bilinear" returned values are interpolated from the values of the four nearest raster cells
                              method="simple",
                              ID = TRUE)
pts81Cov <- left_join(pts81Cov, pts81Aspect, by = "ID")

# continue to repeat for roughness, TRI, TRIriley, TRIrmsd, TPI, flowdir
pts81Rough <- terra::extract(roughness,
                              pts81Vect,
                              na.rm = T,
                              # "simple": values for the cell a point falls in are returned
                              # "bilinear" returned values are interpolated from the values of the four nearest raster cells
                              method="simple",
                              ID = TRUE)
pts81Cov <- left_join(pts81Cov, pts81Rough, by = "ID")

pts81TRI <- terra::extract(tri,
                             pts81Vect,
                             na.rm = T,
                             # "simple": values for the cell a point falls in are returned
                             # "bilinear" returned values are interpolated from the values of the four nearest raster cells
                             method="simple",
                             ID = TRUE)
pts81Cov <- left_join(pts81Cov, pts81TRI, by = "ID")

pts81TRIley <- terra::extract(triRiley,
                           pts81Vect,
                           na.rm = T,
                           # "simple": values for the cell a point falls in are returned
                           # "bilinear" returned values are interpolated from the values of the four nearest raster cells
                           method="simple",
                           ID = TRUE)
pts81Cov <- left_join(pts81Cov, pts81TRIley, by = "ID")

pts81TRmsd <- terra::extract(triRmsd,
                              pts81Vect,
                              na.rm = T,
                              # "simple": values for the cell a point falls in are returned
                              # "bilinear" returned values are interpolated from the values of the four nearest raster cells
                              method="simple",
                              ID = TRUE)
pts81Cov <- left_join(pts81Cov, pts81TRmsd, by = "ID")

pts81TPI <- terra::extract(tpi,
                           pts81Vect,
                           na.rm = T,
                           # "simple": values for the cell a point falls in are returned
                           # "bilinear" returned values are interpolated from the values of the four nearest raster cells
                           method="simple",
                           ID = TRUE)
pts81Cov <- left_join(pts81Cov, pts81TPI, by = "ID")

pts81Flow <- terra::extract(flowDir,
                            pts81Vect,
                            na.rm = T,
                            # "simple": values for the cell a point falls in are returned
                            # "bilinear" returned values are interpolated from the values of the four nearest raster cells
                            method="simple",
                            ID = TRUE)
pts81Cov <- left_join(pts81Cov, pts81Flow, by = "ID")

# repeat with NDWI
# read in raster
ndwi <- terra::rast("./Spatial/NDWI/Autumn_NDWI_Georef.vrt")
#plot(ndwi)
# reprojection takes some time, so get ready to wait
ndwiProj <- terra::project(ndwi, "EPSG:4326") # reproject into consistent project CRS (WGS84)
plot(ndwiProj)

# exract values at points
pts81NDWI <- terra::extract(ndwiProj,
                            pts81Vect,
                            na.rm = T,
                            # "simple": values for the cell a point falls in are returned
                            # "bilinear" returned values are interpolated from the values of the four nearest raster cells
                            method="simple",
                            ID = TRUE) %>% rename(ndwi = "Autumn_NDWI_Georef")

# add it to the final dataset with all spatial variables
pts81Cov <- left_join(pts81Cov, pts81NDWI, by = "ID")

# finally, save the derived data
write_rds(pts81Cov, "./Final Clean/pts81.rds")


##### Vectorise Raster Variables For Polygons #####
# extracting as polygons are more complicated
# left as an exercise for later if necessary
# plotVect <- vect(plotBounds)
# plotElev <- terra::extract(demFull,
#                            plotVect,
#                            na.rm = T,
#                            # function to summarize the extracted data by line or polygon geometry
#                            # fun=table to tabulate raster values for each line or polygon geometry
#                            # If weights=TRUE or exact=TRUE only mean, sum, min, max and table are accepted
#                            fun=mean,
#                            # "simple": values for the cell a point falls in are returned
#                            # "bilinear" returned values are interpolated from the values of the four nearest raster cells
#                            method="bilnear",
#                            # the approximate fraction of each cell that is covered
#                            weights = TRUE,
#                            # the exact fraction of each cell that is covered
#                            #exact = TRUE,
#                            ID = TRUE) %>% rename(elev = "EEA10_INSP_2011_03_31")
# # create a final covariate dataset
# plotBoundsCov <- left_join(plotBounds, plotElev, by = "ID")