##############################
library(tidyverse)
library(Hmisc)
library(ggcorrplot)
library(sf)
library(gganimate)
library(ggnewscale)
library(Polychrome)
library(GGally)


############### NVC & Species Data Merge ###############
# i believe height is in cm for nvcCleaned and in mm for specProc
# as-is, there's no sub-type NVC columns which were purged during creation of nvcCleaned; can go back and re-add
nvcCleaned <- read_rds("./Vegetation/nvcCleaned.rds")
specProc <- read_rds("./Vegetation/speciesProcessed.rds")
specProc <- specProc %>% 
  rename(pinDat = "date") %>% 
  # convert to cm for comparability
  mutate(pinHt = Height/10) %>% 
  select(-c("Treat","Height"))
nvcCleaned <- nvcCleaned %>% rename(Year = "yearS",
                                    nvcDat = "cleanDate") %>% select(-Date)
# haphazardly slap it all together
vegFull <- full_join(specProc, nvcCleaned, by = c("fullID", "Block","Plot", "Point","Year"))

# clean up the variables of interest
vegFull <- vegFull %>% mutate(nvcM1 = ifelse(Year == 2002,
                                             NVC,
                                             nvcM1),
                              # weird syntax to combine POSITx and Date classes
                              nvcDatFull = as.Date(ifelse(Year == 2002, pinDat, as.Date(nvcDat))),
                              # only 1 NVC for 2002 so fill that in
                              nvcSumComs = ifelse(Year == 2002, 1, nvcSumComs),
                              nvcVars = ifelse(Year == 2002, 1, nvcVars),
                              # pre-treatment offtake was always 3 --> "commercial stocking density of 2.7 ewes ha-1
                              # update offtake variables accordingly
                              sheepOff = ifelse(Year == 2002, 3, 
                                                ifelse(Year == 2003,
                                                       # sheep grazing started April 2003 and commercial otherwise
                                                       # 1/3 of normal offtake averaged with 2/3 * 2 = 2
                                                       ((sheepOff/3) + 2),
                                                       sheepOff)),
                              # no pre-treatment cow grazing and left out for the first year
                              # cattle grazing in September 2003, so left out for the second year too
                              cowOff = ifelse(Year <= 2003, 0, cowOff)) %>% 
  select(-c("NVC", "pointID", "plotID", "monthS", 'dayS'))

# c, sC, and v are whole numbers
# v >= sC >= c
# helper function that poorly generalises but suffices here
# numerically equivalent to Faith's phylogenetic diversity assuming that each divergence is 1 branch away
getDistSum <- function(c, sC, v){
  id <- as.numeric(paste0(c,sC,v))
  #dist <- v*(v-1)/2
  dist1 <- v - sC
  dist2 <- ifelse(v <= 2, sC - c,
                  ifelse(id == 244, 3,
                         ifelse(id == 223 | id == 333, 0,
                                ifelse(id == 234 | id == 123, 2,
                                       ifelse(id == 344 | id == 233, 1,
                                              # else is only id = 133
                                              3)))))
  dist3 <- ifelse(v == 0, 0,
                  ifelse(v <= 2, c - 1, 
                         ifelse(c == 1, 0,
                                ifelse(id == 233, c,
                                       ifelse(id == 344, 5, 
                                              # the else here is only id = 234, 333
                                              ifelse(id == 244, 3, sC))))))
  # playing with the weights is very important -- no past work has tried this before
  return(dist1 + 2*dist2 + 3*dist3)
}
# convert for logarithmic form
getDistSumLog <- function(c, sC, v){
  id <- as.numeric(paste0(c,sC,v))
  #dist <- v*(v-1)/2
  dist1 <- v - sC
  dist2 <- ifelse(v <= 2, sC - c,
                  ifelse(id == 244, 3,
                         ifelse(id == 223 | id == 333, 0,
                                ifelse(id == 234 | id == 123, 2,
                                       ifelse(id == 344 | id == 233, 1,
                                              # else is only id = 133
                                              3)))))
  dist3 <- ifelse(v == 0, 0,
                  ifelse(v <= 2, c - 1, 
                         ifelse(c == 1, 0,
                                ifelse(id == 233, c,
                                       ifelse(id == 344, 5, 
                                              # the else here is only id = 234, 333
                                              ifelse(id == 244, 3, sC))))))
  # playing with the weights is very important -- no past work has tried this before
  return(10*dist1 + 100*dist2 + 1000*dist3)
}
# get the mean rather than raw abundance
# equals the Mean Pairwise Distance (MPD)
getDistMean <- function(c, sC, v){
  id <- as.numeric(paste0(c,sC,v))
  #dist <- v*(v-1)/2
  dist1 <- v - sC
  dist2 <- ifelse(v <= 2, sC - c,
                  ifelse(id == 244, 3,
                         ifelse(id == 223 | id == 333, 0,
                                ifelse(id == 234 | id == 123, 2,
                                       ifelse(id == 344 | id == 233, 1,
                                              # else is only id = 133
                                              3)))))
  dist3 <- ifelse(v == 0, 0,
                  ifelse(v <= 2, c - 1, 
                         ifelse(c == 1, 0,
                                ifelse(id == 233, c,
                                       ifelse(id == 344, 5, 
                                              # the else here is only id = 234, 333
                                              ifelse(id == 244, 3, sC))))))
  # playing with the weights is very important -- no past work has tried this before
  return((dist1 + 2*dist2 + 3*dist3)/(dist1 + dist2 + dist3))
}
# finally get the variance
# equals the VPD (variation of pairwise distances)
getDistVar <- function(c, sC, v){
  id <- as.numeric(paste0(c,sC,v))
  #dist <- v*(v-1)/2
  dist1 <- v - sC
  dist2 <- ifelse(v <= 2, sC - c,
                  ifelse(id == 244, 3,
                         ifelse(id == 223 | id == 333, 0,
                                ifelse(id == 234 | id == 123, 2,
                                       ifelse(id == 344 | id == 233, 1,
                                              # else is only id = 133
                                              3)))))
  dist3 <- ifelse(v == 0, 0,
                  ifelse(v <= 2, c - 1, 
                         ifelse(c == 1, 0,
                                ifelse(id == 233, c,
                                       ifelse(id == 344, 5, 
                                              # the else here is only id = 234, 333
                                              ifelse(id == 244, 3, sC))))))
  # playing with the weights is very important -- no past work has tried this before
  return(wtd.var(c(3,2,1), c(dist3,dist2,dist1)))
}

# get NVC community-level presence/absence variables
vegFull <- vegFull %>% rowwise() %>% mutate(isM = ceiling(sum(str_detect(paste0(c_across(nvcM1:nvcM4)),"M"))/4),
                                            isU = ceiling(sum(str_detect(paste0(c_across(nvcM1:nvcM4)),"U"))/4),
                                            isH = ceiling(sum(str_detect(paste0(c_across(nvcM1:nvcM4)),"H"))/4),
                                            isW = ceiling(sum(str_detect(paste0(c_across(nvcM1:nvcM4)),"W"))/4),
                                            nvcComs = sum(c_across(isM:isW), na.rm = TRUE))
# these do not retail NAs so jerryrig an adjustment to reset to NA for cases where it's a problem
# could combine these together in the above call but this is a bit more efficient and readible
vegFull <- vegFull %>% mutate(isM = ifelse(is.na(nvcM1), NA, isM),
                              isU = ifelse(is.na(nvcM1), NA, isU),
                              isH = ifelse(is.na(nvcM1), NA, isH),
                              isW = ifelse(is.na(nvcM1), NA, isW),
                              nvcComs = ifelse(is.na(nvcM1), NA, nvcComs),
                              nvcSumComs = ifelse(is.na(nvcM1), NA, nvcSumComs),
                              nvcVars = ifelse(is.na(nvcM1), NA, nvcVars))
# play around with some novel diversity metrics
# all of these are heavily zero-inflated
# vegFull <- vegFull %>% rowwise() %>% mutate(# old version
#                                             #nvcDiv = ((nvcComs-1)*3 + (nvcSumComs-nvcComs)*2 + (nvcVars-nvcSumComs)),
#                                             nvcDiv0 = getDistSum(nvcComs, nvcSumComs, nvcVars),
#                                             #nvcDivAlt = getDistSumLog(nvcComs, nvcSumComs, nvcVars),
#                                             # the logged version is actually a much nicer predictor
#                                             #nvcDivAltLog = ifelse(nvcDivAlt == 0, 0, log(nvcDivAlt)),
#                                             # both of these struggle w/ NaNs
#                                             nvcDiv1 = getDistMean(nvcComs, nvcSumComs, nvcVars))#,
#                                             # variance looks to be a poor metric
#                                             #nvcDiv2 = getDistVar(nvcComs, nvcSumComs, nvcVars))#,
#                                             #nvcDiv2 = ifelse(is.nan(nvcDiv2Int), 0, nvcDiv2Int))


# check values
# works alright now so it's hidden
#vegTest <- vegFull %>% #filter(nvcVars >= 3) %>% 
#  mutate(uniqueComb = paste(nvcComs, nvcSumComs, nvcVars, sep = ";"))
#table(vegTest$uniqueComb, vegTest$nvcDiv)

# check co-occurence of the main variables
crossprod(as.matrix(vegFull %>% select(c("isM", "isU", "isH", "isW")) %>% filter(!is.na(isM))))
# see how this translates to correlation
cor(vegFull %>% select(c("isM", "isU", "isH", "isW")), use = "complete.obs")
#chisq.test(vegFull %>% select(c("isM", "isU", "isH", "isW")) %>% filter(!is.na(isM)))

# quick corrplot of species data to understand covariates
ggcorrplot(cor(specProc %>% select(where(is.numeric)), use = "complete.obs"))

# afforestation signal
#nrow(vegFull %>% filter(nTree >= 1))
#nrow(vegFull %>% filter(isW == 1))
#nrow(vegFull %>% filter(isW == 1 & nvcComs > 1))
#table(vegFull$isW, vegFull$nTree)
table(vegFull$isW, as.numeric(vegFull$nTree >= 1))

# vegFull <- vegFull %>% mutate(hasTree = ifelse(is.na(isW) & is.na(nTree), NA,
#                                                ifelse(ceiling(sum(isW, (nTree >= 1), na.rm = TRUE)/2) >= 1, 
#                                                1, 0)))

# get binary columns for all NVC sub-communities
vegTest <- vegFull
# a way to get a unique column for each NVC sub-community presence/absence
uniqueNVC <- sort(unique(append(unique(vegFull$nvcM1),
                                c(unique(vegFull$nvcM2),
                                  unique(vegFull$nvcM3),
                                  unique(vegFull$nvcM4)))))

# yes loops are bad...but this feels too trivial to define a whole function for
for (i in uniqueNVC) {
  isPresent <- colSums(rbind(str_detect(vegFull$nvcM1, i),
                             str_detect(vegFull$nvcM2, i),
                             str_detect(vegFull$nvcM3, i),
                             str_detect(vegFull$nvcM4, i)), 
                       na.rm=TRUE)
  # tests for NA to avoid casting NA to 0 and inflating data availability
  vegTest[, i] <-  ifelse(!is.na(vegFull$nvcM1),
                               ceiling(isPresent/4), # hack to turn any non-zero into 1
                               NA)
}



# generalisable helper function for community-level NVC turnover
nvcTurnoverCom <- function(H1, M1, U1, W1, H2, M2, U2, W2){
  # define needed variables
  tot1 <- sum(H1, M1, U1, W1, na.rm = TRUE)
  propH1 <- H1/tot1
  propM1 <- M1/tot1
  propU1 <- U1/tot1
  propW1 <- W1/tot1
  
  tot2 <- sum(H2, M2, U2, W2, na.rm = TRUE)
  propH2 <- H2/tot2
  propM2 <- M2/tot2
  propU2 <- U2/tot2
  propW2 <- W2/tot2
  
  # this simplification is fine when sample sizes are equal (same # NVCs surveyed at each point, which is true)
  hornInd <- 2*(H1*H2 + M1*M2 + U1*U2 + W1*W2)/((H1^2 + H2^2) + (M1^2 + M2^2) + (U1^2 + U2^2) + (W1^2 + W2^2))
  return(hornInd)
}

# manual brute force method to get subcommunity-level NVC turnover 
# NOT generalisable
nvcTurnoverSubCom <- function(H10_1, H12_1, H17_1, H18_1, H21_1, M11_1, M15_1, M16_1, M17_1, M18_1, M19_1, M20_1, M22_1, M23_1, 
                              M25_1, M26_1, M35_1, M6_1, U13_1, U19_1, U20_1, U4_1, U5_1, U6_1, W4_1, W6_1,
                              H10_2, H12_2, H17_2, H18_2, H21_2, M11_2, M15_2, M16_2, M17_2, M18_2, M19_2, M20_2, M22_2, M23_2, 
                              M25_2, M26_2, M35_2, M6_2, U13_2, U19_2, U20_2, U4_2, U5_2, U6_2, W4_2, W6_2){
  # define needed variables
  tot1 <- sum(H10_1, H12_1, H17_1, H18_1, H21_1, M11_1, M15_1, M16_1, M17_1, M18_1, M19_1, M20_1, M22_1, M23_1, 
              M25_1, M26_1, M35_1, M6_1, U13_1, U19_1, U20_1, U4_1, U5_1, U6_1, W4_1, W6_1, 
              na.rm = TRUE)
  propH10_1 = H10_1/tot1
  propH12_1 = H12_1/tot1
  propH17_1 = H17_1/tot1
  propH18_1 = H18_1/tot1
  propH21_1 = H21_1/tot1
  propM11_1 = M11_1/tot1
  propM15_1 = M15_1/tot1
  propM16_1 = M16_1/tot1
  propM17_1 = M17_1/tot1
  propM18_1 = M18_1/tot1
  propM19_1 = M19_1/tot1
  propM20_1 = M20_1/tot1
  propM22_1 = M22_1/tot1
  propM23_1 = M23_1/tot1
  propM25_1 = M25_1/tot1
  propM26_1 = M26_1/tot1
  propM35_1 = M35_1/tot1
  propM6_1 = M6_1/tot1
  propU13_1 = U13_1/tot1
  propU19_1 = U19_1/tot1
  propU20_1 = U20_1/tot1
  propU4_1 = U4_1/tot1
  propU5_1 = U5_1/tot1
  propU6_1 = U6_1/tot1
  propW4_1 = W4_1/tot1
  propW6_1 = W6_1/tot1
  
  tot2 <- sum(H10_2, H12_2, H17_2, H18_2, H21_2, M11_2, M15_2, M16_2, M17_2, M18_2, M19_2, M20_2, M22_2, M23_2, 
              M25_2, M26_2, M35_2, M6_2, U13_2, U19_2, U20_2, U4_2, U5_2, U6_2, W4_2, W6_2, 
              na.rm = TRUE)
  propH10_2 = H10_2/tot2
  propH12_2 = H12_2/tot2
  propH17_2 = H17_2/tot2
  propH18_2 = H18_2/tot2
  propH21_2 = H21_2/tot2
  propM11_2 = M11_2/tot2
  propM15_2 = M15_2/tot2
  propM16_2 = M16_2/tot2
  propM17_2 = M17_2/tot2
  propM18_2 = M18_2/tot2
  propM19_2 = M19_2/tot2
  propM20_2 = M20_2/tot2
  propM22_2 = M22_2/tot2
  propM23_2 = M23_2/tot2
  propM25_2 = M25_2/tot2
  propM26_2 = M26_2/tot2
  propM35_2 = M35_2/tot2
  propM6_2 = M6_2/tot2
  propU13_2 = U13_2/tot2
  propU19_2 = U19_2/tot2
  propU20_2 = U20_2/tot2
  propU4_2 = U4_2/tot2
  propU5_2 = U5_2/tot2
  propU6_2 = U6_2/tot2
  propW4_2 = W4_2/tot2
  propW6_2 = W6_2/tot2
  
  # this simplification is fine when sample sizes are equal (same # NVCs surveyed at each point, which is true)
  numerator = 2*(propH10_1*propH10_2 + propH12_1*propH12_2 + propH17_1*propH17_2 + propH18_1*propH18_2 + propH21_1*propH21_2 + 
                   propM11_1*propM11_2 + propM15_1*propM15_2 + propM16_1*propM16_2 +
                   propM17_1*propM17_2 + propM18_1*propM18_2 + propM19_1*propM19_2 + propM20_1*propM20_2 + propM22_1*propM22_2 +
                   propM23_1*propM23_2 + propM25_1*propM25_2 + propM26_1*propM26_2 + propM35_1*propM35_2 +
                   propM6_1*propM6_2 + propU13_1*propU13_2 + propU19_1*propU19_2 + propU20_1*propU20_2 + propU4_1*propU4_2 + 
                   propU5_1*propU5_2 + propU6_1*propU6_2 + propW4_1*propW4_2 + propW6_1*propW6_2)
  denominator = propH10_1^2 + propH12_1^2 + propH17_1^2 + propH18_1^2 + propH21_1^2 + propM11_1^2 + propM15_1^2 + propM16_1^2 +
    propM17_1^2 + propM18_1^2 + propM19_1^2 + propM20_1^2 + propM22_1^2 + propM23_1^2 + propM25_1^2 + propM26_1^2 + propM35_1^2 +
    propM6_1^2 + propU13_1^2 + propU19_1^2 + propU20_1^2 + propU4_1^2 + propU5_1^2 + propU6_1^2 + propW4_1^2 + propW6_1^2 + 
    propH10_2^2 + propH12_2^2 + propH17_2^2 + propH18_2^2 + propH21_2^2 + propM11_2^2 + propM15_2^2 + propM16_2^2 + propM17_2^2 +
    propM18_2^2 + propM19_2^2 + propM20_2^2 + propM22_2^2 + propM23_2^2 + propM25_2^2 + propM26_2^2 + propM35_2^2 + propM6_2^2 + 
    propU13_2^2 + propU19_2^2 + propU20_2^2 + propU4_2^2 + propU5_2^2 + propU6_2^2 + propW4_2^2 + propW6_2^2
  
  return(numerator/denominator)
}

# add the sub-community turnover in a clunky way
vegFull <- vegTest %>% mutate(ptCode = paste0(Plot, Block, Point)) %>% 
  group_by(ptCode) %>% 
  mutate(nvcTurn = nvcTurnoverCom(isH, isM, isU, isW, lag(isH), lag(isM), lag(isU), lag(isW)),
         nvcTurnSub = nvcTurnoverSubCom(H10, H12, H17, H18, H21, M11, M15, M16, M17, M18, M19, M20, M22, M23, M25, M26, M35, M6, U13,
                                        U19, U20, U4, U5, U6, W4, W6, 
                                        lag(H10), lag(H12), lag(H17), lag(H18), lag(H21), lag(M11), lag(M15), lag(M16), lag(M17), 
                                        lag(M18), lag(M19), lag(M20), lag(M22), lag(M23), lag(M25), lag(M26), lag(M35), lag(M6), 
                                        lag(U13), lag(U19), lag(U20), lag(U4), lag(U5), lag(U6), lag(W4), lag(W6)))
# test alternative parametrisations; not used
vegFull <- vegFull %>% mutate(nvcTurnWt = (nvcTurn*2 + nvcTurnSub)/3,
                               nvcTurnWt2 = (nvcTurn*10 + nvcTurnSub)/11) %>% ungroup()



# export merged data
write_rds(vegFull, "./Vegetation/vegCleaned.rds")



# quickly inspect possible response variables
vegCorr <- vegFull %>% #filter(nvcSumComs > 1) %>% 
  select(c("hill0", "hill1", "hill2", 
           "hillBryo0", "hillBryo1",
           "hillVasc0", "hillVasc1",
           #"nvcDiv0", #"nvcDivAlt", "nvcDivAltLog", 
           #"nvcDiv1", #"nvcDiv2",
           "isM", "isU", "isH", "isW",
           "nTree", #"hasTree",
           "nvcComs", "nvcSumComs", "nvcVars"))
# hill 1 and hill 2 are nearly perfectly correlated, so leaving out taxa-specific hill2 may not be an issue
# our nvcDiv is most related to hill0 and then hill1, hill2, and vascular metrics are tied; bad measure for bryophyte diversity
# played around a lot, and looks like lots of NVC correlation is inflated by values w/ only 1 sub-community (zero-inflation)
ggcorrplot(cor(vegCorr, use = "complete.obs"))
#ggplot(vegCorr, aes(x = hill0, y = nvcDiv1)) + geom_point()

# confirm if these relationships depend at all on the sampling points (1-25 vs 25+)
#vegCorNVC <- vegCorr %>% select(-c("hill0", "hill1", "hill2", 
#                                    "hillBryo0", "hillBryo1",
#                                    "hillVasc0", "hillVasc1"))
# looks like that part is OK
#ggcorrplot(cor(vegCorNVC, use = "complete.obs"))



# finally, add climate data
# rename variables ASAP to facilitate merging
climDat <- read_rds("./Climate/climBasic.rds") %>% rename(blockGroup = "block",
                                                          Year = "sampleYr") %>% 
  # both of these are always the same: every summer day is a grow day and there are no frost days
  select(-c("growDaysMin_Summer", "frostDays_Summer"))

# visualise the relationship
ggcorrplot(cor(climDat %>% ungroup() %>% select(where(is.numeric)), #& ends_with("Spring")), 
               use = "complete.obs"))
# create a block group by which to merge
vegFull <- vegFull %>% mutate(blockGroup = factor(Block, 
                                                  levels = LETTERS[1:6],
                                                  labels = c("AB", "AB", "CD","CD", "EF", "EF")))
# cleanly merge survey and climate data by block group
vegClim <- left_join(vegFull, climDat, by = c("blockGroup", "Year"))


############### Spatial Processing ############### 
# add spatial data (at the end bc sf is a pain)
points25 <- st_read("./Spatial/21points/25PointsCorrected.shp")
pts25 <- st_transform(points25, "EPSG:4326") # reproject to WGS84
pts25$lat <- st_coordinates(pts25)[,1]
pts25$long <- st_coordinates(pts25)[,2]
ggplot(pts25, aes(x = lat, y = long)) + 
  geom_point()

# repeat with 81 points
# read in pre-processed DEM covariates with pre-formatted spatial points data
pts81Cov <- read_rds("./Final Clean/pts81.rds")

# create a unified ID in the NVC
vegSpat <- vegClim %>% mutate(id = paste0(Block, Plot, Point),
                              blockGroup = factor(Block,
                                                  levels = LETTERS[1:6],
                                                  labels = c("A&B", "A&B", "C&D","C&D", "E&F", "E&F")))
# merge together
vegMerge <- left_join(vegSpat, pts81Cov, by = "id") 

# reset as sf object bc merging strips out the spatial formatting
coords <- c("lat", "long")
vegGeo <- st_as_sf(vegMerge %>% filter(id %in% unique(pts81Cov$id)), 
                   coords = coords, crs = "EPSG:4326") %>% select(-c("SAMPLE", "ID"))
# remake coordinates here bc they got wiped in the setting
vegGeo$lat <- st_coordinates(vegGeo)[,1]
vegGeo$long <- st_coordinates(vegGeo)[,2]

write_rds(vegGeo, "./Vegetation/vegGeo.rds")

# some minor issues where points don't have spatial coordinates because they're 81+ in full plots
#sum(is.na(vegGeo$lat))
#ptsMissed <- vegGeo %>% filter(is.na(lat))
#unique(ptsMissed$id)
#table(ptsMissed$id, ptsMissed$Year)
# FIXED: no more data thrown out for any missing points after manual map cleaning
# nvc1, all the heights and densities all missing the same values
#vegGeoIssues <- ptsMissed %>% filter(!is.na(`RIGHT Density`))
#vegGeoIssues <- vegGeoIssues %>% select(-c("White stick found?", "Period", "Year", "monthS", "dayS", "cleanDate", 
#                                           "fullID", "pointID", "nvc3", "nvc4", "nvcM1", "nvcM2", "nvcM3", "nvcM4", "whiteStick"))
#table(vegGeoIssues$id, vegGeoIssues$Year)

############### Spatial Plotting ############### 
# create visuals to understand spatiotemporal dynamics
data(Dark24)
basePlot <- ggplot(data = vegGeo, aes(x = lat, y = long, colour = as.factor(nvcM1))) +
  geom_point(size = 2) +
  theme_bw() + 
  labs(x = "",
       y = "",
       color = 'NVC') +
  facet_wrap(~fct_relevel(blockGroup, 
                          "C&D", "E&F", "A&B"), scales = "free") + 
  scale_colour_manual(values = as.vector(Dark24))


timePlot <- basePlot + transition_manual(Year) + 
  labs(title = "Survey Year: {current_frame}") +
  theme(strip.text.x = element_text(size = 16),
        plot.title = element_text(face = "bold",
                                  size = "20",
                                  hjust = 0.5))

# render as GIF in-window for viewing
animate(timePlot, renderer = gifski_renderer())
# export as GIF
anim_save("./Vegetation/NVCPlots/nvcPrimarySpatial.gif", 
          animation = timePlot,
          renderer = gifski_renderer(),
          height = 800, width = 1600)
# export as mkv (mp4 output not compatible on macOS but works in-browser and Windows media players)
anim_save("./Vegetation/NVCPlots/nvcPrimarySpatial.mkv", 
          animation = timePlot,
          renderer = ffmpeg_renderer(),
          height = 800, width = 1600)

# add plot boundaries to visuals
plotBounds <- st_read("./Spatial/boundaries/newplotsCleaned.shp")
plotBounds <- st_transform(plotBounds, "EPSG:4326") # reproject to WGS84
plotBounds <- plotBounds %>% mutate(plot = str_sub(Plot_ID, 2),
                                    block = str_sub(Plot_ID, 0, 1),
                                    blockGroup = factor(block, 
                                                        levels = LETTERS[1:6],
                                                        labels = c("A&B", "A&B", "C&D","C&D", "E&F", "E&F")))
#plotBounds <- plotBounds %>% rowwise() %>% mutate(lat = list(st_coordinates(geometry)[,1]),
#                                                  long = list(st_coordinates(geometry)[,2]))

# plot the 25 points and plot boundaries
ggplot() +
  geom_sf(data = plotBounds, aes(fill = plot)) + 
  scale_fill_grey() + 
  labs(fill = "Plot") + 
  ## New scale
  new_scale_colour() + new_scale_fill() +
  geom_point(data = pts25, aes(x = lat, y = long, colour = "#2E7D32"), inherit.aes = FALSE, size = 1) +
  scale_color_identity() +
  #scale_colour_manual(values = as.vector(Dark24)) + 
  theme_bw() + 
  labs(x = "",
       y = "",
       color = 'NVC',
       title = "Test Title") + 
  theme(plot.title = element_text(face = "bold",
                                  hjust = 0.5))

# visualise
ggplot() +
  geom_sf(data = plotBounds, aes(fill = plot)) + 
  scale_fill_grey() + 
  labs(fill = "Plot") + 
  ## New scale
  new_scale_colour() + new_scale_fill() +
  geom_point(data = vegGeo %>% filter(Year == 2023), aes(x = lat, y = long, colour = as.factor(nvcM1)), inherit.aes = FALSE, size = 2) +
  scale_colour_manual(values = as.vector(Dark24)) + 
  theme_bw() + 
  labs(x = "",
       y = "",
       color = 'NVC',
       title = "Test Title") + 
  theme(plot.title = element_text(face = "bold",
                                  hjust = 0.5))
# facets don't work with geom_sf()
#          +
#  facet_wrap(~fct_relevel(blockGroup, 
#                          "C&D", "E&F", "A&B")) 


basePlot2 <- ggplot() +
  geom_sf(data = plotBounds, aes(fill = plot)) + 
  scale_fill_grey() + 
  labs(fill = "Plot") + 
  ## New scale
  new_scale_colour() + new_scale_fill() +
  geom_point(data = vegGeo, aes(x = lat, y = long, colour = as.factor(nvcM1)), inherit.aes = FALSE, size = 0.5) +
  scale_colour_manual(values = as.vector(Dark24)) + 
  theme_bw() + 
  labs(x = "",
       y = "",
       color = 'NVC')

timePlot2 <- basePlot2 + transition_manual(Year) + 
  labs(subtitle = "Survey Year: {current_frame}") +
  theme(legend.title = element_text(face = "bold",
                                    size = 18),
        legend.text = element_text(size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = "20",
                                     face = "bold")) + 
  guides(colour = guide_legend(override.aes = list(shape = 15,
                                                   size = 6)))


# render as GIF in-window for viewing
animate(timePlot2, renderer = gifski_renderer())
# export as GIF
anim_save("./Vegetation/NVCPlots/nvcPrimarySpatialBlocks.gif", 
          animation = timePlot2,
          renderer = gifski_renderer(),
          width = 4000,
          height = 3600,
          res = 300)
# export as mkv (mp4 output not compatible on macOS but works in-browser and Windows media players)
anim_save("./Vegetation/NVCPlots/nvcPrimarySpatialBlocks.mkv", 
          animation = timePlot2,
          renderer = ffmpeg_renderer(),
          width = 4000,
          height = 3600,
          res = 300)


# render as GIF in-window for viewing
animate(timePlot, renderer = gifski_renderer())