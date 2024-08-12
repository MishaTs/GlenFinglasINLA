##############################
library(tidyverse)
library(readxl)
library(zoo)
library(lubridate)
library(ggcorrplot)

##############################
# import species-specifc data
specRaw <- read_delim("Vegetation/Plant_species/Pakeman_Data_updated2023.csv", 
                      delim = ";", escape_double = FALSE, trim_ws = TRUE) #%>% select(c-())
# these are the columns containing specific labels
colnames(specRaw)[9:192]
# derive some biodiversity indices? we already have 


##### Clean Covariates #####
specFull2017 <- read_excel("Vegetation/Plant_species/Raw/LIESE-POINTQUADRAT_RAWDATA~_2017.xlsx", skip = 3, 
                           col_types = "text")
# empty: ...39, ...225, ...226
specFull2017 <- specFull2017 %>% select(-c("...39", "...47", "...225", "...226"))
# a bit of sanity checking in Excel is removed here
# there is one row in 2023 that doesn't match but only by 1 observation so can realistically caveat it out for now
specFull2017 <- specFull2017 %>% select(-c("check", "...229"))
# 224 columns is intimidating
#colnames(specFull2017)
# 14th to 29th unused
#summary(specFull2017)
# everything is ok here except the date
#unique(specFull2017$...1)
# this one value doesn't have month/day and messes up cleaning unless manually corrected
specFull2017 <- specFull2017 %>% 
  mutate(across(c("...1"), .fns = ~replace(., . ==  ". .05" , NA))) %>% 
  mutate(rawDate = str_replace(`...1`, "0705","07/05")) %>% relocate(rawDate)
#unique(specFull2017$rawDate)
# now process the rest
specFull2017 <- specFull2017 %>% rowwise() %>% mutate(date = as.Date(ifelse(str_detect(rawDate,"[[:digit:]]{5}"),
                                                                    as.Date(as.numeric(rawDate), origin = "1899-12-30"),
                                                                    as.Date(format(as.Date(rawDate, 
                                                                                           tryFormats = c("%d.%m.%Y",
                                                                                                          "%d/%m/%Y")), 
                                                                                   "20%y-%m-%d"))),
                                                                    origin ="1970-01-01")) %>% relocate(c(date)) %>% 
  # throw out all presence/absence info and get the first pin hit for reference
  select(c("date","Year","Block","Plot","Point","Distance","Max ht", "Ground ht", "height diff.","1st",
           "Dominant", "Dung", "Grazed", "Tussockiness", "Flowers", "Litter", "Comments"))

# impute date based on alphabetic order to compensate for excel formatting
specFull2017$cleanDate <- na.locf(specFull2017$date, na.rm = FALSE)

# finally, fix some incorrect year stuff (DIV Point 21 is set to 2010 for some reason)
specFull2017 <- specFull2017 %>% mutate(cleanDate = as.Date(ifelse(Year == 2010,
                                                           cleanDate %m+% years(1),
                                                           cleanDate)),
                                        Year = ifelse(Year == "2010",
                                                      "2011",
                                                      Year))

# create unique ID
specFull2017 <- specFull2017 %>% mutate(fullID = paste0(Year,
                                                        "_",
                                                        Block,
                                                        "_",
                                                        as.numeric(as.roman(Plot)), 
                                                        "_",
                                                        ifelse(Point < 10, "0", ""), Point,
                                                        "_",
                                                        ifelse(Distance < 10, 
                                                               "00", 
                                                               ifelse(Distance < 100, 
                                                                      "0", "")), Distance)) %>% 
  arrange(fullID)
# clean up columns
specFull2017 <- specFull2017 %>% select(-c("date"))
# cleaned
#View(specFull2017)
#colnames(specFull2017)


# repeat for 2020
# 2020 data has no summary column, would need to import all sheets
#specFull2020 <- read_excel("Vegetation/Plant_species/Raw/Pin_Frame_Data_2020.xls")

# create method to deal with multiple excel sheets
read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}
# read everything in as a list with multiple contained dataframes
specFull2020Raw <- read_excel_allsheets("Vegetation/Plant_species/Raw/Pin_Frame_Data_2020.xls")
# get all unique column names since they differ between years
spec2020Names <- c()
for(i in length(specFull2020Raw):1){
  spec2020Names <- append(spec2020Names, colnames(specFull2020Raw[[i]]))
}
# create empty dataframe with all unique column names for a unified dataframe
specFull2020 <- read_csv(I(paste(unique(spec2020Names), collapse=",")))
# loop through and bind each sheet to the master dataframe
for(i in 1:length(specFull2020Raw)){
  specFull2020 <- bind_rows(specFull2020, specFull2020Raw[[i]] %>%
                            # convert everything to character because type mismatches otherwise
                            mutate(across(everything(), as.character)))
}
# clean columns up-front this time, everything is the same except no "Year" column
specFull2020 <- specFull2020 %>% select(c("Date","Block","Plot","Point","Distance","Max ht", "Ground ht", "height diff.","1st",
                                          "Dominant", "Dung", "Grazed", "Tussockiness", "Flowers", "Litter", "Comments")) %>% 
  # Date is a lot better behaved with one weird value that's clearly an NA
  mutate(rawDate = str_replace(Date, "`", NA_character_)) %>% relocate(rawDate)

# now process dates
specFull2020 <- specFull2020 %>% rowwise() %>% mutate(date = as.Date(ifelse(str_detect(rawDate,"[[:digit:]]{5}"),
                                                                            as.Date(as.numeric(rawDate), origin = "1899-12-30"),
                                                                            as.Date(format(as.Date(rawDate, 
                                                                                                   tryFormats = c("%Y-%m-%d",
                                                                                                                  "%d/%m/%Y",
                                                                                                                  "%d/%m//%Y")), 
                                                                                           "20%y-%m-%d"))),
                                                                     origin ="1970-01-01")) %>% relocate(c(date)) %>% 
  # get rid of the pre-processed stuff since it looks legit
  select(-c("rawDate", "Date"))
# impute date based on alphabetic order to compensate for excel formatting
specFull2020$cleanDate <- na.locf(specFull2020$date, na.rm = FALSE)
# create Year column to match other data
specFull2020 <- specFull2020 %>% mutate(Year = year(cleanDate))
# generate ID
specFull2020 <- specFull2020 %>% mutate(fullID = paste0(Year,
                                                        "_",
                                                        Block,
                                                        "_",
                                                        as.numeric(as.roman(Plot)), 
                                                        "_",
                                                        ifelse(Point < 10, "0", ""), Point,
                                                        "_",
                                                        ifelse(Distance < 10, 
                                                               "00", 
                                                               ifelse(Distance < 100, 
                                                                      "0", "")), Distance)) %>% 
  arrange(fullID) %>% select(-c("date"))

# cleaned
#View(specFull2020)
#colnames(specFull2020)

# 2023 looks broadly fine except for one mismatch in the "check" column, same as prior to 2017
# for our purposes this is fine and good to use
specFull2023 <- read_excel("Vegetation/Plant_species/Raw/Pin_Frame_Data_2023.xlsx", 
                           sheet = "Summary", skip = 3)
# clean columns up-front this time, everything is the same except no "Year" column
specFull2023 <- specFull2023 %>% select(c("Date...1","Block...2","Plot...3","Point...4","Distance...5",
                                          "Max ht", "Ground ht", "height diff.","1st",
                                          "Dominant", "Dung", "Grazed", "Tussockiness", "Flowers", "Litter", "Comments")) %>% 
  # fix column names to match
  rename(Date = "Date...1",
         Block = "Block...2",
         Plot = "Plot...3",
         Point = "Point...4",
         Distance = "Distance...5")
# repeat date cleaning, but date is already cleaned
#class(specFull2023$Date)
# convert from POSIXct to Date for formatting
specFull2023 <- specFull2023 %>% mutate(rawDate = as.Date(Date))
# impute date based on alphabetic order to compensate for excel formatting
specFull2023$cleanDate <- na.locf(specFull2023$rawDate, na.rm = FALSE)
# create Year column to match other data
specFull2023 <- specFull2023 %>% mutate(Year = year(cleanDate))
# generate ID
specFull2023 <- specFull2023 %>% mutate(fullID = paste0(Year,
                                                        "_",
                                                        Block,
                                                        "_",
                                                        as.numeric(as.roman(Plot)), 
                                                        "_",
                                                        ifelse(Point < 10, "0", ""), Point,
                                                        "_",
                                                        ifelse(Distance < 10, 
                                                               "00", 
                                                               ifelse(Distance < 100, 
                                                                      "0", "")), Distance)) %>% 
  arrange(fullID) %>% 
  select(-c("Date", "rawDate"))
# cleaned
#View(specFull2023)
#colnames(specFull2023)

# combine it all together
specFullRaw <- rbind(specFull2017, specFull2020, specFull2023)
# take out species from dung column
specFullRaw <- specFullRaw %>% rowwise() %>% mutate(dungClean = as.numeric(str_extract(Dung,"[[:digit:]]")))
# clean Flowers, Grazed, and Tussockiness,
specFullRaw <- specFullRaw %>% mutate(grazedClean = as.numeric(str_extract(Grazed,"[[:digit:]]")),
                                      # these columns have no strings in them anyways, so that's fine
                                      flowerClean = as.numeric(ifelse(Flowers > 3, 3, Flowers)),
                                      tussockClean = as.numeric(ifelse(Tussockiness > 2, 2, Tussockiness)))

# merge
# i'm going with sums as better than either averages/variance in separate variables
# there's a few missing values that pose a challenge, so I rescale assuming that the missing point is equal to the average
# this applies to 4 for all covariates (2008_A_IV_25; 2020_A_II_25; 2023_F_III_14) 
# grazing has this specifically on 2020_C_IV_1
# trussock on 2005_C_II_20
# flower specifically have this uniquely twice (2005_F_I_10; 2020_A_II_13)
specCovs <- specFullRaw %>% group_by(Year, Block, Plot, Point) %>% 
  summarise(date = min(cleanDate),
            graze = (sum(grazedClean, na.rm = TRUE) * 5/sum(!is.na(grazedClean))),
            flower = (sum(flowerClean, na.rm = TRUE) * 5/sum(!is.na(flowerClean))),
            tussock = (sum(tussockClean, na.rm = TRUE) * 5/sum(!is.na(tussockClean))),
            dung = (sum(dungClean, na.rm = TRUE) * 5/sum(!is.na(dungClean))))
# replace all NaN with NA
# don't run on date bc imputed, but also because it strips out the formatting
specCovs <- specCovs %>% mutate(across(!starts_with("date"),
                                        ~ifelse(is.nan(.), NA, .))) %>% 
  # remove NAs that will mess things up
  filter(!is.na(Block)) %>% 
  # create unique ID for merging
  mutate(fullID = paste0(Block, 
                         as.numeric(as.roman(Plot)), 
                         "_", 
                         # for some reason Point doesn't work as a string
                         ifelse(as.numeric(Point) < 10, "0", ""), 
                         Point, 
                         Year))

# process the original data in parallel 
specRaw <- specRaw %>% mutate(pointTrue = ifelse(Point %% 25 == 0,
                                                 25,
                                                 Point %% 25),
                              fullID = paste0(Block, 
                                              as.numeric(as.roman(Trt)), 
                                              "_", 
                                              ifelse(pointTrue < 10, "0", ""), 
                                              pointTrue, 
                                              Year))
# merge together and separate out crucial columns
specFull <- left_join(specRaw, specCovs, by = "fullID") %>% relocate(c(fullID, 
                                                                       date,
                                                                       Year.x, Year.y, 
                                                                       Block.x, Block.y,
                                                                       Trt,
                                                                       Plot.x, Plot.y,
                                                                       Treat,
                                                                       pointTrue,
                                                                       Point.x, Point.y)) %>% 
  select(-c(Year.y, Block.y, Plot.x, Plot.y, Point.x, Point.y)) %>% 
  rename(Year = "Year.x",
         Block = "Block.x",
         Plot = "Trt",
         Point = "pointTrue")
# dung column still looks a bit weird, so I don't think we'll use the special one
#summary(specFull %>% select(D, dung))
specFull <- specFull %>% select(-c(dung))
# at this point the data is CLEAN
#summary(specFull)



##### Get Target Variables of Interest #####
# derive some basic species richness and hill numbers
# first clean up the data a bit more
specProc <- specFull %>% rename(bareGround = "Bg",
                                litter = "Lt",
                                dung = "D",
                                mud = "M",
                                rock = "Rk",
                                water = "W",
                                wood = "Wd") %>% 
  relocate(c(bareGround, dung, litter, mud, rock, water, wood), .after = "Sum") %>% 
  rowwise() %>% 
  mutate(nTree = sum(c_across(any_of(c("Bpd", "Soa", "Sx", "Sxa"))), na.rm = TRUE))

# function to calculate hill number
# special case of alpha = 1 is treated differently as it is undefined at exactly alpha = 1, but tends towards this value
#hill <- function(pis, alpha = 0){
#  if(alpha==1){
#    H <- exp(-1*sum(pis*log(pis)))
#  }
#  if(alpha!=1){
#    H <- sum(pis^alpha)^(1/(1-alpha))
#  }
#  return(H)
#}

# a function that also converts the counts to pis and then calculates hill numbers
#ni_hill <- function(nis, alpha = 2){
#  pis <- nis / sum(nis)
#  h <- hill(pis, alpha = alpha)
#  return(h)
#}

# richness is just hill number = 0, but our function doesn't work that nicely
# sanity chceck confirms that it does work same as hillR for q = 2; simpson's diversity
#specProc <- specFull %>% rowwise() %>% mutate(nSpec = sum(c_across(Ac:Are)!=0),
#                                              hillInv2 = ni_hill(c_across(Ac:Are), alpha = -0.5))

# columns 10 to 187; Ac to Are are the species data of interest
library(hillR)
# get some basic hill numbers from hillR
# rarefaction, extrapolation, functional trait-based stuff from vegan(), BiodiversityR(), and hillR may come up in the future
# this code is just a template for some basic response variables and its syntax
# can probably be done well in tidy format too but was too finnicky to spend time on
specProc$hill0 <- hill_taxa(comm = specProc %>% select(Ac:Are), q = 0)
specProc$hill1 <- hill_taxa(comm = specProc %>% select(Ac:Are), q = 1)
specProc$hill2 <- hill_taxa(comm = specProc %>% select(Ac:Are), q = 2)


# taxonomy is probably tangential here, but functional traits?
# LEDA database
# relevant traits are 
# Canopy height
# Canopy structure
# Leaf Dry Matter Content
# Leaf Size
# Specific Leaf Area
# Seed mass

# BioFlor variables no longer available
# Canopy Structure
# Vegetative Spread
# Flowering Start (month)

# Junus articulatus and acutiflorus both exist; need to pick one?
specDict <- read_delim("Vegetation/Plant_species/SpeciesDict.csv", 
                       delim = ";", escape_double = FALSE, col_names = FALSE, 
                       trim_ws = TRUE) %>% 
  select(c("X1","X2","X3","X4"))
colnames(specDict) <- c("code", "commonName", "mergeName", "origName")

# 10 "sp" variants as well
# also "Liverwort" as a category is a hard one
# Poacea == grasses is another one that's going to be hard to classify
# Mnium and Sphagnum are the only sp variants that are non-vascular (i.e., bryophytes); other 8 are all vascular/herbaceous
specDict <- specDict %>% mutate(isSp = as.numeric(str_detect(tolower(mergeName), "\\bsp[.p]") | 
                                                    str_detect(tolower(mergeName), "\\bsp$")))
#View(specDict)

# speculative but Thymus polytrichus reclassified as serpyllum
# Carex ovalis not in dict but leporina is; ovalis removed from coding
# Taraxacum officinale not listed but Taraxacum palustre is: accoridng to Kew these may be the same? 
# https://powo.science.kew.org/taxon/254151-1
cehVascular <- read_csv("../CEH Plant Database/Vascular/data/BI_main.csv") %>% 
  rename(mergeName = "taxon_name_binom")

specDat <- left_join(specDict, cehVascular, by = "mergeName") %>% mutate(hasValVasc = as.numeric(!is.na(kew_id)))
sum(!is.na(specDat$kew_id))

cehBryo <- read_excel("~/Desktop/Blass/StAndrews/Classes/Dissertation/CEH Plant Database/Bryophyte/Bryoatt_updated_2017.xls") %>% 
  rename(mergeName = "Taxon name")

specDat2 <- left_join(specDict, cehBryo, by = "mergeName") %>% mutate(hasValBry = as.numeric(!is.na(Name_new)))
sum(!is.na(specDat2$Name_new))

# species without a name are....?
nrow(specDict) - sum(!is.na(specDat2$Name_new)) - sum(!is.na(specDat$kew_id))
# which are missing
# 8 of these aren't plants
# Parmelia saxatilis is a lichen too...which is complicated
# 4 more as "Herb1", "Fungus", "Moss", "Lichen"
# 10 Sps
# so...8 + 4 + 1 + 10 ~ approx 23 are expected to not be identified
validDat <- left_join(specDict, specDat %>% select(mergeName, hasValVasc), by = "mergeName") %>% 
  left_join(., specDat2 %>% select(mergeName, hasValBry), by = "mergeName") %>% rowwise() %>% 
  mutate(hasVal = sum(c_across(hasValVasc:hasValBry), na.rm = TRUE))

#View(validDat %>% filter(hasVal == 0))# & isSp != 1))


# get the codes corresponding to each taxa subgroup
vascCols <- as.vector(validDat %>% filter(hasValVasc == 1 | mergeName %in% c("Carex sp.", "Euphrasia sp", "Hieracium sp.", 
                                                                             "Juncus articulatus / Juncus acutiflorus", "Luzula sp.", 
                                                                             "Poa sp.", "Poacea", "Salix sp.", "Viola spp.", 
                                                                             "Crepis spp", "Herb1")))$code
bryoCols <- as.vector(validDat %>% filter(hasValBry == 1 | code %in% c("Mn", "Sp", "L", "Ms")) %>% select(code))$code

# get hill numbers for specific taxa
specProc$hillBryo0 <- hill_taxa(comm = specProc %>% select(any_of(bryoCols)), q = 0)
specProc$hillBryo1 <- hill_taxa(comm = specProc %>% select(any_of(bryoCols)), q = 1)
specProc$hillBryo2 <- hill_taxa(comm = specProc %>% select(any_of(bryoCols)), q = 2)

specProc$hillVasc0 <- hill_taxa(comm = specProc %>% select(any_of(vascCols)), q = 0)
specProc$hillVasc1 <- hill_taxa(comm = specProc %>% select(any_of(vascCols)), q = 1)
specProc$hillVasc2 <- hill_taxa(comm = specProc %>% select(any_of(vascCols)), q = 2)
# Inverse Simpson's diversity breaks (hill2) breaks when no species are observed
# might be better to convert to 1 - D, so 1 - 1/hill_taxa with a catch exception for "Inf" values (meaning that hill = 0...?)
# but is observing no species really indicative of "perfect" diversity?
# remove taxa-specific hill2 for now
specProc <- specProc %>% select(-c("hillVasc2", "hillBryo2"))

View(specProc)
# remove raw data 
specExport <- specProc %>% select(-c(Ac:Are))


# vascular plant traits from Pakeman et al. 2020
# rescale the traits to eliminate outsize influence of large values in euclidean distance?
# not needed as differences are marginal
vascTraits <- read_table("./Vegetation/Plant_species/Traits.txt")

# for consistency, remove flowering start month, Rhizome, Stolon?


# Excluded: Herb1, Urtica dioica, Galium aparine, Hyacinthoides non-scripta, Carex dioica, Arrhenatherum elatius
# these were all only recorded once except for Galium aparine (3 times), so this is fine
#vascCols[!(vascCols %in% vascTraits$code)]

specVascTraits <- specProc %>% select(any_of(append(vascTraits$code, "fullID"))) %>% 
  mutate(sum = rowSums(across(where(is.numeric)))) %>% 
  filter(sum != 0) %>% 
  select(-sum)

vascFD <- data.frame(
  hill_func(comm = specVascTraits %>% select(-fullID),
            traits = vascTraits %>% column_to_rownames(var="code"),
            traits_as_is = FALSE,
            q = 1))
# Rao's Q is not symmetric, so it's not great to use within the confines of our analysis
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0100014#pone.0100014.e005
vascFD_t <- data.frame(t(vascFD)) %>% select(-Q)
# re-add IDs to resulting metrics
vascFD_t$fullID <- specVascTraits$fullID
# when only 1 species is present, all metrics are 0
# not too unreasonable to assume that for "diversity"'s sake, maybe this is fine to force to 0 as well?
#View(specExport %>% filter(fullID %in% unique(vascFD_t %>% filter(D_q == 0) %>% select(fullID))$fullID))
# missing values had no vascular plants (even including the ones we left out for insufficient observations)
#View(specExport %>% filter(!(fullID %in% vascFD_t$fullID)))


colnames(vascFD_t)[1:4] <- paste0("vasc", colnames(vascFD_t)[1:4])





# bryophyte traits derivation
bryoTraits <- left_join(specDict, cehBryo, by = "mergeName") %>% 
  filter(!is.na(Name_new) | code %in% c("L", "Ms")) %>%
  # remove all columns where everything is NA
  select(where(~!all(is.na(.x)))) %>% 
  # remove Number of 10-km squares (hectads) in GB, IE, and Channel Islands; proxy for commonness?
  select(-ends_with("no")) %>% 
  # remove climate data because we're modelling that anyways; leave maximum altitude for now?
  # include for consistency w/ vascular plant traits?
  #select(-c("TJan", "TJul", "Prec")) %>% 
  # remove unnecessary labels 
  select(-c("commonName", "mergeName", "isSp", "Name_new","Concept", "BRC_num")) %>% 
  # we don't care about taxonomy, so remove Order
  select(-Ord) %>% 
  # all native; 1 non-native (Campylopus introflexus), so this isn't very relevant to us either
  select(-Stat) %>% 
  # little variation in perennial/annual life history. All perennial but Ceratodon purpureum (sometimes annual)
  select(-Per) %>% 
  # too many omissions in Tubers, Gemmae, Branches, and Leaves
  select(-c("Tub", "Gem", "Bra", "Lvs")) %>% 
  # remove geographic indicators 
  select(-c("E", "W", "Sc", "NI", "IR", "CI")) %>% 
  # too many missing values in maturity/fruiting
  select(-c("Spbeg", "Spend", "Spbeg2", "Spend2")) %>% 
  # we can't average categorical traits, so remove life form, sex
  select(-c("LF1", "LF2", "Sex")) %>% 
  # salt and heavy metal tolerance not relevant to this data
  select(-c("S", "HM")) %>% 
  # convert fruiting sporophyte frequency to continuous (Abundant = 4; Frequent = 3; Occasional = 2; Rare = 1)
  mutate(FrNum = as.numeric(factor(Fr,
                                   levels = c("A", "F", "O", "R"),
                                   labels = c(4, 3, 2, 1)))) %>% select(-Fr) %>% 
  # replace all NAs w/ 0 --> Substrate and EUNIS habitat classes NA suggest no evidence of appearing and set to 0
  mutate(across(Len:J4, .fns = ~replace_na(.,0))) 

# missing values
#colSums(is.na(bryoTraits))

# do average computation in a non-tidy way because it's more intuitive for large "moss" and "liverwort" groups
# Liverworts: Pep, Ptc, Lb
bryoTraits[bryoTraits$code == "L", 4:ncol(bryoTraits)] <- as.list(colMeans(bryoTraits %>% 
                                                                             filter(ML == "L") %>% 
                                                                             select(-c("code", "origName", "ML")), na.rm=TRUE))
# Mosses: everything except liverworts
bryoTraits[bryoTraits$code == "Ms", 4:ncol(bryoTraits)] <- as.list(colMeans(bryoTraits %>% 
                                                                             filter(ML == "M") %>% 
                                                                             select(-c("code", "origName", "ML")), na.rm=TRUE))

# remove categorical columns
bryoTraits <- bryoTraits %>% select(-c("ML", "origName"))
# rescaling only changes values at the 11th decimal point and later, so not necessary --> this is beyond our anticipated precision
#bryoTraits_Alt <- bryoTraits %>% mutate_if(is.numeric, scale)

# matrix is not symmetric
#isSymmetric.matrix(as.matrix(bryoTraits_Alt %>% column_to_rownames(var="code")))

# compute functional biodiversity
specBryoTraits <- specProc %>% select(any_of(append(bryoTraits$code, "fullID"))) %>% 
  mutate(sum = rowSums(across(where(is.numeric)))) %>% 
  filter(sum != 0) %>% 
  select(-sum)

bryoFD <- data.frame(
  hill_func(comm = specBryoTraits %>% select(-fullID),
            traits = bryoTraits %>% column_to_rownames(var="code"),
            traits_as_is = FALSE,
            q = 1))
# Rao's Q is not symmetric, so it's not great to use within the confines of our analysis
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0100014#pone.0100014.e005
bryoFD_t <- data.frame(t(bryoFD)) %>% select(-Q)
# re-add IDs to resulting metrics
bryoFD_t$fullID <- specBryoTraits$fullID
colnames(bryoFD_t)[1:4] <- paste0("bryo",colnames(bryoFD_t)[1:4])

specFuncDiv <- left_join(specExport, vascFD_t, by = "fullID") 
specFuncDiv <- left_join(specFuncDiv, bryoFD_t, by = "fullID") %>% 
  select(-c("bryoFDis", "vascFDis")) %>% 
  relocate(hillBryo0:hillBryo1, .after = vascFD_q)

specDivPlot <- specFuncDiv  %>% 
  relocate(hill2, .after = "bryoFD_q") %>% 
  rename(`Full richness` = "hill0",
         `Full diversity` = "hill1",
         `Vasc richness` = "hillVasc0", 
         `Vasc diversity` = "hillVasc1", 
         `Vasc FD` = "vascD_q",   
         `Vasc mean FD` = "vascMD_q",  
         `Vasc total FD` = "vascFD_q", 
         `Bryo richness` = "hillBryo0", 
         `Bryo diversity` = "hillBryo1", 
         `Bryo FD` = "bryoD_q", 
         `Bryo mean FD` = "bryoMD_q",    
         `Bryo total FD` = "bryoFD_q")

# ggcorrplot(cor(specDivPlot %>% select(`Full richness`:`Bryo total FD`), use = "complete.obs"),
#            type = "upper",
#            legend.title = "Correlation",
#            colors = c("#fde725","#21918c", "#440154"),
#            lab = TRUE,
#            tl.col = "#fde725")

library(reshape2)

ggplot(melt(cor(specDivPlot %>% select(`Full richness`:`Bryo total FD`), use = "complete.obs")), aes(Var1, Var2, fill=value)) +
  geom_tile(height=0.8, width=0.8) +
  scale_fill_gradient2(low="#fde725", mid="#21918c", high="#440154", midpoint = 0.55) +
  theme_minimal() +
  coord_equal() +
  labs(x="",y="",fill="Correlation") +
  theme(axis.text.x=element_text(size=13, angle=45, vjust=1, hjust=1, 
                                 margin=margin(-3,0,0,0)),
        axis.text.y=element_text(size=13, margin=margin(0,-3,0,0)),
        panel.grid.major=element_blank()) 


# actual traits across both species --> canopy height, Ellenberg values (L, F, R, N)





# save the final processed data
write_rds(specFuncDiv, "./Vegetation/speciesProcessed.rds")
