##############################
library(tidyverse)
library(zoo)
library(ggcorrplot)
library(sf)
library(gganimate)
library(ggnewscale)
library(Polychrome)
library(GGally)
library(viridis)

##############################
# batch import raw NVC data from annual surveys
file_paths <- paste0("./Vegetation/NVC/", list.files(path='./Vegetation/NVC'))
nvcRaw <- file_paths %>% 
  map_dfr(~ read_csv(.x, col_types = cols(Dung = col_character(),
                                          Other = col_character())))

# fix data entry problems in U4 to make dummy variable
nvcRaw$U4[nvcRaw$U4 > 1] <- 1
# create IDs
nvc <- nvcRaw %>% mutate(pointID = paste0(Block, as.numeric(as.roman(Plot)), "_", ifelse(Point < 10, "0", ""), Point),
                         plotID = paste0(Block, Plot),
                         # nothing ever goes smoothly and these dates are no exception
                         cleanDate = make_datetime(Year, Month, Day),
                         fullID = paste0(pointID,Year))

# main NVC columns
majorNVC <- colnames(nvc)[c(12:22,32:46)]
# get a sum of NVCs
nvc <- nvc %>% mutate(nvcTot = rowSums(across(any_of(majorNVC)), na.rm = TRUE))


# deal with Other and additional NVC columns later


# row numbers of multiple NVC
# both of these have empty Other and additional NVC columns
nvcDouble <- which(grepl("2", nvc$nvcTot))
# loop through the duplicate IDs; another hack-y approach
for(i in nvcDouble) {
  # select columns where the value is 1; could use median, mean, sum, etc for this
  # no other columns have a value of 1 except the ones we're interested in
  nvcDups <- colnames(nvc[i,] %>% select_if(function(col) mean(col) == 1))
  # get the column indices for base R references
  colInd <- grep(nvcDups[2], colnames(nvc))
  # add the second one to additional NVC column for processing later
  nvc$`additional NVC`[i] <- nvcDups[2]
  # replace the 1 with NA for the second match to resolve the issue
  nvc[i,colInd] <- NA
}
# safely remove the diagnostic column
nvc <- nvc %>% select(-nvcTot)

# attempt to 
# melt() from reshape2 deleted 568 rows
# full code for that found @ lines 124-127 as of NVC_density_height_dung.R
# retry with pivot_longer, but same results
# looks like rows with an NVC in "Other" or "additional NVC" may have non-NA 0 values in main NVC columns
# conditionally replacing 0s with NAs when those other columns are non-NA may help?
# left as an exercise for the reader
#nvcTest <- nvc %>% pivot_longer(
#  cols = majorNVC, 
#  names_to = "nvcClean",
#  values_to = "nvcCheck") %>% 
#  filter(nvcCheck != 0)

# combine all rows together in a hack-y way
# blank placeholder column
nvc <- nvc %>% mutate(nvcClean = NA)
# yes loops are bad...but this feels too trivial to define a whole function for
for (i in majorNVC) {
  nvc$nvcClean[nvc[[i]] == 1] <- i
}

###############
# unique ID issues are next
# 40894 rows originally
#nrow(nvcRaw)
# 31 duplicate plot, point, times
#nrow(nvcRaw) - n_distinct(nvc$fullID)
# manual check shows that these represent unsurveyed points with absolutely no data of any kind
#View(nvc %>% filter(duplicated(.[["fullID"]])))
# some manual checks to make sure that points are broadly consistent
# mostly OK except for EIII - Point 81
#nvcIssues <- nvc %>% filter(duplicated(.[["fullID"]]))
#unique(nvcIssues$plotID)
View(nvc %>% filter(plotID %in% unique(nvcIssues$plotID)[11]) %>% arrange(Point))
# track down row numbers of missing dates for one case (turns out it's in 2010 when most of these happened)
#which(grepl("B1_82NA", nvc$fullID))
# losing this isn't the end of the world...but non-unique IDs is an issue
# 2 solutions
# (1) throw out duplicates (no issue in theory)
# (2) impute dates even though surveys didn't happen

# (2) imputation is harder, so implement this first until there's a problem with it
# missing values are only for points 79-83, so there will always be a preceding date
# probably due to variation in sampler walks (intermediate point selection)
# first create some duplicate columns for completeness
nvcImpute <- nvc %>% mutate(yearS = Year,
                            monthS = Month,
                            dayS = Day)
# backfill with preceding non-NAs
nvcImpute$yearS <- na.locf(nvcImpute$Year, na.rm = FALSE)
nvcImpute$monthS <- na.locf(nvcImpute$Month, na.rm = FALSE)
nvcImpute$dayS <- na.locf(nvcImpute$Day, na.rm = FALSE)
# re-compute IDs using the cleaner Y/M/D
nvcFilled <- nvcImpute %>% mutate(cleanDate = make_datetime(yearS, monthS, dayS),
                                  fullID = paste0(pointID,yearS))
# all fixed
#n_distinct(nvcFilled$fullID)


###############
# now returning to the grimey NVCs
# 480 rows, of which 340 are for rows without N/A
#sum(!is.na(nvcRaw$Other))
# 2003 - 2004 these were coded as 0/1 for just some other NVC (also once in 2013)
# frankly, an unlabelled "other" NVC has no real value here
#unique(nvcFilled$Other)

# this is a mess...
#unique(nvcFilled$`additional NVC`)
# pre-cleaning "M25a+b" and "M23a/b" manually to simplify cleaning
nvcFilled$`additional NVC`[nvcFilled$`additional NVC` == "M23a/b"] <- "M23a + M23b"
nvcFilled$`additional NVC`[nvcFilled$`additional NVC` == "M25a+b"] <- "M25a + M25b"

# use regex to extract meaningful NVC categories
# first Other
# [[:upper:]]{1} --> exactly 1 uppercase character, followed by
# [:digit:]{1,2} --> 1-2 digits, followed by
# [:alpha:]{0,1} --> 0-1 letters
# [-1] to ignore the first element, which is NA
#nvcOther <- unique(unlist(str_extract_all(unique(nvcFilled$Other),
#                                          "[[:upper:]]{1}[:digit:]{1,2}[:alpha:]{0,1}")))[-1]
# next, additional NVC
# same as above except all subtypes are helpfully lowercase, so the last expression is changed to
# [:lower:]{0,1} --> 0-1 lowercase letters
#nvcExtra <- unique(unlist(str_extract_all(unique(nvcFilled$`additional NVC`),
#                                          "[[:upper:]]{1}[:digit:]{1,2}[:lower:]{0,1}")))[-1]
# get an NVC master list with subtypes
#nvcSupp <- append(nvcExtra, nvcOther)
# put in alphabetical order for easier comprehension
# convert all uppercase bc of some subtypes in uppercase
# there are 45 total NVCs recorded including subtypes
#nvcAll <- sort(unique(toupper(append(nvcSupp, unique(nvcFilled$nvcClean)))))
# removing the subtypes, there are 25 recorded categories
#nvcMain <- sort(unique(unlist(str_extract_all(nvcAll,"[[:upper:]]{1}[:digit:]{1,2}"))))
# compared against 12 main columns in the base dataset
#sort(unique(unlist(str_extract_all(nvcFilled$nvcClean,"[[:upper:]]{1}[:digit:]{1,2}"))))

# previous iterations considered throwing out Comments
# this is true, except for NVC due to methodology
# official methods state "Communities often grade into each other. Note the one which is most likely, and note any other community in the comments column"
# however, other info like species IDs and etc is still omitted
# high work; low return -- only 8.7% of observations have a comment
# these contain information on some more NVCs, but these are more difficult to deal with (lots of "nearby" to mess up spatial analysis)'
# decent species encoding but not consistent spatiotemporally
#n_distinct(nvcRaw$Comments)
#sum(!is.na(nvcRaw$Comments))
# some scratch work
#unique(nvcRaw %>% filter(!is.na(Comments)) %>% select(Year))
# 20 items
# these include some "maybe" or "?" uncertainty
# also comment on stuff "nearby" or as well as effects on the "edge" and "side"
# one row which said "if further across" begs some questions, but left in given spatial uncertainty
#nvcComments <- unique(unlist(str_extract_all(unique(nvcExtract$Comments), "[HMUW]{1}[:digit:]{1,2}[:alpha:]{0,1}")))#[-1]
# 74 rows of interest which is minor, but doesn't hurt to just get it?
# for future years, this will fail in cases like "M25a or b" which showed up once
# luckily, both of these were already documented in normal columns
# the old regex was a bit weird when "R2D" happened to show up randomly, so changed the regex slightly
#nvcCom <- nvcExtract %>% mutate(extraNVC = str_detect(Comments, "[HMUW]{1}[:digit:]{1,2}[:alpha:]{0,1}"))
#View(nvcCom %>% filter(extraNVC == 1))

###############
# now the major clean-up
# using methods above (same regex) without screening for duplicates via "unique()"
# cram all NVC IDs into a single string based on order of priority (main > Other > additional)
nvcCompiled <- nvcFilled %>% select(-all_of(majorNVC)) %>% rowwise() %>% 
  mutate(nvcFull = paste(nvcClean, # take the melted nvcClean colun as-is
                         paste(unlist(str_extract_all(
                           Other, # process Other column w/ regex
                           "[[:upper:]]{1}[:digit:]{1,2}[:alpha:]{0,1}")),
                           collapse = ";"), # use ; as a separator
                         paste(unlist(str_extract_all(
                           `additional NVC`, # repeat with additional NVC
                           "[[:upper:]]{1}[:digit:]{1,2}[:lower:]{0,1}")),
                           collapse = ";"),
                         paste(unlist(str_extract_all(
                           Comments, # repeat with Comments
                           # change regex slightly bc this column is more messy
                           "[HMUW]{1}[:digit:]{1,2}[:alpha:]{0,1}")),
                           collapse = ";"),
                         sep = ";"), # always use ; separators for consistency
         # repeat the process for NVCs without subtypes
         # easier to just handle in parallel than deal with duplicates after-the-fact
         nvcSimp = paste(unlist(str_extract_all(
                           nvcClean,
                           "[[:upper:]]{1}[:digit:]{1,2}")),
                         paste(unlist(str_extract_all(
                           Other,
                           "[[:upper:]]{1}[:digit:]{1,2}")),
                           collapse = ";"),
                         paste(unlist(str_extract_all(
                           `additional NVC`,
                           "[[:upper:]]{1}[:digit:]{1,2}")),
                           collapse = ";"),
                         paste(unlist(str_extract_all(
                           Comments,
                           "[HMUW]{1}[:digit:]{1,2}")),
                           collapse = ";"),
                         sep = ";"))
# then, remove NAs, duplicates, and extra ";" to get a "clean" list of NVCs in priority order separated by ";"
# finally, extract each NVC from the collapsed string and place into respective column if exists
# remove unneeded columns
# maybe leave Other and additional NVC for now bc they contain some notes too?
nvcExtract <- nvcCompiled %>% mutate(nvcStr = gsub('^\\;|\\;$', '', # remove leading and trailing ;
                                                   gsub(";;", ";", # replace all ";;" with single ; to make process easier
                                                        str_remove_all( # remove all NAs within the string
                                                          # un-compress to remove duplicates (force to upper case bc case mismatches)
                                                          paste(unique(unlist(str_split(toupper(nvcFull), ";"))), 
                                                                collapse = ";"), # re-compress using ; separators
                                                          "NA"))),
                                     # unlist and portion out
                                     nvc1 = unlist(str_split(nvcStr, ";"))[1], 
                                     nvc2 = unlist(str_split(nvcStr, ";"))[2],
                                     nvc3 = unlist(str_split(nvcStr, ";"))[3],
                                     nvc4 = unlist(str_split(nvcStr, ";"))[4],
                                     # repeat with 
                                     nvcStr2 = gsub('^\\;|\\;$', '', gsub(";;", ";", str_remove_all(
                                       paste(unique(unlist(str_split(nvcSimp, ";"))), collapse = ";"),
                                       "NA"))),
                                     nvcM1 = unlist(str_split(nvcStr2, ";"))[1],
                                     nvcM2 = unlist(str_split(nvcStr2, ";"))[2],
                                     nvcM3 = unlist(str_split(nvcStr2, ";"))[3],
                                     nvcM4 = unlist(str_split(nvcStr2, ";"))[4]) %>% 
  select(-c("Other", "additional NVC", "nvcClean","nvcFull", "nvcStr", "nvcSimp", "nvcStr2"))
# replace leftover empty strings with NA
nvcExtract$nvc1[nvcExtract$nvc1 == ""] <- NA
nvcExtract$nvcM1[nvcExtract$nvcM1 == ""] <- NA

# adding Comments gives 49 more NVC2 and 10 more NVC3
# 40894 total rows, of which
#nrow(nvcExtract)
# all but 220 have one NVC class (99.4%), so we're recovered ~346 rows of data tossed by the initial script
#sum(!is.na(nvcExtract$nvcM1))
# 7832 have a second major NVC class (19.2% of those with one NVC)
#sum(!is.na(nvcExtract$nvcM2))
# 427 have a third major NVC class (5.5% of those with a second NVC)
#sum(!is.na(nvcExtract$nvcM3))
# 13 have 4 classifications (3.0% of those with a third NVC)
#sum(!is.na(nvcExtract$nvcM4))

# some more checks
# few of these seem systemic except EIII
# something happened to improve data collection in 2015: 2008-2014 have some notable absences
# almost all points 77 (and below) have height and density measured, so something's just up with NVC
# EIII is only the last points (79-81) with no data of note
nvcMissing <- nvcExtract %>% filter(is.na(nvcM1))
View(nvcMising)
table(nvcMissing$plotID, nvcMissing$yearS)

###############
# what actually are these multi-matches?
#nvcMultiple <- nvcExtract %>% filter(!is.na(nvc2))
#nvcMultiple2 <- nvcMultiple %>% count(nvc1, nvc2, nvc3, nvc4, sort = TRUE)
#nvcMultiple3 <- nvcMultiple2 %>%
#  mutate(uniqueCombo = pmap_chr(list(nvc1, nvc2, nvc3, nvc4), function(x, y, z, a) toString(sort(c(x, y, z, a))))) %>%
#  group_by(uniqueCombo) %>% summarise(sum = sum(n)) %>% 
#  mutate(length = str_length(uniqueCombo),
#         percent = sum/sum(!is.na(nvcExtract$nvc2))) %>%
#  arrange(-sum)
# export for review
#write_csv(nvcMultiple3, file = "./Vegetation/MultiCodedNVCs.csv")

# repeat for no-subtype NVCs
nvcMultipleMain <- nvcExtract %>% filter(!is.na(nvcM2))
nvcMultipleMain2 <- nvcMultipleMain %>% count(nvcM1, nvcM2, nvcM3, nvcM4, sort = TRUE)
nvcMultipleMain3 <- nvcMultipleMain2 %>%
  mutate(uniqueCombo = pmap_chr(list(nvcM1, nvcM2, nvcM3, nvcM4), function(x, y, z, a) toString(sort(c(x, y, z, a))))) %>%
  group_by(uniqueCombo) %>% summarise(occurrencesInData = sum(n)) %>% 
  mutate(distinctNVC = ceiling(str_length(uniqueCombo)/4.9),
         percentOfObs = occurrencesInData/sum(!is.na(nvcExtract$nvc2))) %>%
  arrange(-occurrencesInData)
# export for review
#write_csv(nvcMultipleMain3, file = "./Vegetation/MultiCodedNVC_main.csv")



###############
# other variables

# could clean up date rows (two columns w/ substr before and after the &?)
# first, let's check how much info is missing....1527 rows (4% of the whole dataset)
# actually 1581 if we count the easily-imputed stuff
#sum(!is.na(nvcRaw$Month)) - sum(!is.na(nvcRaw$Date))

# first, density (units?) and height (cm) 
# Right densities are always larger (by no more than 1 unit) relative to front and left (except AII, CIII, EIV, FII)
# IQRs slightly more divergent due to narrower data spread
#nvcExtract %>% select(matches("Density|Block|Plot")) %>% split(.$Plot) %>% map(summary)
# IQR are almost exclusively identical (except for Max), especially with higher-level groupings (plot and block only)
# Difference in means always < 2 cm; likely driven by large outliers
# Front heights are larger overall (not uniformly across all plots), especially A&B blocks
# Right heights slightly larger in most C&D, E&F plots but only by a bit
#nvcExtract %>% select(matches("Height|Block|Plot")) %>% split(.$plotID) %>% map(summary)
# without looking at the YoY breakdown, this looks suitable to average (everything can be chalked up to observer bias or stochasticity)
# so let's just create the mean/var of these
nvcExtract <- nvcExtract %>% rowwise() %>% mutate(denAvg = mean(c(`FRONT Density`, `LEFT Density`, `RIGHT Density`), na.rm=TRUE),
                                                  denVar =  var(c(`FRONT Density`, `LEFT Density`, `RIGHT Density`), na.rm=TRUE),
                                                  htAvg = mean(c(`FRONT Height`, `LEFT Height`, `RIGHT Height`), na.rm=TRUE),
                                                  htVar = var(c(`FRONT Height`, `LEFT Height`, `RIGHT Height`), na.rm=TRUE))

# slightly more compact and safer, but less computationally efficient version
# comes out to the same values but is left for completeness
#nvcExtract <- nvcExtract %>% rowwise() %>% mutate(denAvg = mean(c_across(ends_with("Density")), na.rm=TRUE),
#                                                  denVar =  var(c_across(ends_with("Density")), na.rm=TRUE),
#                                                  htAvg = mean(c_across(ends_with("Height")), na.rm=TRUE),
#                                                  htVar = var(c_across(ends_with("Height")), na.rm=TRUE))

# D/De - deer
# S/Sh - sheep
# carnivore - carnivore
# rd - deer
# C - cow
# V - vole
# 5/55 - sheep
# 25 - sheep
# 1, y - sheep bc “if the dung is not sheep, note the species"
# this coding is supported as all "1" and "Y" are in early years or in recorder DJR's first survey season of 2014
# 0, 1, n, y, are not helpful and set to 0
#unique(nvcRaw$Dung)
# back to the regex
nvcExtract <- nvcExtract %>% mutate(dungC = as.numeric(str_detect(tolower(Dung), "[c](?!a)")),
                                    dungD = as.numeric(str_detect(tolower(Dung), "[d](?!a)")),
                                    dungV = as.numeric(str_detect(tolower(Dung), "[v](?!o)")),
                                    dungS = as.numeric(str_detect(tolower(str_replace_all(Dung, 
                                                                                          pattern=" ",
                                                                                          repl="")), "sh(?![as])|[s](?![:alpha])|[15y]")),
                                    dungCar = as.numeric(str_detect(tolower(Dung), "carnivore")))

# not very common to see dung either way, but I guess this makes sense
#summary(nvcExtract %>% select(contains("Dung") | contains("dung")) %>% mutate_all(as.factor))

# y/Y, S = yes
# n/a, n, N, x = no
nvcExtract <- nvcExtract %>% mutate(whiteStick = factor(`White stick found?`,
                                                        levels = sort(unique(nvcRaw$`White stick found?`)),
                                                        labels = c(0, 0, 0, 1, 0, 1, 1)))

# replace all NaN with NA
# messes up cleanDate if run on the whole df
#summary(sapply(nvcExtract, FUN = is.nan))
# run only on the offending rows
nvcExtract <- nvcExtract %>% mutate(across(c(denAvg, htAvg),
                                           ~ifelse(is.nan(.), NA, .)))

# add columns for sheep/cow offtake
nvcExtract <- nvcExtract %>% mutate(sheepOffRaw = strtoi(factor(Plot,
                                                                levels = unique(nvcRaw$Plot),
                                                                labels = c(9, 3, 2, 0))),
                                    # correct for increased grazing in 2021
                                    sheepOff = ifelse(yearS == 2021, (sheepOffRaw*11 + 12)/12, sheepOffRaw),
                                    cowOff = ifelse(Plot == "III", 1, 0))
nvcExtract <- nvcExtract %>% select(-"sheepOffRaw")

# export for plotting later
write_rds(nvcExtract, "./Vegetation/nvcPlots.rds")

# remove unneeded columns, most of these are self-explanatory
# NVC subtypes have too much error going on, so they're left out
nvcFinal <- nvcExtract %>% 
  mutate(nvcSumComs = sum(!is.na(nvcM1), 
                          !is.na(nvcM2),
                          !is.na(nvcM3),
                          !is.na(nvcM4)),
         nvcVars = sum(!is.na(nvc1), 
                        !is.na(nvc2),
                        !is.na(nvc3),
                        !is.na(nvc4)))
nvcFinal <- nvcFinal %>% select(-c("White stick found?", "Period", "Dung", "nvc1", "nvc2", "nvc3", "nvc4", "FRONT Density", 
                                    "LEFT Density", "RIGHT Density", "FRONT Height", "LEFT Height", "RIGHT Height", "Year", "Month", 
                                    "Day"))
nvcFinalCorr <- nvcFinal

# a way to get a unique column for each NVC sub-community presence/absence
uniqueNVC <- sort(unique(append(unique(nvcFinalCorr$nvcM1),
                                c(unique(nvcFinalCorr$nvcM2),
                                  unique(nvcFinalCorr$nvcM3),
                                  unique(nvcFinalCorr$nvcM4)))))

# yes loops are bad...but this feels too trivial to define a whole function for
for (i in uniqueNVC) {
  isPresent <- colSums(rbind(str_detect(nvcFinalCorr$nvcM1, i),
                             str_detect(nvcFinalCorr$nvcM2, i),
                             str_detect(nvcFinalCorr$nvcM3, i),
                             str_detect(nvcFinalCorr$nvcM4, i)), 
                       na.rm=TRUE)
  # tests for NA to avoid casting NA to 0 and inflating data availability
  nvcFinalCorr[, i] <-  ifelse(!is.na(nvcFinalCorr$nvcM1),
                           ceiling(isPresent/4), # hack to turn any non-zero into 1
                           NA)
}

# types of NVC co-occurence
ggcorrplot(cor(nvcFinalCorr[,30:54], use = "complete.obs"))

# further correlation plot with all numeric variables
nvcCorr <- nvcFinalCorr %>% filter(!is.na(nvcM1)) %>% 
  mutate(#plot = as.numeric(as.roman(Plot)),
         block = as.numeric(chartr("ABCDEF","123456",Block))) %>% 
  select(-c(colnames(nvcFinalCorr)[c(4:10,12:17)])) %>% 
  relocate(c("block", "sheepOff", "cowOff"), .before = "Point")
summary(nvcCorr)
# probably some issues with casting all numeric NAs to 0
nvcCorr <- nvcCorr %>% mutate(across(where(is.numeric), 
                                     ~ifelse(is.na(.), 0, .)))
#block and plot are kind of meaningless here but exist to capture other correlation for inspection later
ggcorrplot(cor(nvcCorr %>% select(where(is.numeric)) %>% select(-contains("dung")), use = "complete.obs"))
# for this plot, it'd be best to cast NAs in dung columns to 0
nvcCorrD <- nvcCorr %>% mutate(across(contains("dung"), 
                              ~ifelse(is.na(.), 0, .)))
ggcorrplot(cor(nvcCorrD %>% select(where(is.numeric)), use = "complete.obs"))

# remove the rest of the unneeded columns
# density is just a worse version of height and not in species data; throw it out
# sheep and cow dung are kinda meaningless here bc we already have treatment offtake info
# leave comments for now because merging with species-specific data may be interesting as a sanity check
nvcCleaned <- nvcFinal %>% select(-c("dungS", "dungC", "denAvg", "denVar"))
write_rds(nvcCleaned, "./Vegetation/nvcCleaned.rds")

###############
# absolute simplest plot
nvcSimple <- nvcExtract %>% mutate(nvcCom1 = str_sub(nvc1, 1, 1),
                                   nvcCom2 = str_sub(nvc2, 1, 1),
                                   nvcCom3 = str_sub(nvc3, 1, 1),
                                   nvcCom4 = str_sub(nvc4, 1, 1),
                                   isW = grepl("W", paste(nvcCom1, nvcCom2, nvcCom3, nvcCom4, sep = "")),
                                   isM = grepl("M", paste(nvcCom1, nvcCom2, nvcCom3, nvcCom4, sep = "")),
                                   isU = grepl("U", paste(nvcCom1, nvcCom2, nvcCom3, nvcCom4, sep = "")),
                                   isH = grepl("H", paste(nvcCom1, nvcCom2, nvcCom3, nvcCom4, sep = "")))


simplePlot <- ggplot(nvcSimple, aes(x = yearS, fill = nvcCom1)) + 
  geom_bar() + 
  facet_wrap(~plotID) + 
  scale_fill_viridis_d(labels=c("Heath", "Mire", "Grassland", "Woodland")) + 
  theme_bw() + 
  labs(x = "Year",
       y = "# Points",
       fill = "NVC")

ggsave(filename = "./Vegetation/NVCPlots/nvcSimplePlot.jpg", simplePlot, 
       width = 10, height = 6, dpi = 300)


# plot starting with simpler categories
# first look at difference between treatments
data(Dark24)
nvcPrimaryMainPlot <- ggplot(nvcExtract, aes(x = yearS, fill = nvcM1)) + 
  geom_bar() + 
  facet_wrap(~Plot) + 
  scale_fill_manual(values = as.vector(Dark24)) + 
  theme_bw() + 
  labs(x = "Year",
       y = "# Points",
       fill = "NVC")

ggsave(filename = "./Vegetation/NVCPlots/nvcPrimaryMainPlot.jpg", nvcPrimaryMainPlot, dpi = 300)

# split out by block to see inherent vegetation characteristics based on location
nvcPrimaryMainBlock <-ggplot(nvcExtract, aes(x = yearS, fill = nvcM1)) + 
  geom_bar() + 
  facet_wrap(~Block) + 
  scale_fill_manual(values = as.vector(Dark24)) + 
  theme_bw() + 
  labs(x = "Year",
       y = "# Points",
       fill = "NVC")

ggsave(filename = "./Vegetation/NVCPlots/nvcPrimaryMainBlock.jpg", nvcPrimaryMainBlock, dpi = 300)


# break it down into the detailed sub-categories to see
# the first 2-3 years clearly use different (broader) NVC than later
# might be worth building a custom colour palette where the sub-types are hues of 20 main colours
data(palette36)
nvcPrimarySubPlot <- ggplot(nvcExtract, aes(x = yearS, fill = nvc1)) + 
  geom_bar() + 
  facet_wrap(~Plot) + 
  scale_fill_manual(values = as.vector(palette36)) + 
  theme_bw() + 
  labs(x = "Year",
       y = "# Points",
       fill = "NVC")

ggsave(filename = "./Vegetation/NVCPlots/nvcPrimarySubPlot.jpg", nvcPrimarySubPlot, dpi = 300)

nvcPrimarySubBlock <- ggplot(nvcExtract, aes(x = yearS, fill = nvc1)) + 
  geom_bar() + 
  facet_wrap(~Block) + 
  scale_fill_manual(values = as.vector(palette36)) + 
  theme_bw() + 
  labs(x = "Year",
       y = "# Points",
       fill = "NVC")

ggsave(filename = "./Vegetation/NVCPlots/nvcPrimarySubBlock.jpg", nvcPrimarySubBlock, dpi = 300)


# repeat plots above but accounting for columns with multiple NVCs
nvcExtra2 <- nvcExtract %>% 
  pivot_longer(
    cols = nvcM1:nvcM4, 
    names_to = "nvcType",
    values_to = "nvc") %>% 
  filter(!is.na(nvc))

data(alphabet)
nvcFullMainPlot <- ggplot(nvcExtra2, aes(x = yearS, fill = nvc)) + 
  geom_bar() + 
  facet_wrap(~Plot) + 
  scale_fill_manual(values = as.vector(alphabet)) + 
  theme_bw() + 
  labs(x = "Year",
       y = "# Points",
       fill = "NVC")

ggsave(filename = "./Vegetation/NVCPlots/nvcFullMainPlot.jpg", nvcFullMainPlot, dpi = 300)


nvcFullMainBlock <- ggplot(nvcExtra2, aes(x = yearS, fill = nvc)) + 
  geom_bar() + 
  facet_wrap(~Block) + 
  scale_fill_manual(values = as.vector(alphabet)) + 
  theme_bw() + 
  labs(x = "Year",
       y = "# Points",
       fill = "NVC")

ggsave(filename = "./Vegetation/NVCPlots/nvcFullMainBlock.jpg", nvcFullMainBlock, dpi = 300)

# re-create dataset but with all subtypes
nvcExtra <- nvcExtract %>% 
  pivot_longer(
    cols = nvc1:nvc4, 
    names_to = "nvcType",
    values_to = "nvc") %>% 
  filter(!is.na(nvc)) 

nvcCols <- as.vector(createPalette(n_distinct(nvcExtra$nvc), c("#2E8445","#FFFFFF", "#E4371C")))

nvcFullSubPlot <- ggplot(nvcExtra, aes(x = yearS, fill = nvc)) + 
  geom_bar() + 
  facet_wrap(~Plot) + 
  scale_fill_manual(values = nvcCols) + 
  theme_bw() + 
  labs(x = "Year",
       y = "# Points",
       fill = "NVC")

ggsave(filename = "./Vegetation/NVCPlots/nvcFullSubPlot.jpg", nvcFullSubPlot, dpi = 300)

nvcFullSubBlock <- ggplot(nvcExtra, aes(x = yearS, fill = nvc)) + 
  geom_bar() + 
  facet_wrap(~Block) + 
  scale_fill_manual(values = nvcCols) + 
  theme_bw() + 
  labs(x = "Year",
       y = "# Points",
       fill = "NVC")

ggsave(filename = "./Vegetation/NVCPlots/nvcFullSubBlock.jpg", nvcFullSubBlock, dpi = 300)

# make some association plots
# first with only the main NVC
nvcAssocMain <- nvcExtract %>% group_by(yearS, nvcM1) %>% summarise(n = n())
  
ggplot(nvcAssocMain, aes(nvcM1, yearS)) +
  geom_tile(aes(fill = n)) + 
  scale_fill_viridis_c(trans = scales::pseudo_log_trans(sigma = 5)) + 
  theme_bw()

# good ideas on weighting that seem to have no impact and are removed for now
# then for all the NVCs
#nvcExtra2 <- nvcExtra2 %>% mutate(type = 1/as.numeric(factor(str_sub(nvcType,-1),
#                                                    levels = c(1,2,3,4),
#                                                    labels = c(4,3,2,1))))
nvcAssocFull <- nvcExtra2 %>% group_by(yearS, nvc) %>% 
  summarise(n = n())#, 
           # weight = sum(type, na.rm = TRUE)) %>% 
  # in n, rows w/ multiple NVCs have their importance inflated
  # this divides additional NVCs by the scale of their contribution
  #mutate(nW = n - (weight - n)/12)

# filter out the NVCs with only one observation to improve readability
ggplot(nvcAssocFull %>% filter(!(nvc %in% c("H17", "M11", "M16", "M22", "M26", "W6"))),
       aes(nvc, yearS)) +
  geom_tile(aes(fill = n)) + 
  scale_fill_viridis_c(trans = scales::pseudo_log_trans(sigma = 5)) + 
  theme_bw()

# plot a heatplot for density average
# quite cluttered, so worth zooming in to see more
#ggplot(nvcExtract %>% 
#         filter(nvcM1 != "W6" & nvcM1 != "U19" & nvcM1 != "M11" & nvcM1 != "M18" & nvcM1 != "U13"), 
#       aes(yearS, plotID, fill = denAvg))+
#  geom_tile(color = "white", linewidth = 0.1) + 
  #set the fill colour to viridis and label the legend
#  scale_fill_viridis_c(name = "Density (cm)") + 
#  facet_grid(~nvcM1) 


# plot a heatplot for height average
# quite cluttered, so worth zooming in to see more
ggplot(nvcExtract %>% 
         filter(nvcM1 != "W6" & nvcM1 != "U19" & nvcM1 != "M11" & nvcM1 != "M18" & nvcM1 != "U13"), 
       aes(yearS, plotID, fill = htAvg))+
  geom_tile(color = "white", linewidth = 0.1) + 
  #set the fill colour to viridis and label the legend
  scale_fill_viridis_c(name = "Height (cm)") + 
  facet_grid(~nvcM1) 

# some more plotting of heights to get a sense of what's going on
nvcFinal %>% 
  # woodlands are notorious outliers and will skew results
  # comment out the line below to see full data: rescaling the legend/colours to log may help
  filter(W4 != 1 & W6 != 1) %>% 
  drop_na(htAvg) %>% 
  ggplot(aes(x = yearS, y = Plot, fill = htAvg)) + 
  #scale_fill_viridis_c(trans = scales::pseudo_log_trans(sigma = 5)) + 
  geom_tile()

nvcFinal %>% ggplot(aes(x = sheepOff, y = htAvg, fill = Plot)) + 
  #scale_fill_viridis_c(trans = scales::pseudo_log_trans(sigma = 5)) + 
  geom_point()


