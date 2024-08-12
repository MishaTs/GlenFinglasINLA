# load packages
library(inlabru)
library(fmesher)
library(INLA)
library(sf)
#library(terra)
library(tidyverse)

# read in data
vegINLA <- read_rds("./Vegetation/vegGeo.rds") %>% mutate(fBlock = as.factor(Block),
                                                          Plot_code = paste0(Plot, Block),
                                                          # as.numeric converts [0, 1] to [1, 2] for some reason
                                                          whiteStick = as.numeric(whiteStick) - 1,
                                                          dungD = as.numeric(dungD),
                                                          treat = factor(Plot,
                                                                         levels = c("III", "I", "II", "IV")),
                                                          sameObserver = as.numeric(Recorder == lag(Recorder))) %>% 
  filter(Year > 2002)
# prepare data to match (just in case?)
squeeze_TurnBeta <- (1 - sort(unique(vegINLA$nvcTurnSub))[n_distinct(vegINLA$nvcTurnSub) - 2])/10
vegINLA <- vegINLA %>% mutate(nvcTurnSubSq = ifelse(nvcTurnSub == 1, 1 - squeeze_TurnBeta, 
                                                    ifelse(nvcTurnSub == 0, squeeze_TurnBeta,
                                                           nvcTurnSub)))

fitSpatTempAR1_TurnBeta <- read_rds("./Models/Objects/fitSpatTempAR1_TurnBeta.rds")

# helper function to avoid typing it out
inv.logit <- function(x){
  return(exp(x)/(1+exp(x)))
}

# first, the mesh
mesh_TurnBeta <- fm_mesh_2d_inla(
  loc = vegINLA$geometry, 
  # doesn't look like any changes where max.edge > 0.05 have much impact
  # much finer mesh possible otherwise, but left as-is for computing benefits for now
  # outer edge set to 1 but doesn't matter
  max.edge = c(0.05,1), # map units inside and outside
  # cutoff is min edge
  #cutoff = 50, offset = c(100, 300),
  crs = fm_crs(vegINLA)
)

# generate spatial predictions to get an "average" intercept incorporating the spatial field
# these predictions give a slightly different realisation with each run, but it's within acceptable margins
ppxl_TurnBeta <- fm_pixels(mesh_TurnBeta) %>% mutate(Year = 2023)
pred_TurnBeta <- predict(
  fitSpatTempAR1_TurnBeta, ppxl_TurnBeta,
  ~ field + Intercept
)

# all variables to include in the model
vars_TurnBeta <- rownames(fitSpatTempAR1_TurnBeta$summary.fixed)
varsPlot_TurnBeta <- vars_TurnBeta[which(vars_TurnBeta != "Intercept")]

margEffects_TurnBeta <- data.frame()


for(varName in varsPlot_TurnBeta){
  # build the needed file structures
  # variable for partial effects plotting hypothetical sequence within realistic values
  #varName <- "meanT"
  minVal <- floor(min(st_drop_geometry(vegINLA[varName]), na.rm = TRUE))
  maxVal <- ceiling(max(st_drop_geometry(vegINLA[varName]), na.rm = TRUE))
  
  var <- seq(minVal, 
             maxVal, 
             by = (maxVal - minVal)/100)
  # mean posterior coefficient estimates
  means <- fitSpatTempAR1_TurnBeta$summary.fixed %>% select(mean) %>% 
    arrange(tolower(rownames(fitSpatTempAR1_TurnBeta$summary.fixed)))
  # lower 95% CI coefficient estimates
  lower <- fitSpatTempAR1_TurnBeta$summary.fixed %>% select(`0.025quant`) %>% 
    arrange(tolower(rownames(fitSpatTempAR1_TurnBeta$summary.fixed)))
  # upper 95% CI coefficient estimates
  upper <- fitSpatTempAR1_TurnBeta$summary.fixed %>% select(`0.975quant`) %>% 
    arrange(tolower(rownames(fitSpatTempAR1_TurnBeta$summary.fixed)))
  # remove our variable of interest here 
  varsFilt <- vars_TurnBeta[! vars_TurnBeta %in% varName]
  
  # get the average of all variables to feed in "all else held equal"
  # this includes binary variables which are treated as proportions here
  colMeans <- st_drop_geometry(vegINLA) %>% ungroup() %>% 
    select(any_of(vars_TurnBeta)) %>%
    summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))
  
  # generate values to "predict" across 
  predDf <- data.frame(var, 
                       colMeans %>% select(!all_of(varName)),
                       Intercept = 1) 
  # rename column to match for sorting
  colnames(predDf)[which(colnames(predDf) %in% "var")] <- varName
  # sort in alphabetical order for easier base R computation later
  predDf <- predDf %>% 
    relocate(sort(names(predDf)))
  
  # replace intercept with spatial predictions
  means["Intercept",] <- mean(pred_TurnBeta$mean)
  # assume all other coefficient estimates are held constant at their mean values
  lower[varsFilt,] <- means[varsFilt,]
  upper[varsFilt,] <- means[varsFilt,]
  
  # get weighted sum (e.g., evaluate the linear addative regression model) on the link scale
  # convert to response scale too
  predOutput <- predDf %>% rowwise() %>% 
    mutate(linPred = sum(as.numeric(c_across(cowOff:windSpeed))*as.vector(means)$mean),
           respPred = inv.logit(linPred),
           linUpper = sum(as.numeric(c_across(cowOff:windSpeed))*as.vector(upper)$`0.975quant`),
           respUpper = inv.logit(linUpper),
           linLower = sum(as.numeric(c_across(cowOff:windSpeed))*as.vector(lower)$`0.025quant`),
           respLower = inv.logit(linLower))
  
  # unify naming for consistency
  predExport <- predOutput %>% select(all_of(c(varName, "linPred","respPred","linUpper","respUpper","linLower","respLower"))) %>% 
    mutate(param = varName)
  colnames(predExport)[which(colnames(predExport) %in% varName)] <- "val"
  # update working object
  # rbind leads to huge headaches; use bind_rows tidy format instead
  margEffects_TurnBeta <- bind_rows(margEffects_TurnBeta, predExport)
  # remove unused objects
  rm(minVal, maxVal, var, lower, upper, varsFilt, colMeans, predDf, predOutput, predExport)
}

# for easier facet wrapping
margEffects_TurnBeta <- margEffects_TurnBeta %>% mutate(param = factor(param,
                                                       levels = varsPlot_TurnBeta))



# plot predictions for marginal effects plot
ggplot(data = margEffects_TurnBeta) + 
  geom_ribbon(aes(x=val, y=respPred, ymin=respLower, ymax=respUpper),
              fill="#5ec962", alpha = 0.25) +
  geom_line(aes(x=val, y=respPred), linewidth=1, color="#3b528b") +
  facet_wrap(~param, scales = "free") +
  xlab("") +
  ylab("Morisita-Horn overlap index") +
  theme_bw()
