library(tidyverse)
library(inlabru)

vegINLA <- read_rds("./Vegetation/vegGeo.rds") %>% mutate(fBlock = as.factor(Block),
                                                          Plot_code = paste0(Plot, Block),
                                                          whiteStick = as.logical(whiteStick),
                                                          dungD = as.logical(dungD),
                                                          treat = factor(Plot,
                                                                         levels = c("III", "I", "II", "IV"))) %>% 
  filter(Year > 2002)

squeeze_TurnBeta <- (1 - sort(unique(vegINLA$nvcTurnSub))[n_distinct(vegINLA$nvcTurnSub) - 2])/10
vegINLA <- vegINLA %>% mutate(nvcTurnSubSq = ifelse(nvcTurnSub == 1, 1 - squeeze_TurnBeta, 
                                                    ifelse(nvcTurnSub == 0, squeeze_TurnBeta,
                                                           nvcTurnSub)))

fitSpatTempAR1_TurnBeta <- read_rds("./Models/Objects/fitSpatTempAR1_TurnComBeta.rds")

inv.logit <- function(x){
  return(exp(x)/(1+exp(x)))
}

logit <- function(x){
  return(x/(1-x))
}


varName <- "precip"
minVal <- floor(min(st_drop_geometry(vegINLA[varName]), na.rm = TRUE))
maxVal <- ceiling(max(st_drop_geometry(vegINLA[varName]), na.rm = TRUE))

var <- seq(minVal, 
           maxVal, 
           by = (maxVal - minVal)/100)
vars <- rownames(fitSpatTempAR1_TurnBeta$summary.fixed)

colMeans <- st_drop_geometry(vegINLA) %>% ungroup() %>% 
  select(any_of(vars)) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))

predDf <- data.frame(precip, 
                     colMeans %>% select(-precip), 
                     Intercept = 1)
predDf <- predDf %>% 
  relocate(sort(names(predDf)))

ppxl_TurnBeta <- fm_centroids(mesh_TurnBeta)#, mask = plotBounds)

avgPreds_TurnBeta <- data.frame(mean = numeric(),
                                sd = numeric(),
                                q0.025 = numeric(),
                                q0.975 = numeric(),
                                iter = numeric(),
                                `.triangle` = numeric())
for(i in 1:nrow(predDf)){
  ppxl_all_TurnBeta <- fm_cprod(ppxl_TurnBeta, data.frame(Year = 2023,
                                                          sheepOff = predDf[i,]$sheepOff,
                                                          cowOff = predDf[i,]$cowOff,
                                                          windSpeed = predDf[i,]$windSpeed,
                                                          meanT = predDf[i,]$meanT,
                                                          solarRad = predDf[i,]$solarRad,
                                                          precip = predDf[i,]$precip,
                                                          litter = predDf[i,]$litter,
                                                          elev = predDf[i,]$elev,
                                                          isW = predDf[i,]$isW,
                                                          isH = predDf[i,]$isH,
                                                          isU = predDf[i,]$isU,
                                                          isM = predDf[i,]$isM,
                                                          whiteStick = predDf[i,]$whiteStick,
                                                          sameObserver = predDf[i,]$sameObserver))
  fullPred_TurnBeta <- predict(
    fitSpatTempAR1_TurnBeta, ppxl_all_TurnBeta,
    # it looks like we need to provide actual values for covariates (and year) for predictions to work
    ~ field + Intercept + sheepOff + cowOff + windSpeed + meanT + solarRad + precip + litter + isW + isH + isU + isM  + whiteStick + sameObserver
  ) %>% select(c("mean", "sd", "q0.025", "q0.975", ".triangle")) %>% 
    mutate(iter = i)

  avgPreds_TurnBeta <- bind_rows(avgPreds_TurnBeta, st_drop_geometry(fullPred_TurnBeta))
}





# test <- predDf %>% rowwise() %>% 
#   mutate(linPred = sum(as.numeric(c_across(everything()))*as.vector(means)$mean),
#          respPred = inv.logit(linPred),
#          linUpper = sum(as.numeric(c_across(cowOff:windSpeed))*as.vector(upper)$`0.975quant`),
#          respUpper = exp(linUpper)/(1 + exp(linUpper)),
#          linLower = sum(as.numeric(c_across(cowOff:windSpeed))*as.vector(lower)$`0.025quant`),
#          respLower = exp(linLower)/(1 + exp(linLower)))

# ggplot(test, aes(x = precip, y = respLower)) + 
#   geom_line() + 2
#   theme_bw()

test <- predDf
test$mean <- inv.logit(avgPreds_TurnBeta$mean)
test$lower <- inv.logit(avgPreds_TurnBeta$lower)
test$upper <- inv.logit(avgPreds_TurnBeta$upper)

test2 <- avgPreds_TurnBeta %>% filter(`.triangle` == 2)
test2$precip <- predDf$precip
test2$mean2 <- inv.logit(test2$mean)
test2$lower <- inv.logit(test2$q0.025)
test2$upper <- inv.logit(test2$q0.975)


ggplot(data = test2) + 
    # geom_ribbon(aes(x=precip, y=mean, ymin=lower, ymax=upper),
    #             fill="cyan", alpha = 0.25) +
  geom_line(aes(x=precip, y=mean2), linewidth=1, color="blue") + 
  theme_bw()
