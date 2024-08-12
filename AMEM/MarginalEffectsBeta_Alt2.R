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

ppxl_TurnBeta <- st_centroid(plotBounds) %>% select(c("Plot_ID", "plot", "block")) %>% mutate(Year = 2023)
pred_TurnBeta <- predict(
  fitSpatTempAR1_TurnBeta, ppxl_TurnBeta,
  ~ field + Intercept
)

repPtField <- st_drop_geometry(pred_TurnBeta) %>% 
  select(where(is.numeric)) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))

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
                     colMeans %>% select(-precip))
predDf <- predDf %>% 
  relocate(sort(names(predDf)))

avgPreds_TurnBeta <- data.frame(mean = numeric(),
                                sd = numeric(),
                                q0.025 = numeric(),
                                q0.975 = numeric(),
                                iter = numeric(),
                                Plot_ID = character())
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
    ~ sheepOff + cowOff + windSpeed + meanT + solarRad + precip + litter + isW + isH + isU + isM  + whiteStick + sameObserver
  ) %>% select(c("mean", "sd", "q0.025", "q0.975", "Plot_ID")) %>% 
    mutate(iter = i)

  avgPreds_TurnBeta <- bind_rows(avgPreds_TurnBeta, st_drop_geometry(fullPred_TurnBeta))
}


test3 <- avgPreds_TurnBeta %>% filter(Plot_ID == "BII") %>% 
  select(-Plot_ID) %>% #group_by(iter) %>%
  #summarise(across(everything(), \(x) mean(x, na.rm = TRUE))) %>% 
  mutate(fieldPredM = mean(repPtField$mean),
         fieldPredL = mean(repPtField$q0.025),
         fieldPredU = mean(repPtField$q0.975),
         linkPredM = fieldPredM + mean,
         linkPredL = fieldPredM + q0.025, 
         linkPredU = fieldPredM + q0.975,
         respPredM = inv.logit(linkPredM),
         respPredL = inv.logit(linkPredL),
         respPredU = inv.logit(linkPredU))

test3$precip <- predDf$precip

test4 <- avgPreds_TurnBeta %>% #filter(Plot_ID == "BII") %>% 
  select(-Plot_ID) %>% 
  group_by(iter) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE))) %>% 
  mutate(fieldPredM = mean(repPtField$mean),
         fieldPredL = mean(repPtField$q0.025),
         fieldPredU = mean(repPtField$q0.975),
         linkPredM = fieldPredM + mean,
         linkPredL = fieldPredM + q0.025, 
         linkPredU = fieldPredM + q0.975,
         respPredM = inv.logit(linkPredM),
         respPredL = inv.logit(linkPredL),
         respPredU = inv.logit(linkPredU))

test4$precip <- predDf$precip


















# test <- predDf
# test$mean <- inv.logit(avgPreds_TurnBeta$mean)
# test$lower <- inv.logit(avgPreds_TurnBeta$lower)
# test$upper <- inv.logit(avgPreds_TurnBeta$upper)
# 
# test2 <- avgPreds_TurnBeta %>% filter(`.triangle` == 2)
# test2$precip <- predDf$precip
# test2$mean2 <- inv.logit(test2$mean)
# test2$lower <- inv.logit(test2$q0.025)
# test2$upper <- inv.logit(test2$q0.975)


ggplot(data = test4) + 
  geom_ribbon(aes(x=precip, y=respPredM, ymin=respPredL, ymax=respPredU),
               fill="cyan", alpha = 0.25) +
  geom_line(aes(x=precip, y=respPredM), linewidth=1, color="blue") + 
  theme_bw()
