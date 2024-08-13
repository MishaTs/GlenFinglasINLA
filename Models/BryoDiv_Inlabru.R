library(tidyverse)
library(INLA)
library(fmesher)
library(inlabru)

ggplot(vegINLA, aes(x = meanT, y = growDays)) + geom_point() + stat_smooth(method = "lm", formula=y ~  poly(x, 2)) +
  scale_y_continuous(limits = c(650, 950)) +
  scale_x_continuous(limits = c(6.5, 9.5))


#### Data Setup ####
vegGeo <- read_rds("./Vegetation/vegGeo.rds")

vegINLA <- vegGeo %>% mutate(fBlock = as.factor(Block),
                             Plot_code = paste0(Plot, Block),
                             treat = factor(Plot,
                                            levels = c("III", "I", "II", "IV"))) #%>% filter(!is.na(nvcM1)) 

#vegINLA <- vegINLA %>% mutate(meanT2 = meanT)

# t(apply(table(vegINLA$isW, vegINLA$flower),1, function(x) x/sum(x)))

# first, the mesh
mesh_BryoDiv <- fm_mesh_2d_inla(
  loc = vegINLA$geometry, 
  # doesn't look like any changes where max.edge > 0.05 have much impact
  # much finer mesh possible otherwise, but left as-is for computing benefits for now
  # outer edge set to 1 but doesn't matter
  max.edge = c(0.05,1), # map units inside and outside
  # cutoff is min edge
  #cutoff = 50, offset = c(100, 300),
  crs = fm_crs(vegINLA)
)
ggplot() + gg(mesh_BryoDiv) + theme_bw()
#mesh_BryoDiv$n # 2044 points

# Making SPDE
# Quantify distance between points and define Matern correlation on the mesh
#inla.spde2.pcmatern allows for penalised complexity priors
aMatSPDE_BryoDiv <- inla.spde2.matern(mesh = mesh_BryoDiv)#, 
#prior.range = c(10, 0.5), 
#prior.sigma = c(.5, .5))#,
#alpha = 2) 
# code from slide 10 of Finn's ppt
# https://www.maths.ed.ac.uk/~flindgre/talks/Lindgren_RSS2023.pdf

# map across years for the temporal component
yearMap_BryoDiv <- bru_mapper(fm_mesh_1d(sort(unique(vegINLA$Year))), indexed = TRUE)

#### Modelling ####
# first define the model components
# high-leverage variables end up being quite impactful (e.g., mostly 0 with a few rows of 1)
inlaNullC_BryoDiv <- ~ #Intercept(1) + 
  # interaction terms don't change sheepOff but kill the signal in cowOff due to construction
  sheepOff + cowOff + 
  #treat +
  # humidity, precipitation, and pressure are too correlated to include in a single model
  # simply by effect size and significance, pressure is best; harder to interpret though
  # frost days have no effect
  windSpeed + meanT + solarRad + precip + 
  # dung, bareGround, and wood have almost no observations
  # water, rock, and mud also fairly low
  # flowers and tussock take away from the NVC community variables (less positive)
  litter + rock + #graze +
  # flower and tussock are mainly M and U
  elev +
  isW + isH + isU + isM

inlaTempC_BryoDiv <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + rock + #graze +
  elev +
  isW + isH + isU + isM +
  field1(Year, model='ar1', replicate = Point) +
  field2(fBlock, model = "iid") + 
  field3(Plot_code, model = "iid")

# tests for alternative AR structures
inlaSpatTempAR1C_BryoDiv <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + rock + #graze +
  elev +
  isW + isH + isU + isM +
  field(geometry, 
        model = aMatSPDE_BryoDiv, 
        group = Year, 
        group_mapper = yearMap_BryoDiv,
        control.group = list(model = "ar",
                             # order <= 10
                             order = 1))

inlaSpatTempAR2C_BryoDiv <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + rock + #graze +
  elev +
  isW + isH + isU + isM +
  field(geometry, 
        model = aMatSPDE_BryoDiv, 
        group = Year, 
        group_mapper = yearMap_BryoDiv,
        control.group = list(model = "ar",
                             # order <= 10
                             order = 2))

# as high as we for AR terms
inlaSpatTempAR3C_BryoDiv <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + rock + #graze +
  elev +
  isW + isH + isU + isM +
  field(geometry, 
        model = aMatSPDE_BryoDiv, 
        group = Year, 
        group_mapper = yearMap_BryoDiv,
        control.group = list(model = "ar",
                             # order <= 10
                             order = 3))

# define linear combinations
lc1 <- inla.make.lincomb(sheepOff = 1, cowOff = 1)
names(lc1) <- "offtake"
lc2 <- inla.make.lincomb(isW = 1, isH = 1, isU = 1, isM = 1)
names(lc2) <- "nvc"

# build the likelihoods
likeNull_BryoDiv <- like(
  formula = hillBryo1 ~ .,
  # possible distributions: gamma and tweedie
  # tweedie is a superset of gaussian, gamma, poisson, and inverse gaussians
  family = "tweedie", 
  # dtweedie: Reached upper limit. Increase TWEEDIE_MAX_IDX or scale your data!
  scale = 0.5,
  data=vegINLA,
  # domain optional for non-spatial model
  domain = list(geometry = mesh_BryoDiv)
)

likeTemp_BryoDiv <- like(
  formula = hillBryo1 ~ .,
  # tweedie is a superset of gaussian, gamma, poisson, and inverse gaussians
  family = "tweedie", 
  scale = 0.5,
  data=vegINLA,
  # domain optional for non-spatial model
  domain = list(geometry = mesh_BryoDiv)
)

likeSpatTempAR_BryoDiv <- like(
  formula = hillBryo1 ~ .,
  # tweedie is a superset of gaussian, gamma, poisson, and inverse gaussians
  family = "tweedie", 
  scale = 0.5,
  data=vegINLA,
  domain = list(geometry = mesh_BryoDiv)
)

# run the actual fit
# null model
fitNull_BryoDiv <- bru(components = inlaNullC_BryoDiv,
                       likeNull_BryoDiv,
                       options = list(
                         # add linear combinations
                         lincomb = c(lc1, lc2),
                         #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                         control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                         control.predictor = list(compute=TRUE),
                         #control.inla = list(int.strategy = "eb"),
                         verbose = FALSE))

# temporal model
fitTemp_BryoDiv <- bru(components = inlaTempC_BryoDiv,
                       likeTemp_BryoDiv,
                       options = list(
                         # add linear combinations
                         lincomb = c(lc1, lc2),
                         #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                         control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                         control.predictor = list(compute=TRUE),
                         #control.inla = list(int.strategy = "eb"),
                         verbose = FALSE))

# spatiotemporal model
# fairly memory-intensive
fitSpatTempAR1_BryoDiv <- bru(components = inlaSpatTempAR1C_BryoDiv,
                           likeSpatTempAR_BryoDiv,
                           options = list(
                             # add linear combinations
                             lincomb = c(lc1, lc2),
                             #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                             control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                             control.predictor = list(compute=TRUE),
                             #control.inla = list(int.strategy = "eb"),
                             verbose = FALSE))


fitSpatTempAR2_BryoDiv <- bru(components = inlaSpatTempAR2C_BryoDiv,
                              likeSpatTempAR_BryoDiv,
                              options = list(
                                # add linear combinations
                                lincomb = c(lc1, lc2),
                                #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                                control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                                control.predictor = list(compute=TRUE),
                                #control.inla = list(int.strategy = "eb"),
                                verbose = FALSE))

fitSpatTempAR3_BryoDiv <- bru(components = inlaSpatTempAR3C_BryoDiv,
                              likeSpatTempAR_BryoDiv,
                              options = list(
                                # add linear combinations
                                lincomb = c(lc1, lc2),
                                #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                                control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                                control.predictor = list(compute=TRUE),
                                #control.inla = list(int.strategy = "eb"),
                                verbose = FALSE))

#### Model Output ####
# WAIC, DIC, CPO
# CPO needs no NAs to work perfectly; na.rm is good otherwise
cbind(fitNull_BryoDiv$waic$waic[1], fitNull_BryoDiv$dic$dic[1], -sum(log(fitNull_BryoDiv$cpo$cpo), na.rm = TRUE))
cbind(fitTemp_BryoDiv$waic$waic[1], fitTemp_BryoDiv$dic$dic[1], -sum(log(fitTemp_BryoDiv$cpo$cpo), na.rm = TRUE))
cbind(fitSpatTempAR1_BryoDiv$waic$waic[1], fitSpatTempAR1_BryoDiv$dic$dic[1], -sum(log(fitSpatTempAR1_BryoDiv$cpo$cpo), na.rm = TRUE))
cbind(fitSpatTempAR2_BryoDiv$waic$waic[1], fitSpatTempAR2_BryoDiv$dic$dic[1], -sum(log(fitSpatTempAR2_BryoDiv$cpo$cpo), na.rm = TRUE))
cbind(fitSpatTempAR3_BryoDiv$waic$waic[1], fitSpatTempAR3_BryoDiv$dic$dic[1], -sum(log(fitSpatTempAR3_BryoDiv$cpo$cpo), na.rm = TRUE))
# ar1
#[1,] 12025.41 11994.13 6041.836
#[1,] 12024.58 11993.68 6041.598
# ar2
#[1,] 12003.54 11987.77 6038.287
#[1,] 12007.56 11974.84 6039.013
# ar3
#[1,] 11993.83 11949.38 6030.014
#[1,] 11989.41 11953.99 6027.315
# ar4 -- worse DIC and CPO and only slightly better WAIC (maybe)
#[1,] 11992.9 11967.33 6032.676

# info on priors
#fitSpatTempAR1_BryoDiv$bru_info$model$effects$field
#inla.set.control.fixed.default()
#names(inla.models()$latent$ar)
#names(fitSpatTempAR1_BryoDiv$all.hyper$predictor$hyper$theta)

# summary
summary(fitSpatTempAR1_BryoDiv)

#### Fixed effects ####
# note that there are no observations without at least one of H, M, U, and W community classifications
# since these are non-exclusive, all are slightly positively linked to diversity
brinla::bri.fixed.plot(fitNull_BryoDiv)
brinla::bri.fixed.plot(fitTemp_BryoDiv)
brinla::bri.fixed.plot(fitSpatTempAR1_BryoDiv)
# get custom quantiles
# doesn't work due to Inf and -Inf values
#inla.pmarginal(c(0.05, 0.95), fitSpatTempAR1_BryoDiv$marginals.fitted.values)
#inla.qmarginal(p = c(0.12, 0.88),
#               marginal= fitSpatTempAR1_BryoDiv$marginals.fixed$sheepOff)
# this tests the 0 level for a variable
inla.pmarginal(q = 0,
               marginal= fitSpatTempAR1_BryoDiv$marginals.fixed$sheepOff)
inla.pmarginal(q = 0,
               marginal= fitSpatTempAR1_BryoDiv$marginals.fixed$cowOff)
inla.pmarginal(q = 0,
               marginal= fitSpatTempAR1_BryoDiv$marginals.lincomb.derived$offtake)


#### Marginal distribution ####

plots_BryoDiv <- list()

# Null model
# define empty dataframe
inlaNullPvalT_BryoDiv<-rep(NA, nrow=(vegINLA))
# go through each observation
for(i in 1:nrow(vegINLA)){
  # get bayesian p-value
  inlaNullPvalT_BryoDiv[i]<-inla.pmarginal(q=vegINLA$hillBryo1[i], # compare against actual observed value
                                           # use full posteriors models to compute it
                                           marginal=fitNull_BryoDiv$marginals.fitted.values[[i]])
}

# generate dataset for plotting with predicted, observed, lower and upper CI values
# add the p-value too
postProb1_BryoDiv <- data.frame(predicted = fitNull_BryoDiv$summary.fitted.values$mean[1:length(vegINLA$hillBryo1)],
                                observed = vegINLA$hillBryo1,
                                lower = fitNull_BryoDiv$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo1)],
                                upper = fitNull_BryoDiv$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo1)],
                                p.value = inlaNullPvalT_BryoDiv) 
# get max and min for plots
min1_BryoDiv <- min(postProb1_BryoDiv[, c('upper', 'observed')], na.rm = TRUE)
max1_BryoDiv <- max(postProb1_BryoDiv[, c('upper', 'observed')], na.rm = TRUE)

# posterior probability plot of bayesian p-values
plots_BryoDiv[[1]] <- ggplot(postProb1_BryoDiv, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

# observed vs fitted
plots_BryoDiv[[2]] <- ggplot(postProb1_BryoDiv, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min1_BryoDiv, max1_BryoDiv), y = c(min1_BryoDiv, max1_BryoDiv)) + 
  theme_bw()

# Temporal model -- repeat as above
inlaTempPvalT_BryoDiv<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaTempPvalT_BryoDiv[i]<-inla.pmarginal(q=vegINLA$hillBryo1[i],
                                           # lots of infinity values, so this will not work
                                           marginal=fitTemp_BryoDiv$marginals.fitted.values[[i]])
}
postProb2_BryoDiv <- data.frame(predicted = fitTemp_BryoDiv$summary.fitted.values$mean[1:length(vegINLA$hillBryo1)],
                                observed = vegINLA$hillBryo1,
                                lower = fitTemp_BryoDiv$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo1)],
                                upper = fitTemp_BryoDiv$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo1)],
                                p.value = inlaTempPvalT_BryoDiv) 

min2_BryoDiv <- min(postProb2_BryoDiv[, c('upper', 'observed')], na.rm = TRUE)
max2_BryoDiv <- max(postProb2_BryoDiv[, c('upper', 'observed')], na.rm = TRUE)

plots_BryoDiv[[3]] <- ggplot(postProb2_BryoDiv, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_BryoDiv[[4]] <- ggplot(postProb2_BryoDiv, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min2_BryoDiv, max2_BryoDiv), y = c(min2_BryoDiv, max2_BryoDiv)) + 
  theme_bw()

# Spatial temporal model -- repeat as above
#AR1
inlaSpatTempAR1PvalT_BryoDiv<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaSpatTempAR1PvalT_BryoDiv[i]<-inla.pmarginal(q=vegINLA$hillBryo1[i],
                                               marginal=fitSpatTempAR1_BryoDiv$marginals.fitted.values[[i]])
}

postProb3_BryoDiv <- data.frame(predicted = fitSpatTempAR1_BryoDiv$summary.fitted.values$mean[1:length(vegINLA$hillBryo1)],
                                observed = vegINLA$hillBryo1,
                                lower = fitSpatTempAR1_BryoDiv$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo1)],
                                upper = fitSpatTempAR1_BryoDiv$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo1)],
                                p.value = inlaSpatTempAR1PvalT_BryoDiv) 

min3_BryoDiv <- min(postProb3_BryoDiv[, c('upper', 'observed')], na.rm = TRUE)
max3_BryoDiv <- max(postProb3_BryoDiv[, c('upper', 'observed')], na.rm = TRUE)

plots_BryoDiv[[5]] <- ggplot(postProb3_BryoDiv, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_BryoDiv[[6]] <- ggplot(postProb3_BryoDiv, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min3_BryoDiv, max3_BryoDiv), y = c(min3_BryoDiv, max3_BryoDiv)) + 
  theme_bw()

inlaSpatTempAR2PvalT_BryoDiv<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaSpatTempAR2PvalT_BryoDiv[i]<-inla.pmarginal(q=vegINLA$hillBryo1[i],
                                                  marginal=fitSpatTempAR2_BryoDiv$marginals.fitted.values[[i]])
}

# originally, we plotted here
#cowplot::plot_grid(plotlist = plots_BryoDiv, nrow = 3)

#AR2
inlaSpatTempAR2PvalT_BryoDiv<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaSpatTempAR2PvalT_BryoDiv[i]<-inla.pmarginal(q=vegINLA$hillBryo1[i],
                                                  marginal=fitSpatTempAR2_BryoDiv$marginals.fitted.values[[i]])
}

postProb4_BryoDiv <- data.frame(predicted = fitSpatTempAR2_BryoDiv$summary.fitted.values$mean[1:length(vegINLA$hillBryo1)],
                                observed = vegINLA$hillBryo1,
                                lower = fitSpatTempAR2_BryoDiv$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo1)],
                                upper = fitSpatTempAR2_BryoDiv$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo1)],
                                p.value = inlaSpatTempAR2PvalT_BryoDiv) 

min4_BryoDiv <- min(postProb4_BryoDiv[, c('upper', 'observed')], na.rm = TRUE)
max4_BryoDiv <- max(postProb4_BryoDiv[, c('upper', 'observed')], na.rm = TRUE)

plots_BryoDiv[[7]] <- ggplot(postProb4_BryoDiv, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_BryoDiv[[8]] <- ggplot(postProb4_BryoDiv, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min4_BryoDiv, max4_BryoDiv), y = c(min4_BryoDiv, max4_BryoDiv)) + 
  theme_bw()

#AR3
inlaSpatTempAR3PvalT_BryoDiv<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaSpatTempAR3PvalT_BryoDiv[i]<-inla.pmarginal(q=vegINLA$hillBryo1[i],
                                                  marginal=fitSpatTempAR3_BryoDiv$marginals.fitted.values[[i]])
}

postProb5_BryoDiv <- data.frame(predicted = fitSpatTempAR3_BryoDiv$summary.fitted.values$mean[1:length(vegINLA$hillBryo1)],
                                observed = vegINLA$hillBryo1,
                                lower = fitSpatTempAR3_BryoDiv$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo1)],
                                upper = fitSpatTempAR3_BryoDiv$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo1)],
                                p.value = inlaSpatTempAR3PvalT_BryoDiv) 

min5_BryoDiv <- min(postProb5_BryoDiv[, c('upper', 'observed')], na.rm = TRUE)
max5_BryoDiv <- max(postProb5_BryoDiv[, c('upper', 'observed')], na.rm = TRUE)

plots_BryoDiv[[9]] <- ggplot(postProb5_BryoDiv, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_BryoDiv[[10]] <- ggplot(postProb5_BryoDiv, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min5_BryoDiv, max5_BryoDiv), y = c(min5_BryoDiv, max5_BryoDiv)) + 
  theme_bw()

# plot everything for comparison
cowplot::plot_grid(plotlist = plots_BryoDiv, nrow = 5)

#### Plotting ####
# hyperparameter plots
hyper_BryoDiv <- data.frame(p = fitSpatTempAR3_BryoDiv$marginals.hyperpar$`p parameter for Tweedie`,
                            dispersion = fitSpatTempAR3_BryoDiv$marginals.hyperpar$`Dispersion parameter for Tweedie`)
# initial conditions don't seem to affect things too much; lots of information about the parametrisation
# mean and CI for hyperparameters

# non-negative continuous if 1>p>2 w/ zero inflation
# p = 1; quasipoisson; p = 2 is a gamma
ggplot(hyper_BryoDiv, aes(x =p.x, y = p.y)) + geom_line() + theme_bw()
ggplot(hyper_BryoDiv, aes(x =dispersion.x, y = dispersion.y)) + geom_line() + theme_bw()

# spatial prediction plots

# create data for 2023 predictions and block group labels for plotting later
vegPred <- st_drop_geometry(vegINLA) %>% ungroup() %>% filter(Year == 2023) %>% 
  mutate(blockGroup = factor(Block,
                             levels = LETTERS[1:6],
                             labels = c("AB", "AB", "CD", "CD", "EF", "EF")))

# create mesh only within the study area and not intermediate mesh
ppxl_BryoDiv <- fm_pixels(mesh_BryoDiv, mask = plotBounds)
# add the raw data for predictions
ppxl_all_BryoDiv <- fm_cprod(ppxl, vegPred %>% select(c(Year, sheepOff, cowOff, windSpeed, meanT, solarRad, precip, litter, rock,
                                                        elev, isW, isH, isU, isM, hillBryo1)))

# run the built-in inlabru function to get predictions
spatPred_BryoDiv <- predict(
  fitSpatTempAR3_BryoDiv, ppxl_all,
  # use the full formula with all covariates
  ~ field + Intercept + sheepOff + cowOff + windSpeed + meanT + solarRad + precip + litter + rock + elev + isW + isH + isU + isM
)

# modify data to get it to the response scale for mean and lower/upper CIs  
spatPred_BryoDiv <- spatPred_BryoDiv %>% mutate(linkMean = exp(mean),
             linkLower = exp(q0.025),
             linkUpper = exp(q0.975)) %>% 
  # roundabout way using spatial checks to categorise each observation to the block group
  mutate(ab = ifelse(st_intersects(geometry, plotBounds %>% filter(blockGroup == "A&B")),"A&B",""),
         cd = ifelse(st_intersects(geometry, plotBounds %>% filter(blockGroup == "C&D")),"C&D",""),
         ef = ifelse(st_intersects(geometry, plotBounds %>% filter(blockGroup == "E&F")),"E&F",""),
         blockGroup = paste0(ab, cd, ef))

# plot each block group separately bc geom_sf and gg don't support facet_wrapping
# specific plot elements are "sausage making" and not really relevant beyond visibility
spatPlots <- list()

spatPlots[[1]] <- ggplot() +
  gg(lambda1 %>% filter(blockGroup == "NAC&DNA") %>% relocate(linkLower), geom = "tile") +
  scale_fill_viridis_c(begin = 0, end = 1, limits=c(min(lambda1$linkLower), max(lambda1$linkUpper))) + 
  theme_bw() + 
  theme(axis.title = element_blank(),
        panel.grid = element_line(colour = "#ABB0B8",
                                  linetype = 2),
        legend.position = 'none') +
  scale_x_continuous(breaks = seq(-4.453, -4.447, by = 0.002),
                     minor_breaks = seq(-4.453, -4.446, by = 0.001)) +
  scale_y_continuous(breaks = seq(56.275, 56.281, by = 0.002),
                     minor_breaks = seq(56.274, 56.281, by = 0.001)) 

spatPlots[[2]] <- ggplot() +
  gg(lambda1 %>% filter(blockGroup == "NANAE&F") %>% relocate(linkLower), geom = "tile") +
  scale_fill_viridis_c(begin = 0, end = 1, limits=c(min(lambda1$linkLower), max(lambda1$linkUpper))) + 
  theme_bw() + 
  theme(axis.title = element_blank(),
        panel.grid = element_line(colour = "#ABB0B8",
                                  linetype = 2),
        legend.position = 'none') +
  scale_x_continuous(minor_breaks = seq(-4.408, -4.396, by = 0.001)) +
  scale_y_continuous(breaks = seq(56.295, 56.301, by = 0.002),
                     minor_breaks = seq(56.295, 56.301, by = 0.001)) 

spatPlots[[3]] <- ggplot() +
  gg(lambda1 %>% filter(blockGroup == "A&BNANA") %>% relocate(linkLower), geom = "tile") +
  scale_fill_viridis_c(begin = 0, end = 1, limits=c(min(lambda1$linkLower), max(lambda1$linkUpper))) + 
  theme_bw() + 
  theme(axis.title = element_blank(),
        panel.grid = element_line(colour = "#ABB0B8",
                                  linetype = 2),
        panel.grid.minor = element_line(colour = "#CCCCCC",
                                        linetype = 2)) +
  labs(fill='Pred Bryo\nDiversity') +
  scale_x_continuous(breaks = seq(-4.380, -4.374, by = 0.002),
                     minor_breaks = seq(-4.381, -4.374, by = 0.001)) +
  scale_y_continuous(minor_breaks = seq(56.262, 56.270, by = 0.001)) 

spatPlots[[4]] <- ggplot() +
  gg(lambda1 %>% filter(blockGroup == "NAC&DNA") %>% relocate(linkMean), geom = "tile") +
  scale_fill_viridis_c(begin = 0, end = 1, limits=c(min(lambda1$linkLower), max(lambda1$linkUpper))) + 
  theme_bw() + 
  theme(axis.title = element_blank(),
        panel.grid = element_line(colour = "#ABB0B8",
                                  linetype = 2),
        legend.position = 'none') +
  scale_x_continuous(breaks = seq(-4.453, -4.447, by = 0.002),
                     minor_breaks = seq(-4.453, -4.446, by = 0.001)) +
  scale_y_continuous(breaks = seq(56.275, 56.281, by = 0.002),
                     minor_breaks = seq(56.274, 56.281, by = 0.001)) 

spatPlots[[5]] <- ggplot() +
  gg(lambda1 %>% filter(blockGroup == "NANAE&F") %>% relocate(linkMean), geom = "tile") +
  scale_fill_viridis_c(begin = 0, end = 1, limits=c(min(lambda1$linkLower), max(lambda1$linkUpper))) + 
  theme_bw() + 
  theme(axis.title = element_blank(),
        panel.grid = element_line(colour = "#ABB0B8",
                                  linetype = 2),
        legend.position = 'none') +
  scale_x_continuous(minor_breaks = seq(-4.408, -4.396, by = 0.001)) +
  scale_y_continuous(breaks = seq(56.295, 56.301, by = 0.002),
                     minor_breaks = seq(56.295, 56.301, by = 0.001)) 

spatPlots[[6]] <- ggplot() +
  gg(lambda1 %>% filter(blockGroup == "A&BNANA") %>% relocate(linkMean), geom = "tile") +
  scale_fill_viridis_c(begin = 0, end = 1, limits=c(min(lambda1$linkLower), max(lambda1$linkUpper))) + 
  theme_bw() + 
  theme(axis.title = element_blank(),
        panel.grid = element_line(colour = "#ABB0B8",
                                  linetype = 2),
        panel.grid.minor = element_line(colour = "#CCCCCC",
                                        linetype = 2)) +
  labs(fill='Pred Bryo\nDiversity') +
  scale_x_continuous(breaks = seq(-4.380, -4.374, by = 0.002),
                     minor_breaks = seq(-4.381, -4.374, by = 0.001)) +
  scale_y_continuous(minor_breaks = seq(56.262, 56.270, by = 0.001)) 

spatPlots[[7]] <- ggplot() +
  gg(lambda1 %>% filter(blockGroup == "NAC&DNA") %>% relocate(linkUpper), geom = "tile") +
  scale_fill_viridis_c(begin = 0, end = 1, limits=c(min(lambda1$linkLower), max(lambda1$linkUpper))) + 
  theme_bw() + 
  theme(axis.title = element_blank(),
        panel.grid = element_line(colour = "#ABB0B8",
                                  linetype = 2),
        legend.position = 'none') +
  scale_x_continuous(breaks = seq(-4.453, -4.447, by = 0.002),
                     minor_breaks = seq(-4.453, -4.446, by = 0.001)) +
  scale_y_continuous(breaks = seq(56.275, 56.281, by = 0.002),
                     minor_breaks = seq(56.274, 56.281, by = 0.001)) 

spatPlots[[8]] <- ggplot() +
  gg(lambda1 %>% filter(blockGroup == "NANAE&F") %>% relocate(linkUpper), geom = "tile") +
  scale_fill_viridis_c(begin = 0, end = 1, limits=c(min(lambda1$linkLower), max(lambda1$linkUpper))) + 
  theme_bw() + 
  theme(axis.title = element_blank(),
        panel.grid = element_line(colour = "#ABB0B8",
                                  linetype = 2),
        legend.position = 'none') +
  scale_x_continuous(minor_breaks = seq(-4.408, -4.396, by = 0.001)) +
  scale_y_continuous(breaks = seq(56.295, 56.301, by = 0.002),
                     minor_breaks = seq(56.295, 56.301, by = 0.001)) 

spatPlots[[9]] <- ggplot() +
  gg(lambda1 %>% filter(blockGroup == "A&BNANA") %>% relocate(linkUpper), geom = "tile") +
  scale_fill_viridis_c(begin = 0, end = 1, limits=c(min(lambda1$linkLower), max(lambda1$linkUpper))) + 
  theme_bw() + 
  theme(axis.title = element_blank(),
        panel.grid = element_line(colour = "#ABB0B8",
                                  linetype = 2),
        panel.grid.minor = element_line(colour = "#CCCCCC",
                                        linetype = 2)) +
  labs(fill='Pred Bryo\nDiversity') +
  scale_x_continuous(breaks = seq(-4.380, -4.374, by = 0.002),
                     minor_breaks = seq(-4.381, -4.374, by = 0.001)) +
  scale_y_continuous(minor_breaks = seq(56.262, 56.270, by = 0.001)) 

# plot everything together using cowplot() to get it in a fancy grid
# 3 rows for comparing all sub-categories
cowplot::plot_grid(plotlist = spatPlots, nrow = 3,
                   # add labels for figure caption too
                   labels = c('A','','','B','','','C','',''), label_size = 10,
                   align = "h", axis = "b")
# can align widths to make things nicer, but not needed for now
#, rel_widths = c(0.9495,1.5,1.0105))

# test the fit with kriged diversity -- quite bad, but not sure which is "correct
# ggplot(lambda1, aes(x = linkMean, y = hillBryo1)) + geom_point() +
#        geom_abline(slope = 1, intercept = 0) +
#        labs(y = "Observed", x = "Fitted") +
#        lims(x = c(0, 7), y = c(0, 7)) + 
#        theme_bw()

#### Exports ####
write_rds(fitNull_BryoDiv, "./Models/Objects/fitNull_BryoDiv.rds")
write_rds(fitTemp_BryoDiv, "./Models/Objects/fitTemp_BryoDiv.rds")
write_rds(fitSpatTempAR1_BryoDiv, "./Models/Objects/fitSpatTempAR1_BryoDiv.rds")
write_rds(fitSpatTempAR2_BryoDiv, "./Models/Objects/fitSpatTempAR2_BryoDiv.rds")
write_rds(fitSpatTempAR3_BryoDiv, "./Models/Objects/fitSpatTempAR3_BryoDiv.rds")

write_rds(spatPred_BryoDiv, "./Models/Objects/bryoSpatPreds.rds")
