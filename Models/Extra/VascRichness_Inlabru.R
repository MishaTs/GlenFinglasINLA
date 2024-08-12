library(tidyverse)
library(INLA)
library(fmesher)
library(inlabru)

#### Data Setup ####
vegGeo <- read_rds("./Vegetation/vegGeo.rds")

vegINLA <- vegGeo %>% mutate(fBlock = as.factor(Block),
                             Plot_code = paste0(Plot, Block),
                             treat = factor(Plot,
                                            levels = c("III", "I", "II", "IV"))) #%>% filter(!is.na(nvcM1)) 

# t(apply(table(vegINLA$isW, vegINLA$flower),1, function(x) x/sum(x)))

# first, the mesh
mesh_VascRich <- fm_mesh_2d_inla(
  loc = vegINLA$geometry, 
  # doesn't look like any changes where max.edge > 0.05 have much impact
  # much finer mesh possible otherwise, but left as-is for computing benefits for now
  # outer edge set to 1 but doesn't matter
  max.edge = c(0.05,1), # map units inside and outside
  # cutoff is min edge
  #cutoff = 50, offset = c(100, 300),
  crs = fm_crs(vegINLA)
)
ggplot() + gg(mesh_VascRich) + theme_bw()
#mesh_VascRich$n # 2044 points

# Making SPDE
# Quantify distance between points and define Matern correlation on the mesh
#inla.spde2.pcmatern allows for penalised complexity priors
aMatSPDE_VascRich <- inla.spde2.matern(mesh = mesh_VascRich)#, 
#prior.range = c(10, 0.5), 
#prior.sigma = c(.5, .5))#,
#alpha = 2) 
# code from slide 10 of Finn's ppt
# https://www.maths.ed.ac.uk/~flindgre/talks/Lindgren_RSS2023.pdf

# map across years for the temporal component
yearMap_VascRich <- bru_mapper(fm_mesh_1d(sort(unique(vegINLA$Year))), indexed = TRUE)

#### Modelling ####
# first define the model components
# high-leverage variables end up being quite impactful (e.g., mostly 0 with a few rows of 1)
inlaNullC_VascRich <- ~ #Intercept(1) + 
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

inlaTempC_VascRich <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + rock + #graze +
  elev +
  isW + isH + isU + isM +
  field1(Year, model='ar1', replicate = Point) +
  field2(fBlock, model = "iid") + 
  field3(Plot_code, model = "iid")

inlaSpatTempAR1C_VascRich <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + rock + #graze +
  elev +
  isW + isH + isU + isM +
  field(geometry, 
        model = aMatSPDE_VascRich, 
        group = Year, 
        group_mapper = yearMap_VascRich,
        control.group = list(model = "ar",
                             # order <= 10
                             order = 1))

inlaSpatTempAR2C_VascRich <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + rock + #graze +
  elev +
  isW + isH + isU + isM +
  field(geometry, 
        model = aMatSPDE_VascRich, 
        group = Year, 
        group_mapper = yearMap_VascRich,
        control.group = list(model = "ar",
                             # order <= 10
                             order = 2))

inlaSpatTempAR3C_VascRich <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + rock + #graze +
  elev +
  isW + isH + isU + isM +
  field(geometry, 
        model = aMatSPDE_VascRich, 
        group = Year, 
        group_mapper = yearMap_VascRich,
        control.group = list(model = "ar",
                             # order <= 10
                             order = 3))

# define linear combinations
lc1 <- inla.make.lincomb(sheepOff = 1, cowOff = 1)
names(lc1) <- "offtake"
lc2 <- inla.make.lincomb(isW = 1, isH = 1, isU = 1, isM = 1)
names(lc2) <- "nvc"

# build the likelihoods
likeNull_VascRich <- like(
  formula = hillVasc0 ~ .,
  # full distribution list: inla.list.models()
  # possible distributions: nbionimal and poisson
  family = "nbinomial", 
  data=vegINLA,
  # domain optional for non-spatial model
  domain = list(geometry = mesh_VascRich)
)

likeTemp_VascRich <- like(
  formula = hillVasc0 ~ .,
  # tweedie is a superset of gaussian, gamma, poisson, and inverse gaussians
  family = "nbinomial", 
  data=vegINLA,
  # domain optional for non-spatial model
  domain = list(geometry = mesh_VascRich)
)

likeSpatTempAR_VascRich <- like(
  formula = hillVasc0 ~ .,
  # tweedie is a superset of gaussian, gamma, poisson, and inverse gaussians
  family = "nbinomial", 
  data=vegINLA,
  domain = list(geometry = mesh_VascRich)
)

# run the actual fit
# null model
fitNull_VascRich <- bru(components = inlaNullC_VascRich,
                       likeNull_VascRich,
                       options = list(
                         # add linear combinations
                         lincomb = c(lc1, lc2),
                         #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                         control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                         # might want to consider to setting ``control.predictor=list(link=...)''
                         # otherwise fitted.values in the result-object computed w/ identity for NAs
                         # https://rdrr.io/github/inbo/INLA/man/control.predictor.html
                         control.predictor = list(compute=TRUE),
                         #control.inla = list(int.strategy = "eb"),
                         verbose = FALSE))

# temporal model
fitTemp_VascRich <- bru(components = inlaTempC_VascRich,
                       likeTemp_VascRich,
                       options = list(
                         # add linear combinations
                         lincomb = c(lc1, lc2),
                         #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                         control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                         # might want to consider to setting ``control.predictor=list(link=...)''
                         # otherwise fitted.values in the result-object computed w/ identity for NAs
                         # https://rdrr.io/github/inbo/INLA/man/control.predictor.html
                         control.predictor = list(compute=TRUE),
                         #control.inla = list(int.strategy = "eb"),
                         verbose = FALSE))

# spatiotemporal model
# fairly memory-intensive
fitSpatTempAR1_VascRich <- bru(components = inlaSpatTempAR1C_VascRich,
                              likeSpatTempAR_VascRich,
                              options = list(
                                # add linear combinations
                                lincomb = c(lc1, lc2),
                                #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                                control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                                # might want to consider to setting ``control.predictor=list(link=...)''
                                # otherwise fitted.values in the result-object computed w/ identity for NAs
                                # https://rdrr.io/github/inbo/INLA/man/control.predictor.html
                                control.predictor = list(compute=TRUE),
                                #control.inla = list(int.strategy = "eb"),
                                verbose = FALSE))

fitSpatTempAR2_VascRich <- bru(components = inlaSpatTempAR2C_VascRich,
                              likeSpatTempAR_VascRich,
                              options = list(
                                # add linear combinations
                                lincomb = c(lc1, lc2),
                                #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                                control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                                control.predictor = list(compute=TRUE),
                                #control.inla = list(int.strategy = "eb"),
                                verbose = FALSE))

fitSpatTempAR3_VascRich <- bru(components = inlaSpatTempAR3C_VascRich,
                              likeSpatTempAR_VascRich,
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
cbind(fitNull_VascRich$waic$waic[1], fitNull_VascRich$dic$dic[1], -sum(log(fitNull_VascRich$cpo$cpo), na.rm = TRUE))
cbind(fitTemp_VascRich$waic$waic[1], fitTemp_VascRich$dic$dic[1], -sum(log(fitTemp_VascRich$cpo$cpo), na.rm = TRUE))
cbind(fitSpatTempAR1_VascRich$waic$waic[1], fitSpatTempAR1_VascRich$dic$dic[1], -sum(log(fitSpatTempAR1_VascRich$cpo$cpo), na.rm = TRUE))
#[1,] 20618.84 20620.62 10309.44
#[1,] 20175.77 20221.57 10088
#[1,] 19917.74 19955.97 9958.928

cbind(fitSpatTempAR2_VascRich$waic$waic[1], fitSpatTempAR2_VascRich$dic$dic[1], -sum(log(fitSpatTempAR2_VascRich$cpo$cpo), na.rm = TRUE))
cbind(fitSpatTempAR3_VascRich$waic$waic[1], fitSpatTempAR3_VascRich$dic$dic[1], -sum(log(fitSpatTempAR3_VascRich$cpo$cpo), na.rm = TRUE))
#[1,] 19916.53 19954.73 9958.321
#[1,] 19916.6 19954.19 9958.353

# summary
summary(fitSpatTempAR1_VascRich)

#### Fixed effects ####
# note that there are no observations without at least one of H, M, U, and W community classifications
# since these are non-exclusive, all are slightly positively linked to diversity
brinla::bri.fixed.plot(fitNull_VascRich)
brinla::bri.fixed.plot(fitTemp_VascRich)
brinla::bri.fixed.plot(fitSpatTempAR1_VascRich)
# get custom quantiles
# doesn't work due to Inf and -Inf values
#inla.pmarginal(c(0.05, 0.95), fitSpatTempAR1_VascRich$marginals.fitted.values)
#inla.qmarginal(p = c(0.12, 0.88),
#               marginal= fitSpatTempAR1_VascRich$marginals.fixed$sheepOff)
# this tests the 0 level for a variable
inla.pmarginal(q = 0,
               marginal= fitSpatTempAR1_VascRich$marginals.fixed$sheepOff)
inla.pmarginal(q = 0,
               marginal= fitSpatTempAR1_VascRich$marginals.fixed$cowOff)
inla.pmarginal(q = 0,
               marginal= fitSpatTempAR1_VascRich$marginals.lincomb.derived$offtake)


#### Marginal distribution ####

plots_VascRich <- list()

# Null model
inlaNullPvalT_VascRich<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaNullPvalT_VascRich[i]<-inla.pmarginal(q=vegINLA$hillVasc0[i],
                                           marginal=fitNull_VascRich$marginals.fitted.values[[i]])
}

postProb1_VascRich <- data.frame(predicted = fitNull_VascRich$summary.fitted.values$mean[1:length(vegINLA$hillVasc0)],
                                observed = vegINLA$hillVasc0,
                                lower = fitNull_VascRich$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillVasc0)],
                                upper = fitNull_VascRich$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillVasc0)],
                                p.value = inlaNullPvalT_VascRich) 

min1_VascRich <- min(postProb1_VascRich[, c('upper', 'observed')], na.rm = TRUE)
max1_VascRich <- max(postProb1_VascRich[, c('upper', 'observed')], na.rm = TRUE)

plots_VascRich[[1]] <- ggplot(postProb1_VascRich, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_VascRich[[2]] <- ggplot(postProb1_VascRich, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min1_VascRich, max1_VascRich), y = c(min1_VascRich, max1_VascRich)) + 
  theme_bw()

# Temporal model
inlaTempPvalT_VascRich<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaTempPvalT_VascRich[i]<-inla.pmarginal(q=vegINLA$hillVasc0[i],
                                           # lots of infinity values, so this will not work
                                           marginal=fitTemp_VascRich$marginals.fitted.values[[i]])
}
postProb2_VascRich <- data.frame(predicted = fitTemp_VascRich$summary.fitted.values$mean[1:length(vegINLA$hillVasc0)],
                                observed = vegINLA$hillVasc0,
                                lower = fitTemp_VascRich$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillVasc0)],
                                upper = fitTemp_VascRich$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillVasc0)],
                                p.value = inlaTempPvalT_VascRich) 

min2_VascRich <- min(postProb2_VascRich[, c('upper', 'observed')], na.rm = TRUE)
max2_VascRich <- max(postProb2_VascRich[, c('upper', 'observed')], na.rm = TRUE)

plots_VascRich[[3]] <- ggplot(postProb2_VascRich, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_VascRich[[4]] <- ggplot(postProb2_VascRich, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min2_VascRich, max2_VascRich), y = c(min2_VascRich, max2_VascRich)) + 
  theme_bw()

# Spatial temporal model
inlaSpatTempAR1PvalT_VascRich<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaSpatTempAR1PvalT_VascRich[i]<-inla.pmarginal(q=vegINLA$hillVasc0[i],
                                                  marginal=fitSpatTempAR1_VascRich$marginals.fitted.values[[i]])
}

postProb3_VascRich <- data.frame(predicted = fitSpatTempAR1_VascRich$summary.fitted.values$mean[1:length(vegINLA$hillVasc0)],
                                observed = vegINLA$hillVasc0,
                                lower = fitSpatTempAR1_VascRich$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillVasc0)],
                                upper = fitSpatTempAR1_VascRich$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillVasc0)],
                                p.value = inlaSpatTempAR1PvalT_VascRich) 

min3_VascRich <- min(postProb3_VascRich[, c('upper', 'observed')], na.rm = TRUE)
max3_VascRich <- max(postProb3_VascRich[, c('upper', 'observed')], na.rm = TRUE)

plots_VascRich[[5]] <- ggplot(postProb3_VascRich, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_VascRich[[6]] <- ggplot(postProb3_VascRich, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min3_VascRich, max3_VascRich), y = c(min3_VascRich, max3_VascRich)) + 
  theme_bw()

cowplot::plot_grid(plotlist = plots_VascRich, nrow = 3)

inlaSpatTempAR2PvalT_VascRich<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaSpatTempAR2PvalT_VascRich[i]<-inla.pmarginal(q=vegINLA$hillVasc0[i],
                                                  marginal=fitSpatTempAR2_VascRich$marginals.fitted.values[[i]])
}

postProb4_VascRich <- data.frame(predicted = fitSpatTempAR2_VascRich$summary.fitted.values$mean[1:length(vegINLA$hillVasc0)],
                                observed = vegINLA$hillVasc0,
                                lower = fitSpatTempAR2_VascRich$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillVasc0)],
                                upper = fitSpatTempAR2_VascRich$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillVasc0)],
                                p.value = inlaSpatTempAR2PvalT_VascRich) 

min4_VascRich <- min(postProb4_VascRich[, c('upper', 'observed')], na.rm = TRUE)
max4_VascRich <- max(postProb4_VascRich[, c('upper', 'observed')], na.rm = TRUE)

plots_VascRich[[7]] <- ggplot(postProb4_VascRich, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_VascRich[[8]] <- ggplot(postProb4_VascRich, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min4_VascRich, max4_VascRich), y = c(min4_VascRich, max4_VascRich)) + 
  theme_bw()

inlaSpatTempAR3PvalT_VascRich<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaSpatTempAR3PvalT_VascRich[i]<-inla.pmarginal(q=vegINLA$hillVasc0[i],
                                                  marginal=fitSpatTempAR3_VascRich$marginals.fitted.values[[i]])
}

postProb5_VascRich <- data.frame(predicted = fitSpatTempAR3_VascRich$summary.fitted.values$mean[1:length(vegINLA$hillVasc0)],
                                observed = vegINLA$hillVasc0,
                                lower = fitSpatTempAR3_VascRich$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillVasc0)],
                                upper = fitSpatTempAR3_VascRich$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillVasc0)],
                                p.value = inlaSpatTempAR3PvalT_VascRich) 

min5_VascRich <- min(postProb5_VascRich[, c('upper', 'observed')], na.rm = TRUE)
max5_VascRich <- max(postProb5_VascRich[, c('upper', 'observed')], na.rm = TRUE)

plots_VascRich[[9]] <- ggplot(postProb5_VascRich, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_VascRich[[10]] <- ggplot(postProb5_VascRich, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min5_VascRich, max5_VascRich), y = c(min5_VascRich, max5_VascRich)) + 
  theme_bw()

cowplot::plot_grid(plotlist = plots_VascRich, nrow = 5)




#### Plotting ####
# currently non-functional
pix_VascRich <- fm_pixels(mesh_VascRich)
#Warning in handle_problems(e_input) :
#The input evaluation 'Year' for 'field.group' failed. Perhaps the data object doesn't contain the needed variables? Falling back to '1'.
pred_VascRich <- predict(
  fitSpatTempAR1_VascRich, pix_VascRich,
  # it looks like we need to provide actual values for covariates (and year) for predictions to work
  ~ field #+ Intercept + sheepOff + cowOff + windSpeed + meanT + solarRad + precip + litter + rock + elev + isW + isH + isU + isM 
)
samp <- generate(fitSpatTempAR1_VascRich, pix_VascRich,
                 ~ field + Intercept + sheepOff + cowOff + windSpeed + meanT + solarRad + precip + litter + rock + elev + isW + isH + isU + isM,
                 n.samples = 1
)
pred_VascRich$sample <- samp[, 1]


#### Exports ####
write_rds(fitNull_VascRich, "./Models/Objects/fitNull_VascRich.rds")
write_rds(fitTemp_VascRich, "./Models/Objects/fitTemp_VascRich.rds")
write_rds(fitSpatTempAR1_VascRich, "./Models/Objects/fitSpatTempAR1_VascRich.rds")
write_rds(fitSpatTempAR2_VascRich, "./Models/Objects/fitSpatTempAR2_VascRich.rds")
write_rds(fitSpatTempAR3_VascRich, "./Models/Objects/fitSpatTempAR3_VascRich.rds")