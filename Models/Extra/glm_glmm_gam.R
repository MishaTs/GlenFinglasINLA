library(nlme)
library(tidyverse)
colnames(vegGeo)


# get a point ID for error grouping and filter out NAs for the quantative variable beause lme4 and nmle both struggle w/ that
vegReg <- vegGeo %>% mutate(pointID = paste0(Block, as.numeric(as.roman(Plot)), "_", ifelse(Point < 10, "0", ""), Point)) %>% 
  filter(!is.na(hill2))

firstLM <- glm(isW ~ sheepOff + cowOff + meanT + elev + precip + lat + long + Year,
               family = binomial(link = "logit"),
               data = vegGeo)
summary(firstLM)

# really basic model for testing biodiversity measurements
# looks quite good with indicator variables, but this is only linear GLMMs...so we'd want a logistic 
# for hill0, hill1, hill2, significance is latitude and sheep offtake mainly; not much otherwise
# adding interaction terms isn't very helpful either
# as expected, a spatially explicit model would fare best here
testMM <- lme(hill2 ~ sheepOff + cowOff + meanT + elev + precip + lat + long + Year, random=list(pointID=pdDiag(~Year)),
              data = vegReg, 
              method="REML", correlation = corAR1())
# high-ish correlation between cow & sheep offtake; quite high correlation between latitude, longitude, and elevation
summary(testMM)


# now some GAMs
# first just TPS spatial smooths
library(mgcv)
library(fields)
library(gratia)
# get the treatment as a factor first since keeping it as continuous isn't worth the struggle
vegGeo <- vegGeo %>% mutate(plot = as.factor(Plot))
# set 'IV' rather than 'I' as reference level 
vegGeo$plot <- relevel(vegGeo$plot, ref = "IV") 
write_rds(vegGeo, "./Vegetation/vegCleaned.rds")

# start with independent spatial smooths
# by = "plot" for now but the temporal component is also still missing....
# parameterisation of covariates leaves much to be desired
# no covariates leads to low EDF for everything except plots III and IV in latitude
fit_1Dsmooth <- bam(isW ~ s(long, by=plot) + s(lat, by=plot) + 
                      # covariates below
                      s(meanT) + s(elev) + s(precip) + s(Year), 
                    family = binomial, data = vegGeo, 
                    # speed it up by a lot
                    discrete = TRUE, nthreads = 8)
# still poor performance
summary(fit_1Dsmooth)

# Write function to produce quilt plots for plot differences
spatial_plots <- function(fitted_model, name){
  
  # Predict
  vegGeo$Estimate <- predict(fitted_model, type="response")
  
  # make plots  
  par(mfrow=c(2,2), oma=c(0,0,0,0))
  # I
  fields::quilt.plot(vegGeo$long[vegGeo$plot=="I"],
                     vegGeo$lat[vegGeo$plot=="I"], 
                     vegGeo$Estimate[vegGeo$plot=="I"], 
                     main=paste0(name, ": Plot I"),
                     zlim=range(vegGeo$Estimate),
                     nx = 20, ny = 20,
                     asp=1, add.legend = TRUE)
  # II
  fields::quilt.plot(vegGeo$long[vegGeo$plot=="II"],
                     vegGeo$lat[vegGeo$plot=="II"], 
                     vegGeo$Estimate[vegGeo$plot=="II"], 
                     main=paste0(name, ": Plot II"),
                     zlim=range(vegGeo$Estimate),
                     nx = 20, ny = 20,
                     asp=1, add.legend = TRUE)
  # III
  fields::quilt.plot(vegGeo$long[vegGeo$plot=="III"],
                     vegGeo$lat[vegGeo$plot=="III"], 
                     vegGeo$Estimate[vegGeo$plot=="III"], 
                     main=paste0(name, ": Plot III"),
                     zlim=range(vegGeo$Estimate),
                     nx = 20, ny = 20,
                     asp=1, add.legend = TRUE)
  # IV
  fields::quilt.plot(vegGeo$long[vegGeo$plot=="IV"],
                     vegGeo$lat[vegGeo$plot=="IV"], 
                     vegGeo$Estimate[vegGeo$plot=="IV"], 
                     main=paste0(name, ": Plot IV"),
                     zlim=range(vegGeo$Estimate),
                     nx = 20, ny = 20,
                     asp=1, add.legend = TRUE)
}
# looks like convergence is a mess
spatial_plots(fitted_model = fit_1Dsmooth, name = "Smooth (1D)")
#yeah....this is quite bad
gratia::draw(fit_1Dsmooth, rug = FALSE) & theme_bw()

# now fit a higher-dimensional smooth
# no covariates first
# by = plot doesn't make too much sense here: why would latitude and longitude have different effects depending on the plot?
# including it gives "fitted probabilities numerically 0 or 1 occurred" and model doesn't fit
fit_2Dsmooth <- bam(isW ~ plot + s(lat, long), family = binomial, data = vegGeo,
                    discrete = TRUE, nthreads = 12)
summary(fit_2Dsmooth)
gratia::draw(fit_2Dsmooth, rug=FALSE) & theme_bw()

fit_3Dsmooth <- bam(isW ~ plot + s(lat, long, Year), family = binomial, data = vegGeo,
                    discrete = TRUE, nthreads = 12)
summary(fit_3Dsmooth)
# it looks like this parametrisation is just not good
gratia::draw(fit_3Dsmooth, rug=FALSE) & theme_bw()

fit_3Dsmooth2 <- bam(isW ~ plot + s(lat, long, elev), family = binomial, data = vegGeo,
                     discrete = TRUE, nthreads = 12)
summary(fit_3Dsmooth2)
# elevation seems a bit better...but harder to interpret what's happening
gratia::draw(fit_3Dsmooth2, rug=FALSE) & theme_bw()

# retry the 2D with covariates
fit_2DCov <- bam(isW ~ plot + nvcSumComs + nvcVars + meanT + precip + s(lat, long) + 
                   s(elev) + s(Year), family = binomial, data = vegGeo,
                 discrete = TRUE, nthreads = 12)
summary(fit_2DCov)
# not anything blatantly nonlinear among elevation and year
gratia::draw(fit_2DCov, rug=FALSE) & theme_bw()

# high-dimensional smooths in GAM are too uninterpretable and subtract from the signal of interest
#fit_4Dsmooth <- bam(isW ~ plot + s(lat, long, elev, Year), family = binomial, data = vegGeo,
#                    discrete = TRUE, nthreads = 12)
#summary(fit_4Dsmooth)

# this looks the best at an initial glance
#fit_3DCov <- bam(isW ~ plot + s(lat, long, Year) + s(meanT) + s(elev) + s(precip), family = binomial, data = vegGeo,
#                    discrete = TRUE, nthreads = 12)
#summary(fit_3DCov)

# this model is likely overparameterised
# no real significance to any of the "plot" terms
#fit_4DCov <- bam(isW ~ plot + s(lat, long, elev, Year) + s(meanT) + s(precip), family = binomial, data = vegGeo,
#                 discrete = TRUE, nthreads = 12)
#summary(fit_4DCov)

# follow it up with CReSS
#library(remotes)
#remotes::install_github("lindesaysh/MRSea")
library(MRSea)
#library(splines)

# define the knots
knotGrid <- getKnotgrid(coordData=cbind(vegGeo$lat, vegGeo$long), # coordinates
                        numKnots=300, # no. knots
                        plot=TRUE) # plot result?

# playing around with other response variables is fruitless
# data does not contain successes column: continuous htAvg
# plot cannot be fitted as a factor variable, due to zero counts for some levels: binary isW
vegGeo$response <- vegGeo$nvcSumComs # function needs variable called 'response' 
vegGeo$x.pos <- vegGeo$long # function needs x coord called 'x.pos' 
vegGeo$y.pos <- vegGeo$lat # function needs y coord called 'y.pos'
# Set initial model without spline-based terms
initialModel <- glm(response ~ plot + as.factor(Year), family=quasipoisson, data=vegGeo)

# Set SALSA arguments
factorList <- c("plot","Year")
varList <- c("elev")
salsa1DList <- list(fitnessMeasure = "QAIC", # criteria for k-fold CV
                    minKnots_1d = c(2), # minimum number of knots
                    maxKnots_1d = c(4), # maximum number of knots
                    startKnots_1d = c(2), # starting number (equidistant)
                    degree = c(2), # degree of B-spline
                    #maxIterations = c(100), # max iterations of algorithm
                    gaps=c(0)) # min gaps between knots

# Run SALSA 1D
# this does not work when sf is loaded in the namespace
#detach("package:sf", unload=TRUE)
# detaching sf doesn't help, so it must be one of class_7.3-22, KernSmooth_2.23-22, classInt_0.4-10, e1071_1.7-14 which are in the namespace
# conversely, dotCall64_1.1-1 is removed from the namespace which may be the issue?
salsa1D <- runSALSA1D(initialModel, # starting model (GLM)
                      salsa1DList, # SALSA parameters from above
                      varList, # variable for 1D knot selection
                      factorList, # any factors in the model
                      datain=vegGeo) # name of the data

# Store "best" model
bestModel1D <- salsa1D$bestModel

# Create data-to-knot distance matrix
dMat <- makeDists(datacoords = cbind(vegGeo$x.pos, # X coord
                                     vegGeo$y.pos), # Y coord
                  knotcoords = knotGrid) # starting knots

# Set SALSA parameters (inc defining interaction)
salsa2DList <- list(fitnessMeasure = "QAIC", # criteria for k-fold CV
                    knotgrid=knotGrid, # starting knot locations
                    startKnots = 6, # starting number (equidistant)
                    minKnots = 2, # minimum number of knots
                    maxKnots = 20, # maximum number of knots
                    gaps=0, # min gaps between knots
                    interactionTerm="plot") # allow factor-smooth interaction

# Run SALSA 2D
salsa2D <- runSALSA2D(bestModel1D, # best 1D model
                      salsa2DList, # SALSA parameters from above
                      d2k=dMat$dataDist, # data-to-knot distance
                      k2k=dMat$knotDist) # knot-to-knot distance

# Store "best" model
bestModel2D <- salsa2D$bestModel
# Summary of results
# 10-16 singularities suggesting that CReSS is a poor approach here --> it cannot differentiate between various predictors
# lack of spread in the target variable is another obvious issue of course: not sure which is the underlying issue here
summary(bestModel2D)
# test significance
stats::anova(bestModel2D, test="F")
# partial plots of elevation on link scale
MRSea::runPartialPlots(bestModel2D, varlist.in="elev",
                       showKnots=TRUE, type="link", data=vegGeo)

# quite poor spatial coverage in general with knots and smooths here though
par(mfrow=c(2,2))
# Loop across all phases
for (treat in c("I", "II", "III", "IV")) {
  bWant <- vegGeo$plot %in% treat # a boolean specifying treatment
  fields::quilt.plot(vegGeo$long[bWant],
                     vegGeo$lat[bWant],
                     fitted(bestModel2D)[bWant],
                     nrow=24, ncol=60,
                     zlim=range(fitted(bestModel2D)),
                     main=paste0("Plot ", treat))
}

# broadly, it looks like higher-dimensional TPS are better
# while distance w/ elevation may be possible with CReSS, it's much pickier about the response variable formulation
# next step w/ GAM is to either make hierarchical and/or pool robust SEs under the GLMM framework

