#install.packages('tidyverse')
#install.packages('devtools')
#install.packages("usmap")
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

library(sf)
library(raster)
library(dplyr)
library(spData)
#Note that you will have to install spDataLarge after installing spData
#install.packages('spDataLarge', repos='https://nowosad.github.io/drat/', type='source')
library(spDataLarge)
library(tidyverse)
library(usmap)
library(ggplot2)
library(INLA)
library(rgdal)
library(sp)
library(tableone)

## Setting up directory path
setwd("C:/Users/rpaul9/Dropbox (UNC Charlotte)/CDC_APHA_IVP/aphaworkshop/APHA_Workshop_Underreprting")

#load the alcohol fatality data
alc.data = read.csv("Alcohol_2019_NC_Underreporting.csv", header = T)
## Preview of the variables
head(alc.data)
## Loading Neighborhood Adjacency Matrix
Mat1 <- read.table("WeightNC.txt")
WeightNC <- as.matrix(Mat1)
## Creating a Rural/Urban Binary Variable. Rural = Rural Proportion Above 50%
alc.data$Urban = as.numeric (alc.data$RuralProp < 0.50) 
alc.data$Rural = as.numeric (alc.data$RuralProp >= 0.50)

### Calculating Expected Counts
alc.data$EC = sum(alc.data$fatalcounts)/sum(alc.data$population2019)*alc.data$population2019
### Assigning an ID Variable, 1 to 100 for n = 100 counties in North Carolina
alc.data$ID <- 1:100


#not adjusting for underreporting
## Specifying Priors for Fixed Effects: Gaussian Prior
prior.fixed1 <- list(mean.intercept = 0, prec.intercept = 0.0001,
                     mean = c(0, 0), prec = c(0.0001, 0.0001))
## Specifying Priors for Precision Variables: Gamma Prior
prec.prior1 <- list(prec = list(param = c(0.001, 0.001)))
## Adjacency Wight Matric
NC.Mat <- as.matrix(WeightNC)
## Formula of the model with BYM Model
formula1 <- fatalcounts ~ AgeBelow18 + Female_Prop + f(ID, model = "bym", graph = NC.Mat)
## Fitting model using INLA
dat1.inla <- inla(formula1, family='nbinomial', E=EC,
                  data=alc.data,
                  control.family=list(link='log'),
                  control.fixed = prior.fixed1,
                  control.predictor=list(link = 1, compute = TRUE),
                  control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE,return.marginals.predictor=TRUE))
## Summary of Regression Coefficients
summary(dat1.inla)

## Saving Fixed Effects
write.csv(data.frame(Variables = c("Intercept", "AgeBelow18", "FemaleProp"), Mean = round(dat1.inla$summary.fixed[,1], 2), SD = round(dat1.inla$summary.fixed[,2], 2) ,CrILB = round(dat1.inla$summary.fixed[,3], 2), Median = round(dat1.inla$summary.fixed[,4], 2), CrIUB = round(dat1.inla$summary.fixed[,4], 2)), file = "Summary1.csv")

## Adjusting for same Reporting rates
## Assuming 32% Reporting Rate, we will multiply it with Expected Counts
alc.data$EC1 = 0.32*alc.data$EC

## Fitting INLA with known same reporting rate
dat2.inla <- inla(formula1, family='nbinomial', E=EC1,
                  data=alc.data,
                  control.family=list(link='log'),
                  control.fixed = prior.fixed1,
                  control.predictor=list(link=1, compute=TRUE),
                  control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE))

## Summary of Regression Coefficients
summary(dat2.inla)

## Saving Fixed Effects
write.csv(data.frame(Variables = c("Intercept", "AgeBelow18", "FemaleProp"), Mean = round(dat2.inla$summary.fixed[,1], 2), SD = round(dat2.inla$summary.fixed[,2], 2) ,CrILB = round(dat2.inla$summary.fixed[,3], 2), Median = round(dat2.inla$summary.fixed[,4], 2), CrIUB = round(dat2.inla$summary.fixed[,4], 2)), file = "Summary2.csv")

## Adjusting and Estimating for Variable Reporting Rates (Rural/Urban)
prior.fixed3 <- list(mean.intercept = 0, prec.intercept = 0.0001,
                    mean = c(-1.14, -1.14, 0, 0), prec = c(20, 20, 0.0001, 0.0001))

## Setting prior for precision: Gamma
prec.prior3 <- list(prec = list(param = c(0.001, 0.001)))


formula3 <- fatalcounts ~ Urban + Rural + AgeBelow18 + Female_Prop + f(ID, model = "bym", graph = NC.Mat)

dat3.inla <- inla(formula3, family='nbinomial', E=EC,
                 data=alc.data,
                 control.family=list(link='log'),
                 control.fixed = prior.fixed3,
                 control.predictor=list(link=1, compute=TRUE),
                 control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=TRUE))

## Summary of Regression Coefficients
summary(dat3.inla)

## Saving Fixed Effects
write.csv(data.frame(Variables = c("Intercept", "AgeBelow18", "FemaleProp"), Mean = round(dat3.inla$summary.fixed[c(1,4,5),1], 2), SD = round(dat3.inla$summary.fixed[c(1,4,5),2], 2) ,CrILB = round(dat3.inla$summary.fixed[c(1,4,5),3], 2), Median = round(dat3.inla$summary.fixed[c(1,4,5),4], 2), CrIUB = round(dat3.inla$summary.fixed[c(1,4,5),4], 2)), file = "Summary3.csv")

## Comparing Relative Rates


###### Change FIPS to fips #####
alc.data$fips <- alc.data$FIPS
## Obtaining Relative Rates Under Model with No Adjustments
Fitted_Values <- dat1.inla$summary.fitted.values
## First Column Represents Posterior Means include this in alc.data
alc.data$Post_MeanRR1 <- Fitted_Values[,1]


### Making Maps ####
plot_usmap(data = alc.data, values = "Post_MeanRR1", include = "NC", color = "black") + 
  scale_fill_continuous(low = "white", high = "red", name = "Relative Rates", label = scales::comma) + 
  labs(title = "Relative Rates", subtitle = " ") +
  theme(legend.position = "left")

## Obtaining Relative Rates Under Model with Same Known Reporting Rates
Fitted_Values <- dat2.inla$summary.fitted.values
## First Column Represents Posterior Means include this in alc.data
alc.data$Post_MeanRR2 <- Fitted_Values[,1]
### Making Maps ####
plot_usmap(data = alc.data, values = "Post_MeanRR2", include = "NC", color = "black") + 
  scale_fill_continuous(low = "white", high = "red", name = "Relative Rates", label = scales::comma) + 
  labs(title = "Relative Rates", subtitle = " ") +
  theme(legend.position = "left")

## Obtaining Relative Rates Under Model with unnown Reporting Rates
Fitted_Values <- dat3.inla$summary.fitted.values
## First Column Represents Posterior Means include this in alc.data
alc.data$Post_MeanRR3 <- Fitted_Values[,1]
### Making Maps ####
plot_usmap(data = alc.data, values = "Post_MeanRR3", include = "NC", color = "black") + 
  scale_fill_continuous(low = "white", high = "red", name = "Relative Rates", label = scales::comma) + 
  labs(title = "Relative Rates", subtitle = " ") +
  theme(legend.position = "left")

### Calculating and Plotting Exceedance Probabilities
exc1 <- sapply(dat1.inla$marginals.fitted.values,
               FUN = function(marg){1 - inla.pmarginal(q = 2, marginal = marg)})

alc.data$exc1 <- exc1

### Making Maps of Exceedance Probabilities####
plot_usmap(data = alc.data, values = "exc1", include = "NC", color = "black") + 
  scale_fill_continuous(low = "white", high = "red", name = "", label = scales::comma) + 
  labs(title = "Exceedance Probs.", subtitle = " ") +
  theme(legend.position = "left")

### Under Model 3
exc3 <- sapply(dat3.inla$marginals.fitted.values,
               FUN = function(marg){1 - inla.pmarginal(q = 2, marginal = marg)})

alc.data$exc3 <- exc3

### Making Maps of Exceedance Probabilities####
plot_usmap(data = alc.data, values = "exc3", include = "NC", color = "black") + 
  scale_fill_continuous(low = "white", high = "red", name = "", label = scales::comma) + 
  labs(title = "Exceedance Probs.", subtitle = " ") +
  theme(legend.position = "left")
### END OF WORKSHOP CODES #####
