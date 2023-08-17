#####################
# Script for generating virtual species
# Supplementary materials for case study 2 of 'Species distribution modelling with Bayes rule: a method'. 
# 6 April 2022, Robbert T. van den Dool, Centre for Crop Systems Analysis, Wageningen University 
#####################

# This script needs downloaded bioclimatic spatial data to work. We leave it up to the user to download these files, and change the workingdirectory as well as file names. 
# Files can be found at (last visited 14 April 2022): https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_elev.zip


# Load packages, set seed and set working directory
if(!require("easypackages")) install.packages("easypackages")
library(easypackages)
packages("dplyr", "giscoR", "sp", "sf", "ks", "geodata", "raster", "ggplot2", "rstudioapi",
         prompt = FALSE)
set.seed(123456789)

current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path)) 
rm(current_path)

packages("here")


# Get shapefiles of european countries using the giscoR package (NUTS regions)
# These will be used to delimit the study area and generate background points
countries <- gisco_get_nuts(epsg = "4326", nuts_level = "2")

# Remove outlying areas
countries <- countries %>% 
          filter(., !(NUTS_ID %in% c("ES70","FRY1","FRY2","FRY3","FRY4","FRY5","PT20", "PT30")))  %>% #exclude far away islands
          filter(., !(CNTR_CODE %in% c("TR", "IS", "CY"))) #exclude turkey, iceland and cyprus(kosovo and Bosnia-Herzegovina are not in the NUTS dataset)


# Sample 100,000 background points in regular grid
countries_sp <- as(countries, Class="Spatial")
backgroundpoints <- spsample(countries_sp, type="regular", n=100000)
backgroundpoints <- st_as_sf(backgroundpoints)


#get bioclimatic variables
env <- worldclim_global(var="bio", res=2.5, path=tempdir())
env <- stack(env)
env <- subset(env, c(1,12))

# Extract bioclimatic variable data at each backgroundpoint.
background <- cbind(backgroundpoints, raster::extract(env,backgroundpoints))

# Remove missing values 
background <- background %>% filter(!is.na(wc2.1_2.5m_bio_1) & !is.na(wc2.1_2.5m_bio_12)) 

colnames(background)[1:2] <- c("MAT", "MAP")


# Determine range of variables and fit a multivariate kernel density for the background distribution using 'kde' from package 'ks'. 
mins <- c(min(background$MAT), min(background$MAP))
maxs <- c(max(background$MAT), max(background$MAP))

backgrounddf <- st_drop_geometry(background)

kde_back <- kde(backgrounddf[,1:2], xmin=mins,xmax=maxs) 



# Define BayeSDM response functions, one for each virtual species
# The probability of a location or area being the site of a presence is defined as the density of variables at presence locations divided by the density of variables at all locations. 
# The presence distribution for temperature is defined as a normal distribution.
# The presence distribution for precipitation is defined as a log-normal distribution.
# The background distribution is defined as a multivariate kernel density 
#
# The difference between the two response functions is in the mean of the presence distribution for temperature, species 1: 12.25, species 2: 10.25

resp1 <- function(x1,x2){
  (dnorm(x1,mean=12.25, sd=1.6) * dlnorm(x2, meanlog=6.8, sdlog=0.3)) / dkde(fhat=kde_back,data.frame(x1=x1,x2=x2))
}

resp2 <- function(x1,x2){
  (dnorm(x1,mean=10.25, sd=1.6) * dlnorm(x2, meanlog=6.8, sdlog=0.3)) / dkde(fhat=kde_back,data.frame(x1=x1,x2=x2))
}


# Create scaled BayeSDM response functions (landscape distribution that sums to 1 over all locations)
resp1scaled <- function(x1,x2){
  response = resp1(x1,x2)
  response / sum(response)
}

resp2scaled <- function(x1,x2){
  response = resp2(x1,x2)
  response / sum(response)
}

# Calculate relative probabilities for all background locations for each species
background$spec1 <- resp1scaled(background$MAT, background$MAP)
background$spec2 <- resp2scaled(background$MAT, background$MAP)


# Sample 500 presence locations per species from the background
species1 <- background[sample(nrow(background), size=500, replace=T, prob=background$spec1),]
species2 <- background[sample(nrow(background), size=500, replace=T,prob=background$spec2),]

# Save the presences and background in a single sf file
species1$label = rep(1,500)
species2$label = rep(2,500)
background$label = rep(0, nrow(background))

virtualdata = rbind(species1, species2, background)

# Transform probabilities to percentiles for map plotting
probpercentiles <- function(data){
  brks = quantile(data, probs = seq(0, 1, 0.01), type=7)
  as.numeric(cut(data, breaks = brks, labels = 1:100, include.lowest = TRUE))
}

virtualdata$spec1perc = rep(NA, nrow(virtualdata))
virtualdata$spec2perc = rep(NA, nrow(virtualdata))

virtualdata[virtualdata$label==0,]$spec1perc <- probpercentiles(virtualdata[virtualdata$label==0,]$spec1)
virtualdata[virtualdata$label==0,]$spec2perc <- probpercentiles(virtualdata[virtualdata$label==0,]$spec2)


saveRDS(virtualdata, file=here("data - intermediate", "virtualdata.Rda"))



