#####################
# Script to reproduce the performance benchmark included in the paper
# 'Species distribution modelling using a generative framework'. 
# Robbert T. van den Dool, Alejandro Morales, Wopke van der Werf & J.C. (Bob) Douma
# 17 April 2023, Robbert T. van den Dool, Centre for Crop Systems Analysis, Wageningen University 
#####################


### Load functions and packages. 
source("Scripts/Functions.R")
if(!require("easypackages")) install.packages("easypackages")
library(easypackages)
packages("raster", "dismo", "maxnet", "naivebayes", "gridExtra", "copula", "optimParallel", "viridis", "cowplot", "ggplotify", "ks", "caret", "ggplot2", "sf", "sp", "SDMtune","PRROC", "Rmisc", "Metrics", prompt = FALSE)
set.seed(9999)


#Load data 
#116 presences #9 variables: bio1, bio5, bio6, bio7, bio8, bio12, bio16, bio17, biome(cat)
occ_file <- system.file('ex/bradypus.csv', package='dismo')
occ <- read.table(occ_file, header=TRUE, sep=',')[,-1]  

pred_files <- list.files(system.file('ex', package='dismo'), '\\.grd$', full.names=TRUE )
predictors <- stack(pred_files) #plot(predictors$bio1)

#set extent
extent <- extent(c(-94.8, -34.2, -56.05, 23.55)) #from the 2006 phillips et al. paper 
predictors2 <- crop(predictors, extent) #plot(predictors2$bio1)
predictors3 <- subset(predictors2, subset=c("bio1", "bio7", "bio12", "bio16", "bio5", "bio6", "bio8")) #bio17 not used

#generate psuedo-absences
studyarea <- rasterToPolygons(predictors3$bio1, n=4, na.rm=TRUE, digits=12)
studyarea <- aggregate(studyarea, dissolve=TRUE) #plot(studyarea)

background <- spsample(studyarea, type="regular", n=25000) #creates 100,000 PA points in regular grid #plot(background)

occ_sf <- st_as_sf(x = occ, 
                   coords = c("lon", "lat"),
                   crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
studyarea_sf <- st_as_sf(studyarea)
background_sf <- st_as_sf(background)

#ggplot()+ 
#  geom_sf(data=studyarea_sf)+
#  geom_sf(data=occ_sf, col="blue", cex=1.0)+
#  geom_sf(data=background_sf, col='red', cex=0.5)

#Extract values 
background_df <- as.data.frame(background)
names(background_df)[1] <- "x"
names(background_df)[2] <- "y"
names(occ)[1] <- "x"
names(occ)[2] <- "y"

predictors4 = as(predictors3, "SpatRaster") #since raster is no longer supported by package SDMtune.

data <- prepareSWD(species = "Bradypus variegatus", p = occ, a = background_df,
                   env = predictors4) #data@pa


#Prepare CV data
cvdata = lapply(1:50, function(x)SDMtune::trainValTest(data, test=0.1, val=0.2, only_presence = T, seed=x))



######
###### Predictions from optimal models for each partition:  
######

## MaxEnt 
##

time = Sys.time() #Tries out all possible settings for tuning, on the tuneset. Returns best settings.
optimalsMaxEnt2 <- findoptimalMaxent(data=cvdata, features=c("l", "lq", "lqh", "lqph"), regularization=seq(0.2, 4, 0.4))
Sys.time() - time #runtime 2.3 hours

time = Sys.time() #Creates a model using best settings, tests on test set. Returns predictions.
MaxEntTestpreds2 <- fitoptimalMaxentEach(data=cvdata, settings=optimalsMaxEnt2) #data, settings, seed, valval, valtest
Sys.time() - time

optimalsMaxEnt2_metrics = calclistmetrics2(MaxEntTestpreds2, seeds=1:50) #Returns metrics.
myCI(optimalsMaxEnt2_metrics$nloglike) 
#upper     mean    lower 
#111.4115 110.5137 109.6159 

myCI(optimalsMaxEnt2_metrics$auc)
#upper      mean     lower 
#0.8466236 0.8347152 0.8228069 




### Create 4 variables data
predictors4v <- subset(predictors2, subset=c("bio1", "bio7", "bio12", "bio16")) 

#generate psuedo-absences
studyarea <- rasterToPolygons(predictors4v$bio1, n=4, na.rm=TRUE, digits=12)
studyarea <- aggregate(studyarea, dissolve=TRUE) #plot(studyarea)

background <- spsample(studyarea, type="regular", n=25000) #creates 100,000 PA points in regular grid #plot(background)

occ_sf <- st_as_sf(x = occ, 
                   coords = c("lon", "lat"),
                   crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

studyarea_sf <- st_as_sf(studyarea)
background_sf <- st_as_sf(background)

#Extract values 
background_df <- as.data.frame(background)
names(background_df)[1] <- "x"
names(background_df)[2] <- "y"
names(occ)[1] <- "x"
names(occ)[2] <- "y"

predictors4v = as(predictors4v, "SpatRaster") #since raster is no longer supported by package SDMtune.

data4v <- prepareSWD(species = "Bradypus variegatus", p = occ, a = background_df,
                   env = predictors4v) #data@pa

#make cv data 
cvdata4v = lapply(1:50, function(x)SDMtune::trainValTest(data4v, test=0.1, val=0.2, only_presence = T, seed=x))


## MaxEnt 4v: bio1, bio7, bio12, bio16
##

time = Sys.time()
optimalsMaxEnt4v <- findoptimalMaxent(data=cvdata4v, features=c("l", "lq", "lqh", "lqph"), regularization=seq(0.2, 4, 0.4))
Sys.time() - time #runtime 1.5 hours

time = Sys.time()
MaxEntTestpreds4v <- fitoptimalMaxentEach(data=cvdata4v, settings=optimalsMaxEnt4v) #data, settings, seed, valval, valtest
Sys.time() - time

MaxEntTestpreds4v_metrics = calclistmetrics2(MaxEntTestpreds4v, seeds=1:50)
myCI(MaxEntTestpreds4v_metrics$nloglike) 
#upper     mean    lower 
#112.5799 111.5891 110.5983 

myCI(MaxEntTestpreds4v_metrics$auc)
#upper     mean    lower 
#0.8482350 0.8367286 0.8252221 



## Gen-4v-KDE: generative model using multivariate kernel density estimation
##

mins <- c(min(data@data[data@pa==0,1])-1,min(data@data[data@pa==0,2])-1,min(data@data[data@pa==0,3])-1,min(data@data[data@pa==0,4])-1,min(data@data[data@pa==0,5])-1, min(data@data[data@pa==0,6])-1)
maxs <- c(max(data@data[data@pa==0,1])+1,max(data@data[data@pa==0,2])+1,max(data@data[data@pa==0,3])+1,max(data@data[data@pa==0,4])+1,max(data@data[data@pa==0,5])+1, max(data@data[data@pa==0,6])+1)

outputlist <- vector(mode="list", length=50)

start_time <- Sys.time()
for (i in 1:50){
  svMisc::progress(i)
  datasets <- cvdata[[i]]
  train <- datasets[[1]]
  test <- datasets[[3]] 
  
  trainback <- as.matrix(train@data[train@pa==0,1:4])
  trainpres <- as.matrix(train@data[train@pa==1,1:4])
  testdata <- as.matrix(test@data[,1:4])
  
  BHpi3 <- Hpi(x = trainback)*3
  PHpi3 <- Hpi(x = trainpres)*3
  
  kde_back <- kde(trainback,xmin=mins[1:4],xmax=maxs[1:4], binned=T, H=BHpi3, bgridsize=rep(10,4)) #,xmin=mins[1:2],xmax=maxs[1:2] #rep(21,4)) #H=B_Hlscv ,  H=BHpi5
  dkde_back <- dkde(fhat=kde_back,testdata)
  
  kde_pres <- kde(trainpres,xmin=mins[1:4],xmax=maxs[1:4], binned=T, H=PHpi3, bgridsize=rep(10,4)) #,xmin=mins[1:2],xmax=maxs[1:2] #, H=PV_Hlscv #, H=PV_Hlscv , H=PHpi5
  dkde_pres <- dkde(fhat=kde_pres,testdata) #sum(dkde_pres)
  dkde_pres[dkde_pres < 0] <- 0.000000000000000000000000000000001 #12598 below 0... 10415 without minmax
  
  probpres <- (dkde_pres)/(dkde_back)
  probpres <- data.frame(prob=probpres, pa=test@pa)
  #probpres$nprob <- probpres$prob/sum(probpres$prob)
  
  outputlist[[i]] <- probpres
}
end_time <- Sys.time()
end_time - start_time #51 minutes


BayeSDMTestpreds4v_metrics = calclistmetrics2(outputlist, seeds=1:50)
myCI(BayeSDMTestpreds4v_metrics$nloglike) 
#upper     mean    lower 
#114.8297 112.8731 110.9166 

myCI(BayeSDMTestpreds4v_metrics$auc)
#upper     mean    lower 
#0.8459273 0.8329547 0.8199821 


## Gen-I: generative model assuming Independence, using Kernel Density Estimation. 
## 

start_time <- Sys.time()
optimalsmethod3b_optim <- findoptimalmethod3b_optim2(data=cvdata, par=rep(log(5),14), maxit=4000)
Sys.time() - start_time #1.57 hours

optimalsmethod3b_optim2 <- round(optimalsmethod3b_optim, digits=2)
names(optimalsmethod3b_optim2) = c("seed","bio1_pres", "bio1_back", "bio7_pres", "bio7_back", "bio12_pres", "bio12_back", "bio16_pres", "bio16_back", "bio5_pres", "bio5_back", "bio6_pres", "bio6_back", "bio8_pres", "bio8_back", "NLLP_tune", "convergence")

#fixes very low adjust fits... 
for (i in (2:(ncol(optimalsmethod3b_optim2)-2))){
  optimalsmethod3b_optim2[which(optimalsmethod3b_optim2[,i] < 0.1),i] <- 0.1
}

optimalsmethod3b_optim2list = vector(mode="list", length=7)

#creates list 
for (i in 1:ncol(cvdata[[1]][[1]]@data)){
  optimalsmethod3b_optim2list[[i]] <- data.frame(seed=optimalsmethod3b_optim2$seed, adjusts0 = optimalsmethod3b_optim2[,(2*i+1)] , adjusts1 = optimalsmethod3b_optim2[,(2*i)] )  
} 

fitmethod3b_optim  <- fitoptimalmodel.method1(data=cvdata, settings=optimalsmethod3b_optim2list)

BayeSDMTestpreds3b_metrics = calclistmetrics2(fitmethod3b_optim, seeds=1:50)
myCI(BayeSDMTestpreds3b_metrics$nloglike) 
#upper     mean    lower 
#118.3975 115.9732 113.5489

myCI(BayeSDMTestpreds3b_metrics$auc)
#upper     mean    lower 
#0.8255910 0.8102962 0.7950015 


## Gen-GC: generative model using Gaussian Copulas, using Kernel Density Estimation for the marginals. 
## 

cl <- makeCluster(8)     # set the number of processor cores
setDefaultCluster(cl=cl) # set 'cl' as default cluster

clusterExport(cl, "updateMBvar") 
clusterExport(cl, "predict.nonparametric_naive_bayes2")
clusterExport(cl, "getunifPKDE2")
clusterExport(cl, "approxPKDE2")
clusterExport(cl, "dCopula")
clusterExport(cl, "normalCopula")

startt = Sys.time()
optimals_3bgauscopParallel = findoptimalmethod3bgauscop_optimparallel(data=cvdata, par=rep(log(1),14), seeds=1:50, maxit=2000)
Sys.time()-startt #2.2 hours 

optimals_3bgauscopParallel_list <- vector(mode = "list", length = 7)

for (i in 1:ncol(cvdata[[1]][[1]]@data)){
  optimals_3bgauscopParallel_list[[i]] <- data.frame(seed=optimals_3bgauscopParallel$V1, adjusts1 = optimals_3bgauscopParallel[,(2*i)], adjusts0 = optimals_3bgauscopParallel[,(2*i+1)] )  
}

start_time <- Sys.time()
gauscopParallel_testpreds <- fitoptimalmodel.mbgauscop_correct(data=cvdata, settings=optimals_3bgauscopParallel_list)
end_time <- Sys.time() 
end_time-start_time

gauscopParallel_testpreds2 = gauscopParallel_testpreds[-48]
#failed for 1 run. 48

gauscopParallel_testpreds_metrics = calclistmetrics2(gauscopParallel_testpreds2, seeds=c(1:47,49,50))
myCI(gauscopParallel_testpreds_metrics$nloglike) 
#upper     mean    lower 
#156.4684 148.7712 141.0740 
#156.4684 148.7712 141.0740 
#127.3207 122.8878 118.4549

myCI(gauscopParallel_testpreds_metrics$auc)
#0.8354274 0.8199814 0.8045354 
stopCluster(cl)



# #Nelder-mead, alternative fitting method
# startt = Sys.time()
# optimals_3bgauscop = findoptimalmethod3bgauscop_optim(data=cvdata, par=rep(log(1),14), seeds=1:50, maxit=2000)
# Sys.time()-startt #5h for whole thing. 
# 
# optimals_3bgauscop_list <- vector(mode = "list", length = 7)
# 
# for (i in 1:ncol(cvdata[[1]][[1]]@data)){
#   optimals_3bgauscop_list[[i]] <- data.frame(seed=optimals_3bgauscop$V1, adjusts1 = optimals_3bgauscop[,(2*i)], adjusts0 = optimals_3bgauscop[,(2*i+1)] )  
# }
# 
# start_time <- Sys.time()
# gauscop_testpreds <- fitoptimalmodel.mbgauscop_correct(data=cvdata, settings=optimals_3bgauscop_list)
# end_time <- Sys.time() #49 seconds for 50, 46.5 seconds for 45
# end_time-start_time
# 
# gauscop_testpreds_metrics = calclistmetrics2(gauscop_testpreds, seeds=1:50) 
# myCI(gauscop_testpreds_metrics$nloglike) 
# #upper     mean    lower 
# #156.4684 148.7712 141.0740 
# 
# myCI(gauscop_testpreds_metrics$auc)
# #0.8154276 0.8011781 0.7869285 
# #END Nelder Mead 






## Figure 3
##

modellist <- list(optimalsMaxEnt2_metrics, #basic maxent
                  MaxEntTestpreds4v_metrics, #4v maxent
                  BayeSDMTestpreds3b_metrics, #3b
                  gauscopParallel_testpreds_metrics, #3b + gauscop
                  BayeSDMTestpreds4v_metrics #4v kde
                  )


plotdataframe <- data.frame(index=c(1,2,3,4,5),
                            title = c("Maxent","Maxent-4v","Gen-I", "Gen-GC", "Gen-KDE-4v"),
                            nll = rep(NA,5),
                            nll_u = rep(NA,5),
                            nll_l = rep(NA,5),
                            psum = rep(NA,5),
                            psum_u = rep(NA,5),
                            psum_l = rep(NA,5),
                            auc = rep(NA,5),
                            auc_u = rep(NA,5),
                            auc_l = rep(NA,5),
                            pmean = rep(NA,5),
                            pmean_u = rep(NA,5),
                            pmean_l = rep(NA,5),
                            psd = rep(NA,5),
                            psd_u = rep(NA,5),
                            psd_l = rep(NA,5))


plotdataframe[,3:17] <- t(sapply(1:5, function(x){
  data = sapply(c(6, 7, 8, 2, 3), function(y){Rmisc::CI(modellist[[x]][,y], ci=0.95)})
  output = rep(NA, 15)
  output[seq(1,15,3)] = data[2,]
  output[seq(2,15,3)] = data[1,]
  output[seq(3,15,3)] = data[3,]
  return(output)
}))


#plot 1: NLL with labels
plot1 <- ggplot(plotdataframe, aes(y = index, x = nll, color=as.factor(index))) +
  geom_point(shape = 16, size = 4) +  
  geom_errorbarh(aes(xmin = nll_l, xmax = nll_u), height = 0.25) +
  scale_y_continuous(name = "", breaks=1:5, labels = plotdataframe$title, trans = "reverse") +
  scale_colour_manual("legend", values = c("1" = "black", "2" = brewer.pal(9,"Set1")[1], "3" = brewer.pal(9,"Set1")[3], "4" = brewer.pal(9,"Set1")[4], "5"=brewer.pal(9,"Set1")[2]))+
  xlab("NLLp") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x.bottom = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        legend.position = "none") +
  geom_vline(xintercept = 109.6159, linetype="dotted", 
               color = "black", linewidth=0.5) + 
  geom_vline(xintercept = 111.4115, linetype="dotted", 
             color = "black", linewidth=0.5)
plot1

#plot 2: psum
plot2 <- ggplot(plotdataframe, aes(y = index, x = psum, color=as.factor(index))) +
  geom_point(shape = 16, size = 4) +  
  geom_errorbarh(aes(xmin = psum_l, xmax = psum_u), height = 0.25) +
  scale_y_continuous(name = "", breaks=1:5, labels = c("","","","", ""), trans = "reverse") +
  scale_colour_manual("legend", values = c("1" = "black", "2" = brewer.pal(9,"Set1")[1], "3" = brewer.pal(9,"Set1")[3], "4" = brewer.pal(9,"Set1")[4], "5"=brewer.pal(9,"Set1")[2]))+
  xlab("Psum") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x.bottom = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.ticks.y = element_blank(),
        axis.line.y = element_line(color = "grey", linetype = "solid"),
        legend.position = "none") + 
  scale_x_continuous("Psum", breaks=c(0.0015, 0.0030,0.0045), labels=c(0.0015, 0.0030, 0.0045), limits=c(0.0015, 0.0050)) + 
  geom_vline(xintercept = 0.002018109, linetype="dotted", 
             color = "black", linewidth=0.5) + 
  geom_vline(xintercept = 0.001722162, linetype="dotted", 
             color = "black", linewidth=0.5)
plot2

#plot 3: auc
plot3 <- ggplot(plotdataframe, aes(y = index, x = auc, color=as.factor(index))) +
  geom_point(shape = 16, size = 4) +  
  geom_errorbarh(aes(xmin = auc_l, xmax = auc_u), height = 0.25) +
  scale_y_continuous(name = "", breaks=1:5, labels = c("","","", "", ""), trans = "reverse") +
  scale_colour_manual("legend", values = c("1" = "black", "2" = brewer.pal(9,"Set1")[1], "3" = brewer.pal(9,"Set1")[3], "4" = brewer.pal(9,"Set1")[4], "5"=brewer.pal(9,"Set1")[2]))+
  xlab("AUC") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x.bottom = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_line(color = "grey", linetype = "solid"),
        legend.position = "none") +
  geom_vline(xintercept = 0.8362268, linetype="dotted", 
             color = "black", linewidth=0.5) + 
  geom_vline(xintercept = 0.8589903, linetype="dotted", 
             color = "black", linewidth=0.5)
plot3

#plot 4: pmean
plot4 <- ggplot(plotdataframe, aes(y = index, x = pmean, color=as.factor(index))) +
  geom_point(shape = 16, size = 4) +  
  geom_errorbarh(aes(xmin = pmean_l, xmax = pmean_u), height = 0.25) +
  scale_y_continuous(name = "", breaks=1:5, labels = c("","","","", ""), trans = "reverse") +
  scale_colour_manual("legend", values = c("1" = "black", "2" = brewer.pal(9,"Set1")[1], "3" = brewer.pal(9,"Set1")[3], "4" = brewer.pal(9,"Set1")[4], "5"=brewer.pal(9,"Set1")[2]))+
  xlab("Pmean") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x.bottom = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_line(color = "grey", linetype = "solid"), #element_blank()
        legend.position = "none")+
  geom_vline(xintercept = 0.0001435135, linetype="dotted", 
           color = "black", linewidth=0.5) + 
  geom_vline(xintercept = 0.0001681758, linetype="dotted", 
             color = "black", linewidth=0.5)
plot4



#plot 5: psd
plot5 <- ggplot(plotdataframe, aes(y = index, x = psd, color=as.factor(index))) +
  geom_point(shape = 16, size = 4) +  
  geom_errorbarh(aes(xmin = psd_l, xmax = psd_u), height = 0.25) +
  scale_y_continuous(name = "", breaks=1:5, labels = c("","","","", ""), trans = "reverse") +
  scale_colour_manual("legend", values = c("1" = "black", "2" = brewer.pal(9,"Set1")[1], "3" = brewer.pal(9,"Set1")[3], "4" = brewer.pal(9,"Set1")[4], "5"=brewer.pal(9,"Set1")[2]))+
  xlab("Psd") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x.bottom = element_text(size = 14, colour = "black"), #axis.text.x = element_text(angle=90, hjust=1)
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_line(color = "grey", linetype = "solid"),
        legend.position = "none")+
  geom_vline(xintercept = 0.0001839848, linetype="dotted", 
           color = "black", linewidth=0.5) + 
  geom_vline(xintercept = 0.0001305452, linetype="dotted", 
             color = "black", linewidth=0.5)
plot5

#combined plot 
lay <-  matrix(c(1,1,1,2,2,3,3,4,4,5,5), nrow = 1)

tiff("Figure3.tiff", width = 14, height = 4, units = 'in', res = 800, compression="lzw")
grid.arrange(plot1, plot3, plot2,plot4,plot5, layout_matrix = lay)
dev.off()







## Gen-KDE-4v final model 
##

data0 <- data 
data0@data <- data0@data[data@pa==0,]
data0@pa <- data0@pa[data@pa==0]
data0@coords <- data0@coords[data@pa==0,]
datatest <- cbind(data0@data, data0@coords)

datatest_sf <- st_as_sf(
  datatest,
  coords = c("X", "Y"), crs = 4326
)

data1 <- data 
data1@data <- data1@data[data@pa==1,]
data1@pa <- data1@pa[data@pa==1]
data1@coords <- data1@coords[data@pa==1,]

datapres <- cbind(data1@data, data1@coords)

datapres_sf <- st_as_sf(
  datapres,
  coords = c("X", "Y"), crs = 4326
)

#final model 
trainback <- as.matrix(data@data[data@pa==0,1:4])
trainpres <- as.matrix(data@data[data@pa==1,1:4])
  
BHpi3 <- Hpi(x = trainback)*3
PHpi3 <- Hpi(x = trainpres)*3
  
kde_back <- kde(trainback,xmin=mins[1:4],xmax=maxs[1:4], binned=T, H=BHpi3, bgridsize=rep(10,4)) #,xmin=mins[1:2],xmax=maxs[1:2] #rep(21,4)) #H=B_Hlscv ,  H=BHpi5
dkde_back <- dkde(fhat=kde_back,trainback)
  
kde_pres <- kde(trainpres,xmin=mins[1:4],xmax=maxs[1:4], binned=T, H=PHpi3, bgridsize=rep(10,4)) #,xmin=mins[1:2],xmax=maxs[1:2] #, H=PV_Hlscv #, H=PV_Hlscv , H=PHpi5
dkde_pres <- dkde(fhat=kde_pres,trainback) #sum(dkde_pres)
dkde_pres[dkde_pres < 0] <- 0.000000000000000000000000000000001 #12598 below 0... 10415 without minmax
  
probpres <- (dkde_pres)/(dkde_back)
probpres <- data.frame(prob=probpres)
probpres$nprob <- probpres$prob/sum(probpres$prob)
  


## Figure 5
##

pressample = rkde(10000, kde_pres)
backsample = rkde(10000, kde_back)

#function to generate density function plots
plotkde = function(variable, text, first=F, legend=F){
  
  pdata = pressample[,variable]
  bdata = backsample[,variable]
  
  presdens = density(pdata, from= min(bdata), to = max(bdata))
  backdens = density(bdata, from= min(bdata), to = max(bdata))
  
  ymax = max (c(presdens$y, backdens$y))*1.2
  xlim = c(min(backdens$x), max(backdens$x))
  
  plotfunct = function(){
    
  plot(backdens, lwd=2, ylim=c(0, ymax), xlab=text, xlim=xlim, ylab="", main="", col=brewer.pal(9,"Set1")[1], yaxt="n", cex.lab=1.2, axes=F,frame.plot=TRUE)
  lines(presdens, lwd=2, lty=2, xlim=xlim, col=brewer.pal(9,"Set1")[3])
  Axis(side=1, labels=T)
  if(first) {
    title(ylab = "Density", cex.lab=1.2,line = 0) 
  }
  
  if(legend) {
    legend("topleft", legend=c("background", "presence"), lwd=2, col=c("red", "green"), cex=1.2, lty=1:2, bty = "n")
  }
  }
  return(plotfunct)
}

# function to generate Response plots
plotkderesponse = function(variable, text, first=F){
  
  pdata = pressample[,variable]
  bdata = backsample[,variable]
  
  presdens = density(pdata, from= min(bdata), to = max(bdata))
  backdens = density(bdata, from= min(bdata), to = max(bdata))
  
  responsefunct = function(x){
    pdens = approx(presdens$x,presdens$y,xout=x,rule = 2, ties = "ordered")$y
    bdens = approx(backdens$x,backdens$y,xout=x,rule = 2, ties = "ordered")$y
    
    pdens[pdens < 0] <- 0.000000000000000000000000000000001
    bdens[bdens < 0] <- 0.000000000000000000000000000000001
    
    prob = pdens/bdens
    return(prob)
  }
  
  
  ymax = max(responsefunct(seq(min(bdata), max(bdata), length.out=512)))*1.2
  xlim = c(min(bdata), max(bdata))
  
  plotfunct = function(){
    
    curve(responsefunct(x), lwd=2, ylim=c(0, ymax), xlab=text, xlim=xlim, ylab="", main="", col=brewer.pal(9,"Set1")[5], cex.lab=1.2) 
    abline(h=1,lty=2)
    #Axis(side=1, labels=T)
    if(first) {
      title(ylab = "Relative probability", line = 0, cex.lab=1.2)
    }
    
  }
  return(plotfunct)
}

plot1 = plotkde("bio1","Mean temperature in d°C", first=T, legend=T)
plot2 = plotkde("bio7","Temperature annual range in d°C")
plot3 = plotkde("bio12","Annual precipitation in mm")
plot4 = plotkde("bio16","Precipitation wettest quarter in mm")

plot1 = as.grob(plot1)
plot2 = as.grob(plot2)
plot3 = as.grob(plot3)
plot4 = as.grob(plot4)

plot5 = plotkderesponse("bio1","Mean temperature in d°C", first=T)
plot6 = plotkderesponse("bio7","Temperature annual range in d°C") #d°C
plot7 = plotkderesponse("bio12","Annual precipitation in mm")
plot8 = plotkderesponse("bio16","Precipitation wettest quarter in mm")

plot5 = as.grob(plot5)
plot6 = as.grob(plot6)
plot7 = as.grob(plot7)
plot8 = as.grob(plot8)

tiff("Figure5.tiff", width = 20, height = 8.25, units = 'in', res = 1200, compression="lzw")
plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, labels = c('(a) bio1', '(b) bio7', '(c) bio12', '(d) bio16', '(e) bio1', '(f) bio7', '(g) bio12', '(h) bio16'), nrow=2, label_size = 14, align = "hv") 
dev.off()









####Figure 4
#### KDE-4v
brks.kde4v <- quantile(probpres$nprob, probs = seq(0, 1, 0.01), type=7)
datatest_sf$kde4vprob4 <- as.numeric(cut(probpres$nprob, breaks = brks.kde4v, labels = 1:100, include.lowest = TRUE))


#### Gen-I 
method3b_testpreds_avpars <- getavpar(optimalsmethod3b_optim2) #average model settings

mb3settings <- as.numeric(method3b_testpreds_avpars[2:15])
map.mb3 <- predictMaxBayes.3bgauscop4(data, settings1=mb3settings[c(1,3,5,7,9,11,13)], settings0=mb3settings[c(2,4,6,8,10,12,14)], indep=T)
map.mb3 <- map.mb3$prob[map.mb3$pa==0] #
map.mb3 <- map.mb3/sum(map.mb3)

brks.mb3 <- quantile(map.mb3, probs = seq(0, 1, 0.01), type=7)
datatest_sf$mb3prob4 <- as.numeric(cut(map.mb3, breaks = brks.mb3, labels = 1:100, include.lowest = TRUE))


#### Gen-GC
optimals_3bgauscopParallel2 = optimals_3bgauscopParallel[-48,]
gauscop_testpreds_avpars <- getavpar(optimals_3bgauscopParallel2)
 
mb3gcsettings <- as.numeric(gauscop_testpreds_avpars[2:15])
map.mb3gc <- predictMaxBayes.3bgauscop4(data, settings1=mb3gcsettings[c(1,3,5,7,9,11,13)], settings0=mb3gcsettings[c(2,4,6,8,10,12,14)], indep=F)
map.mb3gc <- map.mb3gc$prob[map.mb3gc$pa==0] #sum(map.mb3gc)
map.mb3gc <- map.mb3gc/sum(map.mb3gc) ##sum(map.mb3gc)
 



brks.mb3gc <- quantile(map.mb3gc, probs = seq(0, 1, 0.01), type=7)
datatest_sf$mb3gcprob4 <- as.numeric(cut(map.mb3gc, breaks = brks.mb3gc, labels = 1:100, include.lowest = TRUE))


#### Maxent 

#average model settings
table(optimalsMaxEnt2$features) #l 4 lq 19 lqh 18 lqph 9  ###lq 19
mean(optimalsMaxEnt2$regularization) #0.824

optmaxentmodel <- SDMtune::train(method = "Maxnet", data = data, fc="lq", reg=0.824)
map.maxent <- predict(optmaxentmodel, data = data0, type = "exponential", progress="text")

brks.maxent <- quantile(map.maxent, probs = seq(0, 1, 0.01), type=7)
datatest_sf$maxentprob4 <- as.numeric(cut(map.maxent, breaks = brks.maxent, labels = 1:100, include.lowest = TRUE))



## Figure 4
##

#Percentiles maps
#4 maps, last one with legend: MaxEnt, Gen-I, Gen-GC, Gen-KDE-4v


###Create raster formats for plots. 
#maxent
raster_1 = datatest_sf[,"maxentprob4"]
raster_1 <- as(raster_1, "Spatial")
raster_1 = rasterFromXYZ(raster_1, res=c(NA,NA), crs="EPSG:4326", digits=10) #21km2
crs(raster_1) = CRS(SRS_string = "EPSG:4326") #3857

raster_1 = as(raster_1, "SpatialPixelsDataFrame")
raster_1 <- as.data.frame(raster_1)
names(raster_1)[1] = "rast1"

#Gen-I
raster_2 = datatest_sf[,"mb3prob4"]
raster_2 <- as(raster_2, "Spatial")
raster_2 = rasterFromXYZ(raster_2, res=c(NA,NA), crs="EPSG:4326", digits=10) #21km2
crs(raster_2) = CRS(SRS_string = "EPSG:4326") #3857

raster_2 = as(raster_2, "SpatialPixelsDataFrame")
raster_2 <- as.data.frame(raster_2)
names(raster_2)[1] = "rast2"

#Gen-GC
raster_3 = datatest_sf[,"mb3gcprob4"]
raster_3 <- as(raster_3, "Spatial")
raster_3 = rasterFromXYZ(raster_3, res=c(NA,NA), crs="EPSG:4326", digits=10) #21km2
crs(raster_3) = CRS(SRS_string = "EPSG:4326") #3857

raster_3 = as(raster_3, "SpatialPixelsDataFrame")
raster_3 <- as.data.frame(raster_3)
names(raster_3)[1] = "rast3"

#Gen-KDE-4v
raster_4 = datatest_sf[,"kde4vprob4"]
raster_4 <- as(raster_4, "Spatial")
raster_4 = rasterFromXYZ(raster_4, res=c(NA,NA), crs="EPSG:4326", digits=10) #21km2
crs(raster_4) = CRS(SRS_string = "EPSG:4326") #3857

raster_4 = as(raster_4, "SpatialPixelsDataFrame")
raster_4 <- as.data.frame(raster_4)
names(raster_4)[1] = "rast4"

plot0 = ggplot()+
  geom_tile(data=raster_1,aes(x = x, y = y, fill = rast1), lwd = 0)+
  #geom_bar( fill="#FF9999", colour="black")
  geom_sf(data=datapres_sf, colour = "black", alpha = 0.5, shape=3, size=0.5)+
  scale_fill_viridis(name="Percentile", discrete=F, option="viridis", begin=0, end=1, direction=-1)+ #, , breaks=c(0,seq(10,100, by=10)), labels=formatC(brks[c(1, 11,21,31,41,51,61,71,81,91,101)]format = "e", digits = 2
  #scale_colour_gradient2(name="Percentile",low = "blue", mid="green", midpoint=50,high = "red")+
  theme(text = element_text(size=7),axis.title.x = element_blank(),axis.title.y = element_blank(),plot.margin=unit(c(0,0.1,0,0.1), "cm"), legend.margin=margin(t = 0, unit='cm'), legend.key.height = unit(1.2, 'cm'),) #t, r, b, l #legend.margin=margin(t = 0, unit='cm'), legend.key.height = unit(1.5, 'cm'),

plot1 = ggplot()+
  geom_tile(data=raster_1,aes(x = x, y = y, fill = rast1), lwd = 0)+
  #geom_bar( fill="#FF9999", colour="black")
  geom_sf(data=datapres_sf, colour = "black", alpha = 0.5, shape=3, size=0.5)+
  scale_fill_viridis(name="Rel. prob.", discrete=F, option="viridis", begin=0, end=1, direction=-1)+
  theme(text = element_text(size=7), axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position="none", plot.margin=unit(c(0,0,0,0), "cm"))+
  annotate("path",
             x=-78+5*cos(seq(0,2*pi,length.out=100)),
             y=4+5*sin(seq(0,2*pi,length.out=100)), col=brewer.pal(9,"Set1")[1], lwd=1) + 
  annotate("path",
           x=-75+5*cos(seq(0,2*pi,length.out=100)),
           y=-46+5*sin(seq(0,2*pi,length.out=100)), col=brewer.pal(9,"Set1")[1], lwd=1) + 
  annotate(geom="text", x=-85, y=4, label="N",
           color="red") +
  annotate(geom="text", x=-82, y=-46, label="S",
           color="red")
  
plot2 =  ggplot()+
  geom_tile(data=raster_2,aes(x = x, y = y, fill = rast2), lwd = 0)+
  #geom_bar( fill="#FF9999", colour="black")
  geom_sf(data=datapres_sf, colour = "black", alpha = 0.5, shape=3, size=0.5)+
  scale_fill_viridis(name="Rel. prob.", discrete=F, option="viridis", begin=0, end=1, direction=-1)+
  theme(text = element_text(size=7), axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position="none", axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin=unit(c(0,0,0,0), "cm"))

plot3 =  ggplot()+
  geom_tile(data=raster_3,aes(x = x, y = y, fill = rast3), lwd = 0)+
  #geom_bar( fill="#FF9999", colour="black")
  geom_sf(data=datapres_sf, colour = "black", alpha = 0.5, shape=3, size=0.5)+
  scale_fill_viridis(name="Rel. prob.", discrete=F, option="viridis", begin=0, end=1, direction=-1)+
  theme(text = element_text(size=7), axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position="none", axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin=unit(c(0,0,0,0), "cm"))

plot4 = ggplot()+
  geom_tile(data=raster_4,aes(x = x, y = y, fill = rast4), lwd = 0)+
  #geom_bar( fill="#FF9999", colour="black")
  geom_sf(data=datapres_sf, colour = "black", alpha = 0.5, shape=3, size=0.5)+
  scale_fill_viridis(name="Rel. prob.", discrete=F, option="viridis", begin=0, end=1, direction=-1)+
  theme(text = element_text(size=7),axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position="none", axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin=unit(c(0,0,0,0), "cm")) #0.1,0,0.1,0

lay <-  matrix(c(1,1,1,2,2,2,3,3,3,4,4,4), nrow = 1)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(plot0)

#fix widths
plots <- list(plot1, plot2, plot3, plot4)
grobs <- list()
widths <- list()

for (i in 1:length(plots)){
  grobs[[i]] <- ggplotGrob(plots[[i]])
  widths[[i]] <- grobs[[i]]$widths[2:5]
}

maxwidth <- do.call(grid::unit.pmax, widths)

for (i in 1:length(grobs)){
  grobs[[i]]$widths[2:5] <- as.list(maxwidth)
}


tiff("Figure4.tiff", width = 11, height = 3.5, units = 'in', compression = "lzw", res=800) #width = 12, height = 4,  # res = 800
plot_grid(grobs[[1]], grobs[[2]],grobs[[3]],grobs[[4]],mylegend, label_size = 8, labels = c('(a)', '(b)', '(c)', '(d)', ''), label_x = 0.07, label_y = 0.985, align = "h", nrow = 1, axis="b", rel_widths = c(4/17,4/17,4/17,4/17, 1/17), rel_heights = c(1/5,1/5,1/5,1/5,1/5)) #hadjust=0.1
dev.off()




###############
###############    Figure 6
###############

#Define south area
box_south = c(xmin = -77, ymin = -55, xmax = -73, ymax = -40) 
extent_south = extent(-77, -73, -55, -40)

#Define north area
box_north = c(xmin = -80, ymin = 1, xmax = -76, ymax = 7.5) 
extent_north = extent(-80, -76, 1, 7.5) 

### Subset data that is within a certain area
data_south = st_crop(datatest_sf, box_south)
data_north = st_crop(datatest_sf, box_north)

datatest_south = datatest_sf[row.names(datatest_sf)%in%row.names(data_south),]
datatest_north = datatest_sf[row.names(datatest_sf)%in%row.names(data_north),]

#South
pobs1vals = pobs(c(datatest_south$bio1, datatest_sf$bio1))[1:nrow(datatest_south)]
pobs7vals = pobs(c(datatest_south$bio7, datatest_sf$bio7))[1:nrow(datatest_south)]

hpts <- chull(pobs1vals, pobs7vals)
hpts <- c(hpts, hpts[1])

#North
pobs1valsN = pobs(c(datatest_north$bio1, datatest_sf$bio1))[1:nrow(datatest_north)]
pobs7valsN = pobs(c(datatest_north$bio7, datatest_sf$bio7))[1:nrow(datatest_north)]

hptsN <- chull(pobs1valsN, pobs7valsN)
hptsN <- c(hptsN, hptsN[1])

#Bio1
curve1 = getBayeSDMdensityfunct(finalBayeSDMImodel, type="1", variable="bio1")
curve0 = getBayeSDMdensityfunct(finalBayeSDMImodel, type="0", variable="bio1")

curvesbio1 = function(){
  curve(curve1(x), xlim= c(min(data@data$bio1), max(data@data$bio1)), ylim=c(0,0.005), cex=1.5, yaxt="n", ylab="",xlab="Mean temperature in d?C", axes=F,frame.plot=TRUE, col="green", lwd=2, lty=2) #bio1: mean temp.
  Axis(side=1, labels=T)
  curve(curve0(x), cex=1.5, col="red", lwd=2, add=T)
  title(ylab = "Density", cex=1, line = 1) 
  
  test = density(data_south$bio1, adjust=2, from=-23, to=289)
  test$y = test$y*0.1
  lines(test, col="darkgrey", lwd=2, cex=1.5, lty=3, xlim=c(-5,300))
  text(69,0.0022, "S", cex=1.5)
  
  
  test = density(data_north$bio1, adjust=2, from=-23, to=289)
  test$y = test$y*0.1
  lines(test, col="darkgrey", lwd=2, lty=3, cex=1.5, xlim=c(-5,300))
  text(256,0.0015, "N", cex=1.5)
  
  legend("topleft", legend=c("background", "presence", "region"), lwd=2, col=c("red", "green", "darkgrey"), cex=1, lty=1:3, bty = "n")
}

f5_plot1 = as.grob(curvesbio1)

#Bio12
curve1 = getBayeSDMdensityfunct(finalBayeSDMImodel, type="1", variable="bio12")
curve0 = getBayeSDMdensityfunct(finalBayeSDMImodel, type="0", variable="bio12")

curvesbio12 = function(){
  curve(curve1(x), xlim= c(min(data@data$bio12), max(data@data$bio12)), cex=1.5, ylim=c(0,0.00035),ylab="", xlab="Annual precipitation in mm", axes=F,frame.plot=TRUE, col="green", lwd=2, lty=2) #bio12: precipitation
  Axis(side=1, labels=T)
  curve(curve0(x),lwd=2, col="red",cex=1.5, add=T)
  title(ylab = "Density", cex=1, line = 1) 
  
  test = density(data_south$bio12, adjust=2, from=0, to=7682)
  test$y = test$y*0.25
  lines(test, col="darkgrey", lwd=2, lty=3, cex=1.5, xlim=c(0,7682))
  text(2300,0.00013, "S", cex=1.5)
  
  
  test = density(data_north$bio12, adjust=2, from=0, to=7682)
  test$y = test$y*0.25
  lines(test, col="darkgrey", lwd=2, cex=1.5, lty=3, xlim=c(0,7682))
  text(6000,0.00005, "N", cex=1.5)
}

curvesbio12()

f5_plot2 = as.grob(curvesbio12)

### 
# tiff("predictionscomparisonplot_maxent.tiff", width = 6, height = 6, units = 'in', res = 300)
# dev.off() 

plotc = function(){
  heat <- heat.colors(100, rev=T)
  plot(pobs(datatest_sf$bio7),pobs(datatest_sf$bio1), col=alpha(heat[datatest_sf$maxentprob4],0.15), 
       pch=21, bg=alpha(heat[datatest_sf$maxentprob4],0.15),cex=1.5,
       xlab="F(Temperature annual range in d°C)", ylab="F(Mean temperature in d°C)"
  ) 
  points(pobs(c(datapres_sf$bio7, datatest_sf$bio7))[1:116], pobs(c(datapres_sf$bio1, datatest_sf$bio1))[1:116], col="black", bg="black", cex=1.5, pch=21)
  
  lines(pobs7vals[hpts],pobs1vals[hpts], col=alpha("black",0.7), cex=1.5, lwd=2)
  text(0.4,0.1, "S", cex=1.5)
  
  lines(pobs7valsN[hptsN], pobs1valsN[hptsN], col=alpha("black",0.7), cex=1.5, lwd=2)
  text(0.035,0.35, "N", cex=1.5)
}

f5_plot3 = as.grob(plotc)


plotd = function(){
  plot(pobs(datatest_sf$bio7),pobs(datatest_sf$bio1), col=alpha(heat[datatest_sf$kde4vprob4],0.15), pch=21, bg=alpha(heat[datatest_sf$kde4vprob4],0.15),cex=1.5,
       xlab="F(Temperature annual range in d°C)", ylab="F(Mean temperature in d°C)") 
  points(pobs(c(datapres_sf$bio7, datatest_sf$bio7))[1:116],pobs(c(datapres_sf$bio1, datatest_sf$bio1))[1:116], col="black", bg="black", cex=1.5, pch=21)
  
  lines(pobs7vals[hpts],pobs1vals[hpts], col=alpha("black",0.7), cex=1.5, lwd=2)
  text(0.4,0.1, "S", cex=1.5)
  
  lines(pobs7valsN[hptsN], pobs1valsN[hptsN], col=alpha("black",0.7), cex=1.5, lwd=2)
  text(0.035,0.35, "N", cex=1.5)
}

f5_plot4 = as.grob(plotd)


#construct each plot as grob, then use cowplot format. 
tiff("Figure6.tiff", width = 10, height = 4, units = 'in', res = 1000, compression="lzw")
plot_grid(f5_plot3, f5_plot4, labels = c('(a) Maxent', '(b) Gen-KDE-4v'), ncol=2, label_size = 12, align = "hv", byrow=F) 
dev.off()







##########KS test and cor
datpres = data@pa==1

ks.test(data@data$bio1[datpres],data@data$bio1[!datpres])
#D = 0.30848, p-value = 5.723e-10

ks.test(data@data$bio7[datpres],data@data$bio7[!datpres])
#D = 0.52016, p-value < 2.2e-16

ks.test(data@data$bio12[datpres],data@data$bio12[!datpres])
#D = 0.46059, p-value < 2.2e-16


ks.test(data@data$bio16[datpres],data@data$bio16[!datpres])
#D = 0.39494, p-value = 4.441e-16

ks.test(data@data$bio5[datpres],data@data$bio5[!datpres])
#D = 0.13006, p-value = 0.04024


ks.test(data@data$bio6[datpres],data@data$bio6[!datpres])
#D = 0.45097, p-value < 2.2e-16


ks.test(data@data$bio8[datpres],data@data$bio8[!datpres])
#D = 0.22513, p-value = 1.652e-05

cor(data@data[datpres,])
cor(data@data[!datpres,]); colSums(cor(data@data[!datpres,]))

cor(data@data[datpres,], method="spearman")
cor(data@data[!datpres,], method="spearman"); colSums(cor(data@data[!datpres,], method="spearman"))

plot(density(data@data$bio12[!datpres], adjust=30), xlim=c(min(data@data$bio12[!datpres]),max(data@data$bio12[!datpres])))

