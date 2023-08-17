#####################
# Script to reproduce the performance benchmark included in the paper
# 'Species distribution modelling using a generative framework'. 
# Robbert T. van den Dool, Alejandro Morales, Wopke van der Werf & J.C. (Bob) Douma
# 17 April 2023, Robbert T. van den Dool, Centre for Crop Systems Analysis, Wageningen University 
#####################

### For inspection of results or generation of Figure 2 without running the analysis: Load the performance results directly. 
#allresults = readRDS("Intermediate/30June2023_allresults.rds")
## Then continue at line 189

### Load functions and packages. 
source("Scripts/Functions.R")

if(!require("easypackages")) install.packages("easypackages")
library(easypackages)
packages("raster", "disdat", "furrr", "ks", "rvinecopulib","dismo", "maxnet", "naivebayes", "caret", "ggplot2", "RColorBrewer", "tidyr", "dplyr", "sf", "sp", "SDMtune","PRROC", "Metrics", "gridExtra", "copula", prompt = FALSE)
set.seed(9999)

### Load data using package 'disdat'

### Canada 
can_po <- disPo("CAN")
can_bg <- disBg("CAN")
can_pa <- disPa("CAN")
can_env <- disEnv("CAN")
can_varnames = c("alt", "asp2", "ontprec", "ontprec4", "ontprecsd", "ontslp", "onttemp", "onttempsd", "onttmin4", "watdist")

### Australian Wet Tropics: 20 plants, 20 birds
awt_po <- disPo("AWT") #two groups
awt_bg <- disBg("AWT")

awt_po_bird = awt_po[awt_po$group=="bird",]
awt_pa_bird <- disPa("AWT", group="bird")
awt_env_bird <- disEnv("AWT", group="bird")

awt_po_plant = awt_po[awt_po$group=="plant",]
awt_pa_plant <- disPa("AWT", group="plant")
awt_env_plant <- disEnv("AWT", group="plant")

awt_varnames = c("bc01","bc04","bc05","bc06","bc12","bc15","bc17","bc20","bc31","bc33","slope","topo","tri")

### New South Wales: 54 species from 8 groups. ba   db   nb   ot   ou   rt   ru   sr 
nsw_po <- disPo("NSW")
nsw_bg <- disBg("NSW")

nsw_po_ba = nsw_po[nsw_po$group=="ba",]
nsw_pa_ba <- disPa("NSW", group="ba")
nsw_env_ba <- disEnv("NSW", group="ba")

nsw_po_db = nsw_po[nsw_po$group=="db",]
nsw_pa_db <- disPa("NSW", group="db")
nsw_env_db <- disEnv("NSW", group="db")

nsw_po_nb = nsw_po[nsw_po$group=="nb",]
nsw_pa_nb <- disPa("NSW", group="nb")
nsw_env_nb <- disEnv("NSW", group="nb")

nsw_po_ot = nsw_po[nsw_po$group=="ot",]
nsw_pa_ot <- disPa("NSW", group="ot")
nsw_env_ot <- disEnv("NSW", group="ot")

nsw_po_ou = nsw_po[nsw_po$group=="ou",]
nsw_pa_ou <- disPa("NSW", group="ou")
nsw_env_ou <- disEnv("NSW", group="ou")

nsw_po_rt = nsw_po[nsw_po$group=="rt",]
nsw_pa_rt <- disPa("NSW", group="rt")
nsw_env_rt <- disEnv("NSW", group="rt")

nsw_po_ru = nsw_po[nsw_po$group=="ru",]
nsw_pa_ru <- disPa("NSW", group="ru")
nsw_env_ru <- disEnv("NSW", group="ru")

nsw_po_sr = nsw_po[nsw_po$group=="sr",]
nsw_pa_sr <- disPa("NSW", group="sr")
nsw_env_sr <- disEnv("NSW", group="sr")

nsw_varnames = c("cti","mi","rainann","raindq","rugged","soildepth","solrad","tempann","tempmin","topo")

### New Zealand: 52 plants
nz_po <- disPo("NZ")
nz_bg <- disBg("NZ")
nz_pa <- disPa("NZ")
nz_env <- disEnv("NZ")
nz_varnames = c("deficit","dem","hillshade","mas","mat","r2pet","rain","slope","sseas","tseas","vpd"  )

### South America: 30 plants 
sa_po <- disPo("SA")
sa_bg <- disBg("SA")
sa_pa <- disPa("SA")
sa_env <- disEnv("SA")
sa_varnames = c("sabio1","sabio2","sabio4","sabio5","sabio6","sabio7","sabio8","sabio12","sabio15","sabio17","sabio18")

### Switzerland: 30 trees 
swi_po <- disPo("SWI")
swi_bg <- disBg("SWI")
swi_pa <- disPa("SWI")
swi_env <- disEnv("SWI")
swi_varnames = c("bcc","ccc","ddeg","nutri","pday","precyy","slope","sradyy","swb","tavecc","topo")



### Fit models to each species
#plan(sequential)
plan(multisession, workers = 6) #if your computer has more threads, you can set this value higher for faster computations.

runs_canada = future_map(names(table(can_po$spid)), function(x)fitmodels_species_all(PO=can_po, back=can_bg, PA=can_pa, env=can_env, species=x, varnames=can_varnames)) 
runs_awt_bird = future_map(names(table(awt_po_bird$spid)), function(x)fitmodels_species_all(PO=awt_po_bird, back=awt_bg, PA=awt_pa_bird, env=awt_env_bird, species=x, varnames=awt_varnames)) 
runs_awt_plant = future_map(names(table(awt_po_plant$spid)), function(x)fitmodels_species_all(PO=awt_po_plant, back=awt_bg, PA=awt_pa_plant, env=awt_env_plant, species=x, varnames=awt_varnames)) 
runs_nsw_ba = future_map(names(table(nsw_po_ba$spid)), function(x)fitmodels_species_all(PO=nsw_po_ba, back=nsw_bg, PA=nsw_pa_ba, env=nsw_env_ba, species=x, varnames=nsw_varnames)) 
runs_nsw_db = future_map(names(table(nsw_po_db$spid)), function(x)fitmodels_species_all(PO=nsw_po_db, back=nsw_bg, PA=nsw_pa_db, env=nsw_env_db, species=x, varnames=nsw_varnames))  
runs_nsw_nb = future_map(names(table(nsw_po_nb$spid)), function(x)fitmodels_species_all(PO=nsw_po_nb, back=nsw_bg, PA=nsw_pa_nb, env=nsw_env_nb, species=x, varnames=nsw_varnames))  
runs_nsw_ot = future_map(names(table(nsw_po_ot$spid)), function(x)fitmodels_species_all(PO=nsw_po_ot, back=nsw_bg, PA=nsw_pa_ot, env=nsw_env_ot, species=x, varnames=nsw_varnames))  
runs_nsw_ou = future_map(names(table(nsw_po_ou$spid)), function(x)fitmodels_species_all(PO=nsw_po_ou, back=nsw_bg, PA=nsw_pa_ou, env=nsw_env_ou, species=x, varnames=nsw_varnames))  
runs_nsw_rt = future_map(names(table(nsw_po_rt$spid)), function(x)fitmodels_species_all(PO=nsw_po_rt, back=nsw_bg, PA=nsw_pa_rt, env=nsw_env_rt, species=x, varnames=nsw_varnames))  
runs_nsw_ru = future_map(names(table(nsw_po_ru$spid)), function(x)fitmodels_species_all(PO=nsw_po_ru, back=nsw_bg, PA=nsw_pa_ru, env=nsw_env_ru, species=x, varnames=nsw_varnames))  
runs_nsw_sr = future_map(names(table(nsw_po_sr$spid)), function(x)fitmodels_species_all(PO=nsw_po_sr, back=nsw_bg, PA=nsw_pa_sr, env=nsw_env_sr, species=x, varnames=nsw_varnames))  
runs_nz = future_map(names(table(nz_po$spid)), function(x)fitmodels_species_all(PO=nz_po, back=nz_bg, PA=nz_pa, env=nz_env, species=x, varnames=nz_varnames)) 
runs_sa = future_map(names(table(sa_po$spid)), function(x)fitmodels_species_all(PO=sa_po, back=sa_bg, PA=sa_pa, env=sa_env, species=x, varnames=sa_varnames)) 
runs_swi = future_map(names(table(swi_po$spid)), function(x)fitmodels_species_all(PO=swi_po, back=swi_bg, PA=swi_pa, env=swi_env, species=x, varnames=swi_varnames)) 


#make df. 
runs_canada_df = do.call("rbind",runs_canada)
runs_awt_bird_df = do.call("rbind",runs_awt_bird)
runs_awt_plant_df = do.call("rbind",runs_awt_plant)
runs_nsw_ba_df = do.call("rbind",runs_nsw_ba)
runs_nsw_db_df = do.call("rbind",runs_nsw_db)
runs_nsw_nb_df = do.call("rbind",runs_nsw_nb)
runs_nsw_ot_df = do.call("rbind",runs_nsw_ot)
runs_nsw_ou_df = do.call("rbind",runs_nsw_ou)
runs_nsw_rt_df = do.call("rbind",runs_nsw_rt)
runs_nsw_ru_df = do.call("rbind",runs_nsw_ru)
runs_nsw_sr_df = do.call("rbind",runs_nsw_sr)
runs_nz_df = do.call("rbind",runs_nz)
runs_sa_df = do.call("rbind",runs_sa)
runs_swi_df = do.call("rbind",runs_swi)


runs_canada_df = as.data.frame(runs_canada_df)
runs_canada_df$npres = table(can_po$spid)

runs_awt_bird_df = as.data.frame(runs_awt_bird_df)
runs_awt_bird_df$npres = table(awt_po_bird$spid)

runs_awt_plant_df = as.data.frame(runs_awt_plant_df)
runs_awt_plant_df$npres = table(awt_po_plant$spid)

runs_nsw_ba_df = as.data.frame(runs_nsw_ba_df)
runs_nsw_ba_df$npres = table(nsw_po_ba$spid)

runs_nsw_db_df = as.data.frame(runs_nsw_db_df)
runs_nsw_db_df$npres = table(nsw_po_db$spid)

runs_nsw_nb_df = as.data.frame(runs_nsw_nb_df)
runs_nsw_nb_df$npres = table(nsw_po_nb$spid)

runs_nsw_ot_df = as.data.frame(runs_nsw_ot_df)
runs_nsw_ot_df$npres = table(nsw_po_ot$spid)

runs_nsw_ou_df = as.data.frame(runs_nsw_ou_df) 
runs_nsw_ou_df$npres = table(nsw_po_ou$spid)

runs_nsw_rt_df = as.data.frame(runs_nsw_rt_df)
runs_nsw_rt_df$npres = table(nsw_po_rt$spid)

runs_nsw_ru_df = as.data.frame(runs_nsw_ru_df)
runs_nsw_ru_df$npres = table(nsw_po_ru$spid)

runs_nsw_sr_df = as.data.frame(runs_nsw_sr_df)
runs_nsw_sr_df$npres = table(nsw_po_sr$spid)

runs_nz_df = as.data.frame(runs_nz_df) 
runs_nz_df$npres = table(nz_po$spid)

runs_sa_df = as.data.frame(runs_sa_df)
runs_sa_df$npres = table(sa_po$spid)

runs_swi_df = as.data.frame(runs_swi_df)
runs_swi_df$npres = table(swi_po$spid)


allresults = rbind(runs_canada_df, runs_awt_bird_df, runs_awt_plant_df, 
                   runs_nsw_ba_df, runs_nsw_db_df, runs_nsw_nb_df, runs_nsw_ot_df,
                   runs_nsw_ou_df, runs_nsw_rt_df, runs_nsw_ru_df, runs_nsw_sr_df,
                   runs_nz_df, runs_sa_df, runs_swi_df)

allresults_nona = allresults[!is.na(allresults$maxnet),-6] #remove 5 species which had non-converging models


### Inspect results
myCI(allresults$maxnet) #0.7421372 0.7230491 0.7039609
myCI(allresults$maxnet4) #0.7345072 0.7150013 0.6954955 
myCI(allresults$kde)  #0.7270788 0.7080911 0.6891034 
myCI(allresults$vine) #0.6991780 0.6813369 0.6634958 
myCI(allresults$vine4) #0.7000898 0.6823935 0.6646972 

sum(sapply(1:nrow(allresults_nona), function(x) max(allresults_nona$maxnet[x], allresults_nona$maxnet4[x],na.rm=T) < max(allresults_nona$kde[x],na.rm=T)))
#in 64 cases kde performed better than maxnet 

sum(sapply(1:nrow(allresults_nona), function(x) max(allresults_nona$maxnet[x], allresults_nona$maxnet4[x],na.rm=T) < max(allresults_nona$vine[x], allresults_nona$vine4[x],na.rm=T)))
#in 79 cases the vine performed better than maxnet

sum(sapply(1:nrow(allresults_nona), function(x) max(allresults_nona$maxnet[x], allresults_nona$maxnet4[x],na.rm=T) < max(allresults_nona$kde[x], allresults_nona$vine[x], allresults_nona$vine4[x],na.rm=T)))
#in 115 cases a version of BayeSDM performed better 

sum(sapply(1:nrow(allresults_nona), function(x) max(allresults_nona$maxnet[x], allresults_nona$maxnet4[x],na.rm=T) > max(allresults_nona$kde[x], allresults_nona$vine[x], allresults_nona$vine4[x],na.rm=T)))
#in 106 cases a version of Maxnet performed better 

sum(sapply(1:nrow(allresults_nona), function(x) max(allresults_nona$kde[x],na.rm=T) < max(allresults_nona$vine[x], allresults_nona$vine4[x],na.rm=T)))
# in 107 cases the vine performed better than the kde.

sum(sapply(1:nrow(allresults_nona), function(x) max(allresults_nona$kde[x],na.rm=T) > max(allresults_nona$vine[x], allresults_nona$vine4[x],na.rm=T)))
# in 114 cases the kde performed better than the vine.



### Transform results for Figure 2.
#Determine which model had the best performance for each species, group them together, sorted by performance. 
#Also create a column indicating whether the model was a generative or discriminative model. 

genwins = sapply(1:nrow(allresults_nona), function(x) max(allresults_nona$maxnet[x], allresults_nona$maxnet4[x],na.rm=T) < max(allresults_nona$kde[x], allresults_nona$vine[x], allresults_nona$vine4[x],na.rm=T))
allresults_nona$genwin = genwins

allresults_nona$species = 1:nrow(allresults_nona)

allresults_nona2 = allresults_nona
names(allresults_nona2)[1:5] = c("Maxent", "Maxent-4v", "Gen-KDE-4v", "Gen-Vine", "Gen-Vine-4v")

ranks_maxent = sapply(1:221, function(x) allresults_nona2$Maxent[x] >= max(c(allresults_nona2$`Maxent`[x], allresults_nona2$`Maxent-4v`[x], allresults_nona2$`Gen-KDE-4v`[x], allresults_nona2$`Gen-Vine`[x], allresults_nona2$`Gen-Vine-4v`[x]), na.rm=T))
ranks_maxent4v = sapply(1:221, function(x) allresults_nona2$`Maxent-4v`[x] >= max(c(allresults_nona2$`Maxent`[x], allresults_nona2$`Maxent-4v`[x], allresults_nona2$`Gen-KDE-4v`[x], allresults_nona2$`Gen-Vine`[x], allresults_nona2$`Gen-Vine-4v`[x]), na.rm=T))
ranks_kde = sapply(1:221, function(x) allresults_nona2$`Gen-KDE-4v`[x] >= max(c(allresults_nona2$`Maxent`[x], allresults_nona2$`Maxent-4v`[x], allresults_nona2$`Gen-KDE-4v`[x], allresults_nona2$`Gen-Vine`[x], allresults_nona2$`Gen-Vine-4v`[x]), na.rm=T))
ranks_vine = sapply(1:221, function(x) allresults_nona2$`Gen-Vine`[x] >= max(c(allresults_nona2$`Maxent`[x], allresults_nona2$`Maxent-4v`[x], allresults_nona2$`Gen-KDE-4v`[x], allresults_nona2$`Gen-Vine`[x], allresults_nona2$`Gen-Vine-4v`[x]), na.rm=T))
ranks_vine4v = sapply(1:221, function(x) allresults_nona2$`Gen-Vine-4v`[x] >= max(c(allresults_nona2$`Maxent`[x], allresults_nona2$`Maxent-4v`[x], allresults_nona2$`Gen-KDE-4v`[x], allresults_nona2$`Gen-Vine`[x], allresults_nona2$`Gen-Vine-4v`[x]), na.rm=T))

#now turn into three categories
ranks_maxent = ifelse(ranks_maxent == T, "Maxent", "Loss")
ranks_maxent4v = ifelse(ranks_maxent4v == T, "Maxent", "Loss")
ranks_kde = ifelse(ranks_kde == T, "Gen", "Loss")
ranks_vine = ifelse(ranks_vine == T, "Gen", "Loss")
ranks_vine4v = ifelse(ranks_vine4v == T, "Gen", "Loss")
ranks = c(ranks_maxent, ranks_maxent4v, ranks_kde, ranks_vine, ranks_vine4v) 

allresults_long3 = allresults_nona2[,1:5] %>% 
  gather(key="Model", value="AUC")
allresults_long3$species= rep(allresults_nona2$species, 5)
allresults_long3$genwins = ranks


allresults_long3$winner <- NA
for (i in allresults_long3$species){
  models <- allresults_long3[allresults_long3$species == i,"Model"]
  winner <- models[which.max(allresults_long3[allresults_long3$species == i,"AUC"])]
  allresults_long3[allresults_long3$species == i,"winner"] <- winner
}

allresults_long6 <- do.call(rbind,by(allresults_long3,list(allresults_long3$winner),function(x) x[order(x$AUC,decreasing = T),]))
min.species <- 0
new.species = c()
un.species = c()  
for (j in c("Maxent", "Gen-Vine", "Gen-KDE-4v", "Maxent-4v", "Gen-Vine-4v")){ #unique(allresults_long4$winner
  un <- unique(allresults_long6[allresults_long6$winner == j,"species"])
  un.species = c(un.species,un)
  new.species <- c(new.species,rev(c((min.species+1):(min.species+length(un)))))
  min.species = max(new.species)
}

match.table <-data.frame(un.species,new.species)

allresults_long6$species4<- match.table[match(allresults_long6$species, match.table$un.species),"new.species"]

allresults_long6$genwins = factor(allresults_long6$genwins, levels=c("Loss", "Maxent", "Gen"))

allresults_long6 <- allresults_long6[order(allresults_long6$genwins, decreasing=F),]


# wins for each model
sum(allresults_long6$winner=="Maxent")/5 #76 76/226*100 #33.6%
sum(allresults_long6$winner=="Gen-Vine")/5 #47 #20.8%
sum(allresults_long6$winner=="Gen-KDE-4v")/5 #43 #19.0%
sum(allresults_long6$winner=="Maxent-4v")/5 #30 #13.3%
sum(allresults_long6$winner=="Gen-Vine-4v")/5 #25 #11.1%


### Create Figure 2
tiff("Figure2_30June2022.tiff", width = 6, height = 5, units = 'in', res = 800, compression="lzw")
ggplot(allresults_long6, aes(x=species4, y=AUC, col=genwins, alpha=genwins, shape=Model))+ #, alpha = 0.8 #alpha = ifelse(genwins=="Loss",0.95,1) #alpha = ifelse(genwins=="Loss",0.9,1)#alpha=0.3,#factor(species, levels=unique(allresults_long2$species2))
  geom_point()+ 
  scale_colour_manual(breaks=c('Loss', 'Maxent', 'Gen'), labels=c('Loss', 'Maxent', 'Generative'), values = c("grey", brewer.pal(9,"Set1")[1], brewer.pal(9,"Set1")[2]), aesthetics = c("colour", "fill"), guide=FALSE)+
  
  geom_segment(x = 1, y = 0.7230491, xend = 76, yend = 0.7230491, alpha=0.5, col="black")+ #36.65158 #maxent
  geom_segment(x = 77, y = 0.6813369, xend = 123, yend = 0.6813369, alpha=0.5, col="black")+ #23.52941 #vine
  geom_segment(x = 123, y = 0.7080911, xend = 166, yend = 0.7080911, alpha=0.5, col="black")+ #6.78733 #kde
  geom_segment(x = 167, y = 0.7150013, xend = 196, yend = 0.7150013, alpha=0.5, col="black")+ ##12.66968 #m4
  geom_segment(x = 197, y = 0.6823935, xend = 221, yend = 0.6823935, alpha=0.5, col="black")+ #20.36199 #v4
  
  scale_shape_discrete(breaks=c('Maxent', 'Maxent-4v', 'Gen-KDE-4v', 'Gen-Vine', 'Gen-Vine-4v'))+
  
  annotate(geom="text", x=33, y=0.275, label="Maxent",color="black")+
  annotate(geom="text", x=33, y=0.225, label="33.6%",color="black")+
  annotate(geom="text", x=80, y=0.275, label="Gen-Vine",color="black")+
  annotate(geom="text", x=80, y=0.225, label="20.8%",color="black")+
  annotate(geom="text", x=122, y=0.275, label="Gen-KDE-4v",color="black")+
  annotate(geom="text", x=122, y=0.225, label="19.0%",color="black")+
  annotate(geom="text", x=169, y=0.275, label="Maxent-4v",color="black")+
  annotate(geom="text", x=169, y=0.225, label="13.3%",color="black")+
  annotate(geom="text", x=208, y=0.275, label="Gen-Vine-4v",color="black")+
  annotate(geom="text", x=208, y=0.225, label="11.1%",color="black")+
  
  theme_bw()+
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()#,
  ) +
  ylab("AUC") + ylim(0.2,1) + labs(colour = "Best model") + scale_alpha_manual(genwins,values = c(0.75, 1, 1), guide=FALSE)+ #xlab("Species")+ #scale_alpha_ordinal(range = c(0.9, 0.35, 0.9))#scale_alpha(guide=FALSE) #, guide=FALSE
  scale_x_continuous(name = "Species", breaks=c(-30,251), labels=c("","")) #,limits
dev.off()
