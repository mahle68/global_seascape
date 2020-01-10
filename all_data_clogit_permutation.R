#scripts for meta-analysis of sea-crossing behavior using permutation tests
#follows the data preparation in all_data_prep.R
#Jan 10. 2020. Elham Nourani. Radolfzell, Germany.


library(dplyr)
library(lme4)
library(survival)
library(TwoStepCLogit)
library(rstanarm)
library(INLA)
library(TwoStepCLogit)
library(mclogit)
library(jtools) #summ ftn
#library(ggstance)
library(MuMIn) #adjusted r squared
library(uhcplots)

setwd("/home/enourani/ownCloud/Work/Projects/delta_t")
setwd("/home/mahle68/ownCloud/Work/Projects/delta_t")
wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
meters_proj <- CRS("+proj=moll +ellps=WGS84")
mycl <- makeCluster(detectCores() - 2)

setwd("/home/enourani/ownCloud/Work/Projects/delta_t")


load("R_files/all_spp_temp_sp_filtered_15km_alt_ann_14days.RData") #called ann

#scale the variables, zone-specific
ann_z<-ann%>%
  group_by(zone) %>% #cosider season later 
  mutate_at(c(21,17,14),
            list(z=~scale(.)))%>%
  as.data.frame()

formula <- used ~ delta_t_z + u925_z + v925_z + strata(obs_id)
formula2 <- used ~ delta_t_z + u925_z + strata(obs_id)

################################################################ tradewind zone
#Question: is the effect of delta T larger than wind in each zone?

##### STEP 1: calculate the observed statistic #####
model_twz_obs <- clogit(formula, data = ann_z[ann_z$zone == "tradewind",])
summary(model_twz_obs) #3611.451
AIC(model_twz_obs)

model_twz_obs2 <- clogit(formula2, data = ann_z[ann_z$zone == "tradewind",])
summary(model_twz_obs2)
AIC(model_twz_obs2) #3609.498

coef(model_twz_obs)
stat_obs_twz_u <- abs(coef(model_twz_obs)[1])-abs(coef(model_twz_obs)[2])
stat_obs_twz_v <- abs(coef(model_twz_obs)[1])-abs(coef(model_twz_obs)[3])

##### STEP 2: create alternative datasets #####
permutations <- 1000

#create alternative datasets  
new_data_twz_ls <- lapply(1:permutations, function(x){ #for each permutation
  new_obs_id_ls <-lapply(split(ann_z[ann_z$zone == "tradewind",], ann_z[ann_z$zone == "tradewind","obs_id"]), function(y){ #for each obs_id
    rnd_obs <- sample(1:15,1) #draw a random number to be the new index for used row
    y[rnd_obs,"used"] <- 1
    y[-rnd_obs,"used"] <- 0
    y
  })
  new_obs_id <- do.call(rbind, new_obs_id_ls)
  new_obs_id
}) 

##### STEP 3: calculate the random statistics #####

rnd_st_twz_ls <- lapply(new_data_twz_ls, function(x){
  model <-  clogit(formula, data = x[x$zone == "tradewind",])
  stat_twz_u <- abs(coef(model)[1])-abs(coef(model)[2])
  stat_twz_v <- abs(coef(model)[1])-abs(coef(model)[3])
  data.frame(u = stat_twz_u, v= stat_twz_v, row.names = "")
})

rnd_stat_twz <- do.call(rbind, rnd_st_twz_ls)


##### STEP 4: plot the random and observed values #####
plot(1:permutations, rnd_stat_twz$u, type = "l", main = "tradewind - u wind")
abline(h = stat_obs_twz_u, col = "red")

par(mfrow = c(1,2))
#plot the frequency distribution of u vlaues
hist(rnd_stat_twz$u, breaks = 50, col = "lightgrey", xlim = c(-0.3,0.3), main = "delta_t - wind_u")
abline(v = stat_obs_twz_u, col = "red") #random network is homogenous, that is why cv goes down the more randomization we do
  
#plot the frequency distribution of v vlaues
hist(rnd_stat_twz$v, breaks = 50, col = "lightgrey", xlim = c(-0.3,0.3), main = "delta_t - wind_v")
abline(v = stat_obs_twz_v, col = "red") #random network is homogenous, that is why cv goes down the more randomization we do


##### STEP 5: calculate p-values #####
#calculate the p-value. one-tailed
p_value_u <- sum(stat_obs_twz_u <= rnd_stat_twz$u) / permutations # 1: hypothesis is rejected. delta t does not have a higher effect than u wind
p_value_v <- sum(stat_obs_twz_v <= rnd_stat_twz$v) / permutations # 2: hypothesis is accepted. delta t has a higher effect than v wind



################################################################ temperate zone
#Question: is the effect of delta T larger than wind in each zone?

##### STEP 1: calculate the observed statistic #####

model_tmpz_obs <- clogit(formula, data = ann_z[ann_z$zone == "temperate",])
summary(model_tmpz_obs) 
AIC(model_tmpz_obs) #1759.656

model_tmpz_obs2 <- clogit(formula2, data = ann_z[ann_z$zone == "temperate",])
summary(model_tmpz_obs2)
AIC(model_tmpz_obs2) #1758.918

coef(model_tmpz_obs)
stat_obs_tmpz_u <- abs(coef(model_tmpz_obs)[1])-abs(coef(model_tmpz_obs)[2])
stat_obs_tmpz_v <- abs(coef(model_tmpz_obs)[1])-abs(coef(model_tmpz_obs)[3])

##### STEP 2: create alternative datasets #####
permutations <- 1000

#create alternative datasets  
new_data_tmpz_ls <- lapply(1:permutations, function(x){ #for each permutation
  new_obs_id_ls <-lapply(split(ann_z[ann_z$zone == "temperate",], ann_z[ann_z$zone == "temperate","obs_id"]), function(y){ #for each obs_id
    rnd_obs <- sample(1:15,1) #draw a random number to be the new index for used row
    y[rnd_obs,"used"] <- 1
    y[-rnd_obs,"used"] <- 0
    y
  })
  new_obs_id <- do.call(rbind, new_obs_id_ls)
  new_obs_id
}) 

##### STEP 3: calculate the random statistics #####

rnd_st_tmpz_ls <- lapply(new_data_tmpz_ls, function(x){
  model <-  clogit(formula, data = x[x$zone == "temperate",])
  stat_tmpz_u <- abs(coef(model)[1])-abs(coef(model)[2])
  stat_tmpz_v <- abs(coef(model)[1])-abs(coef(model)[3])
  data.frame(u = stat_tmpz_u, v= stat_tmpz_v, row.names = "")
})

rnd_stat_tmpz <- do.call(rbind, rnd_st_tmpz_ls)


##### STEP 4: plot the random and observed values #####
plot(1:permutations, rnd_stat_tmpz$u, type = "l", main = "temperate - u wind")
abline(h = stat_obs_tmpz_u, col = "red")

par(mfrow = c(1,2))
#plot the frequency distribution of u vlaues
hist(rnd_stat_tmpz$u, breaks = 50, col = "lightgrey", xlim = c(-0.3,0.3), main = "delta_t - wind_u")
abline(v = stat_obs_tmpz_u, col = "red") #random network is homogenous, that is why cv goes down the more randomization we do

#plot the frequency distribution of v vlaues
hist(rnd_stat_tmpz$v, breaks = 50, col = "lightgrey", xlim = c(-0.3,0.3), main = "delta_t - wind_v")
abline(v = stat_obs_tmpz_v, col = "red") #random network is homogenous, that is why cv goes down the more randomization we do


##### STEP 5: calculate p-values #####
#calculate the p-value. one-tailed
p_value_u <- sum(stat_obs_tmpz_u <= rnd_stat_tmpz$u) / permutations # 1: hypothesis is rejected. delta t does not have a higher effect than u wind
p_value_v <- sum(stat_obs_tmpz_v <= rnd_stat_tmpz$v) / permutations # 2: hypothesis is accepted. delta t has a higher effect than v wind




