#scripts for meta-analysis of sea-crossing behavior
#follows the data preparation in all_data_prep.R
#Jan 3. 2020. Elham Nourani. Radolfzell, Germany.

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

##### STEP 1: check the data #####

load("R_files/all_spp_temp_sp_filtered_15km_alt_ann_14days.RData") #called ann

lapply(split(ann,ann$zone), function(x){
  X11()
  par(mfrow = c(1,3))
  hist(ann[ann])
})

#correlation


##### STEP 2: modeling... play around with options #####

#scale the variables, zone-specific
ann_z<-ann%>%
  group_by(zone) %>% 
  mutate_at(c(21,17,14),
            list(z=~scale(.)))%>%
  as.data.frame()


#formula1 <- used ~ scale(delta_t) + scale(u925) + scale(v925) + strata(obs_id)
formula1_sc <- used ~ delta_t_z + u925_z + v925_z + strata(obs_id)

m_twz1 <- clogit(formula1_sc, data = ann_z[ann_z$zone == "tradewind",])
summary(m_twz1)

#m_twz2 <- update(m_twz1,- scale(v925))

m_tmpz1 <- clogit(formula1, data = ann[ann$zone == "temperate",])
summary(m_tmpz1)

rsq(m_twz1)

summ(m_twz1)
r.squaredLR(m_twz1)


model <- glm(used ~ delta_t + u925 +  + v925 + species , family = binomial, data = ann_sc[ann_sc$zone == "tradewind",] ) #higher selection on u-wind
r.squaredLR(model)


clogit_ftn <- function(df) clogit(used ~ delta_t + u925 + v925, data = df)

ann_model <- ann %>%
  group_by(zone) %>% 
  mutate_at(c("delta_t","u925","v925"), scale) %>% #scale the data, zone-specific and copy the species var twice to be used as a mixed effect 
  as.data.frame() %>% 
  nest() %>%
  mutate(model = map(data, clogit_ftn))


#ann_sc <- ann %>%
#  group_by(zone) %>% 
#  mutate_at(c("delta_t","u925","v925"), scale) %>% 
#  mutate(species1 = species,
#         species2 = species) %>% 
#  as.data.frame()

clogit_ftn <- function(df) clogit(formula1, data = df)
clogit_mixed_ftn <- function(df) clogit(formula1, data = df)

ann %>% 
  group_by(zone) %>% 
  nest() %>% 
  mutate(model = map(data,clogit_ftn))


#tradewind zone
formula1 <- used ~ scale(delta_t) + scale(u925) + scale(v925)
formula1_mixed <- used ~ scale(delta_t) + scale(u925) + scale(v925) + (1|species)
formula1_mixed_slopes <- used ~ scale(delta_t) + scale(u925) + scale(v925) + (scale(delta_t)|species) + (scale(u925)|species1) + (scale(v925)|species2) #one variable can't be used multiple times
formula1_mixed_unscaled <- used ~ delta_t + u925 + v925 + (1|species)
formula2 <- used ~ scale(delta_t) + scale(wspd)

#clogit for both seasons together
model_twz3 <- clogit(formula1 + strata(obs_id), data = ann[ann$zone == "tradewind",])
model_twz3

#rstanarm::stan_clogit
#http://mc-stan.org/rstanarm/articles/binomial.html

#make sure data is ordered based on the strata
ann <- ann[order(ann$zone, ann$obs_id),]

model_twz4 <- stan_clogit(formula1_mixed_unscaled, strata = obs_id, data = ann[ann$zone == "tradewind",])
summary(model_twz4)


model_twz4_spd <- stan_clogit(formula2, strata = obs_id, data = ann[ann$zone == "tradewind",])
summary(model_twz4_spd)

d <- Sys.time()
model_twz_ms <- stan_clogit(formula1_mixed_slopes, 
                            strata = obs_id, data = ann2[ann$zone == "tradewind",],
                            cores = 10) #how many cores should the computer utilize
Sys.time()-d

summary(model_twz_ms)

model2<- stan_clogit(used ~ scale())


#Model evaluation
round(posterior_interval(model_twz4, prob = 0.5), 2)

# Predicted probability as a function of x
pr_used <- function(x, ests) plogis(ests[1] + ests[2] * x)
# A function to slightly jitter the binary data
jitt <- function(...) {
  geom_point(aes_string(...), position = position_jitter(height = 0.05, width = 0.1), 
             size = 2, shape = 21, stroke = 0.2)
}
ggplot(ann[ann$zone == "tradewind",], aes(x = scale(delta_t), y = used, color = used)) + 
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  jitt(x="delta_t") + 
  stat_function(fun = pr_used, args = list(ests = coef(model_twz4)), 
                size = 2, color = "gray35")

pr_used2 <- function(x, y, ests) plogis(ests[1] + ests[2] * x + ests[3] * y)
grid <- expand.grid(delta_t = seq(-10, 10, length.out = 100), 
                    wind_u = seq(-10, 10, length.out = 100))
grid$prob <- with(grid, pr_used2("scale(delta_t)", "scale(u925)", coef(model_twz4)))
ggplot(grid, aes(x = "scale(delta_t)", y = "scale(u925)")) + 
  geom_tile(aes(fill = prob)) + 
  geom_point(data = ann[ann$zone == "tradewind",], aes(color = factor(used)), size = 2, alpha = 0.85) + 
  scale_fill_gradient() +
  scale_color_manual("used", values = c("white", "black"), labels = c("No", "Yes"))


### rinla https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2656.13087

formula_inla <- used ~ delta_t + u925 + v925 +
  f(obs_id, model = "iid",
    hyper = list(theta = list(initial = log(1e-6), fixed = TRUE))) +
  f(species, delta_t, values = 1:4, model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = FALSE,
                              prior = "pc.prec", param = c(1, 0.05)))) +
  f(species1, u925, values = 1:4, model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = FALSE,
                              prior = "pc.prec", param = c(1, 0.05)))) +
  f(species2, v925, values = 1:4, model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = FALSE,
                              prior = "pc.prec", param = c(1, 0.05))))

m1 <- inla(formula, family = "poisson", data = ann_sc[ann_sc$zone == "tradewind",]) #Error in formula[2] : object of type 'closure' is not subsettable



### two step clogit
m2 <- Ts.estim( used ~ delta_t + u925 + v925 + strata(obs_id) + cluster(species),
                random = ~ delta_t + u925 + v925,
                data = ann_sc[ann_sc$zone == "tradewind",],
                D = "UN")

#temperate zone 

#one glm for both seasons together... for zone-specific models, make sure to scale the data only based on values in that zone :p
model_tmpz <- glm(used ~ scale(delta_t) + scale(u925) + scale(v925) , family = binomial, data = ann[ann$zone == "temperate",]) #higher selection on u-wind

#model_tmpz7 <- glmer(used ~ -1 + delta_t + u925 + v925 + (delta_t | obs_id) + (u925 | obs_id) + (v925 | obs_id), family = poisson, data =  ann[ann$zone == "temperate",]) #singularity error

model_tmpz2 <- glmer(used ~ scale(delta_t) + scale(u925) + scale(v925) + (1 | obs_id), family = binomial, data = ann[ann$zone == "temperate",]  ) #not converge

#clogit for both seasons together
model_tmpz3 <- clogit(used ~ scale(delta_t) + scale(u925) + scale(v925) + strata(obs_id), data = ann[ann$zone == "temperate",] )
model_tmpz3

model_tmpz4 <- stan_clogit(used ~ scale(delta_t) + scale(u925) + scale(v925), strata = obs_id, data = ann[ann$zone == "temperate",])
model_tmpz4_2 <- stan_clogit(formula2, strata = obs_id, data = ann[ann$zone == "temperate",])
summary(model_tmpz4_2)

d <- Sys.time()
model_tmpz_ms <- stan_clogit(formula1_mixed_slopes, 
                             strata = obs_id, data = ann2[ann$zone == "temperate",],
                             cores = 10) #how many cores should the computer utilize
Sys.time()-d
summary(model_tmpz_ms)



##### STEP 3: model building and evalution using clogit and uhc plots #####
#codes for assessing the robustness of ssf using Fortin et al 2009's method

lapply(split(ann_z,ann_z$zone), function(zone){
  
  k_fold_eval<-rep(NA,100)
  
  for (i in 1:length(k_fold_eval)){
    
    #STEP 1: partition the strata 80%-20%
    partition <- as.numeric(sample(levels(factor(zone$obs_id)),0.2*length(levels(factor(zone$obs_id)))))
    testing <- zone %>% 
      filter(obs_id %in% partition)
    
    training <- zone %>% 
      filter(!(obs_id %in% partition))
    
    #STEP 2: clogit using training set
    #m_eval <- mclogit(cbind(used,obs_id) ~ delta_t + u925 + v925, 
    #                       random=~1|species, data = training)
    m_eval <- clogit(formula1_sc, data = training)
    
    #STEP 3: predict values of suitability for testing set (used and available)
    testing$pred <- as.numeric(predict(m_eval,type = "risk",newdata = testing))
    
    #STEP 4: rank data in each obs_id based on their pred value (1 the lowest, length(levels(testing$obs_id)) or 116 the highest)
    testing <- testing[order(testing$obs_id),]
    testing$ranks <- unlist(with(testing, tapply(pred, obs_id, rank)))
    
    #STEP 5: extract ranks of used points
    used_ranks<-subset(testing,used==1)
    
    #STEP 6:create bins for each possible rank
    bins_df<-data.frame(bins=c(1:116))
    
    for (i in 1:nrow(bins_df)){
      tally<-subset(used_ranks,ranks==i)
      bins_df$freq[i]<-nrow(tally)
    }
    
    #STEP 7: Spearman's rank correlation
    k_fold_eval[30]<-cor(bins_df,method="spearman")[2]
    
    rm(testing,training,bins_df)
    
    
  }
  
})


#results for k_fold_eval[1:50] (model: m_6b_eval <- mclogit(cbind(used,obs_id) ~ cos(rel.angle) + dist + BLH +
#dist_coast.sqrt.l.z * wind_support.z +  dist_coast.sqrt.l.z * cross_wind.z,random=~1|id, data = training))
#mean(k_fold_eval[1:50])
#[1] 0.338549
#range(k_fold_eval[1:50])
#[1] 0.2203467 0.4064363



#codes for calculating mean and range of correlation coefficients assuming random pattern of habitat selection


k_fold_random<-rep(NA,100)

for (i in 1:length(k_fold_random){
  
  #STEP 1: partition the strata 80%-20%
  partition<-sample(levels(zone$obs_id),0.2*length(levels(zone$obs_id))) 
  testing<-subset(zone,zone$obs_id %in% partition)
  testing$obs_id<-droplevels(testing$obs_id) #remove empty levels
  
  training<-subset(zone,!(zone$obs_id %in% partition))
  training$obs_id<-droplevels(training$obs_id) #remove empty levels
  
  #STEP 2: calculate SSF using training set
  
  m_6b_eval <- mclogit(cbind(used,obs_id) ~ cos(rel.angle) + dist + 
                       dist_coast.sqrt.l.z * wind_support.z +  dist_coast.sqrt.l.z * cross_wind.z,random=~1|id, data = training)
  
  
  #STEP 3: predict values of suitability for testing set (used and available)
  testing$pred<-as.numeric(predict(m_6b_eval,type="response",newdata=testing))
  
  #STEP 4: rank data in each obs_id based on their pred value (1 the lowest, length(levels(testing$obs_id)) or 116 the highest)
  testing <- testing[order(testing$obs_id),]
  testing$ranks <- unlist(with(testing, tapply(pred, obs_id, rank)))
  
  #STEP 5: extract ranks of used points
  randoms<-subset(testing,used==0)
  random_ranks_df<-randoms %>% group_by(obs_id) %>% sample_n(size = 1) #package dplyr
  #random_ranks<-random_ranks_df$ranks
  
  #STEP 6:create bins for each possible rank
  bins_df<-data.frame(bins=c(1:116))
  
  for (i in 1:nrow(bins_df)){
    tally<-subset(random_ranks_df,ranks==i)
    bins_df$freq[i]<-nrow(tally)
  }
  
  #STEP 7: Spearman's rank correlation
  k_fold_random[30]<-cor(bins_df,method="spearman")[2]
  
  rm(testing,training,bins_df)
  
  
})



#file:///home/enourani/ownCloud/Work/R_source_codes/uhcplots-master/vignettes/UCHplots-Step-Selection-Functions.html

#tradewind zone
partition_twz <- ann_z %>% 
  filter(zone == "tradewind") %>% 
  distinct(obs_id) %>% 
  sample_frac(0.2) %>% 
  pull(obs_id)

testing_twz <- ann_z %>%
  filter(zone == "tradewind" & obs_id %in% partition_twz)

training_twz <- ann_z %>% 
  filter(zone == "tradewind" & !(obs_id %in% partition_twz))

#formula1_sc <- used ~ delta_t + u925 + v925 + strata(obs_id)
textplot1 <- (expression(y %~% delta_t + u925 + v925 + strata(obs_id)))
form1a <- (used ~ delta_t_z + u925_z + v925_z  + strata(obs_id))
form2a <- ~ delta_t_z + u925_z + v925_z -1

#train full model
ssf_train_full <- clogit(form1a, data = ann_z[ann_z$zone == "tradewind",])
summary(ssf_train_full)

#uhc plots
design.mat.test.full <- model.matrix(form2a, data=testing_twz)

z <- model.matrix(~ delta_t_z + u925_z + v925_z -1, 
                  data = testing_twz)

xchoice.full <- uhcsimstrat(nsims = 2000,
                            xmat = design.mat.test.full, 
                            stratum = testing_twz$obs_id, 
                            fit_ssf = ssf_train_full,
                            z = z)    
denshats.delta_t.full <- uhcdenscalc(rand_sims = xchoice.full[,,1], 
                                     dat = z[testing_twz$used==1,1], 
                                     avail = z[testing_twz$used==0,1]) 
denshats.u925.full <- uhcdenscalc(rand_sims=xchoice.full[,,2], 
                                  dat=z[testing_twz$used==1,2], 
                                  avail=z[testing_twz$used==0,2])  
denshats.v925.full <- uhcdenscalc(rand_sims=xchoice.full[,,3], 
                                  dat=z[testing_twz$used==1,3], 
                                  avail=z[testing_twz$used==0,3]) 


X11()
par(mfrow=c(1,3), mar=c(4,2,2,2), oma=c(3, 4, 7, 0), bty="L")

# First, the density plots
uhcdensplot(densdat = denshats.delta_t.full$densdat, 
            densrand = denshats.delta_t.full$densrand, 
            includeAvail = TRUE, 
            densavail = denshats.delta_t.full$densavail) 
mtext(outer=F, side=2, line=3, "Density")
mtext(outer=F, side=3, line=1, "delta_t", cex=1)
uhcdensplot(densdat = denshats.u925.full$densdat,
            densrand = denshats.u925.full$densrand,
            includeAvail = TRUE, 
            densavail = denshats.u925.full$densavail) 
mtext(outer=F, side=3, line=1, "u_wind", cex=1)
uhcdensplot(densdat = denshats.v925.full$densdat, 
            densrand = denshats.v925.full$densrand,
            includeAvail = TRUE, 
            densavail = denshats.v925.full$densavail) 
mtext(outer=F, side=3, line=1, "v_wind", cex=1)
mtext(outer=T, side=3, line=3,  textplot1, cex=1)

#####
#temperate zone
partition_tmpz <- ann_z %>% 
  filter(zone == "temperate") %>% 
  distinct(obs_id) %>% 
  sample_frac(0.2) %>% 
  pull(obs_id)

testing_tmpz <- ann_z %>%
  filter(zone == "temperate" & obs_id %in% partition_tmpz)

training_tmpz <- ann_z %>% 
  filter(zone == "temperate" & !(obs_id %in% partition_tmpz))

#formula1_sc <- used ~ delta_t + u925 + v925 + strata(obs_id)
textplot1 <- (expression(y %~% delta_t + u925 + v925 + strata(obs_id)))
form1a <- (used ~ delta_t_z + u925_z + v925_z  + strata(obs_id))
form2a <- ~ delta_t_z + u925_z + v925_z -1

#train full model
ssf_train_full <- clogit(form1a, data = ann_z[ann_z$zone == "temperate",])
summary(ssf_train_full)

#uhc plots
design.mat.test.full <- model.matrix(form2a, data=testing_tmpz)

z <- model.matrix(~ delta_t_z + u925_z + v925_z -1, 
                  data = testing_tmpz)

xchoice.full <- uhcsimstrat(nsims = 2000,
                            xmat = design.mat.test.full, 
                            stratum = testing_tmpz$obs_id, 
                            fit_ssf = ssf_train_full,
                            z = z)    
denshats.delta_t.full <- uhcdenscalc(rand_sims = xchoice.full[,,1], 
                                     dat = z[testing_tmpz$used==1,1], 
                                     avail = z[testing_tmpz$used==0,1]) 
denshats.u925.full <- uhcdenscalc(rand_sims=xchoice.full[,,2], 
                                  dat=z[testing_tmpz$used==1,2], 
                                  avail=z[testing_tmpz$used==0,2])  
denshats.v925.full <- uhcdenscalc(rand_sims=xchoice.full[,,3], 
                                  dat=z[testing_tmpz$used==1,3], 
                                  avail=z[testing_tmpz$used==0,3]) 


X11()
par(mfrow=c(1,3), mar=c(4,2,2,2), oma=c(3, 4, 7, 0), bty="L")

uhcdensplot(densdat = denshats.delta_t.full$densdat, 
            densrand = denshats.delta_t.full$densrand, 
            includeAvail = TRUE, 
            densavail = denshats.delta_t.full$densavail) 
mtext(outer=F, side=2, line=3, "Density")
mtext(outer=F, side=3, line=1, "delta_t", cex=1)
uhcdensplot(densdat = denshats.u925.full$densdat,
            densrand = denshats.u925.full$densrand,
            includeAvail = TRUE, 
            densavail = denshats.u925.full$densavail) 
mtext(outer=F, side=3, line=1, "u_wind", cex=1)
uhcdensplot(densdat = denshats.v925.full$densdat, 
            densrand = denshats.v925.full$densrand,
            includeAvail = TRUE, 
            densavail = denshats.v925.full$densavail) 
mtext(outer=F, side=3, line=1, "v_wind", cex=1)
mtext(outer=T, side=3, line=3,  textplot1, cex=1)
