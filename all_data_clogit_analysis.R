#scripts for meta-analysis of sea-crossing behavior
#follows the data preparation in all_data_prep.R

library(lme4)
library(survival)
library(TwoStepCLogit)
library(rstanarm)
library(INLA)
library(TwoStepCLogit)
library(jtools) #summ ftn
#library(ggstance)
library(MuMIn) #adjusted r squared


##### STEP 7: analysis #####

load("R_files/all_spp_temp_sp_filtered_15km_alt_ann_14days.RData") #called ann

lapply(split(ann,ann$zone), function(x){
  X11()
  par(mfrow = c(1,3))
  hist(ann[ann])
})

formula1 <- used ~ scale(delta_t) + scale(u925) + scale(v925) + strata(obs_id)

m_twz1 <- clogit(formula1, data = ann[ann$zone == "tradewind",])
summary(m_twz1)

m_twz2 <- update(m_twz1,- scale(v925))

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


ann_sc <- ann %>%
  group_by(zone) %>% 
  mutate_at(c("delta_t","u925","v925"), scale) %>% 
  mutate(species1 = species,
         species2 = species) %>% 
  as.data.frame()

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
