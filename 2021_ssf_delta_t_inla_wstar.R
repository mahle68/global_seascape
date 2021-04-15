#script to estimate the step selection function for water-crossing raptors.
#mostly copied from ssf_delta_t_inla2.R. For previous steps, see 2021_all_data_prep_analyze.R
#Feb 22. 2021. Radolfzell, Germany. Elham Nourani, PhD.
#data prepped in 2021_ssf_delta_t_inla.R

library(tidyverse)
library(move)
library(sf)
library(circular)
library(CircStats)
library(fitdistrplus)
library(RNCEP)
library(lubridate)
library(mapview)
library(parallel)
library(tidyr)
library(corrr)
library(lme4)
library(MuMIn)
library(mgcv)
library(survival)
library(INLA)
library(ggregplot) #devtools::install_github("gfalbery/ggregplot")
library(maptools)
library(brinla)
library(INLAutils) #library(devtools); install_github('timcdlucas/INLAutils')

setwd("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/")

#meters_proj <- CRS("+proj=moll +ellps=WGS84")
wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

NCEP.loxodrome.na <- function (lat1, lat2, lon1, lon2) {
  deg2rad <- pi/180
  acot <- function(x) {
    return(atan(1/x))
  }
  lat1 <- deg2rad * lat1
  lat2 <- deg2rad * lat2
  lon1 <- deg2rad * lon1
  lon2 <- deg2rad * lon2
  deltaLon <- lon2 - lon1
  pi4 <- pi/4
  Sig1 <- log(tan(pi4 + lat1/2))
  Sig2 <- log(tan(pi4 + lat2/2))
  deltaSig <- Sig2 - Sig1
  if (deltaLon == 0 && deltaSig > 0) {
    head <- 0
  }
  else if (deltaLon == 0 && deltaSig < 0) {
    head <- 180
  }
  else if (deltaSig == 0 && deltaLon > 0) {
    head <- 90
  }
  else if (deltaSig == 0 && deltaLon < 0) {
    head <- 270
  }
  else if (deltaSig < 0 && deltaLon < 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi + 180
  }
  else if (deltaSig < 0 && deltaLon > 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi + 180
  }
  else if (deltaSig > 0 && deltaLon > 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi
  }
  else if (deltaSig > 0 && deltaLon < 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi + 360
  }
  else {
    head <-NA}
  return(head)
}

source("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/wind_support_Kami.R")

rsd <- function(x){
  cv <- sd(x, na.rm = T)/abs(mean(x, na.rm = T))
  rsd <- cv*100
  return(rsd)
}

CubeRoot<-function(x){
  sign(x)*abs(x)^(1/3)
}

w_star <- function(g = 9.81, blh, T2m, s_flux, m_flux) {
  
  z <- blh
  T_k_2m <- T2m
  T_c_2m <- T_k_2m - 273.15
  Thetav_k_z <- (T_k_2m) + 0.006 * z
  wT <- (s_flux * -1) / 1013 / 1.2 #reverse the sign. ECMWF upward fluxes are negative
  wq <- (m_flux * -1) *1000 /1.2 #reverse the sign. ECMWF upward fluxes are negative
  
  wthetav <- wT + 0.61 * T_c_2m * wq
  
  #w_star <- as.complex(g*z*(wthetav/Thetav_k_z)) ^ (1/3)
  w_star <- CubeRoot(g*z*(wthetav/Thetav_k_z))
  
  return(w_star)
  
}


# --- 50 alt. data
load("2021/ssf_input_annotated_60_30_50_updated_sample.RData") #ann_50_sample

colnames(ann_50_sample)[c(28,29)] <- c("location-long","location-lat")
ann_50_sample$timestamp <- paste(as.character(ann_50_sample$date_time),"000",sep = ".")


write.csv(ann_50_sample, file = "2021/annotate_for_wstar.csv")


#open annotated data
ann <- read.csv("2021/annotations/ssf_with_wstar/annotate_for_wstar.csv-6423878732846278525/annotate_for_wstar.csv-6423878732846278525.csv",
                stringsAsFactors = F) %>% 
  mutate(blh = ECMWF.Interim.Full.Daily.SFC.FC.Boundary.Layer.Height,
         s_flux = ECMWF.Interim.Full.Daily.SFC.FC.Instantaneous.Surface.Heat.Flux,
         m_flux = ECMWF.Interim.Full.Daily.SFC.FC.Instantaneous.Moisture.Flux) %>% 
  mutate(wstar = w_star(blh = blh, T2m = t2m, s_flux = s_flux, m_flux = m_flux))

# long-term annotation data (40 year data)
#prep a dataframe with 40 rows corresponding to 40 years (1981,2020), for each point. then i can calculate variance of delta t over 41 years for each point
df_40 <- ann_50_sample %>% #make sure this is not grouped!
  dplyr::select(-c(v925,u925,t2m,sst,delta_t)) %>% 
  slice(rep(row_number(),40)) %>% 
  group_by(row_id) %>% 
  mutate(year = c(1981:2020)) %>%
  ungroup() %>%
  mutate(timestamp = paste(as.character(date_time),"000",sep = ".")) %>% 
  as.data.frame()

str_sub(df_40$timestamp,1,4) <- df_40$year #replace original year with years from 1979-2019
colnames(df_40)[c(23,24)] <- c("location-long","location-lat") #rename columns to match movebank format

#break up into two parts. over 1 million rows
df_40_1 <- df_40 %>% 
  slice(1:999999)
write.csv(df_40_1, "2021/ssf_40_all_spp_60_30_1_updated.csv")

df_40_2 <- df_40 %>% 
  slice(1000000:1999999)
write.csv(df_40_2, "2021/ssf_40_all_spp_60_30_2_updated.csv")

df_40_3 <- df_40 %>% 
  slice(2000000 : nrow(df_40))
write.csv(df_40_3, "2021/ssf_40_all_spp_60_30_3_updated.csv")


#---- after movebank

load("2021/ssf_input_annotated_60_30_50_updated_sample.RData") #ann_50_sample

ann_50_sample <- ann_50_sample %>% 
  mutate(abs_cross_wind = abs(cross_wind(u = u925, v = v925, heading = heading)))

#calculate long-term metrics and merge with previously annotated data
ann_40_ls <- list.files("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/2021/annotations/40_yrs_60_30_updated/",pattern = ".csv", recursive = T,full.names = T) 

ann_cmpl <- lapply(ann_40_ls, read.csv, stringsAsFactors = F) %>% 
  reduce(full_join) %>% 
  rename(sst = ECMWF.ERA5.SL.Sea.Surface.Temperature,
         t2m = ECMWF.ERA5.SL.Temperature..2.m.above.Ground.,
         u925 = ECMWF.ERA5.PL.U.Wind,
         v925 = ECMWF.ERA5.PL.V.wind) %>% 
  mutate(delta_t = sst - t2m,
         wind_support= wind_support(u = u925, v = v925, heading = heading),
         cross_wind= cross_wind(u = u925, v = v925, heading = heading),
         abs_cross_wind = abs(cross_wind(u = u925, v = v925, heading = heading)),
         wind_speed = sqrt(u925^2 + v925^2)) %>% 
  group_by(row_id) %>% 
  summarise_at(c("delta_t", "wind_speed", "wind_support", "abs_cross_wind", "u925", "v925"), #before calculating these, investigate why/if we have NAs??
               list(avg = ~mean(., na.rm = T), var = ~var(., na.rm = T), rsd = ~rsd(.))) %>% 
  ungroup() %>% 
  full_join(ann_50_sample, by = "row_id") %>% 
  rename(location.long = coords.x1,
         location.lat = coords.x2) %>% 
  rowwise() %>% 
  mutate(species = strsplit(group, "_")[[1]][1],
         zone = ifelse(between(location.lat, 0, 30) | between(location.lat, 0, -30), "tradewind",
                       ifelse(between(location.lat, 30,60) | between(location.lat, -30,-60), "temperate",
                              ifelse(between(location.lat, -30,30), "tropical",
                                     "arctic")))) %>% 
  ungroup() %>% 
  as.data.frame()

save(ann_cmpl, file = "2021/ssf_input_ann_cmpl_60_30_updated.RData")


#make sure O_A and EF_S are sampled (10 ind each)
ann_cmpl %>% 
  group_by(group) %>% 
  summarise(n_str = n_distinct(stratum),
            n_ind = n_distinct(ind))


# ---------- STEP 6: clogit to get a feel for things #####

load("2021/ssf_input_ann_60_30_updated_z.RData") #all_data

#### exploration. run a quick clogit to see if the results are what I expect

form_original <- formula(used ~  delta_t_z * wind_speed_z + wind_support_z +
                    strata(stratum))

form_1 <- formula(used ~ lat_at_used * delta_t_z + lat_at_used * wind_support_z +
                      strata(stratum))

form_2 <- formula(used ~ lat_at_used * delta_t_z +  lat_at_used * wind_support_z +
                    strata(stratum))

form_3 <- formula(used ~  delta_t_z * wind_speed_z + delta_t_z * wind_support_z +
                           strata(stratum))

form_4 <- formula(used ~  delta_t_z * wind_support_z +
                    strata(stratum))

form_5 <- formula(used ~  delta_t_z * lat_at_used + wind_support_z * lat_at_used +
                    strata(stratum))

form_6 <- formula(used ~  delta_t_z * wind_support_z + 
                    strata(stratum))

form_7 <- formula(used ~  delta_t_z *  wind_support_z +
                    strata(stratum))

form_8 <- formula(used ~  delta_t_z * wind_support_z + delta_t_z * abs_cross_wind +
                             strata(stratum))

form_9 <- formula(used ~ delta_t_z * wind_support_z + s_elev_angle +
                    strata(stratum))

m1 <- clogit(form_1, data = all_data) #when zone is added, delta t coeff becomes positive. still not sig, but positive. but wind support is negative. but interaction of wind support with zones is positive. so, only negative in the arctic
m1a <- clogit(form_original, data = all_data)
m1b <- clogit(form_2, data = all_data)
m1c <- clogit(form_3, data= all_data)
m1d <- clogit(form_4, data= all_data)
m1e <- clogit(form_5, data = all_data)
m1f <- clogit(form_6, data = all_data)
m1g <- clogit(form_7, data = all_data)
m1h <- clogit(form_8, data = all_data)
m1i <- clogit(form_9, data = all_data)

#only during the day:
day <- all_data[all_data$sun_elev != "night",]
m2c <- clogit(form_3, data = day)



###zone_specific... conclusion: the order of importance is pretty much the same. only the direction of delta_t:windspeed changes... but still small. delta t is negative all through (not sig. in tw)

twtmp <- clogit(form_original, data = all_data[all_data$lat_at_used %in% c("temperate", "tradewind"),]) #similar to all_data
tmp <- clogit(form_original, data = all_data[all_data$lat_at_used == "temperate",])
tw <- clogit(form_original, data = all_data[all_data$lat_at_used == "tradewind",])
arctr <- clogit(form_original, data = all_data[all_data$lat_at_used %in% c("tropical", "arctic"),])

#---
#extract number of strata per zone and consider removing the arctic

all_data %>% 
  group_by (zone) %>% 
  summarise(n_s = n_distinct(stratum),
            n_t = n_distinct(track),
            n_i = n_distinct(ind)) #the arctic only has 26 strata from 4 tracks of 3 individuals.... let's remove it!


# ---------- STEP 7: INLA #####

load("2021/ssf_input_ann_cmpl_60_30_updated.RData") #ann_cmpl

#correlation
ann_cmpl %>% 
  dplyr::select(c("delta_t", "wind_speed", "wind_support", "wind_support_var", "abs_cross_wind", "delta_t_var","step_length")) %>% 
  correlate() %>% 
  stretch() %>% 
  filter(abs(r) > 0.6) #correlated: var_cw with location.lat and var_delta_t with location.lat. avg delta_t and delta_t. avg_ws and var_delta_t
#correlated: wind support var & wind speed var and cross wind var, crosswind var and wind speed var.  delta-t var and wind support var!!! wind speed and abs_crosswind.
#with the sampled dataset, variance of delta t is correlated with delta t and variance of wind support


#z-transform
all_data <- ann_cmpl %>% 
  #group_by(species) 
  mutate_at(c("delta_t", "wind_speed", "wind_support", "wind_support_var", "abs_cross_wind", "delta_t_var"),
            #list(z = ~scale(.)))
            list(z = ~as.numeric(scale(.)))) %>%
  arrange(stratum, desc(used)) %>% 
  group_by(stratum) %>%  
  mutate(lat_at_used = head(zone,1)) %>%  #add a variable for latitudinal zone. This will assign the lat zone of the used point to the entire stratum
  ungroup() 


save(all_data, file = "2021/ssf_input_ann_60_30_updated_z.RData")


#check to make sure each stratum has one lat zone value
all_data %>% 
  group_by(stratum) %>% 
  summarize(n = n_distinct(lat_at_used),
            n_z = n_distinct(zone)) %>% 
  filter(n > 1) #this should be zero


load("2021/ssf_input_ann_60_30_updated_z.RData") #all_data


#repeat variabels that will be used as random slopes
all_data <- all_data %>% 
  mutate(species1 = factor(species),
         species2 = factor(species),
         species3 = factor(species),
         species4 = factor(species),
         species5 = factor(species),
         species6 = factor(species),
         ind1 = factor(ind),
         ind2 = factor(ind),
         ind3 = factor(ind),
         ind4 = factor(ind),
         ind5 = factor(ind),
         ind6 = factor(ind),
         zone1 = factor(lat_at_used), #latitude is correlated with variance values. so, try to use this instead of var variables.. messed everything up...
         zone2 = factor(lat_at_used),
         zone3 = factor(lat_at_used),
         zone4 = factor(lat_at_used),
         zone5 = factor(lat_at_used),
         zone6 = factor(lat_at_used),
         stratum = factor(stratum)) #%>% 
  #dplyr::select(c("used","stratum","delta_t_z","wind_speed_z","wind_support_z","wind_support_var_z", "abs_cross_wind_z","delta_t_var_z",
  #                "species1","species2", "species3", "species4","species5","ind1", "ind2", "ind3", "ind4", "ind5", #"zone1", "zone2","zone3","zone4","zone5","zone6",
  #                "location.lat"))


# Set mean and precision for the priors of slope coefficients
mean.beta <- 0
prec.beta <- 1e-4 

formulaM1 <- used ~ -1 + delta_t_z * wind_support_z + wind_speed_z + wind_support_var_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species1, delta_t_z, model = "iid", 
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species3, wind_speed_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, wind_support_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind1, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ind2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind3, wind_speed_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind4, wind_support_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

(b <- Sys.time())
M1 <- inla(formula = formulaM1, family ="Poisson",  
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           control.inla = list(force.diagonal = T),
           data = all_data,
           num.threads = 10,
           control.compute = list(openmp.strategy = "huge", config = TRUE, mlik = T, waic = T, cpo = F))

Sys.time() - b #2.3 min

save(M1, file = "2021/inla_models/m1_60_30.RData")


summary(M1)#


#remove wind speed
formulaM2 <- used ~ -1 + delta_t_z * wind_support_z + wind_support_var_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species1, delta_t_z, model = "iid", 
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, wind_support_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind1, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ind2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind4, wind_support_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

(b <- Sys.time())
M2 <- inla(formulaM2, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           control.inla = list(force.diagonal = T),
           data = all_data,
           num.threads = 10,
           control.compute = list(openmp.strategy = "huge", config = TRUE, mlik = T, waic = T, cpo = F))
Sys.time() - b #1.3 min



save(M2, file = "2021/inla_models/m2_60_30.RData")


  
#remove wind speed and delta t
formulaM3 <- used ~ -1 + delta_t_z * wind_support_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species1, delta_t_z, model = "iid", 
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind1, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ind4, wind_support_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))


(b <- Sys.time())
M3 <- inla(formulaM3, family ="Poisson", 
            control.fixed = list(
              mean = mean.beta,
              prec = list(default = prec.beta)),
           control.inla = list(force.diagonal = T),
            data = all_data,
            num.threads = 10,
            control.compute = list(openmp.strategy = "huge", config = TRUE, mlik = T, waic = T, cpo = F))
Sys.time() - b #37 sec


save(M3, file = "2021/inla_models/m3_60_30.RData")



# ---------- STEP 7: plots #####

#load models
load("2021/inla_models/m1_60_30.RData")
load("2021/inla_models/m2_60_30.RData")
load("2021/inla_models/m3_60_30.RData")

#FIGURE 2: posterior means of fixed effects 
#easy
Efxplot(list(M1,M2,M3))

#sophisticated
ModelList <- list(M1,M2,M3)
graphlist<-list()
for(i in 1:length(ModelList)){
  model<-ModelList[[i]]
  
  graph<-as.data.frame(summary(model)$fixed)
  colnames(graph)[which(colnames(graph)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
  colnames(graph)[which(colnames(graph)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
  colnames(graph)[which(colnames(graph)%in%c("mean"))]<-c("Estimate")
  
  graph$Model<-i
  graph$Factor<-rownames(graph)
  
  graphlist[[i]]<-graph
}

graph <- bind_rows(graphlist)

graph$Sig <- with(graph, ifelse(Lower*Upper>0, "*", ""))

graph$Model <- as.factor(graph$Model)

position <- ifelse(length(unique(graph$Model))  ==  1, "none", "right")

VarOrder <- rev(unique(graph$Factor))
VarNames <- VarOrder

graph$Factor <- factor(graph$Factor, levels = VarOrder)
levels(graph$Factor) <- VarNames

min<-min(graph$Lower,na.rm = T)
max<-max(graph$Upper,na.rm = T)

graph$Factor_n <- as.numeric(graph$Factor)

X11(width = 4.1, height = 2.7)

pdf("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/2021/coefficients_updated.pdf", width = 4.1, height = 2.7)
jpeg("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/2021/coefficients_updated.jpeg", width = 4.1, height = 2.7, units = "in", res = 300)

par(mfrow=c(1,1), bty="n", #no box around the plot
    #cex.axis= 0.75, #x and y labels have 0.75% of the default size
    #font.axis= 0.75, #3: axis labels are in italics
    #cex.lab = 0.75,
    cex = 0.7,
    oma = c(0,3.7,0,0),
    mar = c(3, 4.1, 0.5, 1),
    bty = "l"
)

plot(0, type = "n", labels = FALSE, tck = 0, xlim = c(-3,5.5), ylim = c(0.2,5.3), xlab = "Estimate", ylab = "")

#add vertical line for zero
abline(v = 0, col = "grey30",lty = 2)
#add points and error bars
points(graph[graph$Model == 1, "Estimate"], graph[graph$Model == 1,"Factor_n"] - 0.25, col = "steelblue1", pch = 20, cex = 1.3)
arrows(graph[graph$Model == 1, "Lower"], graph[graph$Model == 1,"Factor_n"] - 0.25,
       graph[graph$Model == 1, "Upper"], graph[graph$Model == 1,"Factor_n"] - 0.25,
       col = "steelblue1", code = 3, length = 0.03, angle = 90) #angle of 90 to make the arrow head as straight as a line

points(graph[graph$Model == 2, c("Estimate","Factor")], col = "royalblue1", pch = 20, cex = 1.3)
arrows(graph[graph$Model == 2, "Lower"], graph[graph$Model == 2,"Factor_n"],
       graph[graph$Model == 2, "Upper"], graph[graph$Model == 2,"Factor_n"],
       col = "royalblue1", code = 3, length = 0.03, angle = 90)

points(graph[graph$Model == 3, "Estimate"], graph[graph$Model == 3,"Factor_n"] + 0.25, col = "mediumblue", pch =20, cex = 1.3)
arrows(graph[graph$Model == 3, "Lower"], graph[graph$Model == 3,"Factor_n"] + 0.25,
       graph[graph$Model == 3, "Upper"], graph[graph$Model == 3,"Factor_n"] + 0.25,
       col = "mediumblue", code = 3, length = 0.03, angle = 90)
#add axes
axis(side= 1, at= c(-2,0,2,4), labels= c("-2", "0", "2", "4"), 
     tick=T ,col = NA, col.ticks = 1, tck=-.015)

axis(side= 2, at= c(1:5), #line=-4.8, 
     labels = c(expression(paste(italic(paste(Delta,"T"))," : Wind support")),
                "Wind support var","Wind speed", "Wind support", expression(italic(paste(Delta,"T")))),
     tick=T ,col = NA, col.ticks = 1, # NULL would mean to use the defult color specified by "fg" in par
     tck=-.015 , #tick marks smaller than default by this proportion
     las=2) # text perpendicular to axis label 

#add legend
legend(x = 3.5, y = 1.3, legend=c("Model 3","Model 2", "Model 1"), col = c("mediumblue","royalblue1","steelblue1"), #coords indicate top-left
       pch = 19, bg="white",bty="n", cex = 0.9)

dev.off()

#SUPPLEMENTARY FIGURE 1: species-specific coefficients 
#for the best model (M3); original code by Virgilio Gomez-Rubio (Bayesian inference with INLA, 2020)
#species
species_names <- unique(all_data$species)

tab_dt <- data.frame(ID = as.factor(M3$summary.random$species1$ID),
                     mean = M3$summary.random$species1$mean,
                     IClower = M3$summary.random$species1[, 4],
                     ICupper = M3$summary.random$species1[, 6])


tab_wspt <- data.frame(ID = as.factor(M3$summary.random$species2$ID),
                       mean = M3$summary.random$species2$mean,
                       IClower = M3$summary.random$species2[, 4],
                       ICupper = M3$summary.random$species2[, 6])


X11(width = 4, height = 3)

pdf("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/2021/species_var_updated.pdf", width = 4, height = 3)

par(mfrow = c(1,1), bty="n", #no box around the plot
    #cex.axis= 0.75, #x and y labels have 0.75% of the default size
    #font.axis= 0.75, #3: axis labels are in italics
    #cex.lab = 0.75,
    cex = 0.7,
    oma = c(0,3.5,0,0),
    mar = c(3, 2, 0.5, 1),
    bty = "l"
)


plot(0, type = "n", labels = FALSE, tck = 0, xlim = c(-3,6), ylim = c(0,4.5), xlab = "", ylab = "")
#add vertical line for zero
abline(v = 0, col = "grey30",lty = 2)

points(tab_dt$mean, as.numeric(tab_dt$ID) - 0.2, col = "darkgoldenrod2", pch = 19, cex = 1.3)
arrows(tab_dt$IClower, as.numeric(tab_dt$ID) - 0.2,
       tab_dt$ICupper, as.numeric(tab_dt$ID) - 0.2,
       col = "darkgoldenrod2", code = 3, length = 0.03, angle = 90) #angle of 90 to make the arrow head as straight as a line

points(tab_wspt$mean, as.numeric(tab_wspt$ID) + 0.2, col = "cornflowerblue", pch = 19, cex = 1.3)
arrows(tab_wspt$IClower, as.numeric(tab_wspt$ID) + 0.2,
       tab_wspt$ICupper, as.numeric(tab_wspt$ID) + 0.2,
       col = "cornflowerblue", code = 3, length = 0.03, angle = 90) #angle of 90 to make the arrow head as straight as a line

axis(side= 1, at= c(-4,-2,0,2,4), labels= c("-4","-2", "0", "2","4"), 
     tick=T ,col = NA, col.ticks = 1, tck=-.015)

axis(side= 2, at= c(1:4), #line = 6, 
     labels = tab_dt$ID,# c( "Eleonora's falcon","Osprey", "Oriental\n honey buzzard", "Peregrine falcon"), #same order as tab_dt$ID
     tick = T ,col = NA, col.ticks = 1, # NULL would mean to use the defult color specified by "fg" in par
     tck = -.015 , #tick marks smaller than default by this proportion
     las = 2) # text perpendicular to axis label 

#add legend
legend(x = 3.6 , y = 0.7, legend = c("Wind support", expression(italic(paste(Delta,"T")))), 
       col = c("cornflowerblue","darkgoldenrod1"), #coords indicate top-left
       pch = 19, bg="white",bty="n", cex = 0.9)

dev.off()

#SUPPLEMENTARY FIGURE 2: boxplots 


X11(width = 9, height = 7)

pdf("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/2021/boxplots_updated.pdf", width = 9, height = 7)

par(mfrow= c(2,3), 
    oma = c(0,0,3,0), 
    las = 1)

labels <- c(expression(italic(paste(Delta,"T"))), "Wind support", "Wind speed")
variables <- c("delta_t", "wind_support", "wind_speed")
v_variables <- c("delta_t_var", "wind_support_var","wind_speed_var")

for(i in 1:length(variables)){
  
  boxplot(ann_cmpl[,variables[i]] ~ ann_cmpl[,"species"], data = ann_cmpl, boxfill = NA, border = NA, main = labels[i], xlab = "", ylab = "")
  if(i == 1){
    legend("topleft", legend = c("used","available"), fill = c("orange","gray"), bty = "n")
  }
  boxplot(ann_cmpl[ann_cmpl$used == 1, variables[i]] ~ ann_cmpl[ann_cmpl$used == 1,"species"], 
          yaxt = "n", xaxt = "n", add = T, boxfill = "orange",
          boxwex = 0.25, at = 1:length(unique(ann_cmpl[ann_cmpl$used == 1, "species"])) - 0.15)
  boxplot(ann_cmpl[ann_cmpl$used == 0, variables[i]] ~ ann_cmpl[ann_cmpl$used == 0, "species"], 
          yaxt = "n", xaxt = "n", add = T, boxfill = "grey",
          boxwex = 0.25, at = 1:length(unique(ann_cmpl[ann_cmpl$used == 1 , "species"])) + 0.15)
  
}
mtext("Instantaneous values at each step", side = 3, outer = T, cex = 1.3)

for(i in 1:length(v_variables)){
  
  boxplot(ann_cmpl[,v_variables[i]] ~ ann_cmpl[,"species"], data = ann_cmpl, boxfill = NA, border = NA, main = labels[i], xlab = "", ylab = "")
  if(i == 1){
    legend("topleft", legend = c("used","available"), fill = c("orange","gray"), bty = "n")
  }
  boxplot(ann_cmpl[ann_cmpl$used == 1,v_variables[i]] ~ ann_cmpl[ann_cmpl$used == 1,"species"], 
          yaxt = "n",xaxt = "n", add = T, boxfill = "orange",
          boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) - 0.15)
  boxplot(ann_cmpl[ann_cmpl$used == 0,v_variables[i]] ~ ann_cmpl[ann_cmpl$used == 0,"species"], 
          yaxt = "n",xaxt = "n", add = T, boxfill = "grey",
          boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) + 0.15)
} 

mtext("40-year variances at each step", side = 3, outer = T, cex = 1.3, line = -25)

dev.off()
