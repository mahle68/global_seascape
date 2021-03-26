#script for calculating the correlation between delta_t and w_star (using Gil's suggestion)

#Mar 22. 2021, Radolfzell, DE
#Elham Nourani

library(tidyverse)
library(sp)
library(sf)
library(move)
library(scales)

setwd("/home/enourani/ownCloud/Work/Projects/delta_t")
source("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/wind_support_Kami.R")

wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

#open annotated sea-crossing points (prepped in 2021_wind_delta_t_maps.R)
load("R_files/2021/raw_sea_points_annotated.RData") #called ann



#with ECMWF ERA-interim data, units for "ECMWF Interim Full Daily SFC-FC Instantaneous Moisture Flux": kg m^-2 s^-1
#units for "ECMWF Interim Full Daily SFC-FC Instantaneous Surface Heat Flux": J m^-2

#s_flux: "ECMWF Interim Full Daily SFC-FC Instantaneous Surface Heat Flux" (Env-data give the wrong units of J m^-. the actual units are w/m^2)
#m_flux: "ECMWF Interim Full Daily SFC-FC Instantaneous Moisture Flux" (kg m^-2 s^-1)
#blh: boundary layer height (m)
#t2m: 2m temperature

#msshf <- ann$s_flux / 3*3600 #this is now in W/m^2

msshf <- ann$s_flux
mslhf <- ann$m_flux
z <- ann$blh #this is in m
T2m <- ann$t2m #in K
g <- 9.81
T2m <- (ann$t2m - 273.15) * 0.006 * ann$blh

wt <- msshf/1013/1.2 #surface kinetic heat flux
wq <- (mslhf * 1000) /1.2 #surface water vapor flux

wthetav <- wt + (0.61*T*(wq))

ann$w_star <- ((g*z*wthetav)/T2m)^(1/3)




###### Gil's code


w_star <- function(g = 9.81, blh, T2m, s_flux, m_flux) {
  
  z <- blh
  T_k_2m <- T2m
  T_c_2m <- T_k_2m - 273.15
  Thetav_k_z <- (T_k_2m) + 0.006 * z
  wT <- s_flux / 1013 / 1.2
  wq <- m_flux *1000 /1.2
  
  wthetav <- wT + 0.61 * T_c_2m * wq
  
  w_star <- (g*z*(wthetav/Thetav_k_z)) ^ 1/3
  
  return(w_star)
  
}


ann$w_star <- w_star(blh = ann$blh, T2m = ann$t2m, 
                     s_flux = ann$s_flux, m_flux = ann$m_flux)






############
z <- ann$blh
T_k_2m <- ann$t2m
T_c_2m <- T_k_2m - 273.15
Thetav_k_z <- T_k_2m + 0.006 * z
wT <- ann$s_flux / 1013 / 1.2
wq <- ann$m_flux *1000 /1.2

wthetav <- wT + 0.61 * T_c_2m * wq

w_star <- (g*z*(wthetav/Thetav_k_z)) ^ (1/3)



#############
g <- 9.81
z <- ann$blh
T_k_2m <- ann$t2m
T_c_2m <- T_k_2m - 273.15
Thetav_k_z <- (T_k_2m) + 0.006*z
wT <- ann$s_flux/1013/1.2
wq <- ann$m_flux*1000 /1.2

wthetav <- wT + 0.61*T_c_2m*wq

ann$w_star <- (g*z*(wthetav/Thetav_k_z)) ^ (1/3)



w_star <- (g*z*wthetav/Thetav_k_z) ^ (1/3)

wthetv <- wT + (0.61*T_c_2m*wq)



wT <- s_flux/1013/1.2
wq <- m_flux * 1000 /1.2
z <- blh
T_c_2m <- T_k_2m -273.15
Thetav_k_z <- (T_k_2m) + 0.006*z

###

W*= [g*z*(w’thetav’)/(Thetav_k|z)] ^ 1/3
(w’thetv’) = (w’T’) + 0.61*(T_c|2m)*(w’q’)
(Thetav_k|z) = (T_k|2) + 0.006*z
(T_c|2m) = (T_k|2m) -273.15

w’T’ [K*m/s] = s_flux [W/m^2] / 1013 [J/kg/K] / 1.2 [kg/m^3]
w’q’ [g/kg*m/s]= m_flux [kg/s/m^2] *1000 [g/kg] /1.2 [kg/m^3]
z=blh [m]


#ann$w_star <- (9.81*z*(msshf/1013/1.2+0.61*(T2m-273.15)* 0.006*z*mslhf/1.2* 1000) / (T2m*0.006*z))^ 1/3
  
  
postive_dt <- ann %>% 
  filter(delta_t > 0)

#plot
fit <- lm(w_star~ log(delta_t+1),data = postive_dt) 

summary(fit)

with(postive_dt,plot(log(delta_t+1), w_star))
abline(fit)


#plot
fit <- lm(w_star~delta_t,data = postive_dt) 

summary(fit)

with(postive_dt,plot(delta_t, w_star))
abline(fit)

