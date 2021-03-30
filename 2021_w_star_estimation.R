#script for calculating the correlation between delta_t and w_star (using Gil's suggestion)

#Mar 22. 2021, Radolfzell, DE
#Elham Nourani

library(tidyverse)
library(sp)
library(sf)
library(move)
library(scales)
library(maptools)

setwd("/home/enourani/ownCloud/Work/Projects/delta_t")
source("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/wind_support_Kami.R")

wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

#annotated data.(prepped in 2021_wind_delta_t_maps.R)
ann <- read.csv("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/2021/annotations/raw_points_for_maps/interim/raw_points_for_maps.csv-3930001190922922558.csv") %>% 
  mutate(timestamp,timestamp = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>%
  rename(sst = ECMWF.Interim.Full.Daily.SFC.FC.Sea.Surface.Temperature,
         t2m = ECMWF.Interim.Full.Daily.SFC.FC.Temperature..2.m.above.Ground.,
         blh = ECMWF.Interim.Full.Daily.SFC.FC.Boundary.Layer.Height,
         s_flux = ECMWF.Interim.Full.Daily.SFC.FC.Instantaneous.Surface.Heat.Flux,
         m_flux = ECMWF.Interim.Full.Daily.SFC.FC.Instantaneous.Moisture.Flux) %>%
  mutate(delta_t = sst - t2m) %>%
  st_as_sf(coords = c("location.long", "location.lat"), crs = wgs) %>% 
  mutate(s_elev_angle = solarpos(st_coordinates(.), timestamp, proj4string=CRS("+proj=longlat +datum=WGS84"))[,2]) %>% #calculate solar elevation angle
  mutate(sun_elev = ifelse(s_elev_angle < -6, "night", #create a categorical variable for teh position of the sun
                           ifelse(s_elev_angle > 40, "high", "low"))) %>% 
  as("Spatial") %>% 
  as.data.frame() %>% 
  rename(lat = coords.x2,
         lon = coords.x1) %>% 
  arrange(track, timestamp)

#remove tracks with one point
more_than_one_point <- ann %>% 
  group_by(track) %>% 
  summarize(n = n()) %>% 
  filter(n > 1)

ann <- ann %>% 
  filter(track %in% more_than_one_point$track)

save(ann, file = "R_files/2021/df_for_w_star.RData")

#with ECMWF ERA-interim data, units for "ECMWF Interim Full Daily SFC-FC Instantaneous Moisture Flux": kg m^-2 s^-1
#units for "ECMWF Interim Full Daily SFC-FC Instantaneous Surface Heat Flux": J m^-2

#s_flux: "ECMWF Interim Full Daily SFC-FC Instantaneous Surface Heat Flux" (Env-data give the wrong units of J m^-. the actual units are w/m^2)
#m_flux: "ECMWF Interim Full Daily SFC-FC Instantaneous Moisture Flux" (kg m^-2 s^-1)
#blh: boundary layer height (m)
#t2m: 2m temperature

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


load("R_files/2021/df_for_w_star.RData")

pos_dt <- ann %>% 
  filter(delta_t >= 0)

pos_dt$w_star <- w_star(blh = pos_dt$blh, T2m = pos_dt$t2m, 
                     s_flux = pos_dt$s_flux, m_flux = pos_dt$m_flux)


  
postive_dt <- ann %>% 
  filter(delta_t > 0)

#plot
fit_1 <- lm(w_star~ log(delta_t+1),data = postive_dt) 

summary(fit_1)

with(postive_dt,plot(log(delta_t+1), w_star, col= as.factor(sun_elev)))
abline(fit_1)


#plot

#color
#color palette
Pal <- colorRampPalette(c("darkgoldenrod1","lightpink1", "mediumblue")) #colors for negative values
Cols <- paste0(Pal(3), "B3") #add transparency. 50% is "80". 70% is "B3". 80% is "CC". 90% is "E6"

#convert complex numbers to numerics
data <- postive_dt
data$w_star <- as.numeric(data$w_star)
data$color <- as.factor(data$sun_elev)
levels(data$color) <- Cols


fit_2 <- lm(w_star ~ delta_t, data = data) 

summary(fit_2)



#plot in base r #############

#new data for predictions
newx <- seq(min(data$delta_t), 9, by=0.05)
conf_interval <- predict(fit_2, newdata = data.frame(delta_t = newx), interval = "confidence",
                         level = 0.95)



X11(width = 4, height = 3)
pdf("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/2021/w_star.pdf", width = 4, height = 3)

par(mfrow=c(1,1), 
    bty = "l",
    cex.axis = 0.7,
    font.lab = 3,
    cex.lab = 0.9,
    mgp=c(2,0.5,0), #margin line for the axis titles
    mar = c(3.5,3.5,0.5,0.5))

plot(0, type = "n", labels = FALSE, tck = 0, xlim = c(0,8.3), ylim = c(0.3,2.5), xlab = expression(italic(paste(Delta,"T", "(°C)"))), ylab = "w* (m/s)")

with(postive_dt,points(delta_t, w_star, col= as.character(data$color), pch = 20, cex = 0.6))
clip(0, max(data$delta_t), 0, 3)
lines(newx, conf_interval[,2], col = alpha(rgb(0,0,0), 0.15), lwd = 5.5)
lines(newx, conf_interval[,3], col = alpha(rgb(0,0,0), 0.15), lty = 1, lwd = 5.5)
abline(fit_2, col = "black")

axis(side = 1, at = c(0,2,4,6,8), labels = c(0,2,4,6,8), 
     tick = T , col.ticks = 1, col = NA, tck = -.015,lwd = 0, lwd.ticks = 1)
axis(side= 2, at= c(0.5,1,1.5,2,2.5), labels= c(0.5,1,1.5,2,2.5),
     tick=T , col.ticks = 1, col = NA, tck=-.015,
     las=2) # text perpendicular to axis label 
legend(x = 5.5, y = 0.9, legend = c("daytime: high sun", "daytime: low sun", "night"), col = Cols, #coords indicate top-left
       cex = 0.7, pt.cex = 0.9, bg = "white", bty = "n", pch = 20)

text(6.5,0.95, "Time of day", cex = 0.7)

dev.off()

