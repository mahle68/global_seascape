#scripts for building a global GAM for delta-t
#data processed in R_files/global_seascape.R
#May 6, 2020. Radolfzell, DE. Elham Nourani, PhD
#update May 7, 2020: calculate solar position to add a meaningful proxy for time of the day to the model (so, basically, the local time calculation was not necessary)
#update Jun 8, 2020: scale lat, lon, and yday.... or not.... lol

#tutorials
#https://www.r-bloggers.com/advanced-sab-r-metrics-parallelization-with-the-mgcv-package/
#df equals the number of parameters needed to produce the curve: df = number of kntos -1. effective degrees of freedom, edf >= 8 means that curve is non-linear, edf= 1 is a straight line
library(mgcv)
library(parallel)
library(dplyr)
library(sf)
library(raster)
library(maps)
library(lubridate)
library(lutz)
library(maptools)
#library(multidplyr)
library(mapview)
library(parallel)
library(ggregplot)

setwd("/home/enourani/ownCloud/Work/Projects/delta_t/R_files")
wgs<-CRS("+proj=longlat +datum=WGS84 +no_defs")


### STEP 1: open data, spatio-termporal filters #####
load("processed_era_interim_data/samples_with_local_time_hour.RData") #called df_lt
X11(); plot(df_lt$lon, df_lt$lat, pch = 16, cex = 0.3, col = "blue")

#ocean layer with no lakes
ocean <- st_read("/home/enourani/ownCloud/Work/GIS_files/ne_110m_ocean/ne_110m_ocean.shp") %>% 
  slice(2) %>% #remove the caspian sea
  st_crop(xmin = -180, ymin = 0, xmax = 180, ymax = 60)

data <- df_lt %>% 
  filter(between(lat,0,60)) %>% #filter for lat zone 0-60
  as.data.frame() %>% 
  st_as_sf(coords = c("lon","lat"), crs = wgs) %>% 
  st_intersection(ocean) %>% #filter out lakes
  mutate(s_elev_angle = solarpos(st_coordinates(.), date_time, proj4string=CRS("+proj=longlat +datum=WGS84"))[,2]) %>% #calculate solar elevation angle
  mutate(sun_elev = ifelse(s_elev_angle < -6, "night", #create a categorical variable for teh position of the sun
                               ifelse(s_elev_angle > 40, "high", "low")),
         month = month(date_time))
 
save(data, file = "t_data_0_60.RData")

### STEP 2: data exploration and prep #####

load("t_data_0_60.RData")

#running this on the entire dataset gives ram errors! so, just visualize the sample! lol
sample <- data %>% 
  st_crop(xmin = 120, xmax = 130, ymin = 25, ymax = 35) %>% 
  as("Spatial") %>% 
  as.data.frame() %>% 
  rename(lon = coords.x1,
         lat = coords.x2) %>% 
  mutate(sun_elev_f = factor(sun_elev)) %>% 
  as.data.frame()

X11(width = 15, height = 10);par(mfrow= c(2,2), oma = c(0,0,3,0))
#for(i in c("sun_elev", "month","local_hour")){
    
  boxplot(sample$delta_t ~ sample$month, boxfill = "skyblue1", ylab = "delta_t")
  boxplot(sample$delta_t ~ sample$sun_elev, boxfill = "wheat1", ylab = "delta_t")
  boxplot(sample$delta_t ~ sample$local_hour, boxfill = "mediumblue", ylab = "delta_t")
  plot(sample$delta_t ~ sample$lat , pch = 16, col = factor(sample$month), ylab = "delta_t")
  
mtext("40-yr averages at each point", side = 3, outer = T, cex = 1.3)

#prep data
data_df <- data %>% 
  as("Spatial") %>% 
  as.data.frame() %>% 
  rename(lon = coords.x1,
         lat = coords.x2) %>% 
  mutate_at(c("lat","lon","yday"),
            list(z = ~scale(.))) %>% #scale lat, long, and yday
  mutate(sun_elev_f = factor(sun_elev)) %>% 
  as.data.frame() 

  

save(data_df, file = "t_data_df_0_60.RData")

### STEP 3: build the model #####
load ("t_data_df_0_60.RData")

data_df <- data_df %>% 
  mutate(delta_t_b = ifelse(delta_t >= 0, 1,0)) %>% 
  as.data.frame()
  
#keep only rows that I need (for copying to the cluster)
data_c <- data_df %>% 
  dplyr::select(c(1,4,5,7,15,17:22))

data_c$lat_z <- as.numeric(data_c$lat_z)
data_c$lon_z <- as.numeric(data_c$lon_z)
data_c$yday_z <- as.numeric(data_c$yday_z)

write.csv(data_c, "//home/enourani/ownCloud/Work/cluster_computing/global_seascape_proj/delta_t_gam/data_c.csv")


#####use the sample 
m0s <- bam(delta_t ~ s(lat, lon) + s(yday, bs ="cc") + sun_elev_f, data = sample)
AIC(m0s) #24905.22

m1s <- bam(delta_t ~ s(lat, lon, by = sun_elev_f) + s(yday, bs ="cc", by = sun_elev_f) + sun_elev_f, data = sample)
AIC(m1s) #24818.1

X11(width = 15, height = 10);par(mfrow= c(2,2), oma = c(0,0,3,0))
gam.check(m1s)

m2s <- gamm(delta_t ~ s(lon,lat, by = sun_elev_f) + s(yday, bs ="cc", by = sun_elev_f) + sun_elev_f,
            random = list(year = ~1), data = sample)
AIC(m2s$lme) #24911.63

X11(width = 15, height = 10);par(mfrow= c(2,2), oma = c(0,0,3,0))
gam.check(m2s$gam)


#separate smooths for time of day
model2 <- bam(delta_t ~ s(lon,lat) + s(yday,by = factor(sun_elev), bs = "cc") + factor(sun_elev) , data = sample)
AIC(model2) #24897.93
plot(model2)
gam.check(model2)

#separate smooths for time of day
model3 <- bam(delta_t ~ s(lon,lat, by = as.numeric(sun_elev == "night")) +
                s(lon,lat, by = as.numeric(sun_elev == "high")) +
                s(lon,lat, by = as.numeric(sun_elev == "low")) + 
                s(yday,bs="cc") + factor(sun_elev) , data = sample)
AIC(model3) #24826.79
gam.check(model3)

#random eff for elev
model4 <- gamm(delta_t ~ s(lon,lat) + 
                s(yday, bs = "cc"), method = "REML", random = list(sun_elev_f = ~1) , data = sample)
AIC(model4$lme) #24977.86

X11(width = 15, height = 10);par(mfrow= c(2,2), oma = c(0,0,3,0))
gam.check(model4)

#model with whole dataset.

mycl <- makeCluster(9) 

clusterExport(mycl, "data_df") 

clusterEvalQ(mycl, {
  library(mgcv)
})

(b <- Sys.time())

#sun_elev as random effect
m1 <- bam(delta_t ~ s(lat,lon, bs = "sos") + s(yday, bs ="cc") + ##spherical spline for lat and lon
            s(sun_elev_f, bs = "re") + 
            s(year, bs = "re"), #re means random effect
          method = "REML", data = data_df, cluster = mycl) 

m2 <- bam(delta_t ~ s(lat,lon, bs = "sos", by = sun_elev_f) +
            s(yday, by = sun_elev_f, bs = "cc") +
            s(year, bs = "re") +
            sun_elev_f , method = "REML", data = data_df, cluster = mycl)

m2_7 <- bam(delta_t_b ~ s(lat,lon, bs = "sos", by = sun_elev_f) +
              s(yday, by = sun_elev_f, bs = "cc") +
              #s(year, bs = "re") +
              sun_elev_f , method = "REML", family = "binomial", data = data_df, cluster = mycl)

m3 <- bam(delta_t ~ s(lat,lon, bs = "sos") +
            s(yday, by = sun_elev_f, bs = "cc") +
            s(year, bs = "re") + sun_elev_f , 
          method = "REML", data = data_df, cluster = mycl)

m4 <- bam(delta_t ~ s(lat,lon, bs = "sos", by = sun_elev_f) +
            s(yday, bs = "cc") +
            s(year, bs = "re") + sun_elev_f , 
          method = "REML", data = data_df, cluster = mycl)

m5 <- bam(delta_t ~ s(lat,lon, bs = "sos", by = as.numeric(sun_elev_f == "night")) +
            s(lat,lon, bs = "sos",by = as.numeric(sun_elev_f == "high")) +
            s(lat,lon, bs = "sos",by = as.numeric(sun_elev_f == "low")) +
            s(yday, bs = "cc") + 
            s(year, bs = "re") + sun_elev_f , 
          method = "REML", data = data_df, cluster = mycl)

m6 <- bam(delta_t ~ s(lon,lat) +
            s(yday, bs = "cc", by = as.numeric(sun_elev_f == "night")) +
            s(yday, bs = "cc", by = as.numeric(sun_elev_f == "high")) +
            s(yday, bs = "cc", by = as.numeric(sun_elev_f == "low")) +
            s(year, bs = "re") +
            sun_elev_f , 
          method = "REML", data = data_df, cluster = mycl)

m7 <- bam(delta_t ~ s(lat, bs = "cr", by = sun_elev_f) + #lat and lon separately. Kami's suggestion
            s(lon, bs = "cr", by = sun_elev_f) +
            s(yday, by = sun_elev_f, bs = "cc") +
            s(year, bs = "re") +
            sun_elev_f , method = "REML", data = data_df, cluster = mycl) #didnt solve the problem of heteroscedasticity

m8 <- bam(delta_t ~ s(lat,lon, bs = "sos", by = sun_elev_f) + #try a poisson distr
            s(yday, by = sun_elev_f, bs = "cc") +
            s(year, bs = "re") +
            sun_elev_f , method = "REML", family = "poisson", data = data_df, cluster = mycl) #error: negative values not allowed for poisson. lol

stopCluster(mycl)

Sys.time() - b

models <- list(m1,m2,m3,m4,m5, m6)
save(models, file = "gam_models.RData")

load("gam_models.RData")
m2 <- models[[2]]

AIC(m1,m2,m3,m4,m5, m6) #m2 has the lowest AIC and highest df

X11(width = 15, height = 10);par(mfrow= c(2,2), oma = c(0,0,3,0))
gam.check(m2) #higher variation at higher values....
X11(width = 15, height = 10);par(mfrow= c(3,3), oma = c(0,0,3,0))
plot(m2)

#plot resids against variables
X11(width = 15, height = 10);par(mfrow= c(2,2), oma = c(0,0,3,0))
plot(data_df$yday, resid(m2), xlab = "yday", ylab = "residuals")
plot(data_df$lat, resid(m2), xlab = "lat", ylab = "residuals") #higher variance for higher values
plot(data_df$lon, resid(m2), xlab = "lon", ylab = "residuals")
plot(data_df$sun_elev_f, resid(m2), xlab = "sun elevation", ylab = "residuals")





#try adding variance structures
mycl <- makeCluster(9) 

clusterExport(mycl, "data_df") 

clusterEvalQ(mycl, {
  library(mgcv)
})

m2_1 <- bam(delta_t ~ s(lat,lon, by = sun_elev_f, bs = "sos") + #try a spherical spline for lat and lon
            s(yday, by = sun_elev_f, bs = "cc") +
            s(year, bs = "re") +
            sun_elev_f , method = "REML", data = data_df, cluster = mycl)

m2_1 <- gamm(delta_t ~ s(lat,lon,bs = "sos") +
              s(yday, bs = "cc") , data = data_df, 
            weights = varPower(form = ~lat))#, cluster = mycl) #when using weights, use gamm instead of bam.

m2_2 <- bam(delta_t ~ s(lon,lat, by = sun_elev_f) +
              s(yday, by = sun_elev_f, bs = "cc") +
              s(year, bs = "re") +
              sun_elev_f  , data = data_df, 
            weights = varFixed(form = ~1 | lat), cluster = mycl)

m2_3 <- gamm(delta_t ~ s(lat,lon, bs = "sos", by = sun_elev_f) + #did not converge
               s(yday, by = sun_elev_f, bs = "cc") +
               s(year, bs = "re") + sun_elev_f, 
             weights = varPower(), method = "REML", data = data_df) #if no arguments are specified, varPower uses the default which is ~fitted(.). source: Piheiro and Bates 2000

#to solve the convergence issue, use scaled variables and make the model simpler
(b <- Sys.time())
m2_4 <- gamm(delta_t ~ s(lat_z,lon_z, bs = "sos", k = 100) +
               s(yday_z, bs = "cc") +
               #s(year, bs = "re") + 
               sun_elev_f, 
             weights = varPower(), method = "REML", data = data_df) 
Sys.time()-b

save(m2_4, file = "model_2_4.RData")

load("model_2_3.RData")

#try the model with higher k. default is 50
#from ?choose.k: So, exact choice of k is not generally critical: it should be chosen to be large enough that you are reasonably 
#sure of having enough degrees of freedom to represent the underlying ‘truth’ reasonably well, but small enough to maintain reasonable computational efficiency. 
#Clearly ‘large’ and ‘small’ are dependent on the particular problem being addressed.

m2_5 <- bam(delta_t ~ s(lat,lon, bs = "sos", by = sun_elev_f, k = 100) + #AIC decreases compared to AIC(m2), when k=70 (3307754). but gam.check indicated that it may still be too low
              s(yday, by = sun_elev_f, bs = "cc") +
              s(year, bs = "re") +
              sun_elev_f , method = "REML", data = data_df, cluster = mycl)
stopCluster(mycl)

AIC(m2_5) #3297569 #AIC is improved compared to m2 and m2_5 with k=70. but heterescedasticity is still tehre
X11(width = 15, height = 10);par(mfrow= c(2,2), oma = c(0,0,3,0))
gam.check(m2_5)


rsd <- residuals(m2_5)
gam(rsd~s(yday,k=40,bs="cc"),gamma=1.4,data=data_df) ## fine
gam(rsd~s(lat,lon,k=100,bs="sos"),gamma=1.4,data=data_df) ## `k' too low
gam(rsd~s(x3,k=40,bs="cs"),gamma=1.4,data=dat) ## fine


f <- formula(delta_t ~ lat + lon + yday + sun_elev_f)
d <- gls(f, data = data_df)
d3 <- gls(delta_t ~ lat * lon + yday + sun_elev_f, data = data_df)
d2 <- gls(f, data = data_df, weights = varFixed(~lat)) #didnt work. too big




#add a random effect for the year to m1
mycl <- makeCluster(9) 
clusterExport(mycl, "data_df") 
clusterEvalQ(mycl, {
  library(mgcv)
})

(b <- Sys.time())
m1b <- gamm(delta_t ~ s(lon,lat, by = sun_elev_f) + #R crashes
            s(yday, by = sun_elev_f, bs = "cc") + 
            sun_elev_f, random = list(year = ~1)
            , data = data_df, cluster = mycl)


Sys.time() - b

AIC(m1b) #3389866
X11(width = 15, height = 10);par(mfrow= c(2,2), oma = c(0,0,3,0))
gam.check(m1b)


#try model 1 with scaled variables
mycl <- makeCluster(9) 
clusterExport(mycl, "data_df") 
clusterEvalQ(mycl, {
  library(mgcv)
})

(b <- Sys.time())
m1c <- bam(delta_t ~ s(lon_z,lat_z, by = sun_elev_f) +
              s(yday_z, by = sun_elev_f, bs = "cc") + 
              sun_elev_f, data = data_df, cluster = mycl)

stopCluster(mycl)
Sys.time() - b

AIC(m1c) #3357867 ...higher than m1 without scaling
X11(width = 15, height = 10);par(mfrow= c(2,2), oma = c(0,0,3,0))
gam.check(m1c)

#plots are the same for scaled vars
plot(data_df$yday_z, resid(m1c), xlab = "yday", ylab = "residuals")
plot(data_df$lat_z, resid(m1c), xlab = "lat", ylab = "residuals")
plot(data_df$lon_z, resid(m1c), xlab = "lon", ylab = "residuals")
plot(data_df$sun_elev_f, resid(m1c), xlab = "sun elevation", ylab = "residuals")



AIC(m3) #3385093
plot(m3)
gam.check(m3)

#add year
mycl <- makeCluster(9) 

clusterExport(mycl, "data_df") 

clusterEvalQ(mycl, {
  library(mgcv)
})

#add year;
(b <- Sys.time())
m4 <- gamm(delta_t ~ s(lon,lat, by = as.numeric(sun_elev == "night")) +
            s(lon,lat, by = as.numeric(sun_elev == "high")) +
            s(lon,lat, by = as.numeric(sun_elev == "low")) + 
            s(yday) + factor(sun_elev) , random = list(year = ~1),
          data = data_df, cluster = mycl) #singularity error

stopCluster(mycl)

Sys.time() - b

#add sun_elev as random effect
mycl <- makeCluster(9) 

clusterExport(mycl, "data_df") 

clusterEvalQ(mycl, {
  library(mgcv)
})

# include sun_elev as random
(b <- Sys.time())
m5 <- gamm(delta_t ~ s(lon,lat) + s(yday, bs = "cc"), 
           random = list(sun_elev_f = ~1),
           data = data_df, cluster = mycl)

stopCluster(mycl)

Sys.time() - b

#model using scaled variables
mycl <- makeCluster(9) 

clusterExport(mycl, "data_df") 

clusterEvalQ(mycl, {
  library(mgcv)
})

(b <- Sys.time())
m1z <- bam(delta_t ~ s(lon_z,lat_z, by = as.numeric(sun_elev == "night")) +
            s(lon_z,lat_z, by = as.numeric(sun_elev == "high")) +
            s(lon_z,lat_z, by = as.numeric(sun_elev == "low")) + 
            s(yday_z, bs = "cc") + sun_elev_f , method = "REML",
           weights = varIdent(form = ~1 | sun_elev_f)
           , data = data_df, cluster = mycl) #higher var explained compared to m2 and m3

stopCluster(mycl)

Sys.time() - b

X11(width = 15, height = 10);par(mfrow= c(2,2), oma = c(0,0,3,0))
gam.check(m1z)
AIC(m1z) #3366801 #larger than m1

####
#following zuur et al p. 412 to deal with heterogeneity of the residuals
#f1 <- formula(delta_t ~ s(lon,lat) + s(yday) + sun_elev_f)

lmc <- lmeControl(niterEM = 5000, msMaxIter =  1000)

f1<-formula(delta_t ~ s(lon,lat, by = as.numeric(sun_elev == "night")) +
                  s(lon,lat, by = as.numeric(sun_elev == "high")) +
                  s(lon,lat, by = as.numeric(sun_elev == "low")) + 
                  s(yday, bs = "cc") + #apply cyclic smoother for yday. the ends meet
              sun_elev_f) 
            
            

m2A <- gamm(f1, random = list(year = ~1),
            method = "REML", control = lmc, data = sample) #singularity
m2B <- gamm(f1, random = list(year = ~1),
            method = "REML", control = lmc, data = sample,
            weights = varIdent(form = ~1 | sun_elev_f))


mycl <- makeCluster(9) 

clusterExport(mycl, list("data_df", "f1", "lmc")) 

clusterEvalQ(mycl, {
  library(mgcv)
})



(b <- Sys.time())
m1A <- gamm(f1,random = list(sun_elev_f = ~1),
            method = "REML", control = lmc, data = data_df) #assumes homogeneity (added as a point of reference)
m1B <- gamm(f1,random = list(sun_elev_f = ~1),
            method = "REML", control = lmc, data = data_df,
            weights = varIdent(form = ~1 | sun_elev_f)) #assumes heterogeneity per sun_elev_f, but homogeneity within a sun_elev_f along lat and lon

m1C<- gamm(f1, random=list(sun_elev_f = ~ 1), 
           method="REML", control = lmc, data = data_df, weights = varPower(form = ~lat)) #assumes homogeneity between sun elevs but heteorgeneity within a sun elev along lat

m1D<- gamm(f1, random=list(sun_elev_f=~1),
           method="REML", control = lmc, data=data_df, weights = varComb(varIdent(form = ~1 | sun_elev_f),varPower(form = ~ lat))) #allows for heterogeneity between stations and within stations along lat

m1E <- gamm(f1, random=list(sun_elev_f=~1),
            method="REML",control = lmc, data=data_df, weights = varComb(varIdent(form = ~1 | sun_elev_f),varPower(form = ~ lat | sun_elev_f))) #heterogeneity within sun elevs along lat is allowed to differ between sun elevs

stopCluster(mycl)

Sys.time() - b

AIC(m1A$lme, m1B$lme, m1C$lme, m1D$lme,m1E$lme)

###
#scale delta_t
b <- Sys.time()
m1b <- bam(scale(delta_t) ~ s(lon_z,lat_z, by = as.numeric(sun_elev == "night")) +
            s(lon_z,lat_z, by = as.numeric(sun_elev == "high")) +
            s(lon_z,lat_z, by = as.numeric(sun_elev == "low")) + 
            s(yday_z) + factor(sun_elev) ,
           data = data_df, cluster = mycl) #AIC improved compared with m1. or are they even comparable!?

stopCluster(mycl)

Sys.time() - b


#scale delta_t and add year
b <- Sys.time()
m1c <- gamm(scale(delta_t) ~ s(lon,lat, by = as.numeric(sun_elev == "night")) + #singularity error
             s(lon,lat, by = as.numeric(sun_elev == "high")) +
             s(lon,lat, by = as.numeric(sun_elev == "low")) + 
             s(yday) + factor(sun_elev) , data = data_df, random = list(year = ~1))#cluster = mycl) #higher var explained compared to m2 and m3

stopCluster(mycl)

Sys.time() - b



b <- Sys.time()
m2 <- bam(delta_t ~ s(lon,lat) + s(yday,by = factor(sun_elev)) + factor(sun_elev), 
          data = data_df, cluster = mycl)

m3 <- bam(delta_t ~ s(lon,lat) + s(yday) + factor(sun_elev), 
          data = data_df, cluster = mycl)

stopCluster(mycl)

Sys.time() - b

X11(width = 15, height = 10);par(mfrow= c(2,2), oma = c(0,0,3,0))
gam.check(m2)

X11(width = 15, height = 10);par(mfrow= c(2,2), oma = c(0,0,3,0))
gam.check(m3)

model<-bam(delta_t ~ s(lon,lat) + s(yday) + factor(local_hour) , data = sample)
AIC(model) 
plot(model)



model3<-bam(delta_t ~ s(lon,lat, by = factor(month)) + factor(sun_elev) , data = sample)
AIC(model3) 
plot(model3)

gam.check(model3)

plot_smooths(
  model = model3,
  series = lon,
  comparison = factor(month)
) +
  theme(legend.position = "top")

model4 <- gamm(delta_t ~ s(lon,lat) + s(yday,by = factor(hour)) + factor(hour),
               random = list(year = ~1), data = sample)


#compare models with Wald test
anova(model, model2)



#######spatial model with inla ###
#convert lat and long to northing and easting
#working through this book: https://becarioprecario.bitbucket.io/spde-gitbook/ch-INLA.html

#global mesh: (https://groups.google.com/forum/#!topic/r-inla-discussion-group/dcrajXPn-eo)
# To build a spherical mesh, you just need to convert the
# longitude-latitude cutoff&max.edge parameters to radians (multiply by
#                                                           pi/180), and supply
# crs = inla.CRS("sphere")
# to inla.mesh.2d(), which will automatically transform inputs with
# proper crs info (so you'll need to specify a crs=inla.CRS("longlat")
# when you construct the nonconvex hull).


#use default priors (gamma)
m.rw1 <- inla(delta_t ~ f(lat, model = "rw1",scale.model = TRUE) +
                f(lon, model = "rw1",scale.model = TRUE) + #scale.model makes the model to be scaled to have an average variance of 1
                f(yday, model = "rw1",scale.model = TRUE),
              data = sample,
              control.compute = list(dic = T, waic = T)) #the marginal likelihood cannot be used for the improper models, as the renormalization constant is not included. See inla.doc("rw1") for details.

#use pc priors. assuming that the probability of the standard deviation being higher than 1 is quite small. so, set u=1 and alpha to 0.01
pcprior <- list(prec = list(prior = "pc.prec", param = c(1,0.01)))

m.rw2 <- inla(delta_t ~ f(lat, model = "rw1",scale.model = TRUE) +
  f(lon, model = "rw1",scale.model = TRUE, hyper = pcprior) + 
  f(yday, model = "ar1", hyper = pcprior), #autoregressive model... higher log likelihoood than rw1 (rw1 can overfit the data)
  data = sample,
  control.compute = list(dic = T, waic = T))

Efxplot(m.rw2) + theme_bw() #not useful. only shows fixed effects

#inspect the effect of the PC prior
post.sigma.s1 <- inla.tmarginal(function (x) sqrt(1 / exp(x)),
                                m.rw2$internal.marginals.hyperpar[[2]])

post.sigma.s2 <- inla.tmarginal(function (x) sqrt(1 / exp(x)),
                                m.rw2$internal.marginals.hyperpar[[3]])

#compare cpos
 sum(log(climate.ar1$cpo$cpo))
## [1] -43.05
# rw1
 sum(log(climate.rw1$cpo$cpo))

 
 
 
### STEP 4: effect plot #####

#extract range of yday of sea-crossing to use on the plot later (data from ssf_delta_t_inla2.R)
load("ssf_input_ann_cmpl_1hr.RData") #ann_cmpl
ann_cmpl$yday <- yday(ann_cmpl$timestamp)
 
 
#extract data for creating the plots
plotdata <- plot(m2,pages = 1)

plot(delta_t  )

library(ggplot2)
theme_set(theme_bw())
library(dplyr)
library(mgcv)
library(tidymv)



b <- getViz(m2_3$gam)

X11()
plot_smooths(model = m2, series = yday, comparison = sun_elev_f) +
  theme(legend.position = "top") # this looks good!



pdf("delta_T_gam/yday_smoother_fig_model1_3.pdf", height=4, width=4)
#x11()
#dev.new(width=4, height=4)

par(mfrow=c(1,1), bty="n", #no box around the plot
    cex.axis= 0.75, #x and y labels have 0.75% of the default size
    font.axis= 3, #axis labels are in italics
    cex.lab=1
)


plot(m2,select=4, xlab="Julian date",ylab="Delta T", ylim= c(-2,2), shade=T, shade.col= "grey75",
     scheme=0,se=12,bty="l",labels = FALSE, tck=0) #plot only the second smooth term

#rect(xleft=268,ybottom=-2.3,xright=296,ytop=2.3, col="#96CDCD30",border=NA) #for juv...#last two digits are for transparency
#rect(xleft=242,ybottom=-2.3,xright=277,ytop=2.3, col="#FA807230",border=NA) #for adults

rect(xleft=268,ybottom=-2.3,xright=296,ytop=1.14, col="#96CDCD30",border=NA) #for juv...#last two digits are for transparency
rect(xleft=213,ybottom=0.35,xright=268,ytop=1.14, col="#96CDCD30",border=NA)

rect(xleft=242,ybottom=-2.3,xright=277,ytop=0.6, col="#FA807230",border=NA) #for adults
rect(xleft=213,ybottom=-0.5,xright=242,ytop=0.6, col="#FA807230",border=NA)

axis(side= 1, at= c(213,240,260,280,305), line=-0, labels= c(213,240,260,280,305), 
     tick=T , col.ticks = 1, col=NA, tck=-.015)
axis(side= 2, at= c(-1.5,0,1.5), line=0, labels= c(-1.5,0,1.5),
     tick=T , col.ticks = 1,col=NA, tck=-.015, 
     las=2) # text perpendicular to axis label 
legend(x=212, y=2, legend=c("juveniles", "adults"), fill=c("#96CDCD40","#FA807240"), #coords indicate top-left
       cex=0.8, bg="white",bty="n")

dev.off()

#----------------------------------------------------------------------------------------
# Predict with the model!
#----------------------------------------------------------------------------------------
#make predictions with the model for the periods of adult and juvenile migration
#juvenile Sep. 25-Oct. 23 (i.e. ydays 268-296); adults Aug 30-Oct. 4 (i.e. ydays 242-277)
#create a set for each..... conclusion: the two prediction maps look very similar. so, just do one map for the entire autumn
seasons<-list(
  juv=data_aut_sea[data_aut_sea$yday %in% c(268:296),],
  adlt=data_aut_sea[data_aut_sea$yday %in% c(242:277),]
)

#predict
preds_ls<-lapply(seasons, function(x){
  pred<-data.frame(pred=as.numeric(predict(model,x)),lon=x$lon,lat=x$lat)
  
  coordinates(pred)<-~lon+lat
  gridded(pred)<-T
  r<-raster(pred)
  
  return(r)
})


save(preds_ls,file="delta_T_gam/predictions_with_model1.RData")
load("delta_T_gam/predictions_with_model1.RData")


###step 6: modleing-full model #####

###step 7: visualization #####
#predict with the full model for autumn migration season

# interpolate to 1 km for visualization purposes



#####
windows()
map("world",fill=TRUE,col="gray",border="gray")
points(sample$lon,sample$lat,cex=0.5,pch=16,col="red")
points(sample2$lon,sample2$lat,cex=0.3,pch=16,col="green")
plot(st_geometry(land_pol))
lines(land)
plot(lakes,col="red",add=T)