#scripts for building a global GAM for delta-t
#data processed in R_files/global_seascape.R
#May 6, 2020. Radolfzell, DE. Elham Nourani, PhD
#update May 7, 2020: calculate solar position to add a meaningful proxy for time of the day to the model (so, basically, the local time calculation was not necessary)

#tutorials
#https://www.r-bloggers.com/advanced-sab-r-metrics-parallelization-with-the-mgcv-package/
  
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

### STEP 2: data exploration #####
#check to see if there are any monthly, hourly, latitudinal differences

#running this on the entire dataset gives ram errors! so, just visualize the sample! lol

X11(width = 15, height = 10);par(mfrow= c(2,2), oma = c(0,0,3,0))
#for(i in c("sun_elev", "month","local_hour")){
    
  boxplot(sample$delta_t ~ sample$month, boxfill = "skyblue1", ylab = "delta_t")
  boxplot(sample$delta_t ~ sample$sun_elev, boxfill = "wheat1", ylab = "delta_t")
  boxplot(sample$delta_t ~ sample$local_hour, boxfill = "mediumblue", ylab = "delta_t")
  plot(sample$delta_t ~ sample$lat , pch = 16, col = factor(sample$month), ylab = "delta_t")
  

mtext("40-yr averages at each point", side = 3, outer = T, cex = 1.3)


### STEP 3: make the model #####

#prep data
data_df <- data %>% 
  as("Spatial") %>% 
  as.data.frame() %>% 
  rename(lon = coords.x1,
         lat = coords.x2)

save(data_df, file = "t_data_df_0_60.RData")

#try INLA also :p ... is probably more appropriate

#include in the model: delta T ~ s(yday)+ s(lon,lat) + sun_elevation + random(year)

#for modleing, use bam() instead of gam(), it breaks the data into chunks and has a lower memory footprint. also has the option to parallelize

#work on a sample first
sample <- data %>% 
  st_crop(xmin = 120, xmax = 130, ymin = 25, ymax = 35) %>% 
  as("Spatial") %>% 
  as.data.frame() %>% 
  rename(lon = coords.x1,
         lat = coords.x2)

model<-bam(delta_t ~ s(lon,lat) + s(yday) + factor(sun_elev) , data = sample)
AIC(model)
gam.check(model)

#separate smooths for time of day
model2 <- bam(delta_t ~ s(lon,lat) + s(yday,by = factor(sun_elev)) + factor(sun_elev) , data = sample)
AIC(model2) 
plot(model2)

#separate smooths for time of day
model3 <- bam(delta_t ~ s(lon,lat, by = as.numeric(sun_elev == "night")) +
                s(lon,lat, by = as.numeric(sun_elev == "high")) +
                s(lon,lat, by = as.numeric(sun_elev == "low")) + 
                s(yday) + factor(sun_elev) , data = sample)
AIC(model3)
gam.check(model3)


b <- Sys.time()
mycl <- makeCluster(9) 

clusterExport(mycl, "data") #define the variable that will be used within the function

clusterEvalQ(mycl, {
  library(mgcv)
  library(mapview)
})

m1 <- bam(delta_t ~ s(lon,lat, by = as.numeric(sun_elev == "night")) +
                s(lon,lat, by = as.numeric(sun_elev == "high")) +
                s(lon,lat, by = as.numeric(sun_elev == "low")) + 
                s(yday) + factor(sun_elev) , data = data, cluster = mycl)

stopCluster(mycl)


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

#save the best model

##### from previous codes

#-----------------------------------------------------------
#effect plots 
#-----------------------------------------------------------
#previous attempts are in the previous version ;)

pdf("delta_T_gam/yday_smoother_fig_model1_3.pdf", height=4, width=4)
#x11()
#dev.new(width=4, height=4)

par(mfrow=c(1,1), bty="n", #no box around the plot
    cex.axis= 0.75, #x and y labels have 0.75% of the default size
    font.axis= 3, #axis labels are in italics
    cex.lab=1
)


plot(model2,select=2, xlab="Julian date",ylab="Delta T", ylim= c(-2,2), shade=T, shade.col= "grey75",
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