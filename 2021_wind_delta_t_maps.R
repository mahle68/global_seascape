#script for creating maps color-coded with wind and delta-t of the raw data
#also look into Gil's suggestion for showing delta t is correlated with uplift

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

#open sea-crossing points, prepared in 2021_all_data_preo_analyze.R
load("R_files/2021/all_2009_2020_overwater_points.RData") #all_oversea

region <- st_read("/home/enourani/ownCloud/Work/GIS_files/continent_shapefile/continent.shp") %>% 
  st_crop(xmin = -130, xmax = 158, ymin = -74, ymax = 71) %>%
  st_union()

#prepare for movebank annotation

mv <- all_oversea %>% 
  as("Spatial") %>% 
  as.data.frame() %>% 
       mutate(timestamp = paste(as.character(date_time),"000",sep = "."))
     
     
colnames(mv)[c(12,13)] <- c("location-long","location-lat")
       
write.csv(mv, "R_files/2021/raw_points_for_maps.csv")
       

#annotated data
ann <- read.csv("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/2021/annotations/raw_points_for_maps/raw_points_for_maps.csv-5184501718572889126.csv") %>% 
  mutate(timestamp,timestamp = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>%
  rename(sst = ECMWF.ERA5.SL.Sea.Surface.Temperature,
         t2m = ECMWF.ERA5.SL.Temperature..2.m.above.Ground.,
         u925 = ECMWF.ERA5.PL.U.Wind,
         v925 = ECMWF.ERA5.PL.V.wind,
         blh = ECMWF.Interim.Full.Daily.SFC.FC.Boundary.Layer.Height,
         s_flux = ECMWF.Interim.Full.Daily.SFC.FC.Instantaneous.Surface.Heat.Flux,
         m_flux = ECMWF.Interim.Full.Daily.SFC.FC.Instantaneous.Moisture.Flux) %>%
  mutate(delta_t = sst - t2m) %>%
  arrange(track, timestamp)

#remove tracks with one point
more_than_one_point <- ann %>% 
  group_by(track) %>% 
  summarize(n = n()) %>% 
  filter(n > 1)

ann <- ann %>% 
  filter(track %in% more_than_one_point$track)

#remove duplicated rows
rows_to_delete <- unlist(sapply(getDuplicatedTimestamps(x = as.factor(ann$track),timestamps = ann$timestamp),"[",-1)) #get all but the first row of each set of duplicate rows
ann <- ann[-rows_to_delete,]


#create a move object, to calculate heading using the angle function
mv <- move(x = ann$location.long, y = ann$location.lat, time = ann$timestamp, data = ann, animal = ann$track, proj = wgs)

#calculate heading
mv$heading <- unlist(lapply(angle(mv), c, NA))
mv2 <- mv[complete.cases(mv$heading),]
mv2$wind_support <- wind_support(u = mv2$u925,v = mv2$v925, heading = mv2$heading)
mv2$cross_wind <- cross_wind(u = mv2$u925,v = mv2$v925, heading = mv2$heading)


#plot depending on wind support values
Pal <- colorRampPalette(c('thistle1',"slateblue1",'slateblue4'))

#This adds a column of color values
# based on the y values
mv2$Col <- rbPal(8)[as.numeric(cut(mv2$wind_support, breaks = 8))]

plot(region, col="#e5e5e5",border="#e5e5e5")
points(mv2,pch = 20,col = mv2$Col)

legend("topleft",title="Decile",legend=c(1:3),col =rbPal(3),pch=20)


#################################

Pal_p <- colorRampPalette(c("plum3",'slateblue4')) #colors for positive values
Pal_n <- colorRampPalette(c('indianred1','tan1')) #colors for negative values
Cols <- c(Pal_n(3),Pal_p(4))

Cols_t <- paste0(Cols, "E6") #add transparency. 50% is "80". 70% is "B3". 80% is "CC". 90% is "E6"
  

#add a categorical variable for wind levels
breaks <- c(-20,-10,-5,0,5,10,15,35)
tags <- c("< -10","-10 to -5","-5 to 0","0 to 5","5 to 10","10 to 15", "> 15")

mv2$binned_wind <- cut(mv2$wind_support,breaks = breaks, include.lowest = T, right = F, labels = tags)
mv2$col <- as.factor(mv2$binned_wind)
levels(mv2$col) <- Cols_t

negatives <- mv2[mv2$col %in% Cols_t[1:3],]
positives <- mv2[mv2$col %in% Cols_t[4:7],]

plot(region, col="#e5e5e5",border="#e5e5e5")

points(negatives,pch = 20, col = as.character(negatives$col), cex = 0.3) #plot negative points first
points(positives,pch = 20, col = as.character(positives$col), cex = 0.3) #plot positive points

points(mv2, pch = 1, col = as.character(mv2$col), cex = 0.2)


legend("bottomleft", legend = levels(mv2$binned_wind), col = Cols_t, pch = 20, bty = "n")


mapview(mv2,zcol = "wind_support")



df <- as.data.frame(mv) %>% 
  drop_na(heading) %>% 
  mutate(wind_support = ,
         cross_wind= cross_wind(u=u925,v=v925,heading=heading),
         wind_speed = sqrt(u925^2 + v925^2),
         abs_cross_wind = abs(cross_wind(u = u925, v = v925, heading = heading)))


# 
