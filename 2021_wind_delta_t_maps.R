#script for creating maps color-coded with wind and delta-t of the raw data
# maybe also plot all the points. load("R_files/2021/all_spp_unfiltered_updated_lc_0_removed_new_track_id.RData") #called dataset ... as a background.
# or as lines. with the same characteristics as the gam map

#Mar 22. 2021, Radolfzell, DE
#Elham Nourani

library(tidyverse)
library(sp)
library(sf)
library(move)
library(scales)

setwd("/home/mahle68/ownCloud/Work/Projects/delta_t")
source("/home/mahle68/ownCloud/Work/Projects/delta_t/R_files/wind_support_Kami.R")

wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

#open sea-crossing points, prepared in 2021_all_data_preo_analyze.R
load("R_files/2021/all_2009_2020_overwater_probl_pts_removed.RData") #all_oversea

region <- st_read("/home/mahle68/ownCloud/Work/GIS_files/continent_shapefile/continent.shp") %>% 
  st_crop(xmin = -99, xmax = 144, ymin = -30, ymax = 71) %>%
  st_union()

#prepare for movebank annotation
mv <- all_oversea %>% 
  as("Spatial") %>% 
  as.data.frame() %>% 
       mutate(timestamp = paste(as.character(date_time),"000",sep = "."))
     
     
colnames(mv)[c(12,13)] <- c("location-long","location-lat")
       
write.csv(mv, "R_files/2021/raw_points_for_maps_updated.csv")
       

#annotated data
ann <- read.csv("/home/mahle68/ownCloud/Work/Projects/delta_t/R_files/2021/annotations/raw_points_for_maps/era5/raw_points_for_maps.csv-8853543643507767873.csv") %>% 
  mutate(timestamp,timestamp = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>%
  rename(sst = ECMWF.ERA5.SL.Sea.Surface.Temperature,
         t2m = ECMWF.ERA5.SL.Temperature..2.m.above.Ground.,
         u925 = ECMWF.ERA5.PL.U.Wind,
         v925 = ECMWF.ERA5.PL.V.wind) %>% 
  mutate(delta_t = sst - t2m) %>%
  arrange(track, timestamp)

#remove tracks with one point
more_than_two_point <- ann %>% 
  group_by(track) %>% 
  summarize(n = n()) %>% 
  filter(n > 2)

ann <- ann %>% 
  filter(track %in% more_than_two_point$track)

save(ann, file = "R_files/2021/raw_sea_points_annotated.RData")

#remove duplicated rows
rows_to_delete <- unlist(sapply(getDuplicatedTimestamps(x = as.factor(ann$track),timestamps = ann$timestamp),"[",-1)) #get all but the first row of each set of duplicate rows
ann <- ann[-rows_to_delete,]


#create a move object, to calculate heading using the angle function
mv <- move(x = ann$location.long, y = ann$location.lat, time = ann$timestamp, data = ann, animal = ann$track, proj = wgs)

#calculate heading
mv$heading <- unlist(lapply(angle(mv), c, NA))

save(mv, file = "R_files/2021/raw_sea_points_for_maps_mv.RData")


#plot

load("R_files/2021/raw_sea_points_for_maps_mv.RData") #mv

#add a categorical variable for wind levels
breaks_w <- c(-20,-10,-5,0,5,10,15,35)
tags_w <- c("< -10","-10 to -5","-5 to 0","0 to 5","5 to 10","10 to 15", "> 15")
#add a categorical variable for wind levels
breaks_dt <- c(-5,-2,0,2,5,10)
tags_dt <- c("< -5","-5 to -2","0 to 2","2 to 5", "> 5")

#create color palettes and select colors for positive and negative values. same colors for delta t and wind?

#Pal_p <- colorRampPalette(c("plum3",'slateblue4')) #colors for positive values
Pal_n <- colorRampPalette(c('indianred1','darkgoldenrod2')) #colors for negative values
Pal_p <- colorRampPalette(c('cornflowerblue',"mediumblue")) #colors for positive values
Cols_w <- paste0(c(Pal_n(3),Pal_p(4)), "E6") #add transparency. 50% is "80". 70% is "B3". 80% is "CC". 90% is "E6"
Cols_dt <- paste0(c(Pal_n(2),Pal_p(3)), "E6")



df <- mv %>% 
  as.data.frame() %>% 
  drop_na(c("heading","delta_t")) %>% 
  mutate(wind_support= wind_support(u = u925,v = v925,heading = heading),
         cross_wind= cross_wind(u = u925,v = v925,heading = heading)) %>% 
  mutate(binned_w = cut(wind_support,breaks = breaks_w, include.lowest = T, right = F, labels = tags_w),
         binned_dt = cut(delta_t,breaks = breaks_dt, include.lowest = T, right = F, labels = tags_dt)) %>% 
  mutate(cols_w = as.factor(binned_w),
         cols_dt = as.factor(binned_dt)) 

levels(df$cols_w) <- Cols_w
levels(df$cols_dt) <- Cols_dt

df_sp <- SpatialPointsDataFrame(coords = df[,c("location.long", "location.lat")], proj4string = wgs, data = df)

#plot
X11(width = 12, height = 11.5) #make the window proportional to region

pdf("/home/mahle68/ownCloud/Work/Projects/delta_t/paper_prep/figures/2021/raw_wind_dt.pdf", width = 12, height = 11.5)

par(mfrow=c(2,1),
    #fig = c(0,1,0,1), #do this if you want to add the small plots as subplots
    bty="n", #no box around the plot
    cex.axis= 0.6, #x and y labels have 0.75% of the default size
    font = 3, #3: axis labels are in italics
    font.axis = 3,
    cex.lab = 0.6,
    #cex = 0.5,
    oma = c(0,0,1.5,0),
    mar = c(0, 0, 0, 0),
    lend = 1  #rectangular line endings (trick for adding the rectangle to the legend)
)


plot(region, col="#e5e5e5",border="#e5e5e5")
points(df_sp, pch = 1, col = as.character(df_sp$cols_w), cex = 0.2)

#add latitudes
clip(-99, 144, -30, 71)
abline(h = 0, col = "grey70",lty = 2)
abline(h = 30, col = "grey70",lty = 2)
abline(h = 60, col = "grey70",lty = 2)
#text(x = -125, y = c(2,32,62), labels = c("0° ", "30° N", "60° N"), cex = 0.6, col = "grey65")
text(x = -95, y = c(32,62), labels = c("30° N", "60° N"), cex = 0.6, col = "grey65")


#add a frame for the sub-plots and legend
rect(xleft = -100,
     xright = -71,
     ybottom =  -30,
     ytop = 2,
     col="white",
     border = NA)

text(x = -85,y = 0, "Wind support (m/s)", cex = 0.8)
legend(x = -100, y = 0, legend = levels(df_sp$binned_w), col = Cols_w, pch = 20, 
       bty = "n", cex = 0.8)
mtext("Bio-logging data annotated with wind support", 3, outer = F, cex = 1.3, line = -1)

plot(region, col="#e5e5e5",border="#e5e5e5")
points(df_sp, pch = 1, col = as.character(df_sp$cols_dt), cex = 0.2)

#add latitudes
clip(-99, 144, -30, 71)
abline(h = 0, col = "grey70",lty = 2)
abline(h = 30, col = "grey70",lty = 2)
abline(h = 60, col = "grey70",lty = 2)
#text(x = -125, y = c(2,32,62), labels = c("0° ", "30° N", "60° N"), cex = 0.6, col = "grey65")
text(x = -95, y = c(32,62), labels = c("30° N", "60° N"), cex = 0.6, col = "grey65")


#add a frame for the sub-plots and legend
rect(xleft = -100,
     xright = -80,
     ybottom =  -30,
     ytop = 2,
     col="white",
     border = NA)

text(x = -93,y = 0,  expression(italic(paste(Delta,"T", "(°C)"))), cex = 0.8)
legend(x = -100, y = -1, legend = levels(df_sp$binned_dt), col =Cols_dt, pch = 20, 
       bty = "n", cex = 0.8)
mtext(bquote(italic('Bio-logging data annotated with' ~ Delta *"T")), 3, outer = F, cex = 1.3, line = -1)

dev.off()

