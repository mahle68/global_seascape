#code for running gams for each region separately
#Jun 19. 2020
library(mgcv)
library(tidyverse)
library(sp)
library(sf)
library(raster)
library(mapview)
library(parallel)
library(lubridate)
library(fields) #for Tps
library(itsadug) #for gam plots
library(jpeg)
library(TeachingDemos) #for subplot
library(readxl)
library(png)

setwd("/home/enourani/ownCloud/Work/Projects/delta_t/R_files")
wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")

### STEP 1: open and filter ####
load("t_data_0_60.RData")

East_asia <- data %>% 
  st_crop(xmin = 99, xmax = 130, ymin = 2, ymax = 35) %>% 
  #st_difference(st_union(outliers)) %>%  #remove outliers
  as("Spatial") %>% 
  as.data.frame() %>% 
  rename(lon = coords.x1,
         lat = coords.x2) %>% 
  mutate(sun_elev_f = factor(sun_elev)) %>% 
  as.data.frame()
  
#for indian ocean, remove the caspian sea
pg <- st_read("/home/enourani/ownCloud/Work/GIS_files/World_Seas_IHO_v3/World_Seas_IHO_v3.shp") %>% 
  filter(NAME == "Persian Gulf")

Indian_ocean <- data %>% 
  st_crop(xmin = 43, xmax = 78, ymin = 9, ymax = 26) %>% 
  st_difference(pg) %>%  #remove points over the persian gulf
  as("Spatial") %>% 
  as.data.frame() %>% 
  rename(lon = coords.x1,
         lat = coords.x2) %>% 
  mutate(sun_elev_f = factor(sun_elev)) %>% 
  as.data.frame()

#for the Americas, get rid of the pacific ocean

np_ocean <- st_read("/home/enourani/ownCloud/Work/GIS_files/World_Seas_IHO_v3/World_Seas_IHO_v3.shp") %>% 
  filter(NAME == "North Pacific Ocean") %>% 
  st_crop(xmin = -118, xmax = -60, ymin = 0, ymax = 29) 

Americas <- data %>% 
  st_crop(xmin = -98, xmax = -54, ymin = 5.7, ymax = 45) %>% 
  st_difference(np_ocean) %>% 
  st_difference(outlier) %>%  #some points remain on the Pacific ocean. remove them. outlier: outlier <- Americas %>% filter(row_number() == which(Americas$yday == 74 & Americas$year == 2017 & Americas$local_hour == 13))
  #dplyr::select(names(data)) %>%
  #st_drop_geometry() %>% 
  #rename(lon = Longitude,
  #       lat = Latitude) %>% 
  as("Spatial") %>% 
  as.data.frame() %>% 
  rename(lon = coords.x1,
         lat = coords.x2) %>% 
  mutate(sun_elev_f = factor(sun_elev)) %>% 
  as.data.frame()

#for Europe, include mediterranean and the baltic
bs <- st_read("/home/enourani/ownCloud/Work/GIS_files/baltic_sea/iho.shp") %>% 
  st_union()
ms <- st_read("/home/enourani/ownCloud/Work/GIS_files/med_sea/iho.shp") %>% 
  st_union()

eur_sea <- bs %>% 
  st_union(ms)

save(ms, file = "med_sea.RData")
save(bs, file = "blt_sea.RData")
save(eur_sea, file = "eur_sea.RData")
load("eur_sea.RData")

#mask europe layer to keep only waterbodies of itnerest
Europe <- data %>% 
  st_crop(xmin = -11, xmax = 37, ymin = 30, ymax = 62) %>% 
  st_intersection(eur_sea) %>% 
  as("Spatial") %>% 
  as.data.frame() %>% 
  rename(lon = coords.x1,
         lat = coords.x2) %>% 
  mutate(sun_elev_f = factor(sun_elev)) %>% 
  as.data.frame()

save(Europe, file = "Europe_seas.RData")
load("Europe_seas.RData") 

#put all regional data in a list
data_ls <- list(East_asia = East_asia, Americas = Americas, Indian_ocean = Indian_ocean, Europe = Europe)

save(data_ls, file = "data_ls_regional_gam.RData")
save(data_ls, file = "data_ls_regional_gam2.RData") #smaller extend for the Americas

### STEP 2: model####
load("data_ls_regional_gam.RData")

mycl <- makeCluster(9) 

clusterExport(mycl, "data_ls") 

clusterEvalQ(mycl, {
  library(mgcv)
})

models_ls <- lapply(data_ls, function(x){
  
  #m1 <- bam(delta_t ~ s(lat,lon, by = sun_elev_f, k = 100) +
  #      s(yday, by = sun_elev_f, bs = "cc") +
  #      #s(year, bs = "re") +
  #      sun_elev_f , method = "REML", data = x, cluster = mycl)
  
  # bam(delta_t ~ s(lat,lon, by = sun_elev_f, k = 100) +
  #       s(yday, by = sun_elev_f, bs = "cc") +
  #       s(year, bs = "re") +
  #       sun_elev_f , method = "REML", data = x, cluster = mycl)
  
  gamm(delta_t ~ s(lat,lon, by = sun_elev_f, k = 100) +
        s(yday, by = sun_elev_f, bs = "cc") +
        s(year, bs = "re") +
        sun_elev_f , method = "REML", data = x, 
      weights = varPower(form = ~lat))
  
  #models <- list(m1,m2)
  #models
})

stopCluster(mycl)

save(models_ls, file = "models_ls_reg_GAMs.RData")
save(models_ls, file = "models_ls_reg_GAMs2.RData")

load("models_ls_reg_GAMs.RData")

#model checking

lapply(models_ls, function(x){
  X11()
  par(mfrow= c(2,2), oma = c(0,0,3,0))
  gam.check(x$gam)
})

#create output tables in Latex
latex_ls <- lapply(names(models_ls), function(x){
  m <- models_ls[[x]]
  gamtabs(m, caption = x)
})

### STEP 3: prediction and maps ####
#extract migration timing OVER THE SEA for each species
load("all_spp_spatial_filtered_updated.RData") #called dataset_sea. mainly use this. created in all_data_prep_analyze.R

dataset_sea %>% 
  filter(season == "autumn") %>%
  as("Spatial") %>% 
  as.data.frame() %>% 
  mutate(continent = ifelse(coords.x1 < -24, "A","not_relevant")) %>% 
  group_by(continent,species) %>% 
  summarise(min = min(yday(date_time)),
            max = max(yday(date_time)))
  
load("segs_EF_dt.RData") #segs_ann_EF; Eleonora's falcon. also prepped in all_data_prep_analyze.R  
range(yday(segs_ann_EF$timestamp))

timing <- list(OHB = c(260:294), #sea-crossing; ref: almost Yamaguchi
               GFB = c(277:299),
               AF = c(319:323), #mid-november. ref: Bernd's poster
               EF = c(288:301),
               PF_EU = c(261:289),
               O_EU = c(222:277),
               PF_A = c(279:305),
               O_A = c(244:299))

timing_areas <- list(East_asia = c(min(c(timing$GFB,timing$OHB)):max(c(timing$GFB,timing$OHB))),
                     Americas = c(min(c(timing$O_A,timing$PF_A)):max(c(timing$O_A,timing$PF_A))),
                     Indian_ocean = timing$AF,
                     Europe = c(min(c(timing$O_EU,timing$PF_EU,timing$EF)):max(c(timing$O_EU,timing$PF_EU,timing$EF)))) #consider separating this into regions


#make predictions
mycl <- makeCluster(4) 

clusterExport(mycl, list("timing_areas","models_ls", "data_ls", "wgs"))

clusterEvalQ(mycl, {
  library(mgcv)
  library(raster)
  library(dplyr)
  library(sf)
  library(fields) #for Tps
})

preds <- parLapply(cl = mycl, c(names(models_ls)),function(x){
  d <- data_ls[[x]] %>% 
    filter(yday %in% timing_areas[[x]])
  m <- models_ls[[x]]
  
  pred <- data.frame(pred = as.numeric(predict(m,d)), lon = d$lon, lat = d$lat)
  
  coordinates(pred) <- ~lon+lat
  gridded(pred) <- T
  r <- raster(pred)
  proj4string(r) <- wgs
  
  #interpolate. for visualization purposes
  surf.1 <- Tps(as.matrix(as.data.frame(r,xy = T)[,c(1,2)],col = 2),as.data.frame(r,xy = T)[,3])
  
  grd <- expand.grid(x = seq(from = extent(r)[1],to = extent(r)[2],by = 0.1),
                   y = seq(from = extent(r)[3],to = extent(r)[4],by = 0.1))

  grd$coords <- matrix(c(grd$x,grd$y),ncol=2)
  
  surf.1.pred <- predict.Krig(surf.1,grd$coords)
  interpdf <- data.frame(grd$coords,surf.1.pred)
  
  colnames(interpdf)<-c("lon","lat","delta_t")
  
  coordinates(interpdf) <- ~lon+lat
  gridded(interpdf) <- TRUE
  interpr <- raster(interpdf)
  proj4string(interpr) <- wgs
  
  return(interpr)
  
})

stopCluster(mycl)

names(preds) <- names(models_ls)
save(preds, file = "regional_gam_preds.RData")

#mask the rasters with the relevant land layer. or simply plot land over the top
#coastlines layer
region <- st_read("/home/enourani/ownCloud/Work/GIS_files/continent_shapefile/continent.shp") %>% 
  st_crop(xmin = -130, xmax = 157, ymin = -50, ymax = 73) %>%
  st_union()

#mask all with the landmass layer
preds_filt <- lapply(preds, function(x){
  x_f <- raster::mask(x,as(region,"Spatial"), inverse = T) 
    #x %>% 
    #st_difference(region)
  x_f
})

#now for each region, further crop as necessary
waters <- st_read("/home/enourani/ownCloud/Work/GIS_files/World_Seas_IHO_v3/World_Seas_IHO_v3.shp")
  
#Americas
np_ocean <- waters %>% 
  filter(NAME == "North Pacific Ocean") %>% 
  st_crop(xmin = -118, xmax = -60, ymin = 0, ymax = 29) 

preds_filt$Americas <- mask(preds_filt$Americas, np_ocean, inverse = T)

#Europe
load("eur_sea.RData") #eur_sea
preds_filt$Europe <- mask(preds_filt$Europe, as(eur_sea,"Spatial"), inverse = F)

#Indian ocean
pg_rs <- waters %>% 
  filter(NAME %in% c("Persian Gulf","Red Sea"))

preds_filt$Indian_ocean <- mask(preds_filt$Indian_ocean, pg_rs, inverse = T)

#East Asia
ea <- waters %>% 
  filter(NAME %in% c("Andaman or Burma Sea", "Malacca Strait")) %>% 
  st_buffer(0.1) #some points on the edges remain after filtering, so make a buffer to make sure they all go away

preds_filt$East_asia <- mask(preds_filt$East_asia, ea, inverse = T)

save(preds_filt, file = "predictions_regional_gam_map.RData")

load("predictions_regional_gam_map.RData")
load("models_ls_reg_GAMs.RData")

# create a plotting function for effect plots
names <- c("South-East Asia", "the Americas", "Indian Ocean", "Europe")

reg_gam_plot <- function(x){
  m <- models_ls[[x]]
  t <- timing_areas[[x]]
  
  plot(0, type = "n", labels = FALSE, tck = 0, xlim =  c(1,366), ylim = c(-2.5,5), xlab = "", ylab = "")#, main = names[[x]]) #expression(paste(Delta,"T")), main = x
  rect(xleft = min(t),ybottom = -2.7,xright = max(t),ytop = 5, col="#99CC0060",border=NA) #water-crossing window
  plot_smooth(m, view="yday", plot_all="sun_elev_f", rm.ranef=F, lwd = 1.5, #ylim=c(-2,5),
              col = "grey60", hide.label = TRUE, 
              legend_plot_all =  F, 
              h0 = NULL, add = T, lty = c(1,5,3))
  
  axis(side = 1, at = c(100, 200, 300), line = 0, labels = c(100,200,300), 
       tick = T , col.ticks = 1, col = NA, lty = NULL, tck = -.015)
  axis(side = 2, at = c(-2,0,2,4), line = 0.15, labels = c(-2,0,2,4),
       tick = T , col.ticks = 1,col = NA, lty = NULL, tck = -.015, 
       las = 2) # text perpendicular to axis label 
  mtext(expression(italic(paste(Delta,"T"))), 2, line = 1.2 ,las = 2.5, cex = 0.5)
  mtext("Day of year", 1, line = 1.2 , cex = 0.5)
  mtext(names[[x]], 3, line = 0 , cex = 0.7)
  #title(ylab = expression(italic(paste(Delta,"T"))), line = 2.5 ,las = 2, cex = 0.7)
  
}

#####
#load sample species tracks
load("tracks_for_global_map.RData")

#imaginary raster for legend. I want the legend to go from -1 to 5 (range of values in the prediction rasters)
imaginary_r<-preds_filt[[1]]
imaginary_r@data@values[3:9]<- rep(5,7)
imaginary_r@data@values[10:16]<- rep(-1,7)


#create a color palette
cuts<-seq(-1,5,0.01) #set breaks
#pal <- colorRampPalette(c("dodgerblue","darkturquoise","goldenrod1","coral","firebrick1"))
pal <- colorRampPalette(c("dodgerblue","darkturquoise", "goldenrod1","coral","firebrick1","firebrick4"))
colpal <- pal(570)
  
#pdf("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/global_plot_8.pdf", width = 11, height = 5.2)

X11(width = 11, height = 5.2) #make the window proportional to region

par(mfrow=c(1,1),
    fig = c(0,1,0,1), #do this if you want to add the small plots as subplots
    bty="n", #no box around the plot
    cex.axis= 0.6, #x and y labels have 0.75% of the default size
    font = 3, #3: axis labels are in italics
    font.axis = 3,
    cex.lab = 0.6,
    #cex = 0.5,
    oma = c(0,0,0,0),
    mar = c(0, 0, 0, 0)
)

#maps::map("world",fill = TRUE, col = "grey30", border = F)
plot(region, col="#e5e5e5",border="#e5e5e5")
plot(preds_filt[[1]],axes = F, box=F, legend=FALSE,zlim=c(-1,5),breaks=cuts, col = colpal, add = T) 
plot(preds_filt[[2]],axes = F, box=F, legend=FALSE,zlim=c(-1,5),breaks=cuts, col = colpal, add = T) 
plot(preds_filt[[3]],axes = F, box=F, legend=FALSE,zlim=c(-1,5),breaks=cuts, col = colpal, add = T) 
plot(preds_filt[[4]],axes = F, box=F, legend=FALSE,zlim=c(-1,5),breaks=cuts, col = colpal, add = T) 

#add tracks
lwd <- 1.2
col <- "grey25"
plot(st_geometry(sp_samples$EF_sample), add= T, lty = 2, lwd = lwd, col = col)#lwd = 1, col = "violet")
plot(st_geometry(sp_samples$PF_sample), add= T, lty = 3, lwd = lwd, col = col)#lwd = 1, col = "darkblue")
plot(st_geometry(sp_samples$O_sample), add= T, lty = 5, lwd = lwd, col = col)#lwd = 1, col = "forestgreen")
plot(st_geometry(sp_samples$OHB_sample), add= T, lty = 4, lwd = lwd, col = col)#lwd = 1, col = "tan3")
plot(st_geometry(sp_samples$GFB_sample), add= T, lty = 1, lwd = lwd, col = col)#lwd = 1, col = "seagreen4")
plot(st_geometry(sp_samples$AF_sample), add= T, lty = 6, lwd = lwd, col = col)#lwd = 1, col = "seagreen4")

#add latitudes
clip(-130, 157, -50, 73)
abline(h = 0, col = "grey70",lty = 2)
abline(h = 30, col = "grey70",lty = 2)
abline(h = 60, col = "grey70",lty = 2)
text(x = -125, y = c(2,32,62), labels = c("0° ", "30° N", "60° N"), cex = 0.6, col = "grey65")

#add subplots...
centers_x <- c(124,-56,64,4) #distance between centers = 60

for(i in 1:length(centers_x)){
  rect(xleft = centers_x[i] - 27,
       xright = centers_x[i] + 27,
       ybottom =  -38,
       ytop = -2,
       col="white")
  
  subplot(reg_gam_plot(i), x = centers_x[i]+4,y = -19, size = c(1.6,0.7),  
          pars = list(mar=c(0,0,0.6,0),cex = 0.6, bty = "l", mgp = c(0,0.2,0),tck = 0.015, cex.main = 0.8, font.main = 3))
}


#add birds
rasterImage(AF, xleft = 85, xright = 110, ybottom = 20, ytop = 45, angle = +25) #make AF smaller than others
rasterImage(OHB, xleft = 145, xright = 175, ybottom = 25, ytop = 55, angle = +40)
rasterImage(GFB, xleft = 125, xright = 155, ybottom = 15, ytop = 45, angle = 0)
#rasterImage(GFB, xleft = 150, xright = 180, ybottom = 10, ytop = 40, angle = +70)
rasterImage(EF, xleft = 8, xright = 33, ybottom = 8, ytop = 33, angle = 0)
#rasterImage(EFc, xleft = 33, xright = 63, ybottom = 20, ytop = 50, angle = +35)
#rasterImage(EFc, xleft = 7, xright = 37, ybottom = 15, ytop = 45, angle = -50)
rasterImage(PF_E, xleft = 30, xright = 60, ybottom = 45, ytop = 75, angle = 15) 
rasterImage(PF_A, xleft = -65, xright = -35, ybottom = 40, ytop = 70, angle = 20) 
#rasterImage(OSc, xleft = 5, xright = 35, ybottom = 30, ytop = 60)
rasterImage(OS_E, xleft = -10, xright = 20, ybottom = 30, ytop = 60, angle = 20)
rasterImage(OS_A, xleft = -110, xright = -80, ybottom = 35, ytop = 65)

#add legend
plot(imaginary_r, legend.only = TRUE, breaks = cuts, col = colpal, #create an imaginary raster that has the high values from juv and low values from adlt
     #legend.width = 0.3, legend.shrink = 0.3,
     smallplot = c(0.06,0.17, 0.36,0.375), #c(min % from left, max % from left, min % from bottom, max % from bottom)
     #smallplot = c(0.06,0.17, 0.25,0.265),
     axis.args = list(at = c(-1,0,2,4), #same arguments as any axis, to determine the length of the bar and tick marks and labels
                    labels = c(-1,0,2,4), 
                    col = NA, #make sure box type in par is set to n, otherwise axes will be drawn on the legend :p
                    col.ticks = NA,
                    line = -1.3),
     horizontal = T,
     legend.args = list(text = expression(italic(paste(Delta,"T", "(°C)"))), side = 3, font = 2, line = 0.1, cex = 0.7))

text(x = -121.5,y = -17, "Map", cex = 0.7)
legend(-126,-17.5, legend = c("Oriental honey buzzard", "Grey-faced buzzard", "Amur falcon", 
                                  "Eleonora's falcon", "Peregrine falcon", "Osprey"),
       lty = c(4,1,6,2,3,5), cex = 0.55, bty = "n", seg.len = 3)
text(x = -118.5,y = -39.5, "Sub-plots", cex = 0.7)
legend(-126,-40, legend = c("High sun elevation", "Low sun elevation", "Night"),
       lty = c(1,2,3), cex = 0.55, bty = "n", seg.len = 3)
#dev.off()


## add routes
#choose one sample trajecotry for species with good data. unfiltered tracks in track_based_prep_analyze_daily.R
load("all_spp_unfiltered_updated_lc_0_removed_new_track_id.RData") #dataset


pf_ad <- read.csv("/home/enourani/ownCloud/Work/Projects/delta_t/data/Osprey_Americas/Peregrines Ivan.csv", stringsAsFactors = F) %>% 
  filter(Age == "ad")
ao_ad <-read.csv("/home/enourani/ownCloud/Work/Projects/delta_t/data/Osprey_Americas/ROB mig data190411.csv", stringsAsFactors = F) %>% 
  filter(Age2 == "a")

GFB_files <- list.files("/home/enourani/ownCloud/Work/Projects/delta_t/data/Grey_faced_buzzard/",pattern = ".csv",full.names = T)
GFB <- lapply(GFB_files,read.csv,stringsAsFactors = F) %>%
  reduce(full_join) %>% #is locdate is in UTC
  mutate(locdate,date_time = as.POSIXct(strptime(locdate,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>%
  mutate(month = month(date_time),
         year = year(date_time),
         season = ifelse(month %in% c(3,4),"spring",ifelse(month %in% c(8:10),"autumn","other"))) %>% 
  mutate(track = paste(platform,year,season,sep = "_"),
         species = "GFB") %>% 
  rename(location.long = lon,
         location.lat = lat) %>% 
  filter(season %in% c("spring","autumn") &
           class %in% c("1","2","3")) #filter for location classes


PF <- read.csv("/home/enourani/ownCloud/Work/Projects/delta_t/data/LifeTrack Peregrine falcon.csv", stringsAsFactors = F) %>% 
  dplyr::select(1,3:5,16,38,39) %>% #remove columns that are not needed
  filter(individual.local.identifier %in% pf_ad$animal.id) %>% 
  mutate(date_time = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  mutate(month = month(date_time),
         year = year(date_time),
         species = "PF",
         season = ifelse (month %in% c(8:11), "autumn", "other")) %>% 
  mutate(track = paste(tag.local.identifier, year, season,sep = "_")) %>% 
  filter(season != "other")

OE <- read.csv("/home/enourani/ownCloud/Work/Projects/delta_t/data/Osprey in Mediterranean (Corsica, Italy, Balearics).csv", stringsAsFactors = F) %>% 
  dplyr::filter(grepl("ad",individual.local.identifier,, ignore.case = T) & !grepl("juv",individual.local.identifier,, ignore.case = T)) %>% #extract adult data
  dplyr::select(1,3:5,16,35:37) %>% #remove columns that are not needed
  mutate(date_time = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  mutate(month = month(date_time),
         year = year(date_time),
         species = "O",
         season = ifelse(month %in% c(2:4),"spring",ifelse (month %in% c(8:11), "autumn","other"))) %>% 
  mutate(track = paste(tag.local.identifier, year,season, sep = "_")) %>% 
  filter(season != "other")

OA <- read.csv("/home/enourani/ownCloud/Work/Projects/delta_t/data/Osprey_Americas/Osprey Bierregaard North and South America.csv", stringsAsFactors = F) %>% 
  dplyr::select(1,3:5,13,48:52) %>% #remove columns that are not needed
  filter(sensor.type == "gps" | sensor.type == "argos-doppler-shift" & argos.lc %in% c("1","2","3"),
         individual.local.identifier %in% ao_ad$Bird) %>%
  mutate(date_time = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  mutate(month = month(date_time),
         year = year(date_time),
         species = "O",
         season = ifelse(month %in% c(3,4),"spring",ifelse (month %in% c(8:11), "autumn", "other"))) %>% 
  mutate(track = paste(tag.local.identifier, year,season,sep = "_")) %>% 
  filter(season != "other")

dataset <- list(GFB, PF, OE, OA) %>% 
  reduce(full_join, by = c("location.long", "location.lat", "date_time", "track", "month", "year" , "season", "species")) %>% 
  dplyr::select(c("location.long", "location.lat", "date_time", "track", "month", "year" , "season", "species")) %>% 
  drop_na() %>% 
  st_as_sf(coords = c("location.long","location.lat"), crs = wgs)

GFB_sample <- dataset[dataset$track == "93918_2009_autumn",]  %>% 
  summarize(do_union = FALSE) %>% 
  st_cast("LINESTRING")

PF_sample <- dataset[dataset$track %in% c("4390_2016_autumn", "6376_2019_autumn"),] %>% 
  group_by(track) %>% 
  summarize(do_union = FALSE) %>% 
  st_cast("LINESTRING")

O_sample <- dataset[dataset$track %in% c("FOSP04_2013_autumn","18_2012_autumn"),] %>% 
  group_by(track) %>% 
  summarize(do_union = FALSE) %>% 
  st_cast("LINESTRING")

OHB_sample <- read.csv("/home/enourani/ownCloud/Work/Projects/delta_t/data/Tracking_of_the_migration_of_Oriental_Honey_Buzzards.csv") %>% 
  mutate(date_time = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  mutate(month = month(date_time),
         year = year(date_time),
         species = "OHB",
         yday = yday(date_time)) %>%
  mutate(season = ifelse(month %in% c(8:11),"autumn",ifelse(month == 5,"spring","other")),  #11 Sep-20 Oct; spring between 1-5 May
         track = paste(tag.local.identifier, year,season,sep = "_")) %>% 
  filter(season == "autumn" & track == "ngsku1701_2019_autumn") %>% 
  st_as_sf(coords = c("location.long","location.lat"), crs = wgs) %>% 
  summarize(do_union = FALSE) %>% 
  st_cast("LINESTRING")

EF_sample <- read_excel("/home/enourani/ownCloud/Work/Projects/delta_t/data/eleonoras_falcon.xlsx", sheet = 2) %>% 
  mutate(date_time = as.POSIXct(strptime(`timestamp (UTC)`,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>%
  rename(timestamp = `timestamp (UTC)`) %>% 
  mutate(month = month(date_time),
         year = year(date_time),
         species = "EF",
         season = ifelse(month %in% c(3,4),"spring",ifelse (month %in% c(8:11), "autumn", "other"))) %>% 
  mutate(track = paste(ID.individual, year,season,sep = "_")) %>% 
  filter(season == "autumn" & track == "ELEF01_2015_autumn") %>% 
  st_as_sf(coords = c("location.long","location.lat"), crs = wgs) %>% 
  summarize(do_union = FALSE) %>% 
  st_cast("LINESTRING")

AF_sample <- data.frame(x = c(124.742,109.572,94.562,80.410,72.987,66.843,59.596,51.425,46.154,37.689,34.899,28.464),
                        y = c(48.682,35.110,26.158,17.127,15.101,13.156,10.490,8.714,2.661,-6.842,-13.042,-24.980)) %>%  #I made this up based on figures of amur faclon migration route
  st_as_sf(coords = c("x","y"),crs = wgs) %>% 
  summarize(do_union = FALSE) %>% 
  st_cast("LINESTRING")
  

sp_samples <- list(AF_sample,EF_sample,OHB_sample,O_sample,PF_sample,GFB_sample)
names(sp_samples) <- c("AF_sample","EF_sample","OHB_sample","O_sample","PF_sample","GFB_sample")

save(sp_samples, file = "tracks_for_global_map.RData")


#read in bird art
AF <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/AF_1.png")
EF <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/Updated/EF-2.png")
GFB <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/Updated/GFB-2.png")
OHB <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/OHB_2.png")
OS_A <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/Updated/EO-A.png")
OS_E <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/Updated/EO-B.png") #mirrored
PF_A <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/Updated/EP-A.png")
PF_E <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/Updated/EP-B.png")












#add the subset plots as images
ea <- readTIFF("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/East_asia.tif")
rasterImage(ea, 120,35,160,50)


#####################################
#add subset plots
#note on fig: this is the coordinates of where you want the small fig. assumes that coordinates are from 0 to 1. so, use locator to find the location
#on the map, then convert it to 0-1 scale.

rect(147,46,127,60,col="white") #add background and border for the plot

par(mfrow=c(1,1), 
    mar=c(0,0,0,0),
    fig = c((abs(st_bbox(region)[1]) + 120)/(st_bbox(region)[3] - st_bbox(region)[1]), 
            (abs(st_bbox(region)[1]) + 145)/(st_bbox(region)[3] - st_bbox(region)[1]), 
            (abs(st_bbox(region)[2]) + 45)/(st_bbox(region)[4] - st_bbox(region)[2]), 
            (abs(st_bbox(region)[2]) + 60)/(st_bbox(region)[4] - st_bbox(region)[2])),
    new = T,
    cex = 0.6,
    bty = "o",
    mgp = c(0,0.15,0), #margin line for axis labels and axis line
    tck = -0.03
)
reg_gam_plot(1)

par(mfrow=c(1,1),  #the problem is that the sencond inset plot will be on the first inset plot.... hmmm
    mar=c(0,0,0,0),
    fig = c((abs(st_bbox(region)[1]) + 100)/(st_bbox(region)[3] - st_bbox(region)[1]), 
            (abs(st_bbox(region)[1]) + 145)/(st_bbox(region)[3] - st_bbox(region)[1]), 
            (abs(st_bbox(region)[2]) + 25)/(st_bbox(region)[4] - st_bbox(region)[2]), 
            (abs(st_bbox(region)[2]) + 60)/(st_bbox(region)[4] - st_bbox(region)[2])),
    new = T,
    cex = 0.6,
    bty = "o",
    mgp = c(0,0.15,0), #margin line for axis labels and axis line
    tck = -0.03
)
reg_gam_plot(1)

#try TeachingDemos:: subplot


### STEP 4: make effect plots ####
#load data and models
load("data_ls_regional_gam.RData")
load("models_ls_reg_GAMs.RData")

#plot delta_t ~ yday

#use the itsadug package
#summed effects (incl. intercept) mgcv default plots are partial.

#save the plots as jpeg images
lapply(c(names(models_ls)), function(x){
  m <- models_ls[[x]]
  t <- timing_areas[[x]]
  
  tiff(paste0("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/", x, ".tif"), 
       res = 600, width = 3.5, height = 2.5, units = "in")
  
  #X11(width = 3.5, height = 2.5)
  par(mfrow=c(1,1), bty="l", #box around the plot
      cex = 0.7
  )
  
  plot(0, type = "l", labels = FALSE, tck = 0, xlim =  c(1,366), ylim = c(-2.5,5), xlab = "Day of the year", ylab = "") #expression(paste(Delta,"T")), main = x
  rect(xleft = min(t),ybottom = -2.7,xright = max(t),ytop = 5, col="#99CC0060",border=NA) #water-crossing window
  plot_smooth(m, view="yday", plot_all="sun_elev_f", rm.ranef=F, lwd = 1.5, #ylim=c(-2,5),
              col = "grey60", hide.label = TRUE, 
              legend_plot_all =  F, 
              h0 = NULL, add = T, lty = c(1,5,3))
  
  axis(side = 1, at = c(1,seq(50,350,50)), line = 0, labels = c(1,seq(50,350,50)), 
       tick = T , col.ticks = 1, col = NA, tck = -.015)
  axis(side = 2, at = c(-2,0,2,4), line = 0, labels = c(-2,0,2,4),
       tick = T , col.ticks = 1,col = NA, tck = -.015, 
       las = 2) # text perpendicular to axis label 
  mtext(expression(italic(paste(Delta,"T"))), 2, line = 2 ,las = 2.5, cex = 0.7)
  
  

  dev.off()
})



##########################
#if you want to use the par(fig ) option, create a plotting function. placement of the plot is super tricky...
reg_gam_plot <- function(x){
  m <- models_ls[[x]]
  t <- timing_areas[[x]]
  
  plot(0, type = "n", labels = FALSE, tck = 0, xlim =  c(1,366), ylim = c(-2.5,5), xlab = "", ylab = "") #expression(paste(Delta,"T")), main = x
  rect(xleft = min(t),ybottom = -2.7,xright = max(t),ytop = 5, col="#99CC0060",border=NA) #water-crossing window
  plot_smooth(m, view="yday", plot_all="sun_elev_f", rm.ranef=F, lwd = 1.5, #ylim=c(-2,5),
              col = "grey60", hide.label = TRUE, 
              legend_plot_all =  F, 
              h0 = NULL, add = T, lty = c(1,5,3))
  
  axis(side = 1, at = c(100, 200, 300), line = 0, labels = c(100,200,300), 
       tick = T , col.ticks = 1, col = NA, lty = NULL, tck = -.015)
  axis(side = 2, at = c(-2,0,2,4), line = 0, labels = c(-2,0,2,4),
       tick = T , col.ticks = 1,col = NA, lty = NULL, tck = -.015, 
       las = 2) # text perpendicular to axis label 
  #mtext(expression(italic(paste(Delta,"T"))), 2, line = 2 ,las = 2.5, cex = 0.7)
  #title(ylab = expression(italic(paste(Delta,"T"))), line = 2.5 ,las = 2, cex = 0.7)
  
}

###############################

#prepare data to make predictions for
preds_raw <- lapply(c(names(models_ls)),function(x){
  d <- data_ls[[x]]
  m <- models_ls[[x]]
  
  #make sure to include all combinations of sun elev and yday
  pdat <- with(d, expand.grid(sun_elev_f = levels(sun_elev_f),
                              yday = seq(min(yday),max(yday), length = 100)))
  #add a random sample of the other variables             
  pdat <- pdat %>% 
    mutate(lat = sample(c(min(d$lat):max(d$lat)), 300, replace = T),
           lon = sample(c(min(d$lon):max(d$lon)), 300, replace = T),
           year = sample(c(min(d$year):max(d$year)), 300, replace = T))
  
  #make predictions
  pdat <- transform(pdat, pred = predict(m, newdata = pdat, type = "response"))
  pdat
  
  #plot
  ylim <- with(d,range(delta_t))
  
  plot(delta_t ~ yday, data = d)
  levs <- levels(d$sun_elev_f)
  cols <- c("red","green","blue")
  
  for (i in seq_along(levs)){
  dd <- subset(pdat, sun_elev_f == levs[i])
  lines(pred ~ yday, data = dd, col = cols[[i]])
  }
    
  
  })

names(preds_raw) <- names(models_ls)



#extract data for creating the plots
plotdata <- lapply(models_ls, plot,pages = 1)

X11()
plot(x = sample(c(1:365), 100), y = sample(c(-1:1),100, replace = T), type = "n", labels = TRUE, tck = 0, xlab = "Yday", ylab = "delta_t")
#add night
lines(plotdata$East_asia[[6]]$x,plotdata$East_asia[[6]]$fit) 
lines(plotdata$East_asia[[6]]$x ,plotdata$East_asia[[6]]$fit + plotdata$East_asia[[6]]$se, lty = 2)
lines(plotdata$East_asia[[6]]$x ,plotdata$East_asia[[6]]$fit - plotdata$East_asia[[6]]$se, lty = 2)
#add low sun
lines(plotdata$East_asia[[5]]$x,plotdata$East_asia[[5]]$fit) 
lines(plotdata$East_asia[[5]]$x ,plotdata$East_asia[[5]]$fit + plotdata$East_asia[[6]]$se, lty = 2)
lines(plotdata$East_asia[[5]]$x ,plotdata$East_asia[[5]]$fit - plotdata$East_asia[[6]]$se, lty = 2)
#add high sun
lines(plotdata$East_asia[[4]]$x,plotdata$East_asia[[4]]$fit) 
lines(plotdata$East_asia[[4]]$x ,plotdata$East_asia[[4]]$fit + plotdata$East_asia[[6]]$se, lty = 2)
lines(plotdata$East_asia[[4]]$x ,plotdata$East_asia[[4]]$fit - plotdata$East_asia[[6]]$se, lty = 2)



#use predicted values to manually plot
#estimate raw predictions (what I have above is rasters)




plot(pred ~ yday, data = preds_raw$East_asia[preds_raw$East_asia$sun_elev_f == "night",])


#old subplot locations
#add subplots
#rect(95,38,130,60,col="white")#rect(122,38,157,60,col="white")
rect(111,-10,147,-33,col="white")
subplot(reg_gam_plot(1), x = 129, y = -21, size = c(1,0.5), #x = 113, y = 50 
        pars = list(mar=c(0,0,0.6,0),cex = 0.6, bty = "l", mgp = c(0,0.2,0),tck = 0.015, cex.main = 0.8, font.main = 3))
rect(50,-18,85,4,col="white")
subplot(reg_gam_plot(3), x = 68, y = -6, size = c(1,0.5), 
        pars = list(mar=c(0,0,0.6,0),cex = 0.6, bty = "l", mgp = c(0,0.2,0),tck = 0.015, cex.main = 0.8, font.main = 3))
rect(-38,46,-2,69,col="white")#rect(35,40,69,62,col="white")
subplot(reg_gam_plot(4), x = -20, y = 58, size = c(1,0.5), 
        pars = list(mar=c(0,0,0.6,0),cex = 0.6, bty = "l", mgp = c(0,0.2,0),tck = 0.015, cex.main = 0.8, font.main = 3))
rect(-50,-1,-14,21,col="white")#rect(-120,34,-87,56,col="white")        
subplot(reg_gam_plot(2), x = -32, y = 10, size = c(1,0.5), 
        pars = list(mar=c(0,0,0.6,0),cex = 0.6, bty = "l", mgp = c(0,0.2,0),tck = 0.015, cex.main = 0.8, font.main = 3))

#second round
rect(98,-17,134,-39,col="white")
subplot(reg_gam_plot(1), x = 116, y = -28, size = c(1,0.5),  
        pars = list(mar=c(0,0,0.6,0),cex = 0.6, bty = "l", mgp = c(0,0.2,0),tck = 0.015, cex.main = 0.8, font.main = 3))
rect(42,-17,78,-39,col="white")
subplot(reg_gam_plot(3), x = 60, y = -28, size = c(1,0.5), 
        pars = list(mar=c(0,0,0.6,0),cex = 0.6, bty = "l", mgp = c(0,0.2,0),tck = 0.015, cex.main = 0.8, font.main = 3))
rect(-8,-17,28,-39,col="white")
subplot(reg_gam_plot(4), x = 10, y = -28, size = c(1,0.5), 
        pars = list(mar=c(0,0,0.6,0),cex = 0.6, bty = "l", mgp = c(0,0.2,0),tck = 0.015, cex.main = 0.8, font.main = 3))
rect(-85,-17,-49,-39,col="white")       
subplot(reg_gam_plot(2), x = -67, y = -28, size = c(1,0.5), 
        pars = list(mar=c(0,0,0.6,0),cex = 0.6, bty = "l", mgp = c(0,0.2,0),tck = 0.015, cex.main = 0.8, font.main = 3))