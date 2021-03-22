#code for running gams for each region separately. mostly copied from regional_gams.R
#Feb. 22, 2021

library(mgcv)
library(tidyverse)
library(sp)
library(sf)
library(raster)
library(maptools)
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
wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"



### --------STEP 1: determine the extents for each region ####

#spatial extent of ssf data
load("2021/ssf_input_all_df_1hr.RData") #used_av_df_1hr
y_extent <- used_av_df_1hr %>%
  filter(used == 1) %>% 
  summarise(max_lat = round(max(y)),
            min_lat = round(min(y)))

#create a dummy layer
dummy <- st_polygon(list(rbind(c(-180,-180), c(180,-180), c(180,180), c(-180,180),c(-180,-180)))) %>% 
  st_sfc(crs = wgs)

# --- East Asia
East_asia <- dummy %>% 
  st_crop(c(xmin = 99, xmax = 130, ymin = 1.7, ymax = 41))

# --- Indian ocean
Indian_ocean <- dummy %>% 
  st_crop(xmin = 43, xmax = 78, ymin = 9, ymax = 26) #%>% 

# --- American continent
Americas <- dummy %>% 
  st_crop(xmin = -98, xmax = -54, ymin = 5.7, ymax = 47) #%>% 
  #st_difference(np_ocean)

# --- Europe
Europe <- dummy %>% 
  st_crop(xmin = -11, xmax = 50.5, ymin = 30, ymax = 69)

# --- Madagascar
Madagascar <- dummy %>% 
  st_crop(xmin = 33, xmax = 49, ymin = -26, ymax = -9)

mapview(East_asia) + mapview(Indian_ocean) + mapview(Americas) + mapview(Europe) + mapview(Madagascar)

#put all regional data in a list
extent_ls <- list(East_asia = st_bbox(East_asia), Americas = st_bbox(Americas), Indian_ocean = st_bbox(Indian_ocean),
                Europe = st_bbox(Europe), Madagascar = st_bbox(Madagascar))

save(extent_ls, file = "2021/extent_ls_regional_gam.RData") #use this later in 2021_Era_interim_download_prep.R


### --------STEP 2: open data, spatio-termporal filters #####

load("2021/ecmwf_regions_5ksample.RData") #called data_df (from 2021_Era_interim_download_prep.R)
load("2021/ocean.RData") #ocean (prepared in 2021_all_data_prep_analyze.R)

regions_ls <- split(data_df, data_df$region)

regions_ls <- lapply(regions_ls, function(x){
  x %>% 
    st_as_sf(coords = c("lon","lat"), crs = wgs)
})

save(regions_ls, file = "2021/ecmwf_regions_sf_ls.RData")

load("2021/ecmwf_regions_sf_ls.RData") #regions_ls

list2env(regions_ls, envir = .GlobalEnv)

#for East Asia, remove outlier (identified through plotting only one year of data using mapview)
outlier <- East_asia[4904,]
East_asia <- East_asia %>% 
  st_difference(outlier)

save(East_asia, file = "2021/East_asia_sf.RData")

#for indian ocean, remove  persian gulf and the red sea
io_to_delete <- st_read("/home/enourani/ownCloud/Work/GIS_files/World_Seas_IHO_v3/World_Seas_IHO_v3.shp") %>% 
  filter(NAME %in% c("Persian Gulf", "Red Sea")) %>%  
  st_union() %>%
  st_transform(meters_proj) %>% 
  st_buffer(dist = units::set_units(50000, 'm')) %>% 
  st_transform(wgs)

(b <- Sys.time())
Indian_ocean <- Indian_ocean %>%  
  st_difference(io_to_delete) #%>%  #remove points over the persian gulf and the red sea.
  #st_intersection(ocean)
Sys.time() - b #14 min
  
save(Indian_ocean, file = "2021/Indian_ocean_sp.RData")


#for the Americas, get rid of the pacific ocean
np_ocean <- st_read("/home/enourani/ownCloud/Work/GIS_files/World_Seas_IHO_v3/World_Seas_IHO_v3.shp") %>% 
  filter(NAME == "North Pacific Ocean") %>% 
  st_crop(xmin = -118, xmax = -60, ymin = 0, ymax = 29) 

Americas <- Americas %>% 
  st_difference(np_ocean) 

#crop ocean to the extent of americas
ocean_am <- ocean %>% 
  st_crop(st_bbox(Americas))

Americas <- Americas %>% 
  st_intersection(ocean_am)

save(Americas, file = "2021/Americas_sf.RData")

#for Europe, include mediterranean and the baltic, white sea, and Barentsz Sea
bzs <- st_read("/home/enourani/ownCloud/Work/GIS_files/barentsz_sea/iho.shp") %>% 
  st_union()

ws <- st_read("/home/enourani/ownCloud/Work/GIS_files/white_sea/iho.shp") %>% 
  st_union()

load("eur_sea.RData") #eur_sea from regional_gams.R

eur_updated <- eur_sea %>% 
  st_union(bzs) %>% 
  st_union(ws)

save(eur_updated, file = "2021/eur_sea.RData")

#mask europe layer to keep only waterbodies of itnerest   
Europe <- Europe %>% 
 # st_crop(xmin = -11, xmax = 37, ymin = 30, ymax = 68.5) %>% 
  st_intersection(eur_updated) #

save(Europe, file = "2021/Europe_sf.RData")

#for Madagascar, only keep the Mozambique channel
mch <- st_read("/home/enourani/ownCloud/Work/GIS_files/Mozambique_Channel/iho.shp") %>% 
  st_union()

(b <- Sys.time())
Madagascar <- Madagascar %>% 
  st_intersection(mch)
Sys.time() - b #1.497251 hours

save(Madagascar, file = "2021/madagascar_sf.RData")


#put all regional data in a list
load("2021/Indian_ocean_sp.RData")
load("2021/Americas_sf.RData")
load("2021/Europe_sf.RData")
load("2021/madagascar_sf.RData")

data_ls <- list(East_asia = East_asia, Americas = Americas, Indian_ocean = Indian_ocean, Europe = Europe, Madagascar = Madagascar)

# calculate sun position

data_ls_sun <- lapply(data_ls, function(x){
  data <- x %>%
    dplyr::select(c("date_time", "sst","t2m","delta_t", "yday", "hour","year","region","geometry")) %>% 
    mutate(s_elev_angle = solarpos(st_coordinates(.), date_time, proj4string=CRS("+proj=longlat +datum=WGS84"))[,2]) %>% #calculate solar elevation angle
    mutate(sun_elev = ifelse(s_elev_angle < -6, "night", #create a categorical variable for teh position of the sun
                             ifelse(s_elev_angle > 40, "high", "low")),
           month = month(date_time)) %>% 
    as("Spatial") %>% 
    as.data.frame() %>% 
    rename(lon = coords.x1,
           lat = coords.x2)
  
})

save(data_ls_sun, file = "2021/data_ls_regional_gam.RData")



### --------STEP 3: model #####
load("2021/data_ls_regional_gam.RData")

mycl <- makeCluster(5) 

clusterExport(mycl, "data_ls_sun") 

clusterEvalQ(mycl, {
  library(mgcv)
})

(b <- Sys.time())
models_ls <- lapply(data_ls_sun, function(x){
  
  x$sun_elev_f <- as.factor(x$sun_elev)
  
  gamm(delta_t ~ s(lat,lon, by = sun_elev_f, k = 100) +
        s(yday, by = sun_elev_f, bs = "cc") +
        s(year, bs = "re") +
        sun_elev_f , method = "REML", data = x, 
      weights = varPower(form = ~lat))

})

Sys.time() - b #8.8 hours
 
stopCluster(mycl)

save(models_ls, file = "2021/models_ls_reg_GAMs.RData")


load("2021/models_ls_reg_GAMs.RData")

#model checking

lapply(models_ls, function(x){
  X11()
  par(mfrow= c(2,2), oma = c(0,0,3,0))
  gam.check(x$gam)
})

#create output tables in Latex
latex_ls <- lapply(names(models_ls), function(x){
  m <- models_ls[[x]]
  gamtabs(m$gam, caption = x)
})

### --------STEP 4: predictions #####
#extract migration timing OVER THE SEA for each species
load("2021/all_2009_2020_overwater_points.RData") #all_oversea. prepped in 2021_all_data_prep_analyze.R

all_oversea %>% 
  filter(season == "autumn") %>%
  as("Spatial") %>% 
  as.data.frame() %>% 
  mutate(continent = ifelse(coords.x1 < -24, "A", #america
                            ifelse(species == "EF" & coords.x2 < 0, "S.Africa","not_relevant"))) %>% #southern hemisphere for EF
  group_by(continent,species) %>% 
  summarise(min = min(yday(date_time)),
            max = max(yday(date_time)))
  

timing <- list(OHB = c(260:294), #sea-crossing; ref: almost Yamaguchi
               GFB = c(277:299),
               AF = c(319:323), #mid-november. ref: Bernd's poster
               EF_E = c(288:301),
               EF_A = c(303:339),
               PF_EU = c(261:289),
               O_EU = c(222:277),
               PF_A = c(279:305),
               O_A = c(244:299))

timing_areas <- list(East_asia = c(min(c(timing$GFB,timing$OHB)):max(c(timing$GFB,timing$OHB))),
                     Americas = c(min(c(timing$O_A,timing$PF_A)):max(c(timing$O_A,timing$PF_A))),
                     Indian_ocean = timing$AF,
                     Europe = c(min(c(timing$O_EU,timing$PF_EU,timing$EF_E)):max(c(timing$O_EU,timing$PF_EU,timing$EF_E))),
                     Madagascar = timing$EF_A) 


save(timing_areas, file = "2021/timing_for_gam_preds.RData")

#make predictions

load("2021/data_ls_regional_gam.RData") #data_ls_sun
load("2021/models_ls_reg_GAMs.RData") #models_ls
load("2021/timing_for_gam_preds.RData") #timing_areas

#these inputs are too large. don't have enough ram to run multi-core. doesn't take that long anyway

(b <- Sys.time())

preds <- lapply(c(names(models_ls)),function(x){
  d <- data_ls_sun[[x]] %>% 
    filter(yday %in% timing_areas[[x]]) %>% 
    mutate(sun_elev_f = as.factor(sun_elev))
  
  m <- models_ls[[x]]
  
  pred <- data.frame(pred = as.numeric(predict(m$gam, d)), lon = d$lon, lat = d$lat)
  
  coordinates(pred) <- ~ lon + lat
  proj4string(pred) <- wgs
  
  pred_sp <- SpatialPixelsDataFrame(pred, tolerance = 0.916421, pred@data)
  r <- raster(pred_sp)

  #interpolate. for visualization purposes
  surf.1 <- Tps(as.matrix(as.data.frame(r,xy = T)[,c(1,2)],col = 2),as.data.frame(r,xy = T)[,3])
  
  grd <- expand.grid(x = seq(from = extent(r)[1],to = extent(r)[2],by = 0.1),
                   y = seq(from = extent(r)[3],to = extent(r)[4],by = 0.1))

  grd$coords <- matrix(c(grd$x,grd$y),ncol=2)
  
  surf.1.pred <- predict.Krig(surf.1,grd$coords)
  interpdf <- data.frame(grd$coords,surf.1.pred)
  
  colnames(interpdf)<-c("lon","lat","delta_t")
  
  coordinates(interpdf) <- ~ lon + lat
  gridded(interpdf) <- TRUE
  interpr <- raster(interpdf)
  proj4string(interpr) <- wgs
  
  return(interpr)
  
})

Sys.time() -b # 58 sec



names(preds) <- names(models_ls)
save(preds, file = "2021/regional_gam_preds.RData")

#mask the rasters with the relevant land layer. or simply plot land over the top
#coastlines layer
region <- st_read("/home/enourani/ownCloud/Work/GIS_files/continent_shapefile/continent.shp") %>% 
  st_crop(xmin = -130, xmax = 157, ymin = -50, ymax = 70) %>%
  st_union()

#mask all with the landmass layer
preds_filt <- lapply(preds, function(x){
  x_f <- raster::mask(x,as(region,"Spatial"), inverse = T) 
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
load("2021/eur_sea.RData") #eur_updated
eur_sea <- eur_updated %>% 
  st_crop(xmin = -5.425, xmax = 50.875 , ymin = 25.825, ymax = 70) 
  
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
preds_filt$East_asia <- crop(preds_filt$East_asia, extent(99, 130.175, 1.275 , 41.075 ))


#Madagascar
mch <- st_read("/home/enourani/ownCloud/Work/GIS_files/Mozambique_Channel/iho.shp") %>% 
  st_union()

preds_filt$Madagascar <- mask(preds_filt$Madagascar, as(mch,"Spatial"), inverse = F)

save(preds_filt, file = "2021/predictions_regional_gam_map.RData")



### --------STEP 5: plot the map#####
load("2021/predictions_regional_gam_map.RData") #preds_filt
load("2021/models_ls_reg_GAMs.RData") #models_ls
load("2021/timing_for_gam_preds.RData") #timing_areas

#elements prepared in regional_gam.R
load("tracks_for_global_map.RData") #sp_samples
AF <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/AF_1.png")
EF <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/Updated/EF-2.png")
GFB <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/Updated/GFB-2.png")
OHB <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/OHB_2.png")
OS_A <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/Updated/EO-A.png")
OS_E <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/Updated/EO-B.png") #mirrored
PF_A <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/Updated/EP-A.png")
PF_E <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/Updated/EP-B.png")



region <- st_read("/home/enourani/ownCloud/Work/GIS_files/continent_shapefile/continent.shp") %>% 
  st_crop(xmin = -130, xmax = 158, ymin = -74, ymax = 71) %>%
  st_union()

# create a plotting function for effect plots
names <- c("South-East Asia", "The Americas", "Indian Ocean", "Europe", "Mozambique Channel")

reg_gam_plot <- function(x){
  m <- models_ls[[x]]$gam
  t <- timing_areas[[x]]
  
  plot(0, type = "n", labels = FALSE, tck = 0, xlim =  c(1,366), ylim = c(-2.5,5), xlab = "", ylab = "")#, main = names[[x]]) #expression(paste(Delta,"T")), main = x
  abline(h = 0, col = "gray60",lty = 1, lwd = 0.5)
  rect(xleft = min(t),ybottom = -2.7,xright = max(t),ytop = 5, col="#99CC0060",border=NA) #water-crossing window
  plot_smooth(m, view = "yday", plot_all = "sun_elev_f", rm.ranef = F, lwd = 1.5, #ylim=c(-2,5),
              col = "grey50", hide.label = TRUE, 
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
###

#imaginary raster for legend. I want the legend to go from -1 to 5 (range of values in the prediction rasters)
imaginary_r<-preds_filt[[1]]
imaginary_r@data@values[3:9]<- rep(5,7)
imaginary_r@data@values[10:16]<- rep(-1,7)


#create a color palette
cuts <- seq(-1,5,0.01) #set breaks
#pal <- colorRampPalette(c("dodgerblue","darkturquoise","goldenrod1","coral","firebrick1"))
pal <- colorRampPalette(c("dodgerblue","darkturquoise", "goldenrod1","coral","firebrick1","firebrick4"))
colpal <- pal(570)
  
#pdf("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/2021/global_plot.pdf", width = 11, height = 6)

X11(width = 11, height = 6) #make the window proportional to region

par(mfrow=c(1,1),
    fig = c(0,1,0,1), #do this if you want to add the small plots as subplots
    bty="n", #no box around the plot
    cex.axis= 0.6, #x and y labels have 0.75% of the default size
    font = 3, #3: axis labels are in italics
    font.axis = 3,
    cex.lab = 0.6,
    #cex = 0.5,
    oma = c(0,0,0,0),
    mar = c(0, 0, 0, 0),
    lend = 1  #rectangular line endings (trick for adding the rectangle to the legend)
)

#maps::map("world",fill = TRUE, col = "grey30", border = F)
plot(region, col="#e5e5e5",border="#e5e5e5")
plot(preds_filt[[1]],axes = F, box=F, legend=FALSE,zlim=c(-1,5),breaks=cuts, col = colpal, add = T) 
plot(preds_filt[[2]],axes = F, box=F, legend=FALSE,zlim=c(-1,5),breaks=cuts, col = colpal, add = T) 
plot(preds_filt[[3]],axes = F, box=F, legend=FALSE,zlim=c(-1,5),breaks=cuts, col = colpal, add = T) 
plot(preds_filt[[4]],axes = F, box=F, legend=FALSE,zlim=c(-1,5),breaks=cuts, col = colpal, add = T) 
plot(preds_filt[[5]],axes = F, box=F, legend=FALSE,zlim=c(-1,5),breaks=cuts, col = colpal, add = T) 

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
#text(x = -125, y = c(2,32,62), labels = c("0° ", "30° N", "60° N"), cex = 0.6, col = "grey65")
text(x = -125, y = c(32,62), labels = c("30° N", "60° N"), cex = 0.6, col = "grey65")

#add a frame for the sub-plots and legend
rect(xleft = -130,
     xright = 159,
     ybottom =  -74.5,
     ytop = -28,
     col="white",
     border = NA)

rect(xleft = -130,
     xright = -87,
     ybottom =  -15,
     ytop = 7,
     col="white",
     border = NA)

#add subplots...
#centers_x <- c(124,-56,64,4) #distance between centers = 60
centers_x <- c(130,-102,72,-44,14) #distance between centers = 58

for(i in 1:length(centers_x)){
  rect(xleft = centers_x[i] - 27.5,
       xright = centers_x[i] + 28,
       ybottom =  -66,#-42
       ytop = -30, #-4
       col="white")
  
  subplot(reg_gam_plot(i), x = centers_x[i] + 4,y = -47, size = c(1.6,0.7),  #-23 was -19
          pars = list(mar=c(0,0,0.6,0),cex = 0.6, bty = "l", mgp = c(0,0.2,0),tck = 0.015, cex.main = 0.8, font.main = 3))
}


#add birds... try to scale them
rasterImage(AF, xleft = 95, xright = 111, ybottom = 29, ytop = 45, angle = +25) #make AF smaller than others
rasterImage(OHB, xleft = 145, xright = 180, ybottom = 20, ytop = 55, angle = +40)
rasterImage(GFB, xleft = 130, xright = 155, ybottom = 20, ytop = 45, angle = 0)
rasterImage(EF, xleft = 11, xright = 31, ybottom = 10, ytop = 30, angle = 0)
rasterImage(PF_E, xleft = 28, xright = 48, ybottom = 50, ytop = 70, angle = 15) 
rasterImage(PF_A, xleft = -65, xright = -45, ybottom = 40, ytop = 60, angle = 20) 
rasterImage(OS_E, xleft = -10, xright = 20, ybottom = 30, ytop = 60, angle = 20)
rasterImage(OS_A, xleft = -110, xright = -80, ybottom = 35, ytop = 65)


#add legend
plot(imaginary_r, legend.only = TRUE, breaks = cuts, col = colpal, 
     #legend.width = 0.3, legend.shrink = 0.3,
     smallplot = c(0.04,0.15,0.37,.385),#c(0.059,0.169,0.145,0.16), #c(min % from left, max % from left, min % from bottom, max % from bottom)
     #smallplot = c(0.06,0.17, 0.28,0.295),#c(0.06,0.17, 0.36,0.375), 
     axis.args = list(at = c(-1,0,2,4), #same arguments as any axis, to determine the length of the bar and tick marks and labels
                    labels = c(-1,0,2,4), 
                    col = NA, #make sure box type in par is set to n, otherwise axes will be drawn on the legend :p
                    col.ticks = NA,
                    line = -1.3),
     horizontal = T,
     legend.args = list(text = expression(italic(paste(Delta,"T", "(°C)"))), side = 3, font = 2, line = 0.1, cex = 0.7)
     )

text(x = -118,y = 10, "Map legend", cex = 0.8)
legend(-130,8, legend = c("Oriental honey buzzard", "Grey-faced buzzard", "Amur falcon", 
                                  "Eleonora's falcon", "Peregrine falcon", "Osprey"),
       lty = c(4,1,6,2,3,5), cex = 0.55, bty = "n", seg.len = 3)

legend(-45,-67, legend = c("sea-crossing period","High sun elevation", "Low sun elevation", "Night"),
       lty = c(1,1,2,3), lwd = c(9,1,1,1), col = c("#99CC0060", rep("black",3)),cex = 0.55, bty = "n", seg.len = 3, horiz = T)

dev.off()

### --------STEP 5: plot the regions#####

#sample: East Asia

#opean annotated steps
load("2021/ssf_input_ann_cmpl_90_60.RData") #ann_cmpl

#filter out for used steps. only keep wind support

breaks <- c(-15,-10,-5,0,5,10,15,20)
tags <- c("goldenrod2","lemonchiffon1","lightcyan","darkslategray2","darkturquoise","deepskyblue2", "deepskyblue4")

used <- ann_cmpl %>% 
  filter(used == 1) %>% 
  dplyr::select(c("location.long","location.lat","timestamp", "stratum","group", "ind","species", "track","wind_support","delta_t","sun_elev")) %>% 
  mutate(cols = cut(wind_support, breaks = breaks, labels = tags)) %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = wgs) 
  




EF_S <- used %>% 
  filter(group == "EF_S") %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = wgs) 

maps ::map("world", xlim = c(40,50),ylim = c(-17,-10))

points(EF_S[between(EF_S$wind_support, min(used$wind_support),-5),"wind_support"], col = "red")
points(EF_S[between(EF_S$wind_support, -5,0),"wind_support"], col = "orange", add = T)

#convert each step into one line
strata_l <- used %>%
  st_as_sf(coords = c("location.long", "location.lat"), crs = wgs) %>% 
  group_by(track) %>% 
  filter(n() > 1) %>%  #remove tracks with only one point
  arrange(date_time) %>% 
  summarise(track = head(track,1),
            species = head(species,1),
            zone = head(zone, 1), do_union = F) %>% 
  st_cast("LINESTRING")

#bin the wind data
breaks <- c(-10,-5,0,5,10,15,20)
tags <- c("-10 to -5","-5 to 0","0 to 5","5 to 10","10 to 15", "15 to 20")

binned_wind <- cut(used$wind_support,breaks = breaks, include.lowest = T, right = F, labels = tags)

used <- used[order(used$timestamp),]

plot(used[used$group == "EF_S",c("location.long","location.lat")], type = "l")

