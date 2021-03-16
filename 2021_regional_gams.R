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



### --------STEP 5: plot #####
load("2021/predictions_regional_gam_map.RData") #preds_filt
load("2021/models_ls_reg_GAMs.RData") #models_ls
load("2021/timing_for_gam_preds.RData") #timing_areas

region <- st_read("/home/enourani/ownCloud/Work/GIS_files/continent_shapefile/continent.shp") %>% 
  st_crop(xmin = -130, xmax = 157, ymin = -50, ymax = 70) %>%
  st_union()

# create a plotting function for effect plots
names <- c("South-East Asia", "The Americas", "Indian Ocean", "Europe", "Mozambique Channel")

reg_gam_plot <- function(x){
  m <- models_ls[[x]]
  t <- timing_areas[[x]]
  
  plot(0, type = "n", labels = FALSE, tck = 0, xlim =  c(1,366), ylim = c(-2.5,5), xlab = "", ylab = "")#, main = names[[x]]) #expression(paste(Delta,"T")), main = x
  abline(h = 0, col = "gray60",lty = 1, lwd = 0.5)
  rect(xleft = min(t),ybottom = -2.7,xright = max(t),ytop = 5, col="#99CC0060",border=NA) #water-crossing window
  plot_smooth(m, view="yday", plot_all="sun_elev_f", rm.ranef=F, lwd = 1.5, #ylim=c(-2,5),
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

#####
#load sample species tracks
load("tracks_for_global_map.RData")

#imaginary raster for legend. I want the legend to go from -1 to 5 (range of values in the prediction rasters)
imaginary_r<-preds_filt[[1]]
imaginary_r@data@values[3:9]<- rep(5,7)
imaginary_r@data@values[10:16]<- rep(-1,7)


#create a color palette
cuts <- seq(-1,5,0.01) #set breaks
#pal <- colorRampPalette(c("dodgerblue","darkturquoise","goldenrod1","coral","firebrick1"))
pal <- colorRampPalette(c("dodgerblue","darkturquoise", "goldenrod1","coral","firebrick1","firebrick4"))
colpal <- pal(570)
  
#pdf("/home/mahle68/ownCloud/Work/Projects/delta_t/paper_prep/figures/global_plot_scaled.pdf", width = 11, height = 5.2)

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
text(x = -125, y = c(2,32,62), labels = c("0° ", "30° N", "60° N"), cex = 0.6, col = "grey65")

#add a frame for the sub-plots
rect(xleft = -85,
     xright = 153,
     ybottom =  -52,
     ytop = -4,
     col="white",
     border = NA)

#add subplots...
centers_x <- c(124,-56,64,4) #distance between centers = 60

for(i in 1:length(centers_x)){
  rect(xleft = centers_x[i] - 27,
       xright = centers_x[i] + 27,
       ybottom =  -42, #-38
       ytop = -6, #-2
       col="white")
  
  subplot(reg_gam_plot(i), x = centers_x[i] + 4,y = -23, size = c(1.6,0.7),  #-23 was -19
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



# rasterImage(AF, xleft = 85, xright = 110, ybottom = 20, ytop = 45, angle = +25) #make AF smaller than others
# rasterImage(OHB, xleft = 145, xright = 175, ybottom = 25, ytop = 55, angle = +40)
# rasterImage(GFB, xleft = 130, xright = 155, ybottom = 20, ytop = 45, angle = 0)
# rasterImage(EF, xleft = 8, xright = 33, ybottom = 8, ytop = 33, angle = 0)
# rasterImage(PF_E, xleft = 30, xright = 60, ybottom = 45, ytop = 75, angle = 15) 
# rasterImage(PF_A, xleft = -65, xright = -35, ybottom = 40, ytop = 70, angle = 20) 
# rasterImage(OS_E, xleft = -10, xright = 20, ybottom = 30, ytop = 60, angle = 20)
# rasterImage(OS_A, xleft = -110, xright = -80, ybottom = 35, ytop = 65)

#add legend
plot(imaginary_r, legend.only = TRUE, breaks = cuts, col = colpal, #create an imaginary raster that has the high values from juv and low values from adlt
     #legend.width = 0.3, legend.shrink = 0.3,
     smallplot = c(0.059,0.169, 0.145,0.16), #c(min % from left, max % from left, min % from bottom, max % from bottom)
     #smallplot = c(0.06,0.17, 0.28,0.295),#c(0.06,0.17, 0.36,0.375), 
     axis.args = list(at = c(-1,0,2,4), #same arguments as any axis, to determine the length of the bar and tick marks and labels
                    labels = c(-1,0,2,4), 
                    col = NA, #make sure box type in par is set to n, otherwise axes will be drawn on the legend :p
                    col.ticks = NA,
                    line = -1.3),
     horizontal = T,
     legend.args = list(text = expression(italic(paste(Delta,"T", "(°C)"))), side = 3, font = 2, line = 0.1, cex = 0.7)
     )

#text(x = -117,y = -35,  expression(italic(paste(Delta,"T", "(°C)"))), cex = 0.7)
#text(x = -115,y = -8, "Map legend", cex = 0.8)
text(x = -105,y = -8, "Map legend", cex = 0.8)
legend(-126,-10, legend = c("Oriental honey buzzard", "Grey-faced buzzard", "Amur falcon", 
                                  "Eleonora's falcon", "Peregrine falcon", "Osprey"),
       lty = c(4,1,6,2,3,5), cex = 0.55, bty = "n", seg.len = 3)
#text(x = -118.5,y = -39.5, "Sub-plots legend", cex = 0.7)

legend(-30,-43, legend = c("sea-crossing period","High sun elevation", "Low sun elevation", "Night"),
       lty = c(1,1,2,3), lwd = c(9,1,1,1), col = c("#99CC0060", rep("black",3)),cex = 0.55, bty = "n", seg.len = 3, horiz = T)

#legend(-27,-43, legend = c("sea-crossing period","High sun elevation", "Low sun elevation", "Night"),
#       lty = c(1,1,2,3), lwd = c(9,1,1,1), col = c("#99CC0060", rep("black",3)),cex = 0.55, bty = "n", seg.len = 3, ncol = 2)


#legend(-126,-40, legend = c("High sun elevation", "Low sun elevation", "Night"),
#       lty = c(1,2,3), cex = 0.55, bty = "n", seg.len = 3)

#legend(-126,-49, legend = c("sea-crossing timing"),
#       lty = 1, lwd = 9, col = "#99CC0060", cex = 0.55, bty = "n", seg.len = 2.9)
#dev.off()


## add routes
#choose one sample trajecotry for species with good data. unfiltered tracks in track_based_prep_analyze_daily.R

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
load("tracks_for_global_map.RData")


# add all tracks to map
load("all_spp_unfiltered_updated_lc_0_removed_new_track_id.RData") #dataset

all_tracks <- dataset %>% 
  filter(season == "autumn") %>% 
  st_as_sf(coords = c("location.long","location.lat"), crs = wgs) %>% 
  group_by(track) %>% 
  summarize(do_union = FALSE) %>% 
  st_cast("LINESTRING")

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

#test