#script for creating maps color-coded with variance of wind support
#follows up on 2021_wind_delta_t_maps.R

#Jul 12. 2021, Radolfzell, DE
#Elham Nourani

library(tidyverse)
library(sp)
library(sf)
library(move)
library(scales)
library(maptools)
library(mapview)

setwd("/home/enourani/ownCloud/Work/Projects/delta_t")
source("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/wind_support_Kami.R")

wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")


#maps ############

#open file prepped in 2021_wind_delta_t_maps.R

load("R_files/2021/raw_sea_points_for_maps_mv.RData") #mv

data <- mv %>% 
  as.data.frame() %>% 
  dplyr::select(-c("u925","v925","sst","t2m","delta_t","coords.x1","coords.x2","optional", "X", "sensor","timestamps")) %>% 
  mutate(row_id = row_number()) %>% 
  slice(rep(row_number(),40)) %>% 
  group_by(row_id) %>% 
  mutate(year = c(1981:2020)) %>%
  ungroup() %>%
  mutate(timestamp = paste(as.character(date_time),"000",sep = ".")) %>% #add miliseconds per Movebank's requirements
  as.data.frame()


str_sub(data$timestamp,1,4) <- data$year #replace original year with years from 1981-2020
colnames(data)[c(7,8)] <- c("location-long","location-lat") #rename columns to match movebank format


n_chunks <- ceiling(nrow(data)/1e6)
nrow_chunks <- round(nrow(data)/n_chunks)

r <- rep(1: n_chunks,  each = nrow_chunks)

chunks <- split(data, r)

#save as csv
lapply(c(1:length(chunks)), function(i){
  
  write.csv(chunks[[i]], paste0("R_files/2021/public/40_yr_map_input_", i, "_.csv"))
  
})

#open annotated files

#calculate long-term metrics and merge with previously annotated data
ann_40_ls <- list.files("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/2021/annotations/40yr_map/",pattern = ".csv", recursive = T, full.names = T) 

map_input <- lapply(ann_40_ls, read.csv, stringsAsFactors = F) %>% 
  reduce(full_join) %>% 
  rename(u925 = ECMWF.ERA5.PL.U.Wind,
         v925 = ECMWF.ERA5.PL.V.Wind) %>% 
  mutate(wind_support = wind_support(u = u925, v = v925, heading = heading)) %>% 
  group_by(row_id) %>% 
  summarise(wind_support_var = var(wind_support, na.rm = T),
            timestamp = head(timestamp,1),
            species = head(species,1),
            track = head(track,1),
            ind = head(ind,1),
            location.long = head(location.long,1),
            location.lat = head(location.lat,1)) %>% 
  ungroup() %>%
  as.data.frame()


save(map_input, file = "/home/enourani/ownCloud/Work/Projects/delta_t/R_files/2021/public/Dryad/input_var_wspt_map.RData")

#plot

load("R_files/2021/raw_sea_points_for_maps_mv.RData") #mv

#add a categorical variable for wind levels
breaks_wv <- c(0,10,20,40,88)
tags_wv <- c("0 to 10","10 to 20","20 to 40", "> 40")

#create color palettes 
Pal <- colorRampPalette(c("mediumblue", "cornflowerblue", "darkgoldenrod2", "indianred1")) #colors for negative values
Cols_wv <- paste0(Pal(4), "80") #add transparency. 50% is "80". 70% is "B3". 80% is "CC". 90% is "E6"



df <- map_input %>% 
  mutate(binned_wv = cut(wind_support_var,breaks = breaks_w, include.lowest = T, right = F, labels = tags_wv)) %>% 
  mutate(cols_wv = as.factor(binned_wv)) 

levels(df$cols_wv) <- Cols_wv

df_sp <- SpatialPointsDataFrame(coords = df[,c("location.long", "location.lat")], proj4string = wgs, data = df)

#plot
X11(width = 12, height = 15) #make the window proportional to region

pdf("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/2021/raw_wind_dt_updated.pdf", width = 12, height = 11.5)

par(mfrow=c(3,1),
    #fig = c(0,1,0,1), #do this if you want to add the small plots as subplots
    bty="n", #no box around the plot
    cex.axis= 0.6, #x and y labels have 0.75% of the default size
    #font = 3, #3: axis labels are in italics
    font.axis = 3,
    cex.lab = 0.6,
    #cex = 0.5,
    oma = c(0,0,1.5,0),
    mar = c(0, 0, 0.3, 0),
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

text(x = -85,y = 0, "Wind support (m/s)", cex = 0.8, font = 3)
legend(x = -100, y = 0, legend = levels(df_sp$binned_w), col = Cols_w, pch = 20, 
       bty = "n", cex = 0.8, text.font = 3)
mtext("Sea-crossing tracks annotated with wind support", 3, outer = F, cex = 1.3, line = -0.5)

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
       bty = "n", cex = 0.8, text.font = 3)
mtext(bquote(paste('Sea-crossing tracks annotated with', italic(~ Delta *"T"))), 3, outer = F, cex = 1.3, line = -0.5)

dev.off()


#initiation plots ############

#color palette
Pal <- colorRampPalette(c("darkgoldenrod1","lightpink1", "mediumblue")) #colors for negative values
Cols <- paste0(Pal(3), "80") #add transparency. 50% is "80". 70% is "B3". 80% is "CC". 90% is "E6"



load("R_files/2021/raw_sea_points_for_maps_mv.RData") #mv from higher up.... this is all the data, before filtering


starts_all <- mv %>% 
  as.data.frame() %>% 
  drop_na(c("heading","delta_t")) %>% 
  mutate(wind_support= wind_support(u = u925,v = v925,heading = heading),
         cross_wind= cross_wind(u = u925,v = v925,heading = heading)) %>% 
  mutate(group = ifelse(species == "OHB", "OHB",
                        ifelse(species == "GFB", "GFB",
                        ifelse(species == "O" & location.long < -30, "O_A",
                               ifelse(species == "O" & location.long > -30, "O_E",
                                      ifelse(species == "EF" & location.lat > 0, "EF_med",
                                             ifelse(species == "EF" & location.lat < 0, "EF_moz",
                                                    ifelse(species == "PF" & location.long < -30, "PF_A",
                                                           "PF_E"))))))))  %>% 
  group_by(track) %>% 
  arrange(timestamp) %>% 
  slice(1) %>%
  ungroup() %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = wgs) %>% 
  mutate(s_elev_angle = solarpos(st_coordinates(.), timestamp, proj4string = CRS("+proj=longlat +datum=WGS84"))[,2]) %>% #calculate solar elevation angle
  mutate(sun_elev = ifelse(s_elev_angle < -6, "night", #create a categorical variable for teh position of the sun
                           ifelse(s_elev_angle > 40, "high", "low"))) %>% 
  as("Spatial") %>% 
  as.data.frame() %>% 
  mutate(group = as.factor(group)) %>% 
  as.data.frame()


labels <- c("EF \n (Mediterranean)", "EF \n (Mozambique)", "GFB \n", "Osprey \n (America)", "Osprey \n (Europe)", "OHB \n ", "PF \n (America)", "PF \n (Europe)")

variables <- c("wind_support", "delta_t")

pdf("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/2021/initiation_boxplots_updated.pdf", width = 12, height = 5.5)

X11(width = 10.2, height = 6)

par(mfrow= c(2,1), 
    oma = c(1.7,0,2.5,0), 
    mar = c(0.5,4,0.2,0.2),
    las = 1,
    bty = "l",
    cex.axis = 0.8,
    font.axis = 3,
    tck = -0.015,
    mgp=c(1,1,0))



for(i in 1:length(variables)){
  
  boxplot(starts_all[,variables[i]] ~ starts_all[,"group"], data = starts_all, boxfill = NA, border = NA, xlab = "", ylab = "", xaxt = "n")
  abline(h = 0, lty = 2, col = alpha("black", 0.6), lwd = 0.3)
  #if(i == 1){
  #  legend("topleft", legend = c("high","low", "night"), fill = c("orange","gray"), bty = "n")
  #}
  points(jitter((as.numeric(starts_all[starts_all$sun_elev == "high","group"]) -0.3),0.3), starts_all[starts_all$sun_elev == "high", variables[i]], 
         yaxt = "n", xaxt = "n", pch = 20, cex = 0.8, col = alpha("black", 0.6))
  
  boxplot(starts_all[starts_all$sun_elev == "high", variables[i]] ~ starts_all[starts_all$sun_elev == "high","group"], 
          yaxt = "n", xaxt = "n", add = T, boxfill = Cols[1], outline=FALSE, lwd = 0.5, 
          boxwex = 0.25, at = 1:length(unique(starts_all$group)) - 0.3)
  
  points(jitter(as.numeric(starts_all[starts_all$sun_elev == "low","group"]),0.3), starts_all[starts_all$sun_elev == "low", variables[i]], 
         yaxt = "n", xaxt = "n", pch = 20, cex = 0.8, col = alpha("black", 0.6))
  
  boxplot(starts_all[starts_all$sun_elev == "low", variables[i]] ~ starts_all[starts_all$sun_elev == "low", "group"], 
          yaxt = "n", xaxt = "n", add = T, boxfill = Cols[2], outline = FALSE,  lwd = 0.5, 
          boxwex = 0.25, at = 1:length(unique(starts_all$group))+ 0)
  
  points(jitter((as.numeric(starts_all[starts_all$sun_elev == "night","group"]) +0.3),0.3), starts_all[starts_all$sun_elev == "night", variables[i]], 
         yaxt = "n", xaxt = "n", pch = 20, cex = 0.8, col = alpha("black", 0.6))
  
  boxplot(starts_all[starts_all$sun_elev == "night", variables[i]] ~ starts_all[starts_all$sun_elev == "night", "group"], 
          yaxt = "n", xaxt = "n", add = T, boxfill = Cols[3], outline=FALSE,  lwd = 0.5, 
          boxwex = 0.25, at = 1:length(unique(starts_all$group))+ 0.3)
  

  if(i == length(variables)){
    axis(side = 1, at = 1:length(levels(starts_all$group)), labels = labels, 
         tick = T , col.ticks = 1, col = NA, tck = -.015, lwd = 0, lwd.ticks = 1, cex = 0.9)
  }
  
  if(i == 1){
    legend(x = 7.3, y = 19, legend = c("daytime: high sun", "daytime: low sun", "night"), col = Cols, #coords indicate top-left
           cex = 0.78, pt.cex = 1.1, bg = "white", bty = "n", pch = 15)
  }
  
  if(variables[i] == "wind_support"){
    mtext("Wind support (m/s)", side = 2, las = 3, line = 2, font = 3)
  }
  
  if(variables[i] == "delta_t"){
    mtext(expression(italic(paste(Delta,"T", "(°C)"))), side = 2, las = 3, line = 2)
  }
  
}
mtext("Atmospheric conditions at the start of sea-crossing tracks", side = 3, outer = T, cex = 1.2, font = 1, line = 0.7)


dev.off()
