# Scripts for producing the figures (main text and supplementary)
# This is script 6 of 6 for reproducing the results of Nourani et al 2021, ProcB.
# Elham Nourani, PhD. Jun.10. 2021
#-----------------------------------------------------------------



# ---------- Fig 1: w_star #####

#new data for predictions
newx <- seq(min(data$delta_t), 9, by = 0.05)
conf_interval <- predict(fit_2, newdata = data.frame(delta_t = newx), interval = "confidence",
                         level = 0.95)

data$shape <- as.factor(data$species)
levels(data$shape) <- c(0,1,2,3,4,5)

X11(width = 5.2, height = 3)
pdf("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/2021/w_star_spp_gam.pdf", width = 5.2, height = 3)

par(mfrow=c(1,1), 
    bty = "l",
    cex.axis = 0.7,
    font.lab = 3,
    cex.lab = 0.9,
    mgp=c(2,0.5,0), #margin line for the axis titles
    mar = c(3.5,3.5,0.5,6.35),
    #oma = c(0,0,0,4),
    xpd=  TRUE) # allows drawing legend in outer margins

plot(0, type = "n", labels = FALSE, tck = 0, xlim = c(0,8.3), ylim = c(0.5,5), xlab = expression(italic(paste(Delta,"T", "(°C)"))), ylab = "w* (m/s)")

with(pos_dt,points(delta_t, w_star, col= as.character(data$color), pch = as.numeric(data$shape), cex = 0.5))

polygon(x.polygon , y.polygon , col = alpha("grey", 0.5) , border = NA) #confidence intervals

lines(data$delta_t[i.for] , fit[i.for], col = "black" , lwd = 1.1)

axis(side = 1, at = c(0,2,4,6,8), labels = c(0,2,4,6,8), 
     tick = T , col.ticks = 1, col = NA, tck = -.015,lwd = 0, lwd.ticks = 1)
axis(side= 2, at= seq(1,5,1), labels= seq(1,5,1),
     tick=T , col.ticks = 1, col = NA, tck=-.015,
     las=2) # text perpendicular to axis label 

text(9.8,4.6, "Time of day", cex = 0.7)
legend(x = 8.7, y = 4.5, legend = c("daytime: high sun", "daytime: low sun", "night"), col = Cols, #coords indicate top-left
       cex = 0.7, pt.cex = 0.9, bg = "white", bty = "n", pch = 20)


text(9.5,3, "Species", cex = 0.7)
legend(x = 8.7, y = 2.9, legend = c("Pernis ptilorhynchus", "Pandion haliaetus","Butastur indicus", "Falco peregrinus", "F. eleonorae"), #species : unique(data$species), #coords indicate top-left
       cex = 0.7, pt.cex = 0.5, bg = "white", bty = "n", pch = as.numeric(unique(data$shape)),
       text.font = 3)


dev.off()

# ---------- Fig 2: global seascape map #####
#Bird figures are not included in this code for copyright reasons.

load("2021/predictions_regional_gam_map.RData") #preds_filt
load("2021/models_ls_reg_GAMs.RData") #models_ls
load("2021/timing_for_gam_preds.RData") #timing_areas

#elements prepared in regional_gam.R
load("tracks_for_global_map.RData") #sp_samples
EF_S <- read.csv("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/2021/EF_B2378_2017_autumn.csv") %>% 
  mutate(dt,dt = as.POSIXct(strptime(dt,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>%
  filter(dt <= "2017-11-15 14:00:00") %>% #remove wintering points
  st_as_sf(coords = c("long","lat"), crs = wgs) %>% 
  arrange(dt) %>% 
  summarise(do_union = F) %>%
  st_cast("LINESTRING")


# AF <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/AF_1.png")
# EF <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/Updated/EF-2.png")
# GFB <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/Updated/GFB-2.png")
# OHB <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/OHB_2.png")
# OS_A <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/Updated/EO-A.png")
# OS_E <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/Updated/EO-B.png") #mirrored
# PF_A <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/Updated/EP-A.png")
# PF_E <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/Updated/EP-B.png")



region <- st_read("/home/enourani/ownCloud/Work/GIS_files/continent_shapefile/continent.shp") %>% 
  st_crop(xmin = -130, xmax = 158, ymin = -74, ymax = 71) %>%
  st_union()

# create a plotting function for effect plots
names <- c("South-East Asia", "The Americas", "Indian Ocean", "Europe", "Mozambique Channel")


###

#imaginary raster for legend. I want the legend to go from -1 to 5 (range of values in the prediction rasters)
imaginary_r<-preds_filt[[1]]
imaginary_r@data@values[3:9]<- rep(5,7)
imaginary_r@data@values[10:16]<- rep(-1,7)


#create a color palette
cuts <- seq(-1,5,0.01) #set breaks
pal <- colorRampPalette(c("dodgerblue","darkturquoise", "goldenrod1","coral","firebrick1","firebrick4"))
colpal <- pal(570)

#pdf("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/2021/global_plot.pdf", width = 11, height = 6) #for paper
#pdf("/home/enourani/ownCloud/Work/conferences/MPI_YALE_2021/delta_t_global_plot.pdf", width = 22, height = 12) #larger image for presentation

X11(width = 11, height = 6) #make the window proportional to region
X11(width = 22, height = 12) #make the window proportional to region

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

plot(region, col="#e5e5e5",border="#e5e5e5")
plot(preds_filt[[1]],axes = F, box=F, legend=FALSE,zlim=c(-1,5),breaks=cuts, col = colpal, add = T) 
plot(preds_filt[[2]],axes = F, box=F, legend=FALSE,zlim=c(-1,5),breaks=cuts, col = colpal, add = T) 
plot(preds_filt[[3]],axes = F, box=F, legend=FALSE,zlim=c(-1,5),breaks=cuts, col = colpal, add = T) 
plot(preds_filt[[4]],axes = F, box=F, legend=FALSE,zlim=c(-1,5),breaks=cuts, col = colpal, add = T) 
plot(preds_filt[[5]],axes = F, box=F, legend=FALSE,zlim=c(-1,5),breaks=cuts, col = colpal, add = T) 

#add tracks
lwd <- 1.2
col <- "grey25"
plot(st_geometry(sp_samples$EF_sample), add= T, lty = 2, lwd = lwd, col = col)
plot(st_geometry(sp_samples$PF_sample), add= T, lty = 3, lwd = lwd, col = col)
plot(st_geometry(sp_samples$O_sample), add= T, lty = 5, lwd = lwd, col = col)
plot(st_geometry(sp_samples$OHB_sample), add= T, lty = 4, lwd = lwd, col = col)
plot(st_geometry(sp_samples$GFB_sample), add= T, lty = 1, lwd = lwd, col = col)
plot(st_geometry(sp_samples$AF_sample), add= T, lty = 6, lwd = lwd, col = col)
plot(st_geometry(sp_samples$EF_sample), add = T, lty = 2, lwd = lwd, col = col)

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

#for the paper:
#text(x = -118,y = 10, "Map legend", cex = 0.8)
#legend(-130,8, legend = c("Oriental honey buzzard", "Grey-faced buzzard", "Amur falcon", 
#                                  "Eleonora's falcon", "Peregrine falcon", "Osprey"),
#       lty = c(4,1,6,2,3,5), cex = 0.55, bty = "n", seg.len = 3)

#for the slide
text(x = -118,y = 10, "Map legend", cex = 1.3)
legend(-130,8, legend = c("Oriental honey buzzard", "Grey-faced buzzard", "Amur falcon", 
                          "Eleonora's falcon", "Peregrine falcon", "Osprey"),
       lty = c(4,1,6,2,3,5), cex = 1.2, bty = "n", seg.len = 3)


legend(-45,-67, legend = c("sea-crossing period","High sun elevation", "Low sun elevation", "Night"),
       lty = c(1,1,2,3), lwd = c(9,1,1,1), col = c("#99CC0060", rep("black",3)),cex = 0.55, bty = "n", seg.len = 3, horiz = T)

dev.off()
# ---------- Fig S1: species-specific coefficients #####

#SUPPLEMENTARY FIGURE 1: species-specific coefficients 
#for the best model (M1); original code by Virgilio Gomez-Rubio (Bayesian inference with INLA, 2020)
#species
species_names <- unique(all_data$species)

tab_dt <- data.frame(ID = as.factor(M1$summary.random$species1$ID),
                     mean = M1$summary.random$species1$mean,
                     IClower = M1$summary.random$species1[, 4],
                     ICupper = M1$summary.random$species1[, 6])


tab_wspt <- data.frame(ID = as.factor(M1$summary.random$species2$ID),
                       mean = M1$summary.random$species2$mean,
                       IClower = M1$summary.random$species2[, 4],
                       ICupper = M1$summary.random$species2[, 6])


X11(width = 4, height = 3)

pdf("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/2021/species_var_updated.pdf", width = 4, height = 3)

par(mfrow = c(1,1), bty="n", #no box around the plot
    #cex.axis= 0.75, #x and y labels have 0.75% of the default size
    #font.axis= 0.75, #3: axis labels are in italics
    #cex.lab = 0.75,
    cex = 0.7,
    oma = c(0,3.5,0,0),
    mar = c(3, 2, 0.5, 1),
    bty = "l"
)


plot(0, type = "n", labels = FALSE, tck = 0, xlim = c(-3,6), ylim = c(0,4.5), xlab = "", ylab = "")
#add vertical line for zero
abline(v = 0, col = "grey30",lty = 2)

points(tab_dt$mean, as.numeric(tab_dt$ID) - 0.2, col = "darkgoldenrod2", pch = 19, cex = 1.3)
arrows(tab_dt$IClower, as.numeric(tab_dt$ID) - 0.2,
       tab_dt$ICupper, as.numeric(tab_dt$ID) - 0.2,
       col = "darkgoldenrod2", code = 3, length = 0.03, angle = 90) #angle of 90 to make the arrow head as straight as a line

points(tab_wspt$mean, as.numeric(tab_wspt$ID) + 0.2, col = "cornflowerblue", pch = 19, cex = 1.3)
arrows(tab_wspt$IClower, as.numeric(tab_wspt$ID) + 0.2,
       tab_wspt$ICupper, as.numeric(tab_wspt$ID) + 0.2,
       col = "cornflowerblue", code = 3, length = 0.03, angle = 90) #angle of 90 to make the arrow head as straight as a line

axis(side= 1, at= c(-4,-2,0,2,4), labels= c("-4","-2", "0", "2","4"), 
     tick=T ,col = NA, col.ticks = 1, tck=-.015)

axis(side= 2, at= c(1:4), #line = 6, 
     labels = tab_dt$ID,# c( "Eleonora's falcon","Osprey", "Oriental\n honey buzzard", "Peregrine falcon"), #same order as tab_dt$ID
     tick = T ,col = NA, col.ticks = 1, # NULL would mean to use the defult color specified by "fg" in par
     tck = -.015 , #tick marks smaller than default by this proportion
     las = 2) # text perpendicular to axis label 

#add legend
legend(x = 3.6 , y = 0.7, legend = c("Wind support", expression(italic(paste(Delta,"T")))), 
       col = c("cornflowerblue","darkgoldenrod1"), #coords indicate top-left
       pch = 19, bg="white",bty="n", cex = 0.9)

dev.off()


# ---------- Fig S2: boxplot comparing used and available steps #####

load("2021/public/ssf_input_ann_cmpl_60_15.RData") #ann_cmpl; This dataframe includes used and alternative steps and can be reproduced using step_generation.R

X11(width = 9, height = 7)

pdf("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/2021/boxplots_60_15.pdf", width = 9, height = 7)

par(mfrow= c(2,3), 
    oma = c(0,0,3,0), 
    las = 1)

labels <- c(expression(italic(paste(Delta,"T"))), "Wind support", "Wind speed")
variables <- c("delta_t", "wind_support", "wind_speed")
v_variables <- c("delta_t_var", "wind_support_var","wind_speed_var")

for(i in 1:length(variables)){
  
  boxplot(ann_cmpl[,variables[i]] ~ ann_cmpl[,"species"], data = ann_cmpl, boxfill = NA, border = NA, main = labels[i], xlab = "", ylab = "")
  if(i == 1){
    legend("topleft", legend = c("used","available"), fill = c("orange","gray"), bty = "n")
  }
  boxplot(ann_cmpl[ann_cmpl$used == 1, variables[i]] ~ ann_cmpl[ann_cmpl$used == 1,"species"], 
          yaxt = "n", xaxt = "n", add = T, boxfill = "orange",
          boxwex = 0.25, at = 1:length(unique(ann_cmpl[ann_cmpl$used == 1, "species"])) - 0.15)
  boxplot(ann_cmpl[ann_cmpl$used == 0, variables[i]] ~ ann_cmpl[ann_cmpl$used == 0, "species"], 
          yaxt = "n", xaxt = "n", add = T, boxfill = "grey",
          boxwex = 0.25, at = 1:length(unique(ann_cmpl[ann_cmpl$used == 1 , "species"])) + 0.15)
  
}
mtext("Instantaneous values at each step", side = 3, outer = T, cex = 1.3)

for(i in 1:length(v_variables)){
  
  boxplot(ann_cmpl[,v_variables[i]] ~ ann_cmpl[,"species"], data = ann_cmpl, boxfill = NA, border = NA, main = labels[i], xlab = "", ylab = "")
  if(i == 1){
    legend("topleft", legend = c("used","available"), fill = c("orange","gray"), bty = "n")
  }
  boxplot(ann_cmpl[ann_cmpl$used == 1,v_variables[i]] ~ ann_cmpl[ann_cmpl$used == 1,"species"], 
          yaxt = "n",xaxt = "n", add = T, boxfill = "orange",
          boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) - 0.15)
  boxplot(ann_cmpl[ann_cmpl$used == 0,v_variables[i]] ~ ann_cmpl[ann_cmpl$used == 0,"species"], 
          yaxt = "n",xaxt = "n", add = T, boxfill = "grey",
          boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) + 0.15)
} 

mtext("40-year variances at each step", side = 3, outer = T, cex = 1.3, line = -25)

dev.off()

# ---------- Fig S3: maps with annotated tracking points #####

# ---------- Fig S4: boxplots for conditions at sea-crossing initiation #####