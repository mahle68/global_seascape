##script for making a map for Nourani et al 2021 (Proc B) press release

setwd("/home/enourani/ownCloud/Work/Projects/delta_t/R_files")

library(plotrix)

#from 2021_regional_gams.R
AF <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/AF_1.png")
EF <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/Updated/EF-2.png")
GFB <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/Updated/GFB-2.png")
OHB <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/OHB_2.png")
OS_A <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/Updated/EO-A.png")
OS_E <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/Updated/EO-B.png") #mirrored
PF_A <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/Updated/EP-A.png")
PF_E <- readPNG("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/from_james/Updated/EP-B.png")

#from all_figures.R
#load the data  
load("predictions_regional_gams.RData") #preds_filt; rasters of predictions using GAMs produced in seascape_GAM_analysis.R
load("models_ls_reg_GAMs.RData") #models_ls; list of regional gams produced in seascape_GAM_analysis.R
load("timing_for_gam_preds.RData") #timing_areas; timing of sea-crossing
load("sample_tracks.RData") #sp_samples; sample tracks for each species to add to the map

#load continents shapefile as the map background
continents <- shapefile("continent_shapefile/continent.shp")
region <- crop(continents, extent(c(-130,158,-60,71.3)))

#names of regions
#names <- c("South-East Asia", "The Americas", "Indian Ocean", "Europe", "Mozambique Channel")


#imaginary raster for legend. The legend should go from -1 to 5 (range of values in the prediction rasters)
imaginary_r <- preds_filt[[1]]
imaginary_r@data@values[3:9] <- rep(5,7)
imaginary_r@data@values[10:16] <- rep(-1,7)

#create a color palette
cuts <- seq(-1,5,0.01) #set breaks
pal <- colorRampPalette(c("dodgerblue","darkturquoise", "goldenrod1","coral","firebrick1","firebrick4"))
colpal <- pal(570)

X11(width = 10, height = 5) #in inches

png("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/procB/outreach/seascape_map2.png", width = 10, height = 5, units = "in", res = 400)

par(mfrow=c(1,1),
    fig = c(0,1,0,1), #do this if you want to add the small plots as subplots
    bty="n", #no box around the plot
    cex.axis= 0.6, #x and y labels have 0.75% of the default size
    font = 3, #3: axis labels are in italics
    font.axis = 3,
    cex.lab = 0.6,
    oma = c(0,0,0,0),
    mar = c(0, 0, 0, 0),
    lend = 1  #rectangular line endings (trick for adding the rectangle to the legend)
)

#plot the background map
plot(region, col="#dcdcdc",border="#dcdcdc") # used to be ##e5e5e5
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
plot(st_geometry(sp_samples$EF_S_sample), add = T, lty = 2, lwd = lwd, col = col)

#add latitude lines
clip(-130, 157, -50, 73)
abline(h = 0, col = "grey70",lty = 2)
abline(h = 30, col = "grey70",lty = 2)
abline(h = 60, col = "grey70",lty = 2)
text(x = -125, y = c(2,32,62), labels = c("Equator","30° N", "60° N"), cex = 0.6, col = "grey65")

usr <- par("usr")
clip(usr[1], usr[2], usr[3], usr[4])

#add birds... try to scale them
rasterImage(AF, xleft = 95, xright = 111, ybottom = 29, ytop = 45, angle = +25) #make AF smaller than others
rasterImage(OHB, xleft = 145, xright = 180, ybottom = 20, ytop = 55, angle = +40)
rasterImage(GFB, xleft = 130, xright = 155, ybottom = 20, ytop = 45, angle = 0)
rasterImage(EF, xleft = 11, xright = 31, ybottom = 10, ytop = 30, angle = 0)
rasterImage(PF_E, xleft = 28, xright = 48, ybottom = 50, ytop = 70, angle = 15) 
rasterImage(PF_A, xleft = -65, xright = -45, ybottom = 40, ytop = 60, angle = 20) 
rasterImage(OS_E, xleft = -10, xright = 20, ybottom = 30, ytop = 60, angle = 20)
rasterImage(OS_A, xleft = -110, xright = -80, ybottom = 35, ytop = 65)

# #add a frame for the sub-plots and legend
# rect(xleft = -132.67526,
#      xright = -87.56709,
#      ybottom =  -43.191125,
#      ytop =  9.471234,
#      col="white",
#      border = NA)

# #add legend
# plot(imaginary_r, legend.only = TRUE, breaks = cuts, col = colpal, 
#      smallplot = c(0.04,0.15,0.20,.215), #c(min % from left, max % from left, min % from bottom, max % from bottom)
#      axis.args = list(at = c(-1,0,2,4), #same arguments as any axis, to determine the length of the bar and tick marks and labels
#                       labels = c(-1,0,2,4), 
#                       col = NA, #make sure box type in par is set to n, otherwise axes will be drawn on the legend :p
#                       col.ticks = NA,
#                       line = -1.3),
#      horizontal = T,
#      legend.args = list(text = "Uplift (°C)", side = 3, font = 3, line = 0.1, cex = 0.7)
# )
# 
# 
# #text(x = -118,y = 10, "Map legend", cex = 0.8)
# legend(-130,8, legend = c("Oriental honey buzzard", "Grey-faced buzzard", "Amur falcon", 
#                           "Eleonora's falcon", "Peregrine falcon", "Osprey"),
#        lty = c(4,1,6,2,3,5), cex = 0.55, bty = "n", seg.len = 3)
# 
# textbox(x = c(-130, -80), y = -26.65, c("Uplift was estimated as the temperature difference between the sea-surface and the air",
#                                         "(Nourani et al. 2021; Proceedings of the Royal Society B)"),
#         justify = "l", box = F, cex = 0.55, font = 3)

#horizontal legend
#add a frame for the sub-plots and legend
rect(xleft = -128,
     xright = 156,
     ybottom =  -30,
     ytop = -60,
     col="white",
     border = NA)

legend(-107,-31, legend = c("Oriental honey buzzard", "Grey-faced buzzard", "Amur falcon", 
                          "Eleonora's falcon", "Peregrine falcon", "Osprey"),
       lty = c(4,1,6,2,3,5), cex = 0.55, bty = "n", seg.len = 3, horiz = T)  

plot(imaginary_r, legend.only = TRUE, breaks = cuts, col = colpal, 
     #smallplot = c(0.14, 0.14 + (0.11), 0.13, 0.13 + (.015)), #c(min % from left, max % from left, min % from bottom, max % from bottom)
     smallplot = c(0.12, 0.12 + (0.11), 0.16, 0.16 + (.015)), #c(min % from left, max % from left, min % from bottom, max % from bottom)
     axis.args = list(at = c(-1,0,2,4), #same arguments as any axis, to determine the length of the bar and tick marks and labels
                      labels = c(-1,0,2,4), 
                      col = NA, #make sure box type in par is set to n, otherwise axes will be drawn on the legend :p
                      col.ticks = NA,
                      line = -1.3),
     horizontal = T,
     legend.args = list(text = "Uplift (°C)", side = 3, font = 3, line = 0.1, cex = 0.6)
)

textbox(x = c(-65.5, 133), y = -43, c("Uplift was estimated as the temperature difference between the sea-surface and the air. Positive values indicate uplift and negative values correspond to subsidence (Nourani et al. 2021; Proceedings of the Royal Society B)"),
        justify = "l", box = F, cex = 0.6, font = 3)

dev.off()
