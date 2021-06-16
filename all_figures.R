# Scripts for producing the figures (main text and supplementary)
# This is script 6 of 6 for reproducing the results of Nourani et al 2021, ProcB.
# Elham Nourani, PhD. Jun.10. 2021
#-----------------------------------------------------------------



# ---------- Fig 1: w_star #####

# ---------- Fig 2: global seascape map #####

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