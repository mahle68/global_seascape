#script for producing plots for delta_t inla ssf 
#Elham Nourani. Radolfzell, Germany. 
#June 2. 2020


#plot the coefficients for fixed effects ####

#from source code of Efxplot. to extract the info that I need from the inla model
ModelList <- list(m1c,m1e,m1d)
graphlist<-list()
for(i in 1:length(ModelList)){
  model<-ModelList[[i]]
  
    graph<-as.data.frame(summary(model)$fixed)
    colnames(graph)[which(colnames(graph)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
    colnames(graph)[which(colnames(graph)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
    colnames(graph)[which(colnames(graph)%in%c("mean"))]<-c("Estimate")
    
    graph$Model<-i
    graph$Factor<-rownames(graph)
    
    graphlist[[i]]<-graph
}

graph <- bind_rows(graphlist)

graph$Sig <- with(graph, ifelse(Lower*Upper>0, "*", ""))

graph$Model <- as.factor(graph$Model)

position <- ifelse(length(unique(graph$Model))  ==  1, "none", "right")

VarOrder <- rev(unique(graph$Factor))
VarNames <- VarOrder

graph$Factor <- factor(graph$Factor, levels = VarOrder)
levels(graph$Factor) <- VarNames

min<-min(graph$Lower,na.rm = T)
max<-max(graph$Upper,na.rm = T)

graph$Factor_n <- as.numeric(graph$Factor)

X11(width = 3.5, height = 3)
pdf("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/fixed_effects.pdf",width = 3.5, height = 3)

par(mfrow=c(1,1), bty="n", #no box around the plot
    #cex.axis= 0.75, #x and y labels have 0.75% of the default size
    #font.axis= 0.75, #3: axis labels are in italics
    #cex.lab = 0.75,
    cex = 0.7,
    oma = c(0,3.5,0,0),
    mar = c(3, 3.5, 0.5, 1),
    bty = "l"
)

#plot(0, type = "n", bty = "l",labels = FALSE, tck = 0, ann = F, xlim = c(-6,8), ylim = c(0,6.3))
plot(0, type = "n", labels = FALSE, tck = 0, xlim = c(-6,8), ylim = c(0,6.3), xlab = "Estimate", ylab = "")

#add vertical line for zero
abline(v = 0, col = "grey30",lty = 2)
#add points and error bars
points(graph[graph$Model == 1, "Estimate"], graph[graph$Model == 1,"Factor_n"] - 0.2, col = "salmon2", pch = 19, cex = 1.3)
arrows(graph[graph$Model == 1, "Lower"], graph[graph$Model == 1,"Factor_n"] - 0.2,
       graph[graph$Model == 1, "Upper"], graph[graph$Model == 1,"Factor_n"] - 0.2,
       col = "salmon2", code = 3, length = 0.03, angle = 90) #angle of 90 to make the arrow head as straight as a line

points(graph[graph$Model == 2, c("Estimate","Factor")], col = "palegreen3", pch = 19, cex = 1.3)
arrows(graph[graph$Model == 2, "Lower"], graph[graph$Model == 2,"Factor_n"],
       graph[graph$Model == 2, "Upper"], graph[graph$Model == 2,"Factor_n"],
       col = "palegreen3", code = 3, length = 0.03, angle = 90)

points(graph[graph$Model == 3, "Estimate"], graph[graph$Model == 3,"Factor_n"] + 0.2, col = "steelblue1", pch = 19, cex = 1.3)
arrows(graph[graph$Model == 3, "Lower"], graph[graph$Model == 3,"Factor_n"] + 0.2,
       graph[graph$Model == 3, "Upper"], graph[graph$Model == 3,"Factor_n"] + 0.2,
       col = "steelblue1", code = 3, length = 0.03, angle = 90)
#add axes
axis(side= 1, at= c(-5,0,5), labels= c("-5", "0", "5"), 
     tick=T ,col = NA, col.ticks = 1, tck=-.015)

axis(side= 2, at= c(1:6), #line=-4.8, 
     labels= c( expression(paste(Delta,"t"," : wind speed")),
                "wind support var","wind support",expression(paste(Delta,"t"," var")),"wind speed", expression(paste(Delta,"t"))),
     tick=T ,col = NA, col.ticks = 1, # NULL would mean to use the defult color specified by "fg" in par
     tck=-.015 , #tick marks smaller than default by this proportion
     las=2 ) # text perpendicular to axis label 

#add legend
legend(x = 5.3, y = 0.8, legend=c("model 3", "model 2", "model 1"), col = c("steelblue1","palegreen3","salmon2"), #coords indicate top-left
        pch = 19, bg="white",bty="n", cex = 0.75)

dev.off()

# plot species differences (for the best model) ####
#fig 4.8 Virgilio's book
species_names <- unique(all_data$species)

tab_dt <- data.frame(ID = m1d$summary.random$species1$ID,
                  mean = m1d$summary.random$species1$mean,
                  IClower = m1d$summary.random$species1[, 4],
                  ICupper = m1d$summary.random$species1[, 6])

tab_wspd <- data.frame(ID = m1d$summary.random$species2$ID,
                  mean = m1d$summary.random$species2$mean,
                  IClower = m1d$summary.random$species2[, 4],
                  ICupper = m1d$summary.random$species2[, 6])

tab_wspt <- data.frame(ID = m1d$summary.random$species3$ID,
                  mean = m1d$summary.random$species3$mean,
                  IClower = m1d$summary.random$species3[, 4],
                  ICupper = m1d$summary.random$species3[, 6])

X11(width = 3.5, height = 3)
pdf("/home/enourani/ownCloud/Work/Projects/delta_t/paper_prep/figures/species_differences.pdf",width = 3.5, height = 3)

par(mfrow = c(1,1), bty="n", #no box around the plot
    #cex.axis= 0.75, #x and y labels have 0.75% of the default size
    #font.axis= 0.75, #3: axis labels are in italics
    #cex.lab = 0.75,
    cex = 0.7,
    oma = c(0,3.5,0,0),
    mar = c(3, 3.5, 0.5, 1),
    bty = "l"
)


plot(0, type = "n", labels = FALSE, tck = 0, xlim = c(-4,4), ylim = c(0,4.3), xlab = "", ylab = "")
#add vertical line for zero
abline(v = 0, col = "grey30",lty = 2)

points(tab_dt$mean, as.numeric(tab_dt$ID) - 0.2, col = "lightcoral", pch = 19, cex = 1.3)
arrows(tab_dt$IClower, as.numeric(tab_dt$ID) - 0.2,
       tab_dt$ICupper, as.numeric(tab_dt$ID) - 0.2,
       col = "lightcoral", code = 3, length = 0.03, angle = 90) #angle of 90 to make the arrow head as straight as a line

points(tab_wspd$mean, as.numeric(tab_wspd$ID), col = "yellowgreen", pch = 19, cex = 1.3)
arrows(tab_wspd$IClower, as.numeric(tab_wspd$ID),
       tab_wspd$ICupper, as.numeric(tab_wspd$ID),
       col = "yellowgreen", code = 3, length = 0.03, angle = 90) #angle of 90 to make the arrow head as straight as a line

points(tab_wspt$mean, as.numeric(tab_wspt$ID) + 0.2, col = "paleturquoise2", pch = 19, cex = 1.3)
arrows(tab_wspt$IClower, as.numeric(tab_wspt$ID) + 0.2,
       tab_wspt$ICupper, as.numeric(tab_wspt$ID) + 0.2,
       col = "paleturquoise2", code = 3, length = 0.03, angle = 90) #angle of 90 to make the arrow head as straight as a line

axis(side= 1, at= c(-2,0,2), labels= c("-2", "0", "2"), 
     tick=T ,col = NA, col.ticks = 1, tck=-.015)

axis(side= 2, at= c(1:4), #line=-4.8, 
     #labels= c("Falco eleonorae", "Pandion haliaetus", "Pernis ptilorhynchus", "Falco peregrinus"),
     labels= c(expression(italic("F. eleonorae")), expression(italic("P. haliaetus")),
                                                              expression(italic("P. ptilorhynchus")),
                                                              expression(italic("F. peregrinus"))),
     tick=T ,col = NA, col.ticks = 1, # NULL would mean to use the defult color specified by "fg" in par
     tck=-.015 , #tick marks smaller than default by this proportion
     las=2 ) # text perpendicular to axis label 

#add legend
legend(x = 1.8, y = 0.6, legend=c("wind support", "wind speed", expression(paste(Delta,"t"))), 
       col = c("paleturquoise2","yellowgreen","lightcoral"), #coords indicate top-left
       pch = 19, bg="white",bty="n", cex = 0.75)

dev.off()



#plotting interaction effects ####

#INLA does not have a predict function. prediction is the same as fitting the model with missing data and need to be defined at the time of model fitting (INLA FAQ).
#(in statistical rethinking, the function "link" is used for making predictions)
#1) add the new data, with NAs for used, to the dataset all_data. 2) run the model with control predictor = list(compute = T). 3) to show predictive 
#distribution, obtain indices for rows containing the NA values (source: bayesian inference with inla)
used_na <- which(is.na(all_data$used))
#extract summary statistics of the fitted values
model$summary.fitted.values[used_na,c("mean","sd")]

#create a new dataset to be used for prediction
n <- 1000
dt <- sample(all_data$delta_t_z, n,replace = T) #random sample of delta t values (recycled from observed values)
wspd <- sample(c(-2,0,2),n, replace = T) #keep wind speed values at -2,0 and 2
wspt <- rep(0,n) #keep wind support values at 0 
species <- sample(all_data$species, n, replace =T) 
ind <- sample(all_data$ind, n, replace = T)
used <- rep(NA,n)


n = 1000
n.pred = 1000
y <- arima.sim(n = n, model = list(ar = 0.9))
N <- n + n.pred
yy <- c(y,rep(NA,n.pred))
i = 1:N
formula <- yy ~ f(i,model = "ar1")
r <- inla(formula,data = data.frame(i,yy),
          control.family = list(initial = 10, fixed = TRUE))
#predictions
r$summary.random$i[(n+1):N,c("mean","sd")]

#with inspo from rethinking R code 7.28
par(mfrow = c(1,3))

dt.seq <- -1:1
for ( wsp in -1:1 ) {
  dt <- all_data[all_data$wind_speed_z == wsp,]
  plot( used ~ delta_t_z , data=dt , col=rangi2 ,
        main=paste("wspeed =",w) , xaxp=c(-1,1,2) , ylim=c(0,362) ,
        xlab="delta_t (centered)" )
  mu <- link( mc1 , data=data.frame(water.c=w,shade.c=shade.seq) ) #dont know the equivalent of this in inla
  mu.mean <- apply( mu , 2 , mean )
  mu.PI <- apply( mu , 2 , PI , prob=0.97 )
  lines( shade.seq , mu.mean )
  lines( shade.seq , mu.PI[1,] , lty=2 )
  lines( shade.seq , mu.PI[2,] , lty=2 )
}











par(mfrow = c(2, 2), mar = c(2.2, 2, 2, 2))

# Fixed effects
plot(m1$marginals.fixed[[2]], type = "l",
     main = "Treatment effect", xlab = "", ylab = "density")
for(i in 3:4) {
  lines(m1$marginals.fixed[[i]], lty = i)
}
legend("topleft", paste("Treat. ", 2:4, sep = ""), lty = 2:4, bty = "n",
       cex = 0.75)

# Random effects
plot(m1$marginals.random$blend[[1]], type = "l", xlim = c(-5, 5),
     ylim = c(0, 0.8),
     main = "Blend effect", xlab = "", ylab = "density")
for(i in 2:5) {
  lines(m1$marginals.random$blend[[i]], lty = i)
}
legend("topleft", paste("Blend ", 1:5, sep = ""), lty = 1:5, bty = "n",
       cex = 0.75)

# Variances
marg.variances <- lapply(m1$marginals.hyperpar,
                         function(m) inla.tmarginal(function(x) 1 / x, m)
)

plot(marg.variances[[1]], type = "l", main = "Error variance",
     xlab = "", ylab = "density")
plot(marg.variances[[2]], type = "l", main = "Blend variance",
     xlab = "", ylab = "density")



plot1 <- image(Z_ind, xlab = "Random effect", ylab = "Observation",
               main = "ind", sub = "")
plot2 <- image(Zlts, xlab = "Random effect", ylab = "Observation",
   main = "Sample", sub = "")
 
 print(plot2, position = c(0.5, 0, 1, 1), more = TRUE)
print(plot1, position = c(0, 0, 0.5, 1), more = FALSE)



#Match subjects to groups
dd <-all_data
idx <- match(dd$species, m1d$summary.random$species1$ID) #for delta_t_z

# Intercept
dd$intercept <- m1d$summary.fixed[1, "mean"]

# Slope
dd$slope <- m1d$summary.random$species1$mean[idx]


xyplot(used ~ delta_t_z | species, data = dd, 
       a.int = dd$intercept,
       b.slo = dd$slope, type = c("g","p", "r"),
       panel = function(x, y, a.int, b.slo, ..., subscripts){
         panel.xyplot(x, y, ...)
         panel.abline(a.int[subscripts][1], b.slo[subscripts][1])
       },
       xlab = "delta_t_z",
       ylab = "prob used", aspect = "xy")



#interpretation: not much variation among species (? )

#### decent plot for interaction terms
dt <- data.frame(ID = m1b$summary.random$lat_zone1$ID,
                  mean = m1b$summary.random$lat_zone1$mean,
                  IClower = m1b$summary.random$lat_zone1[, 4],
                  ICupper = m1b$summary.random$lat_zone1[, 6])
ws <- data.frame(ID = m1b$summary.random$lat_zone2$ID,
                 mean = m1b$summary.random$lat_zone2$mean,
                 IClower = m1b$summary.random$lat_zone2[, 4],
                 ICupper = m1b$summary.random$lat_zone2[, 6])
cw <- data.frame(ID = m1b$summary.random$lat_zone3$ID,
                 mean = m1b$summary.random$lat_zone3$mean,
                 IClower = m1b$summary.random$lat_zone3[, 4],
                 ICupper = m1b$summary.random$lat_zone3[, 6])
ws_v <- data.frame(ID = m1b$summary.random$lat_zone4$ID,
                 mean = m1b$summary.random$lat_zone4$mean,
                 IClower = m1b$summary.random$lat_zone4[, 4],
                 ICupper = m1b$summary.random$lat_zone4[, 6])
dt_v <- data.frame(ID = m1b$summary.random$lat_zone5$ID,
                   mean = m1b$summary.random$lat_zone5$mean,
                   IClower = m1b$summary.random$lat_zone5[, 4],
                   ICupper = m1b$summary.random$lat_zone5[, 6])

X11();par(mfrow = c(2,3))
plot(ID ~ mean, data = dt, pch = 16, xlab = "coefficient for delta_t", ylab = "species")
arrows( mean - IClower, ID, mean + IClower, ID, data = dt, length = 0.5, code = 3) #angle 90 would make it vertical

points(ID ~ mean, data = ws, pch = 16, xlab = "coefficient for wind support", ylab = "species")
arrows( mean - IClower, ID, mean + IClower, ID, data = tab_ws, length = 0.5, code = 3)


##### smooth effects (virgilio's book) #####
tab.rw1 <- data.frame(x = m4$summary.random$ws_grp[, "ID"],
                      y = m4$summary.random$ws_grp[, "mean"]
)
tab.rw2 <- data.frame(x = m.grp.rw2$summary.random$range[, "ID"],
                      y = m.grp.rw2$summary.random$range.grp[, "mean"]
)

ggplot(aes(x = wind_support_z, y = used), data = sample) + 
  geom_point(aes(x = sample$ws_grp, y = sample$used), colour = "grey") +
  geom_point() +
  ggtitle("sample data") +
  geom_line(aes(x = x, y = y, linetype = "solid"), data = tab.rw1) +
  #geom_line(aes(x = x, y = y, linetype = "dashed"), data = tab.rw2) +
  scale_linetype_manual(name = "Smoothing method (grouped range)",
                        values = c("solid", "dashed"),
                        labels = c("rw2", "rw1")) +
  theme(legend.position = "bottom")

###
tab.rw1 <- data.frame(x = vir.rw1$summary.random$dose[, "ID"],
                      y = vir.rw1$summary.random$dose[, "mean"])
tab.rw2 <- data.frame(x = vir.rw2$summary.random$dose[, "ID"],
                      y = vir.rw2$summary.random$dose[, "mean"])
ggplot() + #aes(x = dose, y = numdead / total), data = H.virescens) +
  #geom_point() +
  xlab(expression(paste("Dose (", mu, "g)"))) +
  ylab("f(dose)") +  
  ggtitle("Dose response data (smooth term)") +
  geom_line(aes(x = x, y = y, linetype = "solid"), data = tab.rw1) +
  geom_line(aes(x = x, y = y, linetype = "dashed"), data = tab.rw2) +
  scale_linetype_manual(name = "Smoothing method",
                        values = c("solid", "dashed"),
                        labels = c("rw2", "rw1")) +
  theme(legend.position = "bottom")




# interaction plots (https://groups.google.com/forum/#!topic/r-inla-discussion-group/DwPyw_fOUGI)
# make sure to add compute = TRUE to the control.predictor list