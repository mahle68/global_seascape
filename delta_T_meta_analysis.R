# script for bringing together annotated data prepared previously for different species and analyzing them together.
#Elham Nourani
#Dec. 31, Radolfzell, Germany

library(tidyverse)
library(sf)
library(lme4)
library(raster)
library(parallel)

setwd("/home/enourani/ownCloud/Work/Projects/delta_t")
wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
meters_proj <- CRS("+proj=moll +ellps=WGS84")
mycl <- makeCluster(detectCores() - 1)

##### STEP 1: open annotated data and add unique obs_id#####

data_ls <- lapply( c(1:4),function(i){
  load(list.files(path = "R_files",pattern = "ann_14days", full.names = T)[[i]])
})



data_ls[[1]]$species <- "O"
data_ls[[2]]$species <- "A"

