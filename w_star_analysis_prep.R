# Scripts for estimating w_star and its relationship with delta_t
# This is script 4 of 6 for reproducing the results of Nourani et al 2021, ProcB.
# Elham Nourani, PhD. Jun.10. 2021
#-----------------------------------------------------------------

#annotate a dataset to be used both for w_star (ERA_interim) and raw maps (ERA5)

# ---------- STEP 1: load data #####

source("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/wind_support_Kami.R")

wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

#annotated data.(prepped in 2021_wind_delta_t_maps.R)
ann <- read.csv("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/2021/public/annotation/raw_points_for_maps_updated.csv-3143330181669857741/raw_points_for_maps_updated.csv-3143330181669857741.csv") %>% 
  mutate(timestamp,timestamp = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>%
  rename(sst_i = ECMWF.Interim.Full.Daily.SFC.FC.Sea.Surface.Temperature,
         t2m_i = ECMWF.Interim.Full.Daily.SFC.FC.Temperature..2.m.above.Ground.,
         blh_i = ECMWF.Interim.Full.Daily.SFC.FC.Boundary.Layer.Height,
         s_flux_i = ECMWF.Interim.Full.Daily.SFC.FC.Instantaneous.Surface.Heat.Flux,
         m_flux_i = ECMWF.Interim.Full.Daily.SFC.FC.Instantaneous.Moisture.Flux,
         sst_5 = ECMWF.ERA5.SL.Sea.Surface.Temperature,
         t2m_5 = ECMWF.ERA5.SL.Temperature..2.m.above.Ground.,
         u925_5 = ECMWF.ERA5.PL.U.Wind,
         v925_5 = ECMWF.ERA5.PL.V.wind) %>%
  mutate(delta_t_5 = sst_5 - t2m_5,
         delta_t_i = sst_i - t2m_i) %>%
  st_as_sf(coords = c("location.long", "location.lat"), crs = wgs) %>% 
  mutate(s_elev_angle = solarpos(st_coordinates(.), timestamp, proj4string=CRS("+proj=longlat +datum=WGS84"))[,2]) %>% #calculate solar elevation angle
  mutate(sun_elev = ifelse(s_elev_angle < -6, "night", #create a categorical variable for the position of the sun
                           ifelse(s_elev_angle > 40, "high", "low"))) %>% 
  as("Spatial") %>% 
  as.data.frame() %>% 
  rename(lat = coords.x2,
         lon = coords.x1) %>% 
  arrange(track, timestamp)

#remove tracks with one point
more_than_one_point <- ann %>% 
  group_by(track) %>% 
  summarize(n = n()) %>% 
  filter(n > 2)

ann <- ann %>% 
  filter(track %in% more_than_one_point$track) %>% 
  dplyr::select(-"X")

#remove duplicated timestamps
#remove duplicated rows
rows_to_delete <- unlist(sapply(getDuplicatedTimestamps(x = as.factor(ann$track),timestamps = ann$timestamp),"[",-1)) #get all but the first row of each set of duplicate rows
ann_pts <- ann[-rows_to_delete,]



save(ann_pts, file = "/home/enourani/ownCloud/Work/Projects/delta_t/R_files/2021/public/Dryad/annotated_points.RData")
