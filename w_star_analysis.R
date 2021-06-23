# Scripts for estimating w_star and its relationship with delta_t
# This is script 4 of 6 for reproducing the results of Nourani et al 2021, ProcB.
# Elham Nourani, PhD. Jun.10. 2021
#-----------------------------------------------------------------



# ---------- STEP 1: load data #####

source("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/wind_support_Kami.R")

wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

#annotated data.(prepped in 2021_wind_delta_t_maps.R)
ann <- read.csv("R_files/2021/annotations/raw_points_for_maps/interim_updated/raw_points_for_maps_updated.csv-5066586225337282923/raw_points_for_maps_updated.csv-5066586225337282923.csv") %>% 
  mutate(timestamp,timestamp = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>%
  rename(sst = ECMWF.Interim.Full.Daily.SFC.FC.Sea.Surface.Temperature,
         t2m = ECMWF.Interim.Full.Daily.SFC.FC.Temperature..2.m.above.Ground.,
         blh = ECMWF.Interim.Full.Daily.SFC.FC.Boundary.Layer.Height,
         s_flux = ECMWF.Interim.Full.Daily.SFC.FC.Instantaneous.Surface.Heat.Flux,
         m_flux = ECMWF.Interim.Full.Daily.SFC.FC.Instantaneous.Moisture.Flux) %>%
  mutate(delta_t = sst - t2m) %>%
  st_as_sf(coords = c("location.long", "location.lat"), crs = wgs) %>% 
  mutate(s_elev_angle = solarpos(st_coordinates(.), timestamp, proj4string=CRS("+proj=longlat +datum=WGS84"))[,2]) %>% #calculate solar elevation angle
  mutate(sun_elev = ifelse(s_elev_angle < -6, "night", #create a categorical variable for teh position of the sun
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
  filter(n > 1)

ann <- ann %>% 
  filter(track %in% more_than_one_point$track)

save(ann, file = "R_files/2021/df_for_w_star_updated.RData")