


library(sf)
library(tidyverse)
library(maptools)
library(mapview)

wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")  

test <- read.csv("/home/mahle68/ownCloud/Work/Projects/delta_t/R_files/2021/annotations/test_points_for_flux_sign.csv-2054261352130105090/test_points_for_flux_sign.csv-2054261352130105090.csv")%>% 
  full_join(read.csv("/home/mahle68/ownCloud/Work/Projects/delta_t/R_files/2021/annotations/test_point_for_flux_sign.csv-1375703174786464509/test_point_for_flux_sign.csv-1375703174786464509.csv")) %>% 
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
                           ifelse(s_elev_angle > 40, "high", "low"))) 
