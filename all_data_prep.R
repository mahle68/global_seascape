#script for preparing all input data simultaneously. previously done separately for each species/region
#Elham Nourani,
#Dec. 31. 2019. Radolfzell, Germany.

library(readxl) #read_excel()
library(lubridate)
library(move)



##### STEP 1: read in the data and filter for adult birds and season and databse-specific filters #####

#also assign date_time, year, month, track, species

OHB_files <- list.files("data/Oriental_honey_buzzard",pattern = ".xls",full.names = T)
OHB <- lapply(OHB_files,read_excel,1,col_types = c("numeric","date","numeric","numeric","numeric","skip","text",rep("numeric",8))) %>%
  reduce(full_join) %>%
  rename(date_time = 'date(gmt)',lon = longitud,lat = latitude) %>%
  mutate(yday = yday(date_time)) %>%
  mutate(season = ifelse(between(yday,253,294),"autumn",ifelse(month == 5,"spring","other")),
         track = paste(ptt,year,sep = "_"),
         species = "OHB") %>% #11 Sep-20 Oct; spring between 1-5 May
  filter(season == "autumn" & #no sea-crossing in spring
           class %in% c("0","1","2","3")) #filter for location classes


GFB_files <- list.files("data/Grey_faced_buzzard/",pattern = ".csv",full.names = T)
GFB <- lapply(GFB_files,read.csv,stringsAsFactors = F) %>%
  reduce(full_join) %>% #is locdate is in UTC
  mutate(locdate,date_time = as.POSIXct(strptime(locdate,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>%
  mutate(month = month(date_time),
         year = year(date_time),
         season = ifelse(month %in% c(3,4),"spring",ifelse(month %in% c(10),"autumn","other")),
         track = paste(platform,year,sep = "_"),
         species = "GFB") %>% #(Nourani et al 2017 for autumn; no sea-crossing in spring based on visually exploring the data.)
  filter(season %in% c("spring","autumn") &
           class %in% c("0","1","2","3")) #filter for location classes

###pf have unknown age! ask Ivan....
PF <- read.csv("data/LifeTrack Peregrine falcon.csv", stringsAsFactors = F) %>% 
  dplyr::select(1,3:5,16,38,39) %>% #remove columns that are not needed
  mutate(date_time = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  mutate(month = month(date_time),
         year = year(date_time),
         track = paste(tag.local.identifier, year,sep = "_"),
         species = "PF",
         season = ifelse (month %in% c(9,10), "autumn", "other")) %>% 
  filter(season != "other")

OE <- read.csv("data/Osprey in Mediterranean (Corsica, Italy, Balearics).csv", stringsAsFactors = F) %>% 
  filter(grepl("ad",individual.local.identifier,, ignore.case = T) & !grepl("juv",individual.local.identifier,, ignore.case = T)) %>% #extract adult data
  dplyr::select(1,3:5,16,35:37) %>% #remove columns that are not needed
  mutate(date_time = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  mutate(month = month(date_time),
         year = year(date_time),
         track = paste(tag.local.identifier, year,sep = "_"),
         species = "O",
         season = ifelse(month %in% c(2:4),"spring",ifelse (month %in% c(8:10), "autumn","other"))) %>% 
  filter(season != "other")

OA <- read.csv("data/Osprey_Americas/Osprey Bierregaard North and South America.csv", stringsAsFactors = F) %>% 
  dplyr::select(1,3:5,48:52) %>% #remove columns that are not needed
  mutate(date_time = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  mutate(month = month(date_time),
         year = year(date_time),
         track = paste(tag.local.identifier, year,sep = "_"),
         species = "O",
         season = ifelse(month %in% c(3,4),"spring",ifelse (month %in% c(9,10), "autumn", "other"))) %>% 
  filter(season != "other") %>% 
  filter(sensor.type == "gps",#keep only the gps tag data
         individual.local.identifier %in% c("Holly","Hackett","Daphne","Shanawdithit","Gundersen","Wausau",
                                            "Crabby","Charlie","Roger Tory")) #keep only adults (based on what i found on the website)

##### STEP 2: merge all data, assign zone, etc #####


