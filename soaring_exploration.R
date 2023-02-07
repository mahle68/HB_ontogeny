#quickly view the soaring behavior
#5.10.2022, Konstanz. DE
#Elham Nourani, PhD

library(tidyverse)
library(move)
library(mapview)
library(lubridate)
library(sf)
library(rgdal)
library(rgl)


wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")


creds <- movebankLogin(username = "mahle68", rstudioapi::askForPassword())
EHB_FN_id <- getMovebankID("European Honey Buzzard_Finland", creds)
data <- getMovebankData(EHB_FN_id, animalName = "D225_236",
                        removeDuplicatedTimestamps = T, login = creds,
                        timestamp_start = "2022100312000000")

soaring <- as.data.frame(data) %>% 
  filter(timestamp < as.POSIXct( "2022-10-03 13:00:00", tz = "UTC")) %>% 
  st_as_sf(coords = c("location_long.1", "location_lat.1"), crs = wgs) %>% 
  #convert projection to utm
  st_transform(crs = "+proj=utm +zone=35 +ellps=GRS80 +datum=NAD83")

extent <- st_bbox(c(xmin = 614657.4, xmax = 619376.9, ymax = 3823869, ymin = 3819000), crs = "+proj=utm +zone=35 +ellps=GRS80 +datum=NAD83")
extent2 <- st_bbox(c(xmin = 3825910, xmax = 3829238, ymin = 262.9544, ymax = 314.3205), crs = "+proj=utm +zone=35 +ellps=GRS80 +datum=NAD83")


subset <- soaring %>%
  st_crop(extent)

png("/home/enourani/ownCloud/Work/Collaborations/powered_glide/altitude.png", width = 7, height = 7, units = "in", res = 300)

plot(st_coordinates(soaring)[,2], soaring$height_above_ellipsoid, pch = 20, ylab = "height above ellipsoid", xlab = "latitude")
points(st_coordinates(subset)[,2], subset$height_above_ellipsoid, pch = 20, ylab = "height above ellipsoid", xlab = "latitude", col = "red", add = T)

dev.off()

png("/home/enourani/ownCloud/Work/Collaborations/powered_glide/latlong_subset.png", width = 7, height = 7, units = "in", res = 300)
plot(st_coordinates(subset)[,1], st_coordinates(subset)[,2], pch = 20, cex = 0.4, xlab = "Longitude", ylab = "latitude", col = "red")
dev.off()

png("/home/enourani/ownCloud/Work/Collaborations/powered_glide/latlong_track.png", width = 7, height = 7, units = "in", res = 300)
plot(st_coordinates(soaring)[,1], st_coordinates(soaring)[,2], pch = 20, cex = 0.4, xlab = "Longitude", ylab = "latitude")
points(st_coordinates(subset)[,1], st_coordinates(subset)[,2], pch = 20, cex = 0.4, xlab = "Longitude", ylab = "latitude", col = "red", add = T)
dev.off()



#3D plot
#soaring_df <- data %>% 
#  as.data.frame() %>% 
#  drop_na(height_above_ellipsoid)


plot3d(x = st_coordinates(soaring)[,1], y = st_coordinates(soaring)[,2], z = soaring$height_above_ellipsoid, 
       xlab="Longitude", ylab="Latitude", zlab = "Height above ellipsoid", col = "blue")#,
       #xlim = c(min(soaring_df$location_long),29), ylim = c(min(soaring_df$location_lat),35), forceClipregion = T)
rglwidget()



plot(soaring_df$timestamp,soaring_df$height_above_ellipsoid, pch = 20, ylab = "height above ellipsoid", xlab = "hour")



plot(soaring_df$location_long,soaring_df$height_above_ellipsoid, pch = 20, ylab = "height above ellipsoid", xlab = "longitude")

plot3d(x = soaring_df$location_long, y = soaring_df$location_lat, z = soaring_df$height_above_ellipsoid, 
       xlab="Longitude", ylab="Latitude", zlab = "Height above ellipsoid", col = "blue", forceClipregion = T)
rglwidget()

one_day_sf <- st_as_sf(soaring_df, coords = c("location_long", "location_lat"), crs = wgs)





one_hr <- soaring_df %>% 
  mutate(day = day(timestamp),
         hr = hour(timestamp)) %>% 
  filter(day == 3 & between(hr, 9,10))


one_day <- soaring_df %>% 
  mutate(day = day(timestamp),
         hr = hour(timestamp)) %>% 
  filter(day == 3 & between(hr, 13,14))

plot(one_day$timestamp,one_day$height_above_ellipsoid, pch = 20, ylab = "height above ellipsoid", xlab = "hour")

plot(one_day$location_lat,one_day$height_above_ellipsoid, pch = 20, ylab = "height above ellipsoid", xlab = "latitude")

plot(one_day$location_long,one_day$height_above_ellipsoid, pch = 20, ylab = "height above ellipsoid", xlab = "longitude")

plot3d(x = one_day$location_long, y = one_day$location_lat, z = one_day$height_above_ellipsoid, 
       xlab="Longitude", ylab="Latitude", zlab = "Height above ellipsoid", col = "blue", forceClipregion = T)
rglwidget()

one_day_sf <- st_as_sf(one_day, coords = c("location_long", "location_lat"), crs = wgs)
  

one_glide <- one_day %>% 
  filter(between(timestamp,"2022-10-03 14:44:13.000", ))
#hester's code
library(rgl)
df <- read.csv("golden_eagle_data.csv", header = T, stringsAsFactors = F)

plot3d(df$location_long[df$flight_type == "soar_circular"], df$location_lat[df$flight_type == "soar_circular"], 
       df$height.above.ellipsoid[df$flight_type == "soar_circular"], xlim=range(df$location_long), 
       ylim=range(df$location_lat), ticktype="detailed", xlab="Longitude", ylab="Latitude", 
       zlab="Height above ellipsoid", type="p", col = "red", size = 2)
plot3d(df$location_long[df$flight_type == "glide"], df$location_lat[df$flight_type == "glide"],
       df$height.above.ellipsoid[df$flight_type == "glide"], col = "blue", add = T)
plot3d(df$location_long[df$flight_type == "soar_linear"], df$location_lat[df$flight_type == "soar_linear"],
       df$height.above.ellipsoid[df$flight_type == "soar_linear"], col = "green", add = T)
plot3d(df$location_long[df$flight_type == "flap"], df$location_lat[df$flight_type == "flap"],
       df$height.above.ellipsoid[df$flight_type == "flap"], col = "purple", add = T) 
legend3d("topright", legend = c('Thermal soaring', 'Gliding', 'Linear soaring', 'Flapping'),  pch = 16, col = c("red","blue","green","purple"), magnify=2, inset=c(0.02)) 


