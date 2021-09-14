###############################################################################################
##Map of Oregon Coast and Otter Rock Marine Reserve
## By: Jenn Fields
## Last updated: 9.13.21

#clear dataframe#####
rm(list=ls())

##Install packages
# load packages
library(tidyverse)
library(lubridate)
library(data.table)
library(ggrepel)
library(ggpubr)
library(ggsn)
library(magick)
library(maps)
library(maptools)
library(ggplot2)
library(grid)
library(ggimage)
#devtools::install_github("dkahle/ggmap")
library(ggmap)
library(patchwork)

#zoomed in map version
#make map and set sources
ggmap::register_google(key = 'AIzaSyAB5u0qqKdGPAbV-JgjkujmycKk3ytdr5o', write= TRUE)
OtterRock <- get_googlemap(center = c(-124.065687, 44.752504), zoom = 17, maptype = 'satellite' ) 

#Get the bounding box for the map and reformat it into a data.frame for scalebar
or <- attr(OtterRock, "bb")
or2 <- data.frame(long = unlist(bb[c(2, 4)]), lat = unlist(bb[c(1,3)]))
 
#Add the scalebar to a ggmap, need ggsn package
OtterRockpools <- ggmap(OtterRock) + 
  labs(x = 'Longitude', y = 'Latitude') +
  scalebar(data = or2, dist =10, st.color = "white",transform = TRUE, 
           model  = "WGS84",  dist_unit = "m", #data are bounding box of map, model gives datum from google maps, 
           #transform recognizes as GPS points as decimal units, location sets location on map, anchor sets bounds of location non map
           location = "bottomright", st.dist = 0.1, box.fill = "white", #sets scalebar bottom left, st.dist is distance between box and numbers,box.fill fills box
           anchor = c( x = or$ll.lon, y = or$ll.lat)) #sets scalebar in specific position long x lat
OtterRockpools
#add north symbol to map
north2(OtterRockpools, x=.28, y=.35, symbol=4) 


#zoomed out map version
#make map and set sources
ORCoast <- get_googlemap(center = c(-121.285627, 43.784211), zoom = 5, maptype = 'satellite' ) 


#Get the bounding box for the map and reformat it into a data.frame for scalebar
bb <- attr(ORCoast, "bb")
bb2 <- data.frame(long = unlist(bb[c(2, 4)]), lat = unlist(bb[c(1,3)]))

#Add the scalebar t o a ggmap, need ggsn package
WestCoastmap <- ggmap(ORCoast) + 
  labs(x = 'Longitude', y = 'Latitude') +
  scalebar(data = bb2, dist = 500, st.color = "white",transform = TRUE, 
           model  = "WGS84",  dist_unit = "km", #data are bounding box of map, model gives dataum from google maps, transform recognizes as GPS points as decimal units, location sets location on map, anchor sets bounds of locatio non map
           location = "topright", st.dist = 0.05, box.fill = "white", #sets scalebar bottom left, st.dist is distance bewteen box and numbers,box.fill fills box
           anchor = c( x = bb$ll.lon+13, y = bb$ll.lat+7))+  #sets scalebar in specific position long x lat
  geom_point(aes(x = -124.066605, y = 44.752906),
             alpha = .5, color=c('#fff7bc'), size = 2, shape = 19)  #geom_point adds specific plot points to map
WestCoastmap
#add north symbol to map
north2(WestCoastmap, x=.12, y=.29, symbol=4)#edit the lettered text

Mapsfigure<-WestCoastmap+OtterRockpools+
  plot_annotation(tag_levels = 'a') &         #label each individual plot with letters A-G
  theme(plot.tag = element_text(size = 50)) 
Mapsfigure

