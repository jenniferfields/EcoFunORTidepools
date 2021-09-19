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

#tide pool gps points and comm comp data
TidePoolGPSpts <- read_csv("Data/PoolPhysicalParameters/TidePoolGPSpts.csv")
OregonTidePoolVolumes <- read_csv("Data/PoolPhysicalParameters/OregonTidePoolVolumes.csv")
TidepoolGPS<-left_join(TidePoolGPSpts,OregonTidePoolVolumes)
TidepoolGPS<-TidepoolGPS[-c(4:27)]
#zoomed in map version
#make map and set sources
ggmap::register_google(key = 'AIzaSyAB5u0qqKdGPAbV-JgjkujmycKk3ytdr5o', write= TRUE)
OtterRock <- get_googlemap(center = c(-124.065687, 44.752504), zoom = 17, maptype = 'satellite' ) 

#Get the bounding box for the map and reformat it into a data.frame for scalebar
bb <- attr(OtterRock, "bb")
bb2 <- data.frame(long = unlist(bb[c(2, 4)]), lat = unlist(bb[c(1,3)]))
 
#Add the scalebar to a ggmap, need ggsn package
OtterRockpools <- ggmap(OtterRock) + 
  labs(x = 'Longitude', y = 'Latitude') +
  scalebar(data = bb2, dist =100, st.color = "white",transform = TRUE, 
           model  = "WGS84",  dist_unit = "m", #data are bounding box of map, model gives datum from google maps, 
           #transform recognizes as GPS points as decimal units, location sets location on map, anchor sets bounds of location non map
          location = "bottomright", st.dist = 0.03, box.fill = "white", #sets scalebar bottom left, st.dist is distance between box and numbers,box.fill fills box
           anchor = c(x = -124.0665, y = 44.7505)) +#sets scalebar in specific position long x lat
  geom_point(TidepoolGPS, mapping=aes(x=lon, y=lat, group=Foundation_spp, color=Foundation_spp, shape=Foundation_spp),
             size=10)+
  scale_color_manual(values=c('#c6dbef','#a1d99b'))+
  theme(legend.position = 'none',
        axis.title =element_blank(),
        axis.text=element_text(size=35))
OtterRockpools

#add north symbol to map
ggsave("output/OtterRockmap.pdf",useDingbats = FALSE, width=75, height=60,dpi=600, unit="cm")


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
           anchor = c( x = -124, y = 35))+  #sets scalebar in specific position long x lat
  geom_point(aes(x = -124.066605, y = 44.752906), color=c('#fff7bc'), size = 15, shape = 8,stroke=4) + #geom_point adds specific plot points to map
theme(legend.position = 'none',
      axis.text =element_text(size=35),
      axis.title=element_text(size=40))
WestCoastmap
#add north symbol to map
north2(WestCoastmap, x=.12, y=.29, symbol=4)#edit the lettered text
ggsave("output/Northarrow.pdf",useDingbats = FALSE, width=75, height=60,dpi=600, unit="cm")
Mapsfigure<-WestCoastmap+OtterRockpools+
  plot_annotation(tag_levels = 'a') &         #label each individual plot with letters A-G
  theme(plot.tag = element_text(size = 50)) 
Mapsfigure
ggsave("output/Coastalmap.pdf",useDingbats = FALSE, width=75, height=60,dpi=600, unit="cm")

