##Tide Pool Physical Parameters
##By: Jenn Fields
##Last updated: 10.23.2020

rm(list=ls()) #Clears the environment
#load libraries
library(tidyverse)
library(factoextra)
library(devtools)
library(ggfortify)
library(RColorBrewer)

##load data
poolSA <- read_csv("Data/PoolPhysicalParameters/poolSA.csv")
TidePooldes <-read_csv("Data/PoolPhysicalParameters/TidePoolDes.csv")
OregonTidePoolVolumes <- read_csv("Data/PoolPhysicalParameters/OregonTidePoolVolumes.csv")

####Calculate surface area####
x<-poolSA$squares #create list for function

SurfaceArea<- (x *  0.01) #calculating surface area from number of squares and converting from cm^2 to m^2
#View(SurfaceArea)
TidePooldes$SurfaceArea<-SurfaceArea #create column in data frame for surface area
#View(poolSA)


##Pool Volumes###
#take average of triplicate measures
PoolVolumes <- OregonTidePoolVolumes %>%
   dplyr::group_by(PoolID, Before_After, Vol_used) %>%
   dplyr::select(Absorbance_Raw, Vol_used, PoolID, Before_After,Volume) %>%
   dplyr::summarise(Absorbance = mean(Absorbance_Raw))

#Equations created from seawater curve script in SilbigerLab TidePoolVolume with SS3300 spec
# y~I(b*x^z)
#for 3 mL dye: tide pool volume = (4.64 * absborbance^ -1.07)
# for 6 mL dye: tide pool volume = (9.36 * absborbance ^ -1.07)
# for 9 mL dye: tide pool volume = (12.7 * absborbance ^ -1.02)
#for 15 mL dye: tide pool volume = (22.3 * absorbance ^ -1.00)

#uses equation for each volume of dye
Volume3<-PoolVolumes %>%
   filter(Vol_used== '3') %>%
   mutate(VolumeL = 4.64 * (Absorbance^-1.07))

Volume6<-PoolVolumes %>%
   filter(Vol_used == '6') %>%
   mutate(VolumeL = 9.36 * (Absorbance^-1.07))

Volume9<-PoolVolumes %>%
   filter(Vol_used == '9') %>%
   mutate(VolumeL = 12.7 * (Absorbance^-1.02))

Volume15<-PoolVolumes %>%
   filter(Vol_used == '15') %>%
   mutate(VolumeL = 22.3 * (Absorbance^-1.00))

PoolVolumes<-rbind(Volume3,Volume6,Volume9,Volume15)

PoolVolumes$Vol<- PoolVolumes$VolumeL * 0.001 #convert volume from L to m3 
Volumes<-PoolVolumes[, c("PoolID", "Before_After", "Vol")] #pull out the necessary columns for joining
TidePooldes<-left_join(TidePooldes,Volumes) #joining with rest of physical parameters
TidePooldes$SAtoV <- (TidePooldes$SurfaceArea/TidePooldes$Vol)  
