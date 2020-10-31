##Temperature and Light one month analysis in OR tide pools
##By: Jenn Fields
##Last updated 10.30.2020

rm(list=ls()) #Clears the environment
#load libraries
library(lubridate)
library(plotrix)
library(ggpubr)
library(MASS)
library(car)
library(cowplot)
library(tidyverse)
library(GGally) #for collinearity with mult regressions
library(ggeffects) #effects for multiple regression
library(ggsidekick) #theme sleek
library(patchwork) #patch plots together with ease!



#Bring in TP physical parameters for covariates
#PCATPDesc <- read_csv("Data/PoolPhysicalParameters/PCATPDesc.csv")
source("scripts/tidepoolphysicalparameters.R")
source("scripts/CommunityComp.R")
#View(PCATPDesc)


#########Cleaning Light and Temp logger data#################
#options(na.action = "na.omit") 
#group temp data from all Hobo pendant loggers together
temp.data <-
  list.files(path = 'Data/LightandTemp',pattern = ".csv", full.names = TRUE) %>% 
  # list files in directory following a particular pattern
  set_names(.) %>% # get the column names
  map_dfr(read.csv, .id = "file.ID") %>% # join all files together in one data frame by file ID
  group_by(file.ID) 

temp.data$Date.Time<- mdy_hm(temp.data$Date.Time, quiet=FALSE, tz="America/Los_Angeles", truncated=0) #format date and time 
#View(temp.data)

#convert lux to par
#parameters from Long 2012 Field exp values
#A1 = -4924.7
#t1 = 20992.9
#y0 =  4929.0
#PARLICOR = A1 e^(–HOBO/t1) + y0
A1 <- -4924.7
t1 <- 20992.9
y0 <- 4929.0

#covert lux to numeric and get rids of pesky , in lux values
temp.data$Intensity.lux <- as.numeric(gsub(",","",temp.data$Intensity.lux))
x <- temp.data$Intensity.lux #create vector for just lux column
x[x == 0] <- NA #convert 0 to NA
#x[x >= 	19671] <- NA #get rid of values over 	19671 lux/30000PAR (not natural?)
#View(x)

# enter into equation: PARLICOR = A1e(–HOBO/t1) + y
#assign to new column for Par
temp.data$Par<- (A1* (exp( -x / t1))) + y0

temp.data$Par[is.na(temp.data$Par)]<- 0
temp.data<-as.data.frame(temp.data)


#Control period data before removal of foundation spp
Control.period<-temp.data %>%
  filter(!(PoolID == "1" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "1" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "2" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "2" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "3" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "3" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "4" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "4" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "5" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "5" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "6" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "6" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "7" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "7" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "8" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "8" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "9" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "9" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "10" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "10" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "11" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "11" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "12" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "12" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "13" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "13" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "14" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "14" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "15" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "15" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "16" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "16" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "17" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "17" & Date.Time > "2019-07-15 0:00:00")) %>% 
  filter(!(PoolID == "18" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "18" & Date.Time > "2019-07-15 0:00:00")) %>% 
  filter(!(PoolID == "19" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "19" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "20" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "20" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "21" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "21" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "22" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "22" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "23" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "23" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "24" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "24" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "25" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "25" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "26" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "26" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "27" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "27" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "28" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "28" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "29" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "29" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "30" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "30" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "31" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "31" & Date.Time > "2019-07-15 0:00:00")) %>%
  filter(!(PoolID == "32" & Date.Time < "2019-06-16 00:00:00")) %>%
  filter(!(PoolID == "32" & Date.Time > "2019-07-15 0:00:00"))
#View(Control.period)

#Temp data after removal of foundation spp
Removal.period<-temp.data %>%
  filter(!(PoolID == "1" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "1" & Date.Time > "2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "2" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "2" & Date.Time > "2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "3" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "3" & Date.Time > "2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "4" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "4" & Date.Time > "2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "5" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "5" & Date.Time > "2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "6" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "6" & Date.Time >"2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "7" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "7" & Date.Time > "2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "8" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "8" & Date.Time > "2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "9" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "9" & Date.Time > "2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "10" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "10" & Date.Time > "2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "11" & Date.Time < "2019-07-17 4:30:00")) %>%
  filter(!(PoolID == "11" & Date.Time > "2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "12" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "12" & Date.Time > "2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "13" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "13" & Date.Time > "2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "14" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "14" & Date.Time > "2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "15" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "15" & Date.Time > "2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "16" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "16" & Date.Time > "2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "17" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "17" & Date.Time > "2019-08-17 4:30:00")) %>% 
  filter(!(PoolID == "18" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "18" & Date.Time > "2019-08-17 4:30:00")) %>% 
  filter(!(PoolID == "19" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "19" & Date.Time > "2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "20" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "20" & Date.Time > "2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "21" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "21" & Date.Time > "2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "22" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "22" & Date.Time > "2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "23" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "23" & Date.Time > "2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "24" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "24" & Date.Time > "2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "25" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "25" & Date.Time > "2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "26" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "26" & Date.Time > "2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "27" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "27" & Date.Time > "2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "28" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "28" & Date.Time > "2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "29" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "29" & Date.Time > "2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "30" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "30" & Date.Time > "2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "31" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "31" & Date.Time > "2019-08-17 4:30:00")) %>%
  filter(!(PoolID == "32" & Date.Time < "2019-07-17 10:00:00")) %>%
  filter(!(PoolID == "32" & Date.Time > "2019-08-17 4:30:00"))
#View(Removal.period)

#convert lux to par
#parameters from Long 2012 Field exp values
#A1 = -4924.7
#t1 = 20992.9
#y0 =  4929.0
#PARLICOR = A1 e^(–HOBO/t1) + y0
A1 <- -4924.7
t1 <- 20992.9
y0 <- 4929.0

#covert lux to numeric and get rids of pesky , in lux values
Control.period$Intensity.lux <- as.numeric(gsub(",","",Control.period$Intensity.lux))
x <- Control.period$Intensity.lux #create vector for just lux column
x[x == 0] <- NA #convert 0 to NA
#x[x >= 	19671] <- NA #get rid of values over 	19671 lux/30000PAR (not natural?)
#View(x)

# enter into equation: PARLICOR = A1e(–HOBO/t1) + y
#assign to new column for Par
Control.period$Par<- (A1* (exp( -x / t1))) + y0
#View(Control.period)


#plot(x~Date.Time, data = Control.period, type = "l")

Removal.period$Intensity.lux <- as.numeric(gsub(",","",Removal.period$Intensity.lux))
y <- Removal.period$Intensity.lux
#plot(y~Date.Time, data = Removal.period, type = "l")
y[y == 0] <- NA #convert 0 to NA

Removal.period$Par<- A1* (exp( -y / t1)) + y0
#View(Removal.period)

#plot(y~Date.Time, data = Removal.period, type = "l")



####PAR DATA GRAPHS#####
LightControlDailyMax <- Control.period %>%
  drop_na(Par) %>%
  dplyr::group_by(Date.Time) %>%
  dplyr::group_by(PoolID, Foundation_spp, Removal_Control) %>%
  dplyr::summarise(Par.max = max(Par, na.rm=T),
                   Temp.max = max(Temp.C, na.rm = T)) 

LightControlDailyMax <- as.data.frame(LightControlDailyMax)

#View(LightControlDailyMax)

LightRemovalDailyMax <- Removal.period %>%
  drop_na(Par) %>%
  dplyr::group_by(Date.Time) %>%
  dplyr::group_by(PoolID, Foundation_spp, Removal_Control) %>%
  dplyr::summarise(Par.max = max(Par, na.rm=T),
                   Temp.max = max(Temp.C, na.rm = T))

LightRemovalDailyMax <-as.data.frame(LightRemovalDailyMax)

LightControlDailyMax$Before_After<-"Before"
LightRemovalDailyMax$Before_After<-"After"

LightData<-rbind(LightControlDailyMax,LightRemovalDailyMax)
#View(LightData)

#manipulating temp and light data for mean and variances
Control.periodsum <- Control.period %>%
  dplyr::group_by(PoolID, Foundation_spp, Removal_Control,LoggerDepth) %>%
  filter(Par > 0) %>%
  dplyr::summarise(Temp.mean = mean(Temp.C, na.rm=T), #do this after editing the dates and times
                   Temp.var = var(Temp.C, na.rm=T),
                   Par.mean = mean(Par, na.rm=T))

Control.periodsum<-as.data.frame(Control.periodsum)
Control.periodsum$Before_After<-"Before" #creating column for before since before removal

Removal.periodsum <- Removal.period %>%
  dplyr::group_by(PoolID, Foundation_spp, Removal_Control,LoggerDepth) %>%
  filter(Par > 0) %>%
  dplyr::summarise(Temp.mean = mean(Temp.C, na.rm=T), #do this after editing the dates and times
                   Temp.var = var(Temp.C, na.rm=T),
                   Par.mean = mean(Par, na.rm=T))

Removal.periodsum<-as.data.frame(Removal.periodsum)

Removal.periodsum$Before_After<-"After"#creating column for after since after removal

#View(Removal.periodsum)
TempandLightdata<-rbind(Control.periodsum,Removal.periodsum)
TempandLightdata<-left_join(TempandLightdata, LightData)
#View(TempandLightdata)

TempandLightdata$PoolID <- as.character(TempandLightdata$PoolID)

###Look at differences between daily max, mean, variances of temp and light
DeltaLightandTempdata<- TempandLightdata %>%
  dplyr::group_by(PoolID, Foundation_spp, Removal_Control)  %>%
  dplyr::summarise(DeltaLightMax = Par.max[Before_After == 'After'] - Par.max[Before_After == 'Before'],
                   PercentLightMax = 100 * (Par.max[Before_After == 'After'] - Par.max[Before_After == 'Before']) /  Par.max[Before_After == 'Before'],
                   DeltaTempMax = Temp.max[Before_After == 'After'] - Temp.max[Before_After == 'Before'],
                   PercentTempMax = 100 * (Temp.max[Before_After == 'After'] - Temp.max[Before_After == 'Before']) /  Temp.max[Before_After == 'Before'],
                   deltaTempmean = Temp.mean[Before_After == 'After'] - Temp.mean[Before_After == 'Before'],
                   deltaTempmax = Temp.max[Before_After == 'After'] - Temp.max[Before_After == 'Before'],
                   deltaTempvar = Temp.var[Before_After == 'After'] - Temp.var[Before_After == 'Before'],
                   deltaParmean = Par.mean[Before_After == 'After'] - Par.mean[Before_After == 'Before'])
              
#taking mean because the volume and surface changed in removal pools in after period
#View(DeltaLightandTempdata)
DeltaLightandTempdata <-as.data.frame(DeltaLightandTempdata)

DeltaLightandTempdata<-left_join(DeltaLightandTempdata,Funsppandpp)

#separate datasets for surfgrass and mussels
DeltaLightandTempPhyllo<-DeltaLightandTempdata %>%
  filter(Foundation_spp =="Phyllospadix")

DeltaLightandTempMytilus<-DeltaLightandTempdata %>%
  filter(Foundation_spp =="Mytilus" & PoolID  != '30')

#check collinearity 

#ggpairs(DeltaLightandTempPhyllo[c(28,32:33)]) #good

#ggpairs(DeltaLightandTempMytilus[c(27,32:33)]) #good
#stats for light/temp plots & ggpredict function for regression line
phyllomaxtempmod<-lm(DeltaTempMax ~Phyllodelta+SAVav +THav,data = DeltaLightandTempPhyllo) 
#model of temp with foundation spp loss and tide pool size and tide height as covariates
#plot(phyllomaxtempmod)#good
qqp(resid(phyllomaxtempmod),"norm") #good
summary(phyllomaxtempmod) #sum of mod

pmaxtempgg<-ggpredict(phyllomaxtempmod, c("Phyllodelta")) #predict marginal effects from model for foundation spp. loss
plot(pmaxtempgg) #plot output 
pmaxtemp<-as.data.frame(pmaxtempgg) #create dataframe 

pmaxtemp<-pmaxtemp %>% #output for values gives you an x for variable. rename variable to match
  rename(Phyllodelta=x) #rename to join to rest of dataframe

pmaxtemp<-left_join(pmaxtemp,DeltaLightandTempPhyllo) #rejoin with main dataframe for ggplot

#display raw data but prediction line and confidence intervals are from ggpredict model
phyllotemp<-ggplot(pmaxtemp, aes(x =Phyllodelta, y=DeltaTempMax)) +
  geom_point(size=4,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=Phyllodelta, y=predicted), color="#006d2c",size =1.5)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.3) +
  theme_classic()+
  theme(axis.title.x=element_text(face="italic", color="black", size=26), 
        axis.title.y=element_text(color="black", size=26),
        axis.text.x =element_blank(),
        axis.text.y =element_text(color="black", size=18)) +
  theme(legend.position="none")+
  labs( x= '', y='Change in daily max temp (°C)') 
phyllotemp 

phyllopercentlightmod<-lm(PercentLightMax ~Phyllodelta +SAVav +THav,data = DeltaLightandTempPhyllo)
#plot(phyllopercentlightmod)
qqp(resid(phyllopercentlightmod),"norm")#good
summary(phyllopercentlightmod)

plightgg<-ggpredict(phyllopercentlightmod, c("Phyllodelta"))
plot(plightgg)
plight<-as.data.frame(plightgg)

plight<-plight %>%
  rename(Phyllodelta=x) #rename to join to rest of dataframe

plight<-left_join(plight,DeltaLightandTempPhyllo)

phyllolight<-ggplot(plight, aes(x =Phyllodelta, y=PercentLightMax)) +
  geom_point(size=4,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=Phyllodelta, y=predicted), color="#006d2c",size =1.5)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.3) +
  theme_classic() +
  theme(axis.title.x=element_text(face="italic", color="black", size=26), 
        axis.title.y=element_text(color="black", size=26),
        axis.text.x =element_text(color="black", size=18),
        axis.text.y =element_text(color="black", size=18)) +
  theme(legend.position="none")+ 
  labs(x ='Surfgrass loss \n (Phyllospadix spp.)', y = 'Percent change in daily max light') 
phyllolight

mytilusmaxtempmod<-lm(DeltaTempMax ~Mytilusdelta +SAVav +THav,data = DeltaLightandTempMytilus)
#plot(mytilusmaxtempmod)#good
qqp(resid(mytilusmaxtempmod),"norm")#good
summary(mytilusmaxtempmod)

mtempgg<-ggpredict(mytilusmaxtempmod, c("Mytilusdelta"))
plot(mtempgg)
mmaxtemp<-as.data.frame(mtempgg)

mmaxtemp<-mmaxtemp %>%
  rename(Mytilusdelta=x) #rename to join to rest of dataframe

mmaxtemp<-left_join(mmaxtemp,DeltaLightandTempMytilus)

mytilustemp<-ggplot(mmaxtemp, aes(x =Mytilusdelta, y=DeltaTempMax)) +
  geom_point(size=4,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=Mytilusdelta, y=predicted), color="#045a8d",size =1.5)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.3) +
  theme_classic() +
  theme(axis.title.x =element_blank(), 
        axis.title.y=element_text(color="black", size=26),
        axis.text.x =element_blank(),
        axis.text.y =element_text(color="black", size=18)) +
  ylim(-5,10)+
  theme(legend.position="none")+ 
  xlab('M. californianus loss')+ ylab("")


DeltaLightandTempMytilusNoTP25<-DeltaLightandTempMytilus%>%
  filter(PercentLightMax >= -50)  #filter out pt less than 50 (Tidepool 25)
mytiluspercentlightmod<-lm(PercentLightMax~Mytilusdelta +SAVav +THav,data = DeltaLightandTempMytilusNoTP25)
#plot(mytiluspercentlightmod) #good
qqp(resid(mytiluspercentlightmod),"norm")#good with quadroot light value
summary(mytiluspercentlightmod)

mlightgg<-ggpredict(mytiluspercentlightmod, c("Mytilusdelta"))
plot(mlightgg)
mlight<-as.data.frame(mlightgg)

mlight<-mlight %>%
  rename(Mytilusdelta=x)

mlight<-left_join(mlight,DeltaLightandTempMytilus)

mytiluslight<-ggplot(mlight, aes(x =Mytilusdelta, y=PercentLightMax)) +
  geom_point(size=4,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=Mytilusdelta, y=predicted), color="#045a8d",size =1.5,linetype=2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.3) +
  theme_classic() +
  theme(legend.position="none")+ 
  theme(axis.title.x=element_text(face="italic", color="black", size=26), 
        axis.title.y=element_text(color="black", size=26),
        axis.text.x =element_text(color="black", size=18),
        axis.text.y =element_text(color="black", size=18)) +
  xlab('CA mussel loss \n (Mytilus californianus)')+ ylab('') 
mytiluslight

#patchwork figs together
figure <-phyllotemp + mytilustemp + phyllolight+ mytiluslight +      #patchwork to combine plots
  plot_annotation(tag_levels = 'a') &         #label each individual plot with letters A-G
  theme(plot.tag = element_text(size = 26, face = "bold"))   #edit the lettered text

figure
ggsave(filename = "Output/LightandTempgraphs.pdf", useDingbats =FALSE,dpi=600,device = "pdf", width = 17, height = 19)
