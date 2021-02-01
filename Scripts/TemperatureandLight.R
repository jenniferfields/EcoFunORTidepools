##Temperature and Light one month analysis in OR tide pools
##By: Jenn Fields
##Last updated 2.01.2021
###################

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
source("scripts/tidepoolphysicalparameters.R")
source("scripts/CommunityComp.R")
#load ocean temp data from ODFW Marine Reserves Mooring
SMURFOceanTemp<-read_csv("Data/LightandTemp/Fields_OR_1mSMURFTemp_2019.csv")

#########Cleaning Light and Temp logger data#################
#options(na.action = "na.omit") 
#group temp data from all Hobo pendant loggers together
temp.data <-
  list.files(path = 'Data/LightandTemp',pattern = ".csv", full.names = TRUE) %>% 
  # list files in directory following a particular pattern
  set_names(.) %>% # get the column names
  map_dfr(read.csv, .id = "file.ID") %>% # join all files together in one data frame by file ID
  group_by(file.ID) 
#assign column to data time
temp.data$Date.Time<- mdy_hm(temp.data$Date.Time, quiet=FALSE, tz="America/Los_Angeles", truncated=0) #format date and time 
SMURFOceanTemp$Date.Time<-mdy_hm(SMURFOceanTemp$Date.Time,quiet=FALSE, tz="America/Los_Angeles", truncated=0) #format date and time

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

#Control period data before removal of foundation spp n=29 days
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


#Temp data after removal of foundation spp n= 29 days
Removal.period<-temp.data %>%
  filter(!(PoolID == "1" & Date.Time < "2019-07-18 0:00:00")) %>%
  filter(!(PoolID == "1" & Date.Time > "2019-08-16 0:00:00")) %>%
  filter(!(PoolID == "2" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "2" & Date.Time > "2019-08-16 00:00:00")) %>%
  filter(!(PoolID == "3" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "3" & Date.Time > "2019-08-16 00:00:00")) %>%
  filter(!(PoolID == "4" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "4" & Date.Time > "2019-08-16 00:00:00")) %>%
  filter(!(PoolID == "5" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "5" & Date.Time > "2019-08-16 00:00:00")) %>%
  filter(!(PoolID == "6" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "6" & Date.Time > "2019-08-16 00:00:00")) %>%
  filter(!(PoolID == "7" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "7" & Date.Time > "2019-08-16 00:00:00")) %>%
  filter(!(PoolID == "8" & Date.Time < "2019-07-18 10:00:00")) %>%
  filter(!(PoolID == "8" & Date.Time > "2019-08-16 00:00:00")) %>%
  filter(!(PoolID == "9" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "9" & Date.Time > "2019-08-16 00:00:000")) %>%
  filter(!(PoolID == "10" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "10" & Date.Time > "2019-08-16 00:00:00")) %>%
  filter(!(PoolID == "11" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "11" & Date.Time > "2019-08-16 00:00:00")) %>%
  filter(!(PoolID == "12" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "12" & Date.Time > "2019-08-16 00:00:00")) %>%
  filter(!(PoolID == "13" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "13" & Date.Time > "2019-08-16 00:00:00")) %>%
  filter(!(PoolID == "14" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "14" & Date.Time > "2019-08-16 00:00:00")) %>%
  filter(!(PoolID == "15" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "15" & Date.Time > "2019-08-16 00:00:00")) %>%
  filter(!(PoolID == "16" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "16" & Date.Time > "2019-08-16 00:00:00")) %>%
  filter(!(PoolID == "17" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "17" & Date.Time > "2019-08-16 00:00:00")) %>% 
  filter(!(PoolID == "18" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "18" & Date.Time > "2019-08-16 00:00:00")) %>% 
  filter(!(PoolID == "19" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "19" & Date.Time > "2019-08-16 00:00:00")) %>%
  filter(!(PoolID == "20" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "20" & Date.Time > "2019-08-16 00:00:00")) %>%
  filter(!(PoolID == "21" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "21" & Date.Time > "2019-08-16 00:00:00")) %>%
  filter(!(PoolID == "22" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "22" & Date.Time > "2019-08-16 00:00:00")) %>%
  filter(!(PoolID == "23" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "23" & Date.Time > "2019-08-16 00:00:00")) %>%
  filter(!(PoolID == "24" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "24" & Date.Time > "2019-08-16 00:00:00")) %>%
  filter(!(PoolID == "25" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "25" & Date.Time > "2019-08-16 00:00:00")) %>%
  filter(!(PoolID == "26" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "26" & Date.Time > "2019-08-16 00:00:00")) %>%
  filter(!(PoolID == "27" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "27" & Date.Time > "2019-08-16 00:00:00")) %>%
  filter(!(PoolID == "28" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "28" & Date.Time > "2019-08-16 00:00:00")) %>%
  filter(!(PoolID == "29" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "29" & Date.Time > "2019-08-16 00:00:00")) %>%
  filter(!(PoolID == "30" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "30" & Date.Time > "2019-08-16 00:00:00")) %>%
  filter(!(PoolID == "31" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "31" & Date.Time > "2019-08-16 00:00:00")) %>%
  filter(!(PoolID == "32" & Date.Time < "2019-07-18 00:00:00")) %>%
  filter(!(PoolID == "32" & Date.Time > "2019-08-16 00:00:00"))


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

# enter into equation: PARLICOR = A1e(–HOBO/t1) + y
#assign to new column for Par
Control.period$Par<- (A1* (exp( -x / t1))) + y0

Removal.period$Intensity.lux <- as.numeric(gsub(",","",Removal.period$Intensity.lux))
y <- Removal.period$Intensity.lux
y[y == 0] <- NA #convert 0 to NA

Removal.period$Par<- A1* (exp( -y / t1)) + y0

####PAR DATA GRAPHS#####
Control.period$Day<-as.factor(as.Date(Control.period$Date.Time, quiet=FALSE, tz="America/Los_Angeles", truncated=0))

ControlDailyMax <- Control.period %>%
  drop_na(Par) %>%#remove NA light values 
  dplyr::group_by(PoolID, Foundation_spp, Removal_Control,Day) %>% #cut by day to get daily max 
  dplyr::summarise(Par.max = max(Par, na.rm=T),
                   Temp.max=max(Temp.C,na.rm=T),Par.mean=mean(Par,na.rm = T))

ControlDailyMax <- as.data.frame(ControlDailyMax)

#View(LightControlDailyMax)
Removal.period$Day<-as.factor(as.Date(Removal.period$Date.Time, quiet=FALSE, tz="America/Los_Angeles", truncated=0))

RemovalDailyMax <- Removal.period %>%
  drop_na(Par) %>% #remove NA light values
  dplyr::group_by(PoolID, Foundation_spp, Removal_Control,Day) %>%
  dplyr::summarise(Par.max = max(Par, na.rm=T),
                   Temp.max=max(Temp.C,na.rm=T),
                   Par.mean=mean(Par,na.rm = T))

RemovalDailyMax <-as.data.frame(RemovalDailyMax)

ControlDailyMax$Before_After<-"Before"
RemovalDailyMax$Before_After<-"After"

MaxTempLightData<-rbind(ControlDailyMax,RemovalDailyMax)

###Look at differences between daily max, mean, variances of temp and light####
DeltaLightandTempdata<- MaxTempLightData %>%
  dplyr::group_by(PoolID, Foundation_spp, Removal_Control,Before_After)  %>%
  dplyr::summarise(Par.maxav = mean(Par.max),Temp.maxav=mean(Temp.max),Parmeanav= mean(Par.mean)) %>%
  dplyr::group_by(PoolID, Foundation_spp, Removal_Control) %>%
  dplyr::summarise(DeltaLightMax = Par.maxav[Before_After == 'After'] - Par.maxav[Before_After == 'Before'],
                   PercentLightMax = 100 * ((Par.maxav[Before_After == 'After'] - Par.maxav[Before_After == 'Before']) /  Par.maxav[Before_After == 'Before']),
                   Percentlightmean =100* ((Parmeanav[Before_After == 'After'] - Parmeanav[Before_After == 'Before']) /  Parmeanav[Before_After == 'Before']),
                   DeltaTempMax = Temp.maxav[Before_After == 'After'] - Temp.maxav[Before_After == 'Before'],
                   PercentTempMax = 100 * (Temp.maxav[Before_After == 'After'] - Temp.maxav[Before_After == 'Before']) /  Temp.maxav[Before_After == 'Before'])

Changesummary<-DeltaLightandTempdata %>%
  group_by(Foundation_spp,Removal_Control)%>%
  summarise(PLMaxmean=mean(PercentLightMax),PLMaxse= std.error(PercentLightMax),
            Maxtempmean=mean(DeltaTempMax),Maxtempse=std.error(DeltaTempMax)) 

#taking mean because the volume and surface changed in removal pools in after period
#View(DeltaLightandTempdata)
DeltaLightandTempdata <-as.data.frame(DeltaLightandTempdata)
DeltaLightandTempdata$PoolID <- as.factor(DeltaLightandTempdata$PoolID)

DeltaLightandTempdata<-left_join(DeltaLightandTempdata,Funsppandpp)

####temp and light analyses####
#separate datasets for surfgrass and mussels
DeltaLightandTempPhyllo<-DeltaLightandTempdata %>%
  filter(Foundation_spp =="Phyllospadix")

DeltaLightandTempMytilus<-DeltaLightandTempdata %>%
  filter(Foundation_spp =="Mytilus" & PoolID  != '30')

#check collinearity between f spp. loss and size of pool and tide height
#use volume and tide height
#ggpairs(DeltaLightandTempPhyllo[c(10:14)]) 
#ggpairs(DeltaLightandTempMytilus[c(9:15)]) 

#Phyllospadix light and temp models and ggplots
#stats for light/temp plots & ggpredict function for regression line
DeltaLightandTempPhyllo$logtemp<-sign(DeltaLightandTempPhyllo$DeltaTempMax)*log(abs(DeltaLightandTempPhyllo$DeltaTempMax))

phyllomaxtempmod<-lm(logtemp~Phyllodelta+Vav +THav,data = DeltaLightandTempPhyllo) 
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
pmaxtemp<-pmaxtemp%>%
  mutate(BTpredict=exp(predicted),BTcl=exp(conf.low),BTch=exp(conf.high),BTtemp=exp(logtemp))

#display raw data but prediction line and confidence intervals are from ggpredict model
phyllotemp<-ggplot(pmaxtemp, aes(x =Phyllodelta, y=logtemp))+ #aes(y = ) using the non-logged data
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=Phyllodelta, y=predicted), color="#006d2c",size =2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  #coord_trans(y='log') +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=45), 
        axis.title.y=element_text(color="black", size=45),
        axis.text.x =element_text(color="black", size=35),
        axis.text.y =element_text(color="black", size=35)) +
  theme(legend.position="none")+
  #scale_x_continuous(limits = c(-28, 100), breaks = seq(-30, 100, by = 25))+
  #ylim(-2,6)+
  labs(x ='', y = 'Log change in average daily max temperature (°C)') 
phyllotemp 


DeltaLightandTempPhyllo$loglight<-sign(DeltaLightandTempPhyllo$PercentLightMax)*log(abs(DeltaLightandTempPhyllo$PercentLightMax))

phyllopercentlightmod<-lm(loglight~Phyllodelta +Vav +THav,data = DeltaLightandTempPhyllo)

#plot(phyllopercentlightmod)
qqp(resid(phyllopercentlightmod),"norm")#good
#plot(phyllopercentlightmod)
summary(phyllopercentlightmod)

plightgg<-ggpredict(phyllopercentlightmod, c("Phyllodelta"))
plot(plightgg)
plight<-as.data.frame(plightgg)

plight<-plight %>%
  rename(Phyllodelta=x) #rename to join to rest of dataframe

plight<-left_join(plight,DeltaLightandTempPhyllo)

phyllolight<-ggplot(plight, aes(x =Phyllodelta, y=loglight))+
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=Phyllodelta, y=predicted), color="#006d2c",size =2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic() +
  theme(axis.title.x=element_text(color="black", size=45), 
        axis.title.y=element_text(color="black", size=45),
        axis.text.x =element_text(color="black", size=35),
        axis.text.y =element_text(color="black", size=35)) +
  theme(legend.position="none")+ 
  #scale_x_continuous(limits = c(-30, 100), breaks = seq(-30, 100, by = 25))+
  #coord_trans(y="log")+
  #scale_y_continuous(breaks = scales::log_breaks())+
  labs(x ='Surfgrass percent loss \n (Phyllospadix spp.)', y = 'Log change in average daily percent max light') 
phyllolight

#Mytilus light and temp models and ggplots
mytilusmaxtempmod<-lm(DeltaTempMax ~Mytilusdelta +Vav +THav,data = DeltaLightandTempMytilus)
#plot(mytilusmaxtempmod)
qqp(resid(mytilusmaxtempmod),"norm")#good
#plot(mytilusmaxtempmod)
summary(mytilusmaxtempmod)

mtempgg<-ggpredict(mytilusmaxtempmod, c("Mytilusdelta"))
plot(mtempgg)
mmaxtemp<-as.data.frame(mtempgg)

mmaxtemp<-mmaxtemp %>%
  rename(Mytilusdelta=x) #rename to join to rest of dataframe

mmaxtemp<-left_join(mmaxtemp,DeltaLightandTempMytilus)

mytilustemp<-ggplot(mmaxtemp, aes(x =Mytilusdelta, y=DeltaTempMax)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=Mytilusdelta, y=predicted), color="#045a8d",size =2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic() +
  theme(axis.title.x =element_text(color="black", size=45),
        axis.title.y=element_text(color="black", size=45),
        axis.text.x =element_text(color="black", size=35),
        axis.text.y =element_text(color="black", size=35)) +
  ylim(-5,10)+
  theme(legend.position="none") + 
  labs(x='', y="Change in average daily max temperature (°C)")
mytilustemp


DeltaLightandTempMytilus$loglight<-sign(DeltaLightandTempMytilus$PercentLightMax)*log(abs(DeltaLightandTempMytilus$PercentLightMax))
mytiluspercentlightmod<-lm(loglight~Mytilusdelta +Vav +THav,data =DeltaLightandTempMytilus)
#plot(mytiluspercentlightmod) #good
qqp(resid(mytiluspercentlightmod),"norm")#good 
summary(mytiluspercentlightmod)

mlightgg<-ggpredict(mytiluspercentlightmod, c("Mytilusdelta"))
plot(mlightgg)
mlight<-as.data.frame(mlightgg)

mlight<-mlight %>%
  rename(Mytilusdelta=x)

mlight<-left_join(mlight,DeltaLightandTempMytilus)

mytiluslight<-ggplot(mlight, aes(x =Mytilusdelta, y=loglight)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=Mytilusdelta, y=predicted), color="#045a8d",size =2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic() +
  theme(legend.position="none")+ 
  theme(axis.title.x=element_text(color="black", size=45), 
        axis.title.y=element_text(color="black", size=45),
        axis.text.x =element_text(color="black", size=35),
        axis.text.y =element_text(color="black", size=35)) +
  xlab('CA mussel percent loss \n (Mytilus californianus)')+ ylab('') 
mytiluslight

#patchwork figs together
figure <-(phyllotemp | mytilustemp) /
  (phyllolight| mytiluslight) +      #patchwork to combine plots
  plot_annotation(tag_levels = 'a') &         #label each individual plot with letters A-G
  theme(plot.tag = element_text(size =40))   #edit the lettered text
figure
#ggsave(filename = "Output/LightandTempgraphsavg.pdf", useDingbats =FALSE,dpi=600,device = "pdf", width = 35, height = 40)

####Supplemental figure--Time series temperature compared to nearshore ocean####
TemptimeseriesC<-Control.period %>%
  dplyr::group_by(PoolID, Foundation_spp, Removal_Control)
TemptimeseriesC$Before_After<-"Before"

TemptimeseriesR<-Removal.period %>%
  dplyr::group_by(PoolID, Foundation_spp, Removal_Control)
TemptimeseriesR$Before_After<-"After"
Temptimeseries<-rbind(TemptimeseriesC,TemptimeseriesR)

Temptimeseries$PoolID<-as.character(Temptimeseries$PoolID)

Temptimeseries<-left_join(Temptimeseries,Funsppandpp)
#filter for time siries graph
SMURFOceanTemp<-SMURFOceanTemp%>%
  filter(Date.Time < "2019-08-16 00:00:00")

#summarizing ocean day by before and after periods
SumOceanB<-SMURFOceanTemp%>%
  dplyr::filter(Date.Time >"2019-06-15 0:00:00" & Date.Time <"2019-07-15 00:00:00")

SumOceanA<-SMURFOceanTemp%>%
  dplyr::filter( Date.Time<"2019-08-16 00:00:00" & Date.Time >"2019-07-18 0:00:00")
                               
SumOceanB$Before_After<-"Before"
SumOceanA$Before_After<-"After" 

Oceantempseries<-rbind(SumOceanB,SumOceanA)
#make day a factor for column
Oceantempseries$Day<-as.factor(as.Date(Oceantempseries$Date.Time, quiet=FALSE, tz="America/Los_Angeles", truncated=0))

Oceansum<-Oceantempseries%>%
  group_by(Day,Before_After) %>%
  summarize(meantemp = mean(OceanTemp.C)) #gets average daily temp

oceansumB<-Oceansum%>%
  filter(Before_After=='Before')

oceansumA<-Oceansum%>%
  filter(Before_After=='After')

oceansumB$row_num <- seq.int(nrow(oceansumB)) #makes column for day #s
oceansumA$row_num <- seq.int(nrow(oceansumA)) 

Oceantempsum<-rbind(oceansumB,oceansumA)

#subtract each daily av between after and before to get daily change between time
Oceantempsum<-Oceantempsum%>%
  dplyr::group_by(row_num) %>% 
  dplyr::summarise(Deltadaily=meantemp[Before_After=='After']-meantemp[Before_After=='Before']) %>%
  dplyr::summarise(meandelta=mean(Deltadaily),se=std.error(Deltadaily))

#summarize ave daily for each time period
OceansumBandA<-Oceantempseries%>%
  group_by(Day,Before_After) %>%
  summarize(meantemp = mean(OceanTemp.C))%>%
  group_by(Before_After)%>%
  summarize(tempBA=mean(meantemp),se=std.error(meantemp))

#Ocean mean:
#before
#11.34917 +/-0.4054801
#after
#10.79524 +/-	0.2405178
#average daily avg change between periods: -0.5494291 +/- 0.4315244

Phyllotimeseries<-Temptimeseries%>%
  filter(Foundation_spp =="Phyllospadix") 

#to calculate daily max values
#Phyllotimeseriesmean<-Temptimeseries%>%
# filter(Foundation_spp =="Phyllospadix") %>%
# dplyr::group_by(PoolID,Removal_Control,Before_After,by60=cut(Date.Time, "1440 min")) %>% #daily max per tide pool
#dplyr::summarise(meantempbypool= mean(Temp.C, na.rm = TRUE)) %>%
#dplyr::group_by(Removal_Control, Before_After,by60) %>% #daily mean max across all tide pools
#dplyr::summarise(meantemp=mean(meantempbypool, na.rm=TRUE))
#Phyllotimeseriesmean$by60<-as_datetime(Phyllotimeseriesmean$by60)

Phyllotimeseries$PoolID<-as.factor(Phyllotimeseries$PoolID)

Temptimeseries$Phyllodelta[Temptimeseries$Before_After=="Before"]<-0 #change phyllo loss before period to 0

#Phyllotimeseries$Day<-as.factor(as.Date(Phyllotimeseries$by60, quiet=FALSE, tz="America/Los_Angeles", truncated=0))
Phyllotempseries<-Temptimeseries %>%
  filter(Foundation_spp =="Phyllospadix") %>%
  ggplot(aes(x=Date.Time,y=Temp.C,colour=Phyllodelta)) +
  geom_line(size=2)+
  scale_colour_distiller(palette = "Greys",direction =1)+
  geom_vline(xintercept = as.POSIXct("2019-07-17"), linetype=2, 
             color = "black", size=3)+
  geom_line(data=SMURFOceanTemp,aes(x=Date.Time,y=OceanTemp.C),colour="#a50f15",size=2)+
  #geom_line(data=Phyllotimeseriesmean,aes(x=by60,y=meantemp),colour="#08519c",size=2)+
  theme_classic()+
  facet_wrap(~Removal_Control) +
  theme(axis.text.x=element_text(size = 30, color = "black"), 
        axis.text.y=element_text(size = 30, color = "black"),
        axis.title.y = element_text(size = 40, color = "black"),
        axis.title.x = element_text(size = 40, color = "black"),
        legend.text =element_text(size = 30, color = "black"),
        legend.title = element_text(size = 35, color = "black"),
        strip.text.x = element_text(size=40, color="black"), #change facet labels
        panel.spacing = unit(4, "lines"), #change facet spacing
        legend.key.size = unit(1.5, "cm"),
        legend.key.width = unit(1.0,"cm")) +
  scale_x_datetime(limits = as.POSIXct(c("2019-06-16 00:00:00", "2019-08-16 00:00:00 ")))+
  annotate(geom="text", y=34, x=as.POSIXct("2019-07-01"), label="Before",
           color="black",size=15)+
  annotate(geom="text", y=34, x=as.POSIXct("2019-08-02"), label="After",
           color="black",size=15)+
  labs(y="Temperature (°C)", x="Date", color ="Surfgrass Loss")
Phyllotempseries

#ggsave(filename = "Output/Surfgrasstemp.pdf", useDingbats =FALSE,dpi=600,device = "pdf", width = 35, height = 10)
Temptimeseries$Mytilusdelta[Temptimeseries$Before_After=="Before"]<-0
Musseltempseries<-Temptimeseries %>%
  filter(Foundation_spp =="Mytilus") %>%
  ggplot(aes(x=Date.Time,y=Temp.C,color=Mytilusdelta)) +
  geom_line(size =2)+
  scale_colour_distiller(palette = "Greys",direction =1)+
  geom_vline(xintercept = as.POSIXct("2019-07-17"), linetype=2, 
             color = "black", size=3)+
  geom_line(data=SMURFOceanTemp,aes(x=Date.Time,y=OceanTemp.C),colour="#a50f15",size=2)+
  theme_classic()+
  facet_wrap(~Removal_Control) +
  theme(axis.text.x=element_text(size = 30, color = "black"), 
        axis.text.y=element_text(size = 30, color = "black"),
        axis.title.y = element_text(size = 40, color = "black"),
        axis.title.x = element_text(size = 40, color = "black"),
        legend.text =element_text(size = 30, color = "black"),
        legend.title = element_text(size = 35, color = "black"),
        strip.text.x = element_text(size=40, color="black"), #change facet labels
        panel.spacing = unit(4, "lines"), #change facet spacing
        legend.key.size = unit(1.5, "cm"),
        legend.key.width = unit(1.0,"cm")) +
  scale_x_datetime(limits = as.POSIXct(c("2019-06-16 00:00:00", "2019-08-16 00:00:00 ")))+
  annotate(geom="text", y=34, x=as.POSIXct("2019-07-01"), label="Before",
           color="black",size=15)+
  annotate(geom="text", y=34, x=as.POSIXct("2019-08-02"), label="After",
           color="black",size=15)+
  labs(y="Temperature (°C)", x="Date", color ="Mussel Loss")
Musseltempseries
#ggsave(filename = "Output/musseltemp.pdf", useDingbats =FALSE,dpi=600,device = "pdf", width = 25, height = 20)
Temptimeseriespm<-(Phyllotempseries/Musseltempseries)+
  plot_annotation(tag_levels = 'a') &         #label each individual plot with letters A-G
  theme(plot.tag = element_text(size =50))   #edit the lettered text
Temptimeseriespm
#ggsave(filename = "Output/temptimeseries.pdf", useDingbats =FALSE,dpi=600,device = "pdf", width = 25, height = 25)
