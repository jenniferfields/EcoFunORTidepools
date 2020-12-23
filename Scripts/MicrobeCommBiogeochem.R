#Biogeochem for microbial samples data processing
## By: Jenn Fields
## Last updated: 12.22.2020
###########################
## clear workspace
rm(list=ls())
####load libraries####
library(seacarb) # for carbonate chem
library(oce) # for carbonate chem
library(vegan) #for multi variate analysis
library(devtools)
library(AUC)
library(cowplot)
library(factoextra)
library(plotrix) #for SE
library(MASS) #for stat analysis assumptions
library(car)#for stat analysis assumptions
library(lubridate) #for data/time
library(hms) #date/time
library(ggrepel) #for pcas
library(heplots)
library(zoo) #rolling averages
library(broom) #for tidy function
library(fitdistrplus) #distribution of data
devtools::install_github("seananderson/ggsidekick")
library(GGally) #for ggpairs function
library(ggsidekick) #for theme sleek
library(tidyverse) #for all things %>%

# load data
#source scripts
source("scripts/tidepoolphysicalparameters.R")#load physical parameters
source("scripts/CommunityComp.R") #load community comp 
#Load Data
Nutrients<-read_csv("Data/Biogeochem/RawNutrientData.csv")
CarbChem<-read_csv("Data/Biogeochem/ChemData.csv")
Salinity<-read_csv("Data/Biogeochem/SampleSalinityData_Adjusted.csv")
TA<- read_csv("Data/Biogeochem/TASamples_Adjusted.csv")
Communitymetrics <- read_csv("Data/CommunityComposition/Communitymetrics.csv")


## bring in pH calibration files
pHcalib<-read.csv('Data/Biogeochem/SummerTrispHcalibration.csv')

# First calculate salinity from conductivity
CarbChem$SalCal<-swSCTp(as.numeric(CarbChem$Conductivity)/1000, CarbChem$Temp.pool, pressure=rep(0, each=nrow(CarbChem)),
                        conductivityUnit="mS/cm")

#plot salinity straight from instrument versus calculated from conuctivity and temp
#plot(CarbChem$SalCal, CarbChem$Salinity, xlab='Salinity from conductivity', ylab="Salinity")
#looks very similar...good

#########Cleaning Light and Temp logger data#################
#Clean up time code:
#change date/time data for lubridate package
#first column date and time in same column

CarbChem$Sampling_Day<-mdy(CarbChem$Sampling_Day, quiet=FALSE, tz="America/Los_Angeles", truncated=0)
CarbChem$Date_Time <- paste(CarbChem$Sampling_Day, CarbChem$Sampling_time)
CarbChem$Date_Time<-ymd_hms(CarbChem$Date_Time, quiet=FALSE, tz="America/Los_Angeles", truncated=0) #format date and time

CarbChemAll<-CarbChem


####Microbe data sets (all time points)#####
#remove time 5 from after period since no microbes taken
CarbChemAll<-CarbChemAll %>%
  filter(!(Before_After =="After" & Time_Point == '5')) #filter out time 5 from after period because no microbes

#create data set with just the tp id and start and stop times
Start<- CarbChemAll %>%
  filter(Foundation_spp != 'Ocean') %>%
  dplyr::group_by(PoolID,Before_After,Day_Night) %>%
  dplyr::filter(Time_Point == 1 ) %>%
  dplyr::mutate(StartTime = Date_Time[Time_Point== 1]) %>%
  dplyr::select(PoolID, StartTime)

Stopbefore<- CarbChemAll %>%
  filter(Foundation_spp != 'Ocean') %>%
  filter(Before_After =='Before') %>%
  dplyr::group_by(PoolID,Before_After,Day_Night) %>%
  dplyr::filter(Time_Point == 5) %>%
  dplyr::mutate(StopTime = Date_Time[Time_Point== 5]) %>%
  dplyr::select(PoolID, StopTime)

stopnightbefore<-CarbChemAll%>%
  filter(Foundation_spp != 'Ocean') %>%
  filter(Before_After =='Before' & Day_Night =="Night" & Sampling_Day == "2019-07-12") %>%
  dplyr::group_by(PoolID,Before_After,Day_Night) %>%
  dplyr::filter(Time_Point == 4) %>%
  dplyr::mutate(StopTime = Date_Time[Time_Point== 4]) %>%
  dplyr::select(PoolID, StopTime)

stopnightbefore<-rbind(Stopbefore,stopnightbefore)

Stopafter<-CarbChemAll %>%
  filter(Foundation_spp != 'Ocean') %>%
  filter(Before_After =="After")%>%
  dplyr::group_by(PoolID,Before_After,Day_Night) %>%
  dplyr::filter(Time_Point == 4) %>%
  dplyr::mutate(StopTime = Date_Time[Time_Point== 4]) %>%
  dplyr::select(PoolID, StopTime)

Stopall<-rbind(stopnightbefore,Stopafter)

StartStopall<-left_join(Start,Stopall) #combine both start and stop

temp.data <-
  list.files(path = 'Data/LightandTemp',pattern = ".csv", full.names = TRUE) %>% #all temp and light files in the is folder
  # list files in directory following a particular pattern
  set_names(.) %>% # get the column names
  map_dfr(read.csv, .id = "file.ID") %>% # join all files together in one data frame by file ID
  group_by(file.ID) 

temp.data$Date_Time<- mdy_hm(temp.data$Date.Time, quiet=FALSE, tz="America/Los_Angeles", truncated=0) #format date and time 
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
temp.data$PoolID<-as.character(temp.data$PoolID)

StartandStopall<-left_join(StartStopall,temp.data) #join with temp data

TempandLightSumall<-StartandStopall %>% #create new data frame finding temp and light at start and stop time of water sampling
  dplyr::filter(!(Date_Time >= StopTime | Date_Time <= StartTime)) %>%
  dplyr::group_by(PoolID,Before_After,Day_Night,by60=cut(Date_Time, "60 min")) %>% #hourly means of pools per ~time point
  dplyr::summarise(Temp.max = max(Temp.C, na.rm=T), 
                   Temp.mean = mean(Temp.C, na.rm = T),
                   Par.mean = mean(Par, na.rm = T),
                   Par.max = max(Par,na.rm = T),
                   Par.sum = sum(Par,na.rm = T))

###########Nutrient Data Cleaning##########
#Missing 8 samples from nutritent analysis out of 612 samples
#N_17_5_BC was lost before taking nutrients \
# D_24_1_BC, N_24_1_BC, N2_Ocean_5_BC, D_1_2_AI, D_16_2_AI, D_4_3_AI, N_20_3_AI were missing
#from lab analyses
#We removed time points that nutrients were missing
#missing data was somewhat randomly distributed in dataframe

CarbChemAll<-left_join(CarbChemAll,Nutrients) #join with carbchem data

#Remove rows missing nutrient values from further analysis
CarbChemAll <- CarbChemAll %>%
  filter(Id_code != 'N_17_5_BC' & Id_code != 'D_24_1_BC' & Id_code != 'N_24_1_BC' &
           Id_code != 'N2_Ocean_5_BC' & Id_code != 'D_1_2_AI' & Id_code != 'D_16_2_AI'&
           Id_code != 'D_4_3_AI' & Id_code != 'N_20_3_AI') 

#Conversions from ug/l to µmol/l
#1 µg P/l = 1/MW P = 0.032285 µmol/l
#1 µg/l NO3 = 1/ MW NO3 µg/l = 0.016128 µmol/l
#1 µg/l NO2 = 1/ MW NO2 µg/l = 0.021736 µmol/l
# NO3 and NO2 = 0.016128 + 0.021736
#1 µg/l NH4 = 1/ MW NH4 µg/l = 0.055437 µmol/l

CarbChemAll$PO_ug_L<-as.numeric(CarbChemAll$PO_ug_L)
CarbChemAll$NN_ug_L<-as.numeric(CarbChemAll$NN_ug_L)
CarbChemAll$NH4_ug_L<-as.numeric(CarbChemAll$NH4_ug_L)
#Convert ug/L to umol/L for phosphate, nitrate, ammonium 
CarbChemAll$PO_umol_L <- CarbChemAll$PO_ug_L * 0.032285
CarbChemAll$NN_umol_L <- CarbChemAll$NN_ug_L * (0.016128 + 0.021736)
CarbChemAll$NH4_umol_L <- CarbChemAll$NH4_ug_L * 0.055437

#Converting negative values that were below detection to 0 values
CarbChemAll$PO_umol_L[which(CarbChemAll$PO_umol_L < 0)] <- 0
CarbChemAll$NN_umol_L [which(CarbChemAll$NN_umol_L  < 0)] <- 0
CarbChemAll$NH4_umol_L[which(CarbChemAll$NH4_umol_L < 0)] <- 0

####Adding TA data and lab salinity data####

CarbChemAll<-left_join(CarbChemAll,TA) #join with TA samples
CarbChemAll<-left_join(CarbChemAll,Salinity) #join with salinity in lab

###########CO2 Calculations##############
##correct the TA for the calibration error
CarbChemAll$TA_CRM<-CarbChemAll$TA_Raw-(2207.03*(CarbChemAll$CRM_off/100)) #compared to TA of batch #169 of CRM used for all samples

## take the mV calibration files by each date and use them to calculate pH and pH insitu
CarbChemAll<-pHcalib %>%
  dplyr::nest_by(date_triscal)%>%
  dplyr::mutate(fitpH = list(lm(mVTris~TTris, data = data))) %>% # linear regression of mV and temp of the tris
  dplyr::summarise(broom::tidy(fitpH)) %>% # make the output tidy
  select(date_triscal, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%# put slope and intercept in their own column
  left_join(CarbChemAll,.) %>% # join with the pH sample data
  mutate(mVTris = Temp.in*TTris + `(Intercept)`) %>% # calculate the mV of the tris at temperature in which the pH of samples were measured
  mutate(pH = pH(Ex=mVolts,Etris=mVTris,S=Salinity,T=Temp.in)) %>% # calculate pH of the samples using the pH seacarb function
  mutate(pH_insitu = pHinsi(pH = pH, ALK = TA_CRM/1000000, Tinsi = Temp.pool, Tlab = Temp.in, S=Salinity, 
                            Pt=PO_umol_L/1000000, Sit=0, k1k2 = "x",  kf = "x", ks="d", pHscale = "T", b="u74")) 

###########DO regression#################
# The incorrect multiparameter sensor was used during Time point 1 of Day TP 1-16 in After period
# We were missing Dissolved Oxygen (DO mg/L and %) values for Ocean, TP 1-15
# We ran a regression of Dissolved Oxygen with pH and Temperature from the Before period to replace DO values
# Due to the strong correlation between DO, Temp, and pH
# We used the regression values to back calculate DO for the missing tide pools
# For the Ocean sample, we took the global average of the DO values of that day from time points 2-5 to replace
# the ocean value from time point 1

#Time point 1 for days (TP 1-32) from the Before period without ocean sample
DayTime1Chem<- CarbChemAll %>%
  filter(PoolID != "Ocean") %>%
  filter(Day_Night == "Day") %>%
  filter(Time_Point == "1") %>%
  filter(Before_After == "Before")
#View(DayTime1Chem)

DOTime1regression<-lm(DO_mg_L~pH_insitu + Temp.pool, data = DayTime1Chem) #DO mg/L,pH, and temp for time point 1 of days
summary(DOTime1regression)
#r2 of 0.9083
#B0 = -36.4751
#BpH =  5.4970  
#B temp  = 0.1504
DOpercentTime1regression<-lm(DO_p~pH_insitu + Temp.pool, data = DayTime1Chem) #%DO,pH, and temp for time point 1 of days
summary(DOpercentTime1regression)
#r2 of 0.9039
#B0 = -470.450
# BpH = 66.208 
# BTemp = 3.935

resid(DOTime1regression) #better to name these as something
DOregressionres<-resid(DOTime1regression)
#Now we can test the normality of the residuals
qqp(DOregressionres, "norm")
#normality
#We can also call the fitted y values as:
fitted(DOTime1regression)

#To test for homogeneity of variance, we want to plot the fitted (predicted) values
#against the residuals
plot(DOregressionres~fitted(DOTime1regression))
#good so can use model for DO values

#plot(DO~pH + Temperature, data=DayTime1Chem)
#abline(DOTime1regression, col="blue")

#select from dataset tide pools 1-15 for time point one on 08052019 sampling day during after period since missing DO values
tpday1AI<-which(CarbChemAll$Time_Point == "1" & CarbChemAll$Before_After == "After" & CarbChemAll$Sampling_Day == "2019-08-05" & 
                  CarbChemAll$PoolID != "Ocean" & CarbChemAll$PoolID != "16")
#View(tpday1AI)

#DO mg/L regression equation
#r2 of 0.9083
#B0 = -36.4751
#BpH =  5.4970  
#B temp  = 0.1504
CarbChemAll$DO_mg_L[tpday1AI]<- -36.4751 + (5.4970 * CarbChemAll$pH_insitu[tpday1AI]) + (0.1504 * CarbChemAll$Temp.pool[tpday1AI])
#View(CarbChem)

#DO % regression equation
#r2 of 0.9039
#B0 = -470.450
# BpH = 66.208 
# BTemp = 3.935
CarbChemAll$DO_p[tpday1AI]<- -470.450 + (66.208  *CarbChemAll$pH_insitu[tpday1AI]) + (3.935 * CarbChemAll$Temp.pool[tpday1AI])

#Now for ocean sample--select from data set time one on 0805209 sampling day
Oceanday1AI<-which(CarbChemAll$Time_Point == "1" & CarbChemAll$Before_After == "After" & CarbChemAll$Sampling_Day == "2019-08-05" & 
                     CarbChemAll$PoolID == "Ocean")
#View(Oceanday1AI)
CarbChemAll$DO_mg_L[Oceanday1AI]<- -36.4751 + (5.4970 * CarbChemAll$pH_insitu[Oceanday1AI]) + (0.1504 * CarbChemAll$Temp.pool[Oceanday1AI])
CarbChemAll$DO_p[Oceanday1AI]<- -470.450 + (66.208  * CarbChemAll$pH_insitu[Oceanday1AI]) + (3.935 * CarbChemAll$Temp.pool[Oceanday1AI])

########NEC and NEP rates##########
#calculate all the CO2Sys params
PoolCO2all<-carb(flag=8, CarbChemAll$pH_insitu, CarbChemAll$TA_CRM/1000000, S=CarbChemAll$Salinity, T=CarbChemAll$Temp.pool, Patm=1, P=0, 
                 Pt=CarbChemAll$PO_umol_L/1000000, Sit=0,
                 k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential")
#TA is divided by 1000 because all calculations are in mol/kg in the seacarb package

# calculate error propagation
er<-errors(flag=8, CarbChemAll$pH_insitu, CarbChemAll$TA_CRM/1000000, 
           S=CarbChemAll$Salinity, T=CarbChemAll$Temp.pool, 
           Patm=1, P=0,Pt=CarbChemAll$PO_umol_L/1000000,
           Sit=0,evar1 = 0.01, evar2 = 5e-6) 

#average error for DIC based on pH and TA
mean(er$DIC*1000000)
#7.091473
sd(er$DIC*1000000)/sqrt(nrow(er))
#0.05899179

#convert CO2, HCO3, CO3, DIC, and Alk back to micromol for easier interpretation
PoolCO2all[,c("CO2","HCO3","CO3","DIC","ALK")]<-PoolCO2all[,c("CO2","HCO3","CO3","DIC","ALK")]*1000000

#combine all the data
CarbChemAll<-cbind(CarbChemAll,PoolCO2all[-c(1:6,10:13)])

#Normalize the TA and DIC to salinity

#TA normalized to constant salinity, but not nutrients here
CarbChemAll$TA_NormSal<-CarbChemAll$TA_CRM*(CarbChemAll$Salinity_Lab/34)

#Normalize DIC constant salinity
CarbChemAll$DIC_Norm<-CarbChemAll$DIC*(CarbChemAll$Salinity_Lab/34)

MicrobeSamplingchem<-CarbChemAll%>%
  filter(Time_Point == 1 | Time_Point == 4 | Time_Point == 5)
write.csv(MicrobeSamplingchem,file="Output/MicrobeCarbChem.csv")