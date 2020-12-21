#C#lean biogeochem and carb chem script
## By: Jenn Fields
## Last updated: 10.23.2020
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

#Calculate delta TA between Time Points
CarbChem<-CarbChem %>%
  filter(Time_Point != '5') %>% #take out time point 5 from further analysis
  filter(PoolID != '30') #removing pool 30 because nutrients values in After period during 
#day were orders of magnitude higher than others (ex:Nh4 avg 4540 ug/L, PO avg 324ug/L, NO avg 670ug/L)
#Something funky happened.

#create data set with just the tp id and start and stop times
Start<- CarbChem %>%
  filter(Foundation_spp != 'Ocean') %>%
  dplyr::group_by(PoolID,Before_After,Day_Night) %>%
  dplyr::filter(Time_Point == 1 ) %>%
  dplyr::mutate(StartTime = Date_Time[Time_Point== 1]) %>%
  dplyr::select(PoolID, StartTime)

Stop<- CarbChem %>%
  filter(Foundation_spp != 'Ocean') %>%
  dplyr::group_by(PoolID,Before_After,Day_Night) %>%
  dplyr::filter(Time_Point == 4) %>%
  dplyr::mutate(StopTime = Date_Time[Time_Point== 4]) %>%
  dplyr::select(PoolID, StopTime)

StartandStop<-left_join(Start,Stop)

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

StartandStop<-left_join(StartandStop,temp.data) #join with temp data

TempandLightSum<-StartandStop %>% #create new data frame finding temp and light at start and stop time of water sampling
  dplyr::filter(!(Date_Time >= StopTime | Date_Time <= StartTime)) %>%
  dplyr::group_by(PoolID,Before_After,Day_Night)%>% 
  dplyr::summarise(Temp.max = max(Temp.C, na.rm=T), 
                   Temp.mean = mean(Temp.C, na.rm = T),
                   Par.mean = mean(Par, na.rm = T),
                   Par.max = max(Par,na.rm = T),
                   Par.sum = sum(Par,na.rm = T))

TempandLightSumall<-StartandStop %>% #create new data frame finding temp and light at start and stop time of water sampling
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

CarbChem<-left_join(CarbChem,Nutrients) #join with carbchem data

#Remove rows missing nutrient values from further analysis
CarbChem <- CarbChem %>%
  filter(Id_code != 'N_17_5_BC' & Id_code != 'D_24_1_BC' & Id_code != 'N_24_1_BC' &
           Id_code != 'N2_Ocean_5_BC' & Id_code != 'D_1_2_AI' & Id_code != 'D_16_2_AI'&
           Id_code != 'D_4_3_AI' & Id_code != 'N_20_3_AI') 

#Conversions from ug/l to µmol/l
#1 µg P/l = 1/MW P = 0.032285 µmol/l
#1 µg/l NO3 = 1/ MW NO3 µg/l = 0.016128 µmol/l
#1 µg/l NO2 = 1/ MW NO2 µg/l = 0.021736 µmol/l
# NO3 and NO2 = 0.016128 + 0.021736
#1 µg/l NH4 = 1/ MW NH4 µg/l = 0.055437 µmol/l

CarbChem$PO_ug_L<-as.numeric(CarbChem$PO_ug_L)
CarbChem$NN_ug_L<-as.numeric(CarbChem$NN_ug_L)
CarbChem$NH4_ug_L<-as.numeric(CarbChem$NH4_ug_L)
#Convert ug/L to umol/L for phosphate, nitrate, ammonium 
CarbChem$PO_umol_L <- CarbChem$PO_ug_L * 0.032285
CarbChem$NN_umol_L <- CarbChem$NN_ug_L * (0.016128 + 0.021736)
CarbChem$NH4_umol_L <- CarbChem$NH4_ug_L * 0.055437

#Converting negative values that were below detection to 0 values
CarbChem$PO_umol_L[which(CarbChem$PO_umol_L < 0)] <- 0
CarbChem$NN_umol_L [which(CarbChem$NN_umol_L  < 0)] <- 0
CarbChem$NH4_umol_L[which(CarbChem$NH4_umol_L < 0)] <- 0

####Adding TA data and lab salinity data####

CarbChem<-left_join(CarbChem,TA) #join with TA samples
CarbChem<-left_join(CarbChem,Salinity) #join with salinity in lab

###########CO2 Calculations##############
##correct the TA for the calibration error
CarbChem$TA_CRM<-CarbChem$TA_Raw-(2207.03*(CarbChem$CRM_off/100)) #compared to TA of batch #169 of CRM used for all samples

## take the mV calibration files by each date and use them to calculate pH and pH insitu
CarbChem<-pHcalib %>%
  dplyr::nest_by(date_triscal)%>%
  dplyr::mutate(fitpH = list(lm(mVTris~TTris, data = data))) %>% # linear regression of mV and temp of the tris
  dplyr::summarise(broom::tidy(fitpH)) %>% # make the output tidy
  select(date_triscal, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%# put slope and intercept in their own column
  left_join(CarbChem,.) %>% # join with the pH sample data
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
DayTime1Chem<- CarbChem %>%
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
tpday1AI<-which(CarbChem$Time_Point == "1" & CarbChem$Before_After == "After" & CarbChem$Sampling_Day == "2019-08-05" & 
                  CarbChem$PoolID != "Ocean" & CarbChem$PoolID != "16")
#View(tpday1AI)

#DO mg/L regression equation
#r2 of 0.9083
#B0 = -36.4751
#BpH =  5.4970  
#B temp  = 0.1504
CarbChem$DO_mg_L[tpday1AI]<- -36.4751 + (5.4970 * CarbChem$pH_insitu[tpday1AI]) + (0.1504 * CarbChem$Temp.pool[tpday1AI])
#View(CarbChem)

#DO % regression equation
#r2 of 0.9039
#B0 = -470.450
# BpH = 66.208 
# BTemp = 3.935
CarbChem$DO_p[tpday1AI]<- -470.450 + (66.208  * CarbChem$pH_insitu[tpday1AI]) + (3.935 * CarbChem$Temp.pool[tpday1AI])

#Now for ocean sample--select from data set time one on 0805209 sampling day
Oceanday1AI<-which(CarbChem$Time_Point == "1" & CarbChem$Before_After == "After" & CarbChem$Sampling_Day == "2019-08-05" & 
                     CarbChem$PoolID == "Ocean")
#View(Oceanday1AI)
CarbChem$DO_mg_L[Oceanday1AI]<- -36.4751 + (5.4970 * CarbChem$pH_insitu[Oceanday1AI]) + (0.1504 * CarbChem$Temp.pool[Oceanday1AI])
CarbChem$DO_p[Oceanday1AI]<- -470.450 + (66.208  * CarbChem$pH_insitu[Oceanday1AI]) + (3.935 * CarbChem$Temp.pool[Oceanday1AI])

########NEC and NEP rates##########
#calculate all the CO2Sys params
PoolCO2<-carb(flag=8, CarbChem$pH_insitu, CarbChem$TA_CRM/1000000, S=CarbChem$Salinity, T=CarbChem$Temp.pool, Patm=1, P=0, 
              Pt=CarbChem$PO_umol_L/1000000, Sit=0,
              k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential")
#TA is divided by 1000 because all calculations are in mol/kg in the seacarb package

# calculate error propagation
er<-errors(flag=8, CarbChem$pH_insitu, CarbChem$TA_CRM/1000000, 
           S=CarbChem$Salinity, T=CarbChem$Temp.pool, 
           Patm=1, P=0,Pt=CarbChem$PO_umol_L/1000000,
           Sit=0,evar1 = 0.01, evar2 = 5e-6) 

#average error for DIC based on pH and TA
mean(er$DIC*1000000)
#7.091473
sd(er$DIC*1000000)/sqrt(nrow(er))
#0.06065831

#convert CO2, HCO3, CO3, DIC, and Alk back to micromol for easier interpretation
PoolCO2[,c("CO2","HCO3","CO3","DIC","ALK")]<-PoolCO2[,c("CO2","HCO3","CO3","DIC","ALK")]*1000000

#combine all the data
CarbChem<-cbind(CarbChem,PoolCO2[-c(1:6,10:13)])

#Normalize the TA and DIC to salinity

#TA normalized to constant salinity, but not nutrients here
CarbChem$TA_NormSal<-CarbChem$TA_CRM*(CarbChem$Salinity_Lab/34)

#Normalize DIC constant salinity
CarbChem$DIC_Norm<-CarbChem$DIC*(CarbChem$Salinity_Lab/34)

#Check Relationship of DIC and TA
#TATime<- CarbChem %>%
#ggplot(aes(x= Time_Point, y = TA_NormSal, shape = Day_Night, color =Removal_Control)) +
#geom_point() +
#geom_smooth(method = "lm") +
#facet_wrap(~Foundation_spp * Before_After*PoolID, scales = "free_y") +
#theme_classic() +
#ggsave("Output/EcoMetabolism/TATime.png", width=55, height=35,dpi=300, unit="cm")
#TATime

#for TA, DIC, Sampling Time to find change over time
DeltaSamples<- CarbChem %>%
  #filter(Day_Night == 'Day') %>% #just do day and night separately since having issues
  filter(Foundation_spp != 'Ocean') %>%
  dplyr::group_by(PoolID,Group,Day_Night) %>% #grouping by pool id and beforeafter/controlremoval
  dplyr::arrange(Time_Point, .by_group = TRUE) %>% #arranges by time point and group so all time points from each tidepool
  dplyr::filter(Id_code != 'D_24_2_BC',Id_code != 'D_24_3_BC', Id_code != 'D_24_4_BC',
                Id_code != 'N_24_2_BC',Id_code != 'N_24_3_BC', Id_code != 'N_24_4_BC') %>%
  dplyr::filter(Time_Point == 1 | Time_Point == 4) %>% #taking the integrated value over low tide
  #and each grouping are in order (D_1_1_BC--D_1_4_BC then D_1_1_AI -- D_1_4_AI etc)
  dplyr::mutate(TADeltaTime = TA_NormSal - lag(TA_NormSal, default = first(TA_NormSal)), # subtracts time point 4 from 1
                DICDeltaTime = DIC_Norm-lag(DIC_Norm, default = first(DIC_Norm)),
                DeltaTime = (difftime(Date_Time,lag(Date_Time,default = first(Date_Time)))/3600), #sampling diff btwn time pts into hrs
                DeltaTime = ifelse(DeltaTime < 0 , DeltaTime + 24, DeltaTime), #if Delta Time is less than 0 then +24 hrs
                #converted to account from time went from 23/24 to 0:00 and added 24 hrs %>%
                DeltaNN = NN_umol_L - lag(NN_umol_L, default = first(NN_umol_L)),
                DeltaNH4 = NH4_umol_L - lag(NH4_umol_L, default = first (NH4_umol_L)),
                DeltaPO = PO_umol_L - lag(PO_umol_L, default =first(PO_umol_L)),
                DeltapH = pH_insitu -lag(pH_insitu, default = first (pH_insitu)),
                DeltapCO2 = pCO2 - lag(pCO2, default = first(pCO2)),
                DeltaDO = DO_mg_L - lag(DO_mg_L, default = first(DO_mg_L))) %>%
  dplyr::filter(TADeltaTime != 0) %>% #filters out row where deltas = 0
  dplyr::select(PoolID, Id_code, Time_Point,Foundation_spp,Before_After,Removal_Control,Day_Night,Group,TADeltaTime, DICDeltaTime,
                DeltaTime,DeltaNN,DeltaNH4,DeltaPO,DeltapH,DeltapCO2,DeltaDO)

#because tp 24 in Before period did not have time pt 1 so took differences btw time point 2 and 4
TP24 <- CarbChem %>%
  #filter(Day_Night == 'Day') %>% #just do day and night separately since having issues
  filter(Foundation_spp != 'Ocean') %>%
  dplyr::group_by(PoolID,Group,Day_Night) %>% #grouping by pool id and beforeafter/controlremoval
  dplyr::arrange(Time_Point, .by_group = TRUE) %>% #arranges by time point and group so all time points from each tidepool
  dplyr::filter(Id_code == 'D_24_2_BC' | Id_code == 'D_24_3_BC'| Id_code == 'D_24_4_BC'|
                  Id_code == 'N_24_2_BC'|Id_code == 'N_24_3_BC'|Id_code == 'N_24_4_BC') %>%
  dplyr::filter(Time_Point == 2 | Time_Point == 4) %>%
  #and each grouping are in order (D_1_1_BC--D_1_4_BC then D_1_1_AI -- D_1_4_AI etc)
  dplyr::mutate(TADeltaTime = TA_NormSal - lag(TA_NormSal, default = first(TA_NormSal)), # subtracts time point 1&2,2&3,3&4
                DICDeltaTime = DIC_Norm-lag(DIC_Norm, default = first(DIC_Norm)),
                DeltaTime = (difftime(Date_Time,lag(Date_Time,default = first(Date_Time)))/3600), #sampling diff btwn time pts into hrs
                DeltaTime = ifelse(DeltaTime < 0 , DeltaTime + 24, DeltaTime), #if Delta Time is less than 0 then +24 hrs
                #converted to account from time went from 23/24 to 0:00 and added 24 hrs %>%
                DeltaNN = NN_umol_L - lag(NN_umol_L, default = first(NN_umol_L)),
                DeltaNH4 = NH4_umol_L - lag(NH4_umol_L, default = first (NH4_umol_L)),
                DeltaPO = PO_umol_L - lag(PO_umol_L, default =first(PO_umol_L)),
                DeltapH = pH_insitu -lag(pH_insitu, default = first (pH_insitu)),
                DeltapCO2 = pCO2 - lag(pCO2, default = first(pCO2)),
                DeltaDO = DO_mg_L - lag(DO_mg_L, default = first(DO_mg_L))) %>%
  filter(TADeltaTime != 0) %>% #filters out row where deltas = 0
  select(PoolID, Id_code, Time_Point,Foundation_spp,Before_After,Removal_Control,Day_Night,Group,TADeltaTime, DICDeltaTime,
         DeltaTime,DeltaNN,DeltaNH4,DeltaPO,DeltapH,DeltapCO2,DeltaDO)

#combine dataframes
DeltaSamples<-rbind(DeltaSamples,TP24)

#bring in physical parameters for nec and nep calculations
PhysicalParameters<-TidePooldes[, c("PoolID", "Before_After", "SurfaceArea", "Vol", "SAtoV", "TideHeight")] #pull out the necessary columns and treatment 
PhysicalParameters$PoolID<-as.character(PhysicalParameters$PoolID)
PhysicalParameters$Before_After<-as.character(PhysicalParameters$Before_After)

#join with delta samples
DeltaSamples<-left_join(DeltaSamples,PhysicalParameters)


#create an ifelse statement for sampling volumes
#Each time point more water was taken out of each pool and need to account for change in volume
#for NEC calculation
DeltaSamples <- DeltaSamples %>%
  mutate(SamplingVolume = #ifelse(Time_Point == '2',Vol - 0.00055, #took out 550mL-->0.00055m3 
           #at time point 1 so volume at time pt 2 is -0.00055 m3
           #ifelse(Time_Point == '3', Vol - 0.00095,
           #took out 400mL at time point 2 so now volume @time3
           #is - 950mL-->0.00095m3
           ifelse(Time_Point == '4', Vol - 0.00135, Vol)) %>%
  #took out 400 mL at time point 3 so volume @time4
  #is 550+400+400-->0.00135m3
  #Some sampling times volume wasn't exact these are take from the note section in the data set 
  #note: these are for the time point after the note was taken
  mutate(AdjSamplingVolume = ifelse(Id_code == 'D_4_4_BC', SamplingVolume - 0.00055,
                                    ifelse(Id_code == 'D_18_4_BC', SamplingVolume - 0.00004,
                                           ifelse(Id_code == 'D_20_4_BC', SamplingVolume - 0.00008,
                                                  ifelse(Id_code == 'D_28_4_BC', SamplingVolume -0.0001,
                                                         ifelse(Id_code == 'D_5_4_BC', SamplingVolume -0.00005,
                                                                ifelse(Id_code == 'D_6_4_BC', SamplingVolume-0.00001,
                                                                       ifelse(Id_code == 'D_17_4_BC', SamplingVolume-0.00015,
                                                                              ifelse(Id_code == 'D_7_4_AI', SamplingVolume-0.00004,
                                                                                     ifelse(Id_code == 'D_9_4_AI', SamplingVolume-0.00004,
                                                                                            ifelse(Id_code == 'N_2_4_AI', SamplingVolume-0.0005,
                                                                                                   ifelse(Id_code == 'N_9_4_AI', SamplingVolume + 0.00001,
                                                                                                          ifelse(Id_code == 'N_1_4_AI', SamplingVolume - 0.00002,
                                                                                                                 ifelse(Id_code == 'N_2_4_AI', SamplingVolume -0.000018,
                                                                                                                        ifelse(Id_code == 'N_4_4_AI', SamplingVolume -0.000015,
                                                                                                                               ifelse(Id_code == 'N_5_4_AI',SamplingVolume -0.00001,
                                                                                                                                      ifelse(Id_code == 'N_6_4_AI',SamplingVolume -0.00004,
                                                                                                                                             ifelse(Id_code == 'N_12_4_AI', SamplingVolume -0.000035,
                                                                                                                                                    ifelse(Id_code == 'N_5_4_AI', SamplingVolume -0.00002,
                                                                                                                                                           ifelse(Id_code == 'N_9_4_AI', SamplingVolume -0.00001, SamplingVolume)))))))))))))))))))) 



#Normalizing change in TA to change in nutrients based on Wolf-Gladrow et al. 2007
DeltaSamples$DeltaTA_N_Norm<-DeltaSamples$TADeltaTime - (DeltaSamples$DeltaNN) - (2*DeltaSamples$DeltaPO) + (DeltaSamples$DeltaNH4)
#for every mol Nitrate and phosophate--increases TA by one and two moles respectively
#so subtract to normalise TA where as TA decreases with every mole of Nh4 so + Nh4 back

DeltaSamples$DeltaTA_N_Norm<- -1 * DeltaSamples$DeltaTA_N_Norm #to change to positive calcification and negative dissolution
DeltaSamples$DICDeltaTime <- -1 * (DeltaSamples$DICDeltaTime) #to change to positive photosynthesis and neg respiration 


DeltaSamples$DeltaTime<-as.numeric(DeltaSamples$DeltaTime) #as numeric since attributes from difftime function were interferring 
#with NEC/NEP calcs

#Equation for NEC rates
#NEC delta TA/2 * seawater density * (volume/ surface area) / time  Divided by 1000 so mmols m2 hr 
DeltaSamples$NEC.mmol.m2.hr<- ((DeltaSamples$DeltaTA_N_Norm)/2) * (1023) * (DeltaSamples$AdjSamplingVolume/DeltaSamples$SurfaceArea) * (1/DeltaSamples$DeltaTime) * (1/1000)
#without nutrients:
#DeltaSamples$NECnotNNorm<- ((DeltaSamples$TADeltaTime)/2) * (1023) * (DeltaSamples$AdjSamplingVolume/DeltaSamples$SurfaceArea) * (1/DeltaSamples$DeltaTime) * (1/1000)
#Function for calculation of air-sea CO2 flux
#adapted from MatLab Code by Cecilia Chapa Balcorta copyright @ 2015

#input:
#pCO2_agua= seawater pCO2 (uatm)
#pCO2_atm=  atmospheric pCO2 (uatm)
#T=  Temperature (Celsius)
#S=  Salinity 
#u = Wind speed (m/s)

#get average salinity, temp and pco2 between time points
AirSeaFlux<- CarbChem %>%
  #filter(Day_Night == 'Day') %>% #just do day and night separately since having issues
  filter(Foundation_spp != 'Ocean') %>%
  dplyr::group_by(PoolID,Group,Day_Night) %>% #grouping by pool id and beforeafter/controlremoval
  dplyr::arrange(Time_Point, .by_group = TRUE) %>% #arranges by time point and group so all time points from each tidepool
  #and each grouping are in order (D_1_1_BC--D_1_4_BC then D_1_1_AI -- D_1_4_AI etc)
  dplyr::filter(Id_code != 'D_24_2_BC',Id_code != 'D_24_3_BC', Id_code != 'D_24_4_BC',
                Id_code != 'N_24_2_BC',Id_code != 'N_24_3_BC', Id_code != 'N_24_4_BC') %>% #take out these ids since not following 4-1 rule
  dplyr::filter(Time_Point == 1 | Time_Point == 4) %>%
  dplyr::mutate(Tempmean = rollmean(Temp.pool,2, na.pad=TRUE, align="right"), #takes rolling means every two pts
                pCO2mean = rollmean(pCO2,2, na.pad=TRUE, align="right"), #avg btw time 1 & 2, time 2&3 etc. 
                Salinitymean = rollmean(Salinity,2, na.pad=TRUE, align="right"),
                Windmean = rollmean(Wind_Speed.m.s,2, na.pad=TRUE, align="right"),
                pHmean = rollmean(pH_insitu, 2, na.pad=TRUE, align="right"),
                DOmean = rollmean(DO_mg_L, 2, na.pad=TRUE, align="right"),
                NNmean = rollmean(NN_umol_L, 2, na.pad=TRUE, align="right"),
                NH4mean = rollmean(NH4_umol_L, 2, na.pad=TRUE, align="right"),
                POmean = rollmean(PO_umol_L, 2, na.pad=TRUE, align="right")) %>% 
  filter(Tempmean != 'NA') %>% #removes row that values are NA (time 1:time1) 
  select(PoolID, Id_code, Time_Point,Foundation_spp,Before_After,Removal_Control,Day_Night,Group,
         Tempmean,pCO2mean,Salinitymean,Windmean,pHmean,DOmean, NNmean, NH4mean,POmean) #selects the columns needed for air-sea flux equation

#TP24
TP24Airseaflux<- CarbChem %>%
  #filter(Day_Night == 'Day') %>% #just do day and night separately since having issues
  filter(Foundation_spp != 'Ocean') %>%
  dplyr::group_by(PoolID,Group,Day_Night) %>% #grouping by pool id and beforeafter/controlremoval
  dplyr::arrange(Time_Point, .by_group = TRUE) %>% #arranges by time point and group so all time points from each tidepool
  #and each grouping are in order (D_1_1_BC--D_1_4_BC then D_1_1_AI -- D_1_4_AI etc)
  dplyr::filter(Id_code == 'D_24_2_BC' | Id_code == 'D_24_3_BC'| Id_code == 'D_24_4_BC'|
                  Id_code == 'N_24_2_BC'|Id_code == 'N_24_3_BC'|Id_code == 'N_24_4_BC') %>%
  dplyr::filter(Time_Point == 2 | Time_Point == 4) %>%
  dplyr::mutate(Tempmean = rollmean(Temp.pool,2, na.pad=TRUE, align="right"), #takes rolling means every two pts
                pCO2mean = rollmean(pCO2,2, na.pad=TRUE, align="right"), #avg btw time 1 & 2, time 2&3 etc. 
                Salinitymean = rollmean(Salinity,2, na.pad=TRUE, align="right"),
                Windmean = rollmean(Wind_Speed.m.s,2, na.pad=TRUE, align="right"),
                pHmean = rollmean(pH_insitu, 2, na.pad=TRUE, align="right"),
                DOmean = rollmean(DO_mg_L, 2, na.pad=TRUE, align="right"),
                NNmean = rollmean(NN_umol_L, 2, na.pad=TRUE, align="right"),
                NH4mean = rollmean(NH4_umol_L, 2, na.pad=TRUE, align="right"),
                POmean = rollmean(PO_umol_L, 2, na.pad=TRUE, align="right")) %>% 
  filter(Tempmean != 'NA') %>% #removes row that values are NA (time 1:time1) 
  select(PoolID, Id_code, Time_Point,Foundation_spp,Before_After,Removal_Control,Day_Night,Group,
         Tempmean,pCO2mean,Salinitymean,Windmean,pHmean,DOmean,NNmean, NH4mean,POmean) #selects the columns needed for air-sea flux equation


#combine dataframes
AirSeaFlux<-rbind(AirSeaFlux,TP24Airseaflux)

#combine with delta samples dataframe
DeltaSamples<-left_join(DeltaSamples,AirSeaFlux)

#Setting parameters to make it easier in functions
T<-DeltaSamples$Tempmean
S<-DeltaSamples$Salinitymean
u<-DeltaSamples$Windmean

pCO2_water<-DeltaSamples$pCO2mean
pCO2_atm<- 410.86 #average uatm CO2 in atms btween July (411.77) and August 2019 (409.95) from Mauna Loa Station


#Air-sea CO2 is calculated as follows:

#FCO2 = K * a * (dpCO2) 

#Where K=is the transfer velocity according to Wanninkhof (1992).
#a = CO2 solibility constant according to Weiss (1974)
#dpCO2 is the difference of air and seawater pCO2 


#Schmidt Number function
#For water of salinity=35 and temperature range 0-30deg C
# A,B,C,D are constants
Sc<-function(T) {
  A <- 2073.1
  B <- 125.62
  C <- 3.6276
  D <- 0.043219
  Sc <- A - (B * T) + (C *(T^2)) - (D * (T^3))
  return(Sc)
}

#Solibuility constant (Weiss, 1974) 
Ko_weiss<- function(T,S) {
  A <- c(-60.2409, 93.4517, 23.3585) #mol/kg.atm
  B <- c(0.023517, -0.023656, 0.0047036)  #mol/kg.atm
  T <- T + 273.15 #Conversion from Celsius degrees to Kelvins
  Ln_Ko = A[1] + (A[2] * (100/T)) + (A[3] * (log(T/100))) + S * (B[1] + (B[2] * (T/100)) + (B[3] * (T/100)^2))
  Ko <- exp(Ln_Ko)
  return(Ko)
}

#CO2 Transfer velocity calculation 
slowwind<-0.31*(u^2) * ((Sc(T)/660)^-0.5) #for wind <=6 m/s
highwind<-0.39*(u^2) * ((Sc(T)/660)^-0.5) #for wind > 6 m/s

DeltaSamples$K<- ifelse(u <= 6, slowwind, highwind) #if else statement for K [transfer velocity according to Wanninkhof (1992)]

DeltaSamples$dpCO2 <- pCO2_water - pCO2_atm #difference of air and seawater pCO2

DeltaSamples$a <- Ko_weiss(T,S) #CO2 solibility constant according to Weiss (1974) solubility in mmol L^-1 atm^-1 or mmol m^-3 uatm^-1

DeltaSamples$F_CO2 <- ((0.24 * DeltaSamples$K * DeltaSamples$a * DeltaSamples$dpCO2)/ 24) #CO2 flux (mmol m^-2 hr^-1) (rate is by day so divide by 24 to convert per hour)


#NEP with FCO2
# NEP = (delta DIC) * SW density * (volume/SA) / time  -  NEC - Fugosity of CO2) *1000 to convert to mmol 
DeltaSamples$NEP.mmol.m2.hr<-((DeltaSamples$DICDeltaTime) * (1023) * (DeltaSamples$AdjSamplingVolume/DeltaSamples$SurfaceArea) * (1/DeltaSamples$DeltaTime) * (1/1000)) - (DeltaSamples$NEC.mmol.m2.hr) - (DeltaSamples$F_CO2)

DeltaSamples$NEC.mmol.m2.hr<-as.numeric(DeltaSamples$NEC.mmol.m2.hr) 
DeltaSamples$NEP.mmol.m2.hr<-as.numeric(DeltaSamples$NEP.mmol.m2.hr)

DeltaSamples<-left_join(DeltaSamples,TempandLightSum)

####All time points#####

#for TA, DIC, Sampling Time to find change over time
Allsamples<- CarbChem %>%
  filter(Foundation_spp != 'Ocean') %>%
  dplyr::group_by(PoolID,Group,Day_Night) %>% #grouping by pool id and beforeafter/controlremoval
  dplyr::arrange(Time_Point, .by_group = TRUE) %>% #arranges by time point and group so all time points from each tidepool
  dplyr::filter(Id_code != 'D_24_2_BC'|Id_code != 'D_24_3_BC'| Id_code != 'D_24_4_BC'|
                Id_code != 'N_24_2_BC'|Id_code != 'N_24_3_BC'| Id_code != 'N_24_4_BC'|
                Id_code != 'D_1_2_AI'|Id_code != 'D_16_2_AI'|Id_code != 'D_4_3_AI'| 
                Id_code != 'N_20_3_AI') %>% #remove samples that don't have nutrients
  #and each grouping are in order (D_1_1_BC--D_1_4_BC then D_1_1_AI -- D_1_4_AI etc)
  dplyr::mutate(TADeltaTime = TA_NormSal - lag(TA_NormSal, default = first(TA_NormSal)), # subtracts time point 1&2,2&3,3&4
                DICDeltaTime = DIC_Norm-lag(DIC_Norm, default = first(DIC_Norm)),
                DeltaTime = (difftime(Date_Time,lag(Date_Time,default = first(Date_Time)))/3600), #sampling diff btwn time pts into hrs
                DeltaTime = ifelse(DeltaTime < 0 , DeltaTime + 24, DeltaTime), #if Delta Time is less than 0 then +24 hrs
                #converted to account from time went from 23/24 to 0:00 and added 24 hrs %>%
                DeltaNN = NN_umol_L - lag(NN_umol_L, default = first(NN_umol_L)),
                DeltaNH4 = NH4_umol_L - lag(NH4_umol_L, default = first (NH4_umol_L)),
                DeltaPO = PO_umol_L - lag(PO_umol_L, default =first(PO_umol_L)),
                DeltapH = pH_insitu -lag(pH_insitu, default = first (pH_insitu)),
                DeltapCO2 = pCO2 - lag(pCO2, default = first(pCO2)),
                DeltaDO = DO_mg_L - lag(DO_mg_L, default = first(DO_mg_L))) %>%
  dplyr::filter(TADeltaTime != 0) %>% #filters out row where deltas = 0
  dplyr::select(PoolID, Id_code, Time_Point,Foundation_spp,Before_After,Removal_Control,Day_Night,Group,TADeltaTime, DICDeltaTime,
                DeltaTime,DeltaNN,DeltaNH4,DeltaPO,DeltapH,DeltapCO2,DeltaDO)

#because tp 24 in Before period did not have time pt 1 so took differences btw time point 2,3 and 4
TP24all <- CarbChem %>%
  #filter(Day_Night == 'Day') %>% #just do day and night separately since having issues
  filter(Foundation_spp != 'Ocean') %>%
  dplyr::group_by(PoolID,Group,Day_Night) %>% #grouping by pool id and beforeafter/controlremoval
  dplyr::arrange(Time_Point, .by_group = TRUE) %>% #arranges by time point and group so all time points from each tidepool
  dplyr::filter(Id_code == 'D_24_2_BC' | Id_code == 'D_24_3_BC'| Id_code == 'D_24_4_BC'|
                  Id_code == 'N_24_2_BC'|Id_code == 'N_24_3_BC'|Id_code == 'N_24_4_BC') %>%
  dplyr::filter(Time_Point == 2 | Time_Point == 3 | Time_Point == 4) %>%
  #and each grouping are in order (D_1_1_BC--D_1_4_BC then D_1_1_AI -- D_1_4_AI etc)
  dplyr::mutate(TADeltaTime = TA_NormSal - lag(TA_NormSal, default = first(TA_NormSal)), # subtracts time point 2&3,3&4
                DICDeltaTime = DIC_Norm-lag(DIC_Norm, default = first(DIC_Norm)),
                DeltaTime = (difftime(Date_Time,lag(Date_Time,default = first(Date_Time)))/3600), #sampling diff btwn time pts into hrs
                DeltaTime = ifelse(DeltaTime < 0 , DeltaTime + 24, DeltaTime), #if Delta Time is less than 0 then +24 hrs
                #converted to account from time went from 23/24 to 0:00 and added 24 hrs %>%
                DeltaNN = NN_umol_L - lag(NN_umol_L, default = first(NN_umol_L)),
                DeltaNH4 = NH4_umol_L - lag(NH4_umol_L, default = first (NH4_umol_L)),
                DeltaPO = PO_umol_L - lag(PO_umol_L, default =first(PO_umol_L)),
                DeltapH = pH_insitu -lag(pH_insitu, default = first (pH_insitu)),
                DeltapCO2 = pCO2 - lag(pCO2, default = first(pCO2)),
                DeltaDO = DO_mg_L - lag(DO_mg_L, default = first(DO_mg_L))) %>%
  filter(TADeltaTime != 0) %>% #filters out row where deltas = 0
  select(PoolID, Id_code, Time_Point,Foundation_spp,Before_After,Removal_Control,Day_Night,Group,TADeltaTime, DICDeltaTime,
         DeltaTime,DeltaNN,DeltaNH4,DeltaPO,DeltapH,DeltapCO2,DeltaDO)

#combine dataframes
Allsamples<-rbind(Allsamples,TP24all)

#bring in physical parameters for nec and nep calculations
PhysicalParameters<-TidePooldes[, c("PoolID", "Before_After", "SurfaceArea", "Vol", "SAtoV", "TideHeight")] #pull out the necessary columns and treatment 
PhysicalParameters$PoolID<-as.character(PhysicalParameters$PoolID)
PhysicalParameters$Before_After<-as.character(PhysicalParameters$Before_After)

#join with delta samples
Allsamples<-left_join(Allsamples,PhysicalParameters)


#create an ifelse statement for sampling volumes
#Each time point more water was taken out of each pool and need to account for change in volume
#for NEC calculation
Allsamples <- Allsamples %>%
  mutate(SamplingVolume = ifelse(Time_Point == '2',Vol - 0.00055, #took out 550mL-->0.00055m3 
           #at time point 1 so volume at time pt 2 is -0.00055 m3
           ifelse(Time_Point == '3', Vol - 0.00095,
           #took out 400mL at time point 2 so now volume @time3
           #is - 950mL-->0.00095m3
           ifelse(Time_Point == '4', Vol - 0.00135, Vol)))) %>%
  #took out 400 mL at time point 3 so volume @time4
  #is 550+400+400-->0.00135m3
  #Some sampling times volume wasn't exact these are take from the note section in the data set 
  #note: these are for the time point after the note was taken
  mutate(AdjSamplingVolume = ifelse(Id_code == 'D_4_2_BC', SamplingVolume - 0.00055,
ifelse(Id_code == 'D_18_2_BC', SamplingVolume - 0.00004,
ifelse(Id_code == 'D_20_2_BC', SamplingVolume - 0.00008,
ifelse(Id_code == 'D_28_2_BC', SamplingVolume -0.0001,
ifelse(Id_code == 'D_5_3_BC', SamplingVolume -0.00005,
ifelse(Id_code == 'D_6_3_BC', SamplingVolume-0.00001,
ifelse(Id_code == 'D_17_3_BC', SamplingVolume-0.00015,
ifelse(Id_code == 'D_7_3_AI', SamplingVolume-0.00004,
ifelse(Id_code == 'D_9_3_AI', SamplingVolume-0.00004,
ifelse(Id_code == 'N_2_3_AI', SamplingVolume-0.0005,
ifelse(Id_code == 'N_9_2_AI', SamplingVolume + 0.00001,
ifelse(Id_code == 'N_1_3_AI', SamplingVolume - 0.00002,
ifelse(Id_code == 'N_2_3_AI', SamplingVolume -0.000018,
ifelse(Id_code == 'N_4_3_AI', SamplingVolume -0.000015,
ifelse(Id_code == 'N_5_3_AI',SamplingVolume -0.00001,
ifelse(Id_code == 'N_6_3_AI',SamplingVolume -0.00004,
ifelse(Id_code == 'N_12_3_AI', SamplingVolume -0.000035,
ifelse(Id_code == 'D_4_4_BC', SamplingVolume - 0.00055,
ifelse(Id_code == 'D_18_4_BC', SamplingVolume - 0.00004,
ifelse(Id_code == 'D_20_4_BC', SamplingVolume - 0.00008,
ifelse(Id_code == 'D_28_4_BC', SamplingVolume -0.0001,
ifelse(Id_code == 'D_5_4_BC', SamplingVolume -0.00005,
ifelse(Id_code == 'D_6_4_BC', SamplingVolume-0.00001,
ifelse(Id_code == 'D_17_4_BC', SamplingVolume-0.00015,
ifelse(Id_code == 'D_7_4_AI', SamplingVolume-0.00004,
ifelse(Id_code == 'D_9_4_AI', SamplingVolume-0.00004,
ifelse(Id_code == 'N_2_4_AI', SamplingVolume-0.0005,
ifelse(Id_code == 'N_9_4_AI', SamplingVolume + 0.00001,
ifelse(Id_code == 'N_1_4_AI', SamplingVolume - 0.00002,
ifelse(Id_code == 'N_2_4_AI', SamplingVolume -0.000018,
ifelse(Id_code == 'N_4_4_AI', SamplingVolume -0.000015,
ifelse(Id_code == 'N_5_4_AI',SamplingVolume -0.00001,
ifelse(Id_code == 'N_6_4_AI',SamplingVolume -0.00004,
ifelse(Id_code == 'N_12_4_AI', SamplingVolume -0.000035,
ifelse(Id_code == 'N_5_4_AI', SamplingVolume -0.00002,
ifelse(Id_code == 'N_9_4_AI', SamplingVolume -0.00001, SamplingVolume)))))))))))))))))))))))))))))))))))))



#Normalizing change in TA to change in nutrients based on Wolf-Gladrow et al. 2007
Allsamples$DeltaTA_N_Norm<-Allsamples$TADeltaTime - (Allsamples$DeltaNN) - (2*Allsamples$DeltaPO) + (Allsamples$DeltaNH4)
#for every mol Nitrate and phosophate--increases TA by one and two moles respectively
#so subtract to normalise TA where as TA decreases with every mole of Nh4 so + Nh4 back

Allsamples$DeltaTA_N_Norm<- -1 * Allsamples$DeltaTA_N_Norm #to change to positive calcification and negative dissolution
Allsamples$DICDeltaTime <- -1 * (Allsamples$DICDeltaTime) #to change to positive photosynthesis and neg respiration 


Allsamples$DeltaTime<-as.numeric(Allsamples$DeltaTime) #as numeric since attributes from difftime function were interferring 
#with NEC/NEP calcs

#Equation for NEC rates
#NEC delta TA/2 * seawater density * (volume/ surface area) / time  Divided by 1000 so mmols m2 hr 
Allsamples$NEC.mmol.m2.hr<- ((Allsamples$DeltaTA_N_Norm)/2) * (1023) * (Allsamples$AdjSamplingVolume/Allsamples$SurfaceArea) * (1/Allsamples$DeltaTime) * (1/1000)
#without nutrients:
#DeltaSamples$NECnotNNorm<- ((DeltaSamples$TADeltaTime)/2) * (1023) * (DeltaSamples$AdjSamplingVolume/DeltaSamples$SurfaceArea) * (1/DeltaSamples$DeltaTime) * (1/1000)
#Function for calculation of air-sea CO2 flux
#adapted from MatLab Code by Cecilia Chapa Balcorta copyright @ 2015

#input:
#pCO2_agua= seawater pCO2 (uatm)
#pCO2_atm=  atmospheric pCO2 (uatm)
#T=  Temperature (Celsius)
#S=  Salinity 
#u = Wind speed (m/s)

#get average salinity, temp and pco2 between time points
AirSeaFluxall<- CarbChem %>%
  #filter(Day_Night == 'Day') %>% #just do day and night separately since having issues
  filter(Foundation_spp != 'Ocean') %>%
  dplyr::group_by(PoolID,Group,Day_Night) %>% #grouping by pool id and beforeafter/controlremoval
  dplyr::arrange(Time_Point, .by_group = TRUE) %>% #arranges by time point and group so all time points from each tidepool
  #and each grouping are in order (D_1_1_BC--D_1_4_BC then D_1_1_AI -- D_1_4_AI etc)
  dplyr::filter(Id_code != 'D_24_2_BC'|Id_code != 'D_24_3_BC'| Id_code != 'D_24_4_BC'|
                Id_code != 'N_24_2_BC'|Id_code != 'N_24_3_BC'| Id_code != 'N_24_4_BC'|
                Id_code != 'D_1_2_AI' |Id_code != 'D_16_2_AI'|Id_code != 'D_4_3_AI'| 
                Id_code != 'N_20_3_AI') %>% #take out these ids since not following 4-1 rule
  #dplyr::filter(Time_Point == 1 | Time_Point == 4) %>%
  dplyr::mutate(Tempmean = rollmean(Temp.pool,2, na.pad=TRUE, align="right"), #takes rolling means every two pts
                pCO2mean = rollmean(pCO2,2, na.pad=TRUE, align="right"), #avg btw time 1 & 2, time 2&3 etc. 
                Salinitymean = rollmean(Salinity,2, na.pad=TRUE, align="right"),
                Windmean = rollmean(Wind_Speed.m.s,2, na.pad=TRUE, align="right"),
                pHmean = rollmean(pH_insitu, 2, na.pad=TRUE, align="right"),
                DOmean = rollmean(DO_mg_L, 2, na.pad=TRUE, align="right"),
                NNmean = rollmean(NN_umol_L, 2, na.pad=TRUE, align="right"),
                NH4mean = rollmean(NH4_umol_L, 2, na.pad=TRUE, align="right"),
                POmean = rollmean(PO_umol_L, 2, na.pad=TRUE, align="right")) %>% 
  filter(Tempmean != 'NA') %>% #removes row that values are NA (time 1:time1) 
  select(PoolID, Id_code, Time_Point,Foundation_spp,Before_After,Removal_Control,Day_Night,Group,
         Tempmean,pCO2mean,Salinitymean,Windmean,pHmean,DOmean, NNmean, NH4mean,POmean) #selects the columns needed for air-sea flux equation

#TP24
TP24Airseafluxall<- CarbChem %>%
  #filter(Day_Night == 'Day') %>% #just do day and night separately since having issues
  filter(Foundation_spp != 'Ocean') %>%
  dplyr::group_by(PoolID,Group,Day_Night) %>% #grouping by pool id and beforeafter/controlremoval
  dplyr::arrange(Time_Point, .by_group = TRUE) %>% #arranges by time point and group so all time points from each tidepool
  #and each grouping are in order (D_1_1_BC--D_1_4_BC then D_1_1_AI -- D_1_4_AI etc)
  dplyr::filter(Id_code == 'D_24_2_BC' | Id_code == 'D_24_3_BC'| Id_code == 'D_24_4_BC'|
                  Id_code == 'N_24_2_BC'|Id_code == 'N_24_3_BC'|Id_code == 'N_24_4_BC') %>%
  dplyr::filter(Time_Point == 2 |Time_Point == 3| Time_Point == 4) %>% #no time point 1
  dplyr::mutate(Tempmean = rollmean(Temp.pool,2, na.pad=TRUE, align="right"), #takes rolling means every two pts
                pCO2mean = rollmean(pCO2,2, na.pad=TRUE, align="right"), #avg btw time 2& 3, time 3&4 etc. 
                Salinitymean = rollmean(Salinity,2, na.pad=TRUE, align="right"),
                Windmean = rollmean(Wind_Speed.m.s,2, na.pad=TRUE, align="right"),
                pHmean = rollmean(pH_insitu, 2, na.pad=TRUE, align="right"),
                DOmean = rollmean(DO_mg_L, 2, na.pad=TRUE, align="right"),
                NNmean = rollmean(NN_umol_L, 2, na.pad=TRUE, align="right"),
                NH4mean = rollmean(NH4_umol_L, 2, na.pad=TRUE, align="right"),
                POmean = rollmean(PO_umol_L, 2, na.pad=TRUE, align="right")) %>% 
  filter(Tempmean != 'NA') %>% #removes row that values are NA (time 1:time1) 
  select(PoolID, Id_code, Time_Point,Foundation_spp,Before_After,Removal_Control,Day_Night,Group,
         Tempmean,pCO2mean,Salinitymean,Windmean,pHmean,DOmean,NNmean, NH4mean,POmean) #selects the columns needed for air-sea flux equation


#combine dataframes
AirSeaFluxall<-rbind(AirSeaFluxall,TP24Airseafluxall)

#combine with delta samples dataframe
Allsamples<-left_join(Allsamples,AirSeaFluxall)

#Setting parameters to make it easier in functions
T<-Allsamples$Tempmean
S<-Allsamples$Salinitymean
u<-Allsamples$Windmean

pCO2_water<-Allsamples$pCO2mean
pCO2_atm<- 410.86 #average uatm CO2 in atms btween July (411.77) and August 2019 (409.95) from Mauna Loa Station


#Air-sea CO2 is calculated as follows:

#FCO2 = K * a * (dpCO2) 

#Where K=is the transfer velocity according to Wanninkhof (1992).
#a = CO2 solibility constant according to Weiss (1974)
#dpCO2 is the difference of air and seawater pCO2 


#Schmidt Number function
#For water of salinity=35 and temperature range 0-30deg C
# A,B,C,D are constants
Sc<-function(T) {
  A <- 2073.1
  B <- 125.62
  C <- 3.6276
  D <- 0.043219
  Sc <- A - (B * T) + (C *(T^2)) - (D * (T^3))
  return(Sc)
}

#Solibuility constant (Weiss, 1974) 
Ko_weiss<- function(T,S) {
  A <- c(-60.2409, 93.4517, 23.3585) #mol/kg.atm
  B <- c(0.023517, -0.023656, 0.0047036)  #mol/kg.atm
  T <- T + 273.15 #Conversion from Celsius degrees to Kelvins
  Ln_Ko = A[1] + (A[2] * (100/T)) + (A[3] * (log(T/100))) + S * (B[1] + (B[2] * (T/100)) + (B[3] * (T/100)^2))
  Ko <- exp(Ln_Ko)
  return(Ko)
}

#CO2 Transfer velocity calculation 
slowwind<-0.31*(u^2) * ((Sc(T)/660)^-0.5) #for wind <=6 m/s
highwind<-0.39*(u^2) * ((Sc(T)/660)^-0.5) #for wind > 6 m/s

Allsamples$K<- ifelse(u <= 6, slowwind, highwind) #if else statement for K [transfer velocity according to Wanninkhof (1992)]

Allsamples$dpCO2 <- pCO2_water - pCO2_atm #difference of air and seawater pCO2

Allsamples$a <- Ko_weiss(T,S) #CO2 solibility constant according to Weiss (1974) solubility in mmol L^-1 atm^-1 or mmol m^-3 uatm^-1

Allsamples$F_CO2 <- ((0.24 * Allsamples$K * Allsamples$a * Allsamples$dpCO2)/ 24) #CO2 flux (mmol m^-2 hr^-1) (rate is by day so divide by 24 to convert per hour)


#NEP with FCO2
# NEP = (delta DIC) * SW density * (volume/SA) / time  -  NEC - Fugosity of CO2) *1000 to convert to mmol 
Allsamples$NEP.mmol.m2.hr<-((Allsamples$DICDeltaTime) * (1023) * (Allsamples$AdjSamplingVolume/Allsamples$SurfaceArea) * (1/Allsamples$DeltaTime) * (1/1000)) - (Allsamples$NEC.mmol.m2.hr) - (Allsamples$F_CO2)

Allsamples$NEC.mmol.m2.hr<-as.numeric(Allsamples$NEC.mmol.m2.hr) 
Allsamples$NEP.mmol.m2.hr<-as.numeric(Allsamples$NEP.mmol.m2.hr)


######PCA of biogeochem########
CarbChem$Sampling_Day<-as.factor(CarbChem$Sampling_Day)


Biogeochem.sum<-CarbChem %>%
  dplyr::group_by(PoolID,Foundation_spp,Before_After,Removal_Control)  %>% #group PoolID
  dplyr::select(PoolID,Foundation_spp,Before_After,Removal_Control,Temp.pool,DO_mg_L,pH_insitu,NN_umol_L,NH4_umol_L,PO_umol_L) %>% #selects columns of data 
  dplyr::summarise_if(is.numeric,.funs=c(mean="mean",var="var",max="max")) #take everything thats numeric and find mean, var,   min and max
#View(Biogeochem.sum)

#add foundation species cover and physical parameters to individual datasets
#create dataframe that has both before (baseline so low number in this case 1) and after
#surfgrass and mussel loss that can be used in plot later as color/size variable 
PFunsppcover<-Funsppcover%>%
  filter(Foundation_spp =="Phyllospadix") 
PFunsppcover<-PFunsppcover[c(1,5)] #select 
PFunsppcover$Before_After<- "After"
PBeforespploss<-PhyllocommunitynMDS%>%
  select(PoolID,Before_After) %>%
  filter(Before_After =='Before')
PBeforespploss$Phyllodelta<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) #create column of baseline of phyllodelta
Phylloloss<-rbind(PBeforespploss,PFunsppcover) #combine before and after together

MFunsppcover<-Funsppcover%>%
  filter(Foundation_spp =="Mytilus") 
MFunsppcover<-MFunsppcover[c(1,4)] #select poolid and mytilus delta
MFunsppcover$Before_After<- "After"
MBeforespploss<-MytiluscommunitynMDS%>%
  select(PoolID,Before_After) %>%
  filter(Before_After =='Before')
MBeforespploss$Mytilusdelta<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) #create column of baseline of phyllodelta
Musselloss<-rbind(MBeforespploss,MFunsppcover) #combine before and after together


PhylloBiogeochem <-Biogeochem.sum %>%
  filter(Foundation_spp !="Mytilus")
MytilusBiogeochem <-Biogeochem.sum %>%
  filter(Foundation_spp !="Phyllospadix")

PhylloBiogeochem<-left_join(PhylloBiogeochem,Phylloloss)
PhysicalParameters$PoolID<-as.character(PhysicalParameters$PoolID)
PhylloBiogeochem<-left_join(PhylloBiogeochem,PhysicalParameters) 
MytilusBiogeochem <-left_join(MytilusBiogeochem,Musselloss)
MytilusBiogeochem<-left_join(MytilusBiogeochem,PhysicalParameters) 

PhylloBiogeochemPCA<-PhylloBiogeochem[,c("DO_mg_L_mean","DO_mg_L_var","DO_mg_L_max",
                                         "pH_insitu_mean","pH_insitu_var","pH_insitu_max",
                                         "NN_umol_L_mean","NN_umol_L_var","NN_umol_L_max",
                                         "NH4_umol_L_mean","NH4_umol_L_var","NH4_umol_L_max",
                                         "PO_umol_L_var","PO_umol_L_mean","PO_umol_L_max","Temp.pool_mean",
                                         "Temp.pool_var","Temp.pool_max")] #Take out columns that will not be in PCA

MytilusBiogeochemPCA<-MytilusBiogeochem[,c("DO_mg_L_mean","DO_mg_L_var","DO_mg_L_max",
                                           "pH_insitu_mean","pH_insitu_var","pH_insitu_max",
                                           "NN_umol_L_mean","NN_umol_L_var","NN_umol_L_max",
                                           "NH4_umol_L_mean","NH4_umol_L_var","NH4_umol_L_max",
                                           "PO_umol_L_var","PO_umol_L_mean","PO_umol_L_max","Temp.pool_mean",
                                           "Temp.pool_var","Temp.pool_max")] #Take out columns that will not be in PCA
#Run the PCA with z-scored data for both phyllo and mytilus model
PhylloBiogeochemPCAmodel<-prcomp(PhylloBiogeochemPCA, center = TRUE, scale. = TRUE)
#PCAmodel<-princomp(BCPCA1.scale)
MytilusBiogeochemPCAmodel<-prcomp(MytilusBiogeochemPCA, center = TRUE, scale. = TRUE)

summary(PhylloBiogeochemPCAmodel)#shows amount of variance explained by each axis
#This is actually a big wide table, with all 38 Principle Components
#It tell us the amount of variance explained by each axis, and the cumulative proportion of variance explained
#PC1  0.3536 PC2 0.2669total 0.6205
summary(MytilusBiogeochemPCAmodel)
#pc1 0.3734 pc2  0.2440 0.6174


PhylloPCAscores<- data.frame(PhylloBiogeochemPCAmodel$x[,c(1,2)])#asks for first rows (PC1 and PC2) of table
MytilusPCAscores<- data.frame(MytilusBiogeochemPCAmodel$x[,c(1,2)])
#View(PCAscores)
biplot(PhylloBiogeochemPCAmodel, xlab="PC1", ylab="PC2")
PhylloBiogeochemPCAgraph<-data.frame(bind_cols(PhylloBiogeochem[,c("PoolID","Foundation_spp","Before_After","Removal_Control","Phyllodelta","SAtoV","TideHeight")], PhylloPCAscores)) #column binds scores back with Tp pool Id, foundation species, removal/control, before/after
MytilusBiogeochemPCAgraph<-data.frame(bind_cols(MytilusBiogeochem[,c("PoolID","Foundation_spp","Before_After","Removal_Control","Mytilusdelta","SAtoV","TideHeight")], MytilusPCAscores)) #column binds scores back with Tp pool Id, foundation species, removal/control, before/after



#Now to do the MANOVA
#remove ocean
PhylloBiogeochemmanova<-PhylloBiogeochemPCAgraph%>%
  filter(Foundation_spp !="Ocean")
phylloManova<-manova(cbind(PC1,PC2)~Phyllodelta  + SAtoV +TideHeight , data=PhylloBiogeochemmanova)
summary(phylloManova, test="Pillai")

mytilusBiogeochemmanova<-MytilusBiogeochemPCAgraph%>%
  filter(Foundation_spp !="Ocean")
mytilusmanova<-manova(cbind(PC1,PC2)~Mytilusdelta + SAtoV +TideHeight, data=mytilusBiogeochemmanova)
summary(mytilusmanova, test="Pillai")

########Combined data with Before/After##########
#Plot with just point data

## add a column combining after before and foundation species
PhylloBiogeochemPCAgraph$AB_F<-factor(paste(PhylloBiogeochemPCAgraph$Before_After, PhylloBiogeochemPCAgraph$Foundation_spp))
MytilusBiogeochemPCAgraph$AB_F<-factor(paste(MytilusBiogeochemPCAgraph$Before_After, MytilusBiogeochemPCAgraph$Foundation_spp))
#create dataframe for centroids with median from x and y axes
pcentroids <- aggregate(cbind(PC1,PC2)~AB_F*Before_After*Removal_Control,PhylloBiogeochemPCAgraph,mean)
mcentroids<-aggregate(cbind(PC1,PC2)~AB_F*Before_After*Removal_Control,MytilusBiogeochemPCAgraph,mean)


#create groupings for shape labels
mgroupings<-c("After Mytilus" = 2,  "Before Mytilus" = 17, 
              "Before Ocean" = 15, "After Ocean" = 0)
pgroupings<-c("After Phyllospadix" = 2,"Before Phyllospadix" = 17,"Before Ocean" = 15, "After Ocean" = 0)

#for arrows fucnction:
x0 <- pcentroids %>%
  filter(Before_After == "Before") %>%
  select(PC1)
x0<-as.matrix(x0)
y0 <- pcentroids %>%
  filter(Before_After == "Before") %>%
  select(PC2)
y0<-as.matrix(y0)
x1<-pcentroids %>%
  filter(Before_After == "After") %>%
  select(PC1)
x1<-as.matrix(x1)
y1<-pcentroids %>%
  filter(Before_After == "After") %>%
  select(PC2)
y1<-as.matrix(y1)


Surfgrassplot<-ggplot(PhylloBiogeochemPCAgraph, aes(x = PC1 , y= PC2,shape = AB_F)) + #basic plot
  geom_point(aes(color =Phyllodelta, size =Phyllodelta, stroke=2), shape=16) +
  scale_color_distiller(palette = "Greens",guide = "legend")+
  scale_size(range = c(1,15)) +
  geom_point(data=pcentroids, size=10, stroke = 2.75) +
  theme_classic() +
  scale_shape_manual(values = c(pgroupings))+
  geom_segment(aes(x = x0[1], y = y0[1], xend = (x1[1]), yend = (y1[1])),size = 1,#segment with arrow for surfgrass before/after control
               colour = "#3182bd", arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  geom_segment(aes(x = x0[3], y = y0[3], xend = (x1[3]), yend = (y1[3])),linetype = 2,size = 1, #segment with arrow for surfgrass before/after removal
               colour = "#bdbdbd",arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  geom_segment(aes(x = x0[2], y = y0[2], xend = (x1[2]), yend = (y1[2])),size = .5, #segment with arrow for ocean before and after
               colour ="#de2d26", arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  labs(x ='PC1 (32.88%)', y = 'PC2 (27.9%)', shape='', color='Surfgrass loss',size='Surfgrass loss', linetype ='Before or after') +
  #PC1  0.3288PC2 0.2790total 0.6079
  theme(axis.text = element_text(color = "black", size = 35), 
        axis.title.x = element_text(color="black", size=40), 
        axis.title.y = element_text(color="black", size=40), 
        legend.title = element_text(color="black", size=40), 
        legend.text = element_text(color = "black", size = 35), 
        legend.position= "top",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  guides(shape = "none")+
  guides(colour = guide_legend(nrow = 1))#makes legend only one row
Surfgrassplot
#ggsave(filename = "Output/Phyllopcagraph.pdf", useDingbats =FALSE,dpi=600,device = "pdf", width = 15, height = 10)
#plot with loadings and point
surfgrass<-autoplot(PhylloBiogeochemPCAmodel, 
                                    loadings = TRUE, loadings.colour = 'black',
                                    loadings.label = TRUE, loadings.label.size = 12 , loadings.label.colour = '#2b8cbe',loadings.label.repel=TRUE, loadings.label.vjust = 1.2) +
  theme_classic() + 
  theme(legend.text=element_text(size=24)) +
  theme(axis.text = element_text(color = "black", size = 35), 
        axis.title.x = element_text(color="black", size=40), 
        axis.title.y = element_text(color="black", size=40), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) 
surfgrass
#ggsave(filename = "Output/Phyllopcaloadings.pdf", useDingbats =FALSE,dpi=600,device = "pdf", width = 25, height = 20)

library(patchwork)

surfgrasspca<-Surfgrassplot+surfgrass+
  plot_annotation(tag_levels = 'A') &         #label each individual plot with letters A-G
  theme(plot.tag = element_text(size =35, face = "bold"))   #edit the lettered text
surfgrasspca
ggsave(filename = "Output/combinedphyllopca.pdf", useDingbats =FALSE,dpi=600,device = "pdf", width = 35, height = 20)

#for arrows fucnction:
v0 <- mcentroids %>%
  filter(Before_After == "Before") %>%
  select(PC1)
v0<-as.matrix(v0)
z0 <- mcentroids %>%
  filter(Before_After == "Before") %>%
  select(PC2)
z0<-as.matrix(z0)
v1<-mcentroids %>%
  filter(Before_After == "After") %>%
  select(PC1)
v1<-as.matrix(v1)
z1<-mcentroids %>%
  filter(Before_After == "After") %>%
  select(PC2)
z1<-as.matrix(z1)


musselplot<-ggplot(MytilusBiogeochemPCAgraph, aes(x = PC1 , y= PC2,shape = AB_F)) + #basic plot
  geom_point(aes(color =Mytilusdelta, size =Mytilusdelta, stroke=2), shape=16) +
  scale_color_distiller(palette = "Blues",guide = "legend")+
  scale_size(range = c(1,15)) +
  geom_point(data=mcentroids, size=10, stroke = 2.75) +
  theme_classic() +
  scale_shape_manual(values = c(mgroupings))+
  geom_segment(aes(x = v0[1], y = z0[1], xend = (v1[1]), yend = (z1[1])),size = 1,#segment with arrow for surfgrass before/after control
               colour = "#3182bd", arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  geom_segment(aes(x = v0[3], y = z0[3], xend = (v1[3]), yend = (z1[3])),linetype = 2,size = 1, #segment with arrow for surfgrass before/after removal
               colour = "#bdbdbd",arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  geom_segment(aes(x = v0[2], y = z0[2], xend = (v1[2]), yend = (z1[2])),size = .5, #segment with arrow for ocean before and after
               colour ="#de2d26", arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  labs(x ='PC1 (38.08%)', y = 'PC2 (23.83%)', shape='', color='Mussel loss',size='Mussel loss', linetype ='Before or after') +
  #pc1 0.3808 pc2 0.2383 0.6191
  theme(axis.text = element_text(color = "black", size = 35), 
        axis.title.x = element_text(color="black", size=40), 
        axis.title.y = element_text(color="black", size=40), 
        legend.title = element_text(color="black", size=40), 
        legend.text = element_text(color = "black", size = 35), 
        legend.position= "top",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  guides(shape = "none")+
  guides(colour = guide_legend(nrow = 1))#makes legend only one row
musselplot
#ggsave(filename = "Output/Phyllopcagraph.pdf", useDingbats =FALSE,dpi=600,device = "pdf", width = 15, height = 10)
#plot with loadings and point
musselloadings<-autoplot(MytilusBiogeochemPCAmodel, 
                    loadings = TRUE, loadings.colour = 'black',
                    loadings.label = TRUE, loadings.label.size = 12 ,
                    loadings.label.colour = '#2b8cbe',loadings.label.repel=TRUE, loadings.label.vjust = 1.2) +
  theme_classic() + 
  theme(legend.text=element_text(size=24)) +
  theme(axis.text = element_text(color = "black", size = 35), 
        axis.title.x = element_text(color="black", size=40), 
        axis.title.y = element_text(color="black", size=40), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) 
musselloadings
#ggsave(filename = "Output/Phyllopcaloadings.pdf", useDingbats =FALSE,dpi=600,device = "pdf", width = 25, height = 20)

library(patchwork)

musselpca<-musselplot+musselloadings+
  plot_annotation(tag_levels = 'A') &         #label each individual plot with letters A-G
  theme(plot.tag = element_text(size =35, face = "bold"))   #edit the lettered text
musselpca
ggsave(filename = "Output/combinedmusselpca.pdf", useDingbats =FALSE,dpi=600,device = "pdf", width = 35, height = 20)


#####for supplemental summary table######

CarbChem$daysample<-as.factor(CarbChem$Sampling_Day)

Sumcarbchem<-CarbChem%>%
  select(PoolID, Foundation_spp, Before_After, Day_Night,DO_mg_L,
         PO_umol_L,NN_umol_L,NH4_umol_L,TA_NormSal,pH_insitu,pCO2,DIC_Norm,Temp.pool)
#ocean samples
Oceansamples<-CarbChem%>%
  select(Foundation_spp, Before_After, Day_Night,DO_mg_L,
         PO_umol_L,NN_umol_L,NH4_umol_L,TA_NormSal,pH_insitu,pCO2,DIC_Norm,Temp.pool,daysample) %>%
  filter(Foundation_spp=="Ocean") %>% 
  group_by(Before_After,Day_Night,daysample) %>%
  summarise_all(.funs=c("min","max")) #min and max values

#write.csv(Oceansamples,file="Output/oceanchemsamp.csv")

TidePooldes$PoolID<-as.character(TidePooldes$PoolID)
Sumcarbchem<-left_join(Sumcarbchem,TidePooldes) #combine with the physical parameters

#other carb chem parameters
Sumcarbchem<-Sumcarbchem %>%
  filter(PoolID !='Ocean') %>%
  group_by(Foundation_spp, Before_After, Day_Night,Removal_Control) %>%
  summarise_all(.funs=c("min","max"))
#write.csv(Sumcarbchem, file="Output/ChemdatanoNEPNEC.csv")

#NEC and NEP
necnepminmax<-Allsamples %>%
  select(Foundation_spp,Removal_Control, Before_After,Day_Night,NEC.mmol.m2.hr, NEP.mmol.m2.hr)%>%
  group_by(Foundation_spp,Removal_Control, Before_After) %>%
  summarise_all(.funs=c("min","max")) 
write.csv(necnepminmax, file="Output/NEPNEC.csv")
