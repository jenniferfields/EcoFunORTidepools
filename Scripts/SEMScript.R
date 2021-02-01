## Data processing and Structural equation model for both surfgrass and mussel tide pools
## By: Jenn Fields
## Last updated: 1.16.2021
###################
## clear workspace
rm(list=ls())

#load libraries
library(piecewiseSEM)
library(ggeffects)
library(patchwork)
library(MASS)
library(car)
library(lavaan)
devtools::install_github("seananderson/ggsidekick")
library(ggsidekick) #for theme sleek
library(tidyverse)
library(lmerTest)
library(lme4)
library(nlme)
library(GGally)

#source code scripts
source("scripts/tidepoolphysicalparameters.R")
source("scripts/CleanCarbChem.R") #sources both community comp and carb chem code
source("scripts/plot.psemfun.R") #plot draft SEM


#change poolid to character
SEMcommunitydata$PoolID<-as.character(SEMcommunitydata$PoolID)

#combine community comp with NEC/NEP
SEMall<-left_join(Allsamples,SEMcommunitydata)
  
#carb chem--averages of pH and nutrients
#All samples--averages of nec and nep
#TempandLightsumall --averages of temp and light

#summariseaverage over all time points for pH and nutrients (usually n=4 day/n=4 night)
#Missing values from one time point in before for tp 24 (both night/day) 
#and one time point in after from tp 1,4,16 (day),& 20(night)
Summarybiogeochem<-CarbChem%>%
  filter(Foundation_spp != "Ocean") %>% 
  dplyr::group_by(PoolID, Foundation_spp, Removal_Control,Day_Night,Before_After,Time_Point)  %>%
  mutate(NtoP = (NH4_umol_L+NN_umol_L)/PO_umol_L) %>% 
  #group by time point to get n:p ratio per tide pool per time point
  dplyr::group_by(PoolID,Foundation_spp, Removal_Control,Day_Night,Before_After)  %>%
  summarise(AvgpH =mean(pH_insitu),
            AvgNtoP =mean(NtoP),
            AvgNN=mean(NN_umol_L),
            AvgP=mean(PO_umol_L),
            AvgNh4=mean(NH4_umol_L))

  
#summarise average max temp over the low tide sample (n=3-4 values per tide pool/time period)
SummaryTemp<-TempandLightSumall%>%
  dplyr::group_by(PoolID, Day_Night,Before_After)  %>%
  summarise(AvgTempmax = mean(Temp.max),
            Avgtemp=mean(Temp.mean),
            AvgPar = mean(Par.mean)) # to show light and temp correlation, not included in SEM model

#leftjoin biogeochem and temp together
Sumbiogeochemtemp<-left_join(Summarybiogeochem,SummaryTemp)


#summarise average NEC and NEP values over the low tide period (n=2-3 samples per pool/time period) &
#community metrics and physical parameters 
SummaryNEPNECcommpp<-SEMall%>%
  dplyr::group_by(PoolID, Foundation_spp, Removal_Control,Day_Night,Before_After)  %>%
  summarise(AvgNEP = mean(NEP.mmol.m2.hr),
            AvgNEC = mean(NEC.mmol.m2.hr),
            AvgMytilus = mean(MusselCover),
            AvgPhyllo = mean(SurfgrassCover),
            AvgFleshyalgae = mean(macroalgae),
            AvgCCA = mean(allCCA),
            Avgproddom = mean(prodphyllodom))

SEMdata<-left_join(DeltaSamples,SEMcommunitydata)
SEMdata<-as.data.frame(SEMdata)

IntegratedNECNEP<- SEMdata %>%
  dplyr::group_by(PoolID, Foundation_spp, Removal_Control, Before_After,Day_Night)  %>%
  dplyr::summarise(IntegratedNEP = mean(NEP.mmol.m2.hr),
                   IntegratedNEC = mean(NEC.mmol.m2.hr))

PP<-SEMall%>%
  dplyr::group_by(Foundation_spp,PoolID,Removal_Control,Day_Night) %>%
  summarise(SAVav = mean(SAtoV),THav = mean(TideHeight),SAav=mean(SurfaceArea),Vav=mean(Vol)) 


SEMcombined<-left_join(SummaryNEPNECcommpp,Sumbiogeochemtemp)
SEMcombined<-left_join(SEMcombined,IntegratedNECNEP)

SEMcombined<-as.data.frame(SEMcombined)
  
#Takes the delta between after-before period
SEMallavg<-SEMcombined%>%
  dplyr::group_by(PoolID, Foundation_spp, Removal_Control,Day_Night)  %>%
  dplyr::summarise(NEPdelta = AvgNEP[Before_After == 'After'] - AvgNEP[Before_After == 'Before'],
                   NECdelta = AvgNEC[Before_After == 'After'] - AvgNEC[Before_After == 'Before'],
                   IntNEP = IntegratedNEP[Before_After == 'After'] - IntegratedNEP[Before_After == 'Before'],
                   IntNEC = IntegratedNEC[Before_After == 'After'] - IntegratedNEC[Before_After == 'Before'],
                   pHdelta= AvgpH[Before_After == 'After'] - AvgpH[Before_After == 'Before'],
                   NtoPdelta = AvgNtoP[Before_After == 'After'] - AvgNtoP[Before_After == 'Before'],
                   NNdelta = AvgNN[Before_After == 'After'] - AvgNN[Before_After == 'Before'],
                   Tempmaxdelta = AvgTempmax[Before_After == 'After'] - AvgTempmax[Before_After == 'Before'],
                   Tempmeandelta = Avgtemp[Before_After == 'After'] - Avgtemp[Before_After == 'Before'],
                   Parmeandelta = AvgPar[Before_After == 'After'] - AvgPar[Before_After == 'Before'],
                   Mytiluscover =AvgMytilus[Before_After == 'After'] - AvgMytilus[Before_After == 'Before'],
                   Mytilusdelta = -1*(AvgMytilus[Before_After == 'After'] - AvgMytilus[Before_After == 'Before']),
                   Phyllodelta= -1*(AvgPhyllo[Before_After == 'After'] - AvgPhyllo[Before_After == 'Before']),
                   #multiply delta of mytilus and phyllospadix by -1 to change to Percent loss
                   micromacroalgaedelta = (AvgFleshyalgae[Before_After == 'After'] - AvgFleshyalgae[Before_After == 'Before']), 
                   #only micro/acroalgae cover
                   CCAdelta = (AvgCCA[Before_After == 'After'] - AvgCCA[Before_After == 'Before']), 
                   #crustose CCA and articulated CA
                   prodomdelta= (Avgproddom[Before_After == 'After'] -  Avgproddom[Before_After == 'Before']))
                   #macroalgae cover - sessile consumers (including Mytilus)
SEMallavg<-left_join(SEMallavg,PP)                  

####Surfgrass SEM####
PhylloDayNightall<-SEMallavg %>%
  filter(Foundation_spp == 'Phyllospadix') %>%
  dplyr::rename(NEP =NEPdelta, NEC=NECdelta, pH=pHdelta, NtoPRatio=NtoPdelta,
                MaxTemp = Tempmaxdelta,Light=Parmeandelta,MytilusLoss=Mytilusdelta,PhyllospadixLoss=Phyllodelta, MicroMacroAlgaeCover=micromacroalgaedelta,
                SAtoVRatio=SAVav,TideHeight=THav,SurfaceArea=SAav,Volume=Vav, NN=NNdelta) #rename cols to match sem

#ggpairs(PhylloDayNightall[c(5:12,16:17,20:23)])

#make day/night a factor for multigroup analysis
PhylloDayNightall$Day_Night<-as.factor(PhylloDayNightall$Day_Night)

#models based off hypotheses
PDNMMAlgaeall<-lm(MicroMacroAlgaeCover~ PhyllospadixLoss+Volume+TideHeight,  data = PhylloDayNightall)
PDNTempall<-lm(MaxTemp~ PhyllospadixLoss+Volume+TideHeight, data =PhylloDayNightall)
PDNNtoPall<-lm(NtoPRatio ~PhyllospadixLoss+Volume+TideHeight, data =PhylloDayNightall)
PDNNECall<- lm(NEC~pH +MaxTemp+ TideHeight,data =PhylloDayNightall) #removed surfgrass loss since very correlated with pH and maxTemp
PDNpHall<- lm(pH ~NEP+PhyllospadixLoss+Volume+TideHeight,data = PhylloDayNightall)
PDNNEPall<-lm(NEP ~MaxTemp +MicroMacroAlgaeCover+NtoPRatio+TideHeight, data =PhylloDayNightall) 


qqp(resid(PDNMMAlgaeall),"norm")
#plot(PDNMMAlgaeall)
#qqp(resid(PDNLightall),"norm")
qqp(resid(PDNpHall),"norm") 
#plot(PDNpHall))
qqp(resid(PDNTempall),"norm") 
#plot(PDNTempall)
qqp(resid(PDNNtoPall),"norm") 
#plot(PDNNtoPall)
qqp(resid(PDNNECall),"norm") 
#plot(PDNNECall)
qqp(resid(PDNNEPall),"norm") 
#plot(PDNNEPall)

PhylloDNSEMall<-psem(
  PDNMMAlgaeall,
  PDNTempall,
  PDNNtoPall,
  PDNNEPall,
  PDNpHall,
  PDNNECall)

OutputPhylloMGSEM<-multigroup(PhylloDNSEMall,standardize = 'scale', test.type = "III",group ="Day_Night")

#summary(PhylloDNSEMall,standardize = 'scale',center = "TRUE") #scale data in sum

PhylloModels<-as.data.frame(OutputPhylloMGSEM$global)
Panovaoutput<-as.data.frame(OutputPhylloMGSEM$anovaInts)
PDayNightcoeffs<-as.data.frame(OutputPhylloMGSEM$group.coefs)
PFisher<-as.data.frame(OutputPhylloMGSEM$Cstat)

write.csv(PhylloModels, 'Output/DaynightphylloglobalUPDATED.csv' )
write.csv(Panovaoutput, 'Output/phylloanovaUPDATED.csv' )
write.csv(PDayNightcoeffs, 'Output/phylloDNcoeffUPDATED.csv' )
write.csv(PFisher, 'Output/phyllodaynightfishersUPDATED.csv' )

######mussels#####
Mytilusdaynightall<-SEMallavg%>%
  filter(Foundation_spp == 'Mytilus') %>%
  dplyr::rename(NEP =NEPdelta, NEC=NECdelta, pH=pHdelta,NtoPRatio=NtoPdelta,
                MaxTemp = Tempmaxdelta, Light=Parmeandelta, 
                MytilusLoss=Mytilusdelta,PhyllospadixLoss=Phyllodelta, MicroMacroAlgaeCover=micromacroalgaedelta,
                SAtoVRatio=SAVav,TideHeight=THav, SurfaceArea=SAav,Volume=Vav, NN=NNdelta) #rename cols to match sem

#ggpairs(Mytilusdaynightall[c(5:11,13,15,18:21)])
#tide height and mussel loss are not correlated with day/night separately 
#must be an artefact of having double values within the dataset

#models based on hypotheses
MDNMMalgaeall<-lm(MicroMacroAlgaeCover ~ MytilusLoss + Volume+TideHeight, data = Mytilusdaynightall)
MDNTempall<-lm(MaxTemp~MytilusLoss +Volume+TideHeight, data = Mytilusdaynightall)
#Mytilusdaynightall$logLight<-sign(Mytilusdaynightall$Light)*log(abs(Mytilusdaynightall$Light+1))
#MDNLightall<-lm(Light ~ MytilusLoss +SAtoVRatio+TideHeight, data = Mytilusdaynightall)
MDNtoPall<-lm(NtoPRatio~ MytilusLoss+Volume+TideHeight,  data = Mytilusdaynightall)
MDNNECall<- lm(NEC~ pH+MaxTemp+TideHeight, data = Mytilusdaynightall) #removed mussel loss since v related to pH&maxtemp
MDNpHall<- lm(pH ~ MytilusLoss+NEP+Volume+TideHeight, data =Mytilusdaynightall)
MDNNEPall<-lm(NEP ~NtoPRatio+MicroMacroAlgaeCover+TideHeight,data = Mytilusdaynightall) 
#used ntoPratio instead of max temp since more background knowledge of how mussels affect nutrient envrionment

qqp(resid(MDNMMalgaeall),"norm") 
#plot(MDNMMalgaeall)
qqp(resid(MDNTempall),"norm") 
#plot(MDNTempall)
qqp(resid(MDNtoPall),"norm") 
#plot(MDNtoPall)
qqp(resid(MDNNECall),"norm") 
#plot(MDNNECall)
qqp(resid(MDNpHall),"norm") 
#plot(MDNpHall)
qqp(resid(MDNNEPall),"norm") 
#plot(MDNNEPall)


ggplot(Mytilusdaynightall,aes(y=NEC,x=MaxTemp))+
  geom_point()+
  geom_smooth(method='lm')

MytilusDNSEMall<-psem(MDNMMalgaeall,
                   MDNTempall,
                   MDNtoPall,
                   MDNNEPall,
                   MDNpHall,
                   MDNNECall)

OutputMytilusSEM<-multigroup(MytilusDNSEMall, group ="Day_Night",standardize = 'scale',test.type = "III")
MytilusModels<-as.data.frame(OutputMytilusSEM$global)
Manovaoutput<-as.data.frame(OutputMytilusSEM$anovaInts)
MDayNightcoeffs<-as.data.frame(OutputMytilusSEM$group.coefs)
MFisher<-as.data.frame(OutputMytilusSEM$Cstat)

write.csv(MytilusModels, 'Output/Daynightmytilusglobalupdated.csv' )
write.csv(Manovaoutput, 'Output/mytilusanovaupdated.csv' )
write.csv(MDayNightcoeffs, 'Output/mytilusDNcoeffupdated.csv' )
write.csv(MFisher, 'Output/mytilusdaynightfishersupdated.csv' )

#####Supplemental Info Figs#####
#significant regression ggplot plots with ggeffects package
####Phyllospadix surfgrass model#####
#1
#Micro/macroalgae and phyllo
PDNMMAlgaeall<-lm(MicroMacroAlgaeCover~ PhyllospadixLoss+Volume+TideHeight,  data = PhylloDayNightall)
pmmagg<-ggpredict(PDNMMAlgaeall, c("PhyllospadixLoss")) #predict marginal effects from model for foundation spp. loss
plot(pmmagg) #plot output 
pmmagg<-as.data.frame(pmmagg) #create dataframe 

pmmagg<-pmmagg %>% #output for values gives you an x for variable. rename variable to match
  rename(PhyllospadixLoss=x) #rename to join to rest of dataframe

pmmagg<-left_join(pmmagg,PhylloDayNightall) #rejoin with main dataframe for ggplot

#display raw data but prediction line and confidence intervals are from ggpredict model
phylloMMA<-ggplot(pmmagg, aes(x =PhyllospadixLoss, y=MicroMacroAlgaeCover)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=PhyllospadixLoss, y=predicted), color="#006d2c",size =2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs( x= '', y='Change in micro/macroalgae % cover') 
phylloMMA

#2
#pH and phyllo
PDNpHall<- lm(pH ~NEP+PhyllospadixLoss+Volume+TideHeight,data = PhylloDayNightall)
ppHphyllogg<-ggpredict(PDNpHall, c("PhyllospadixLoss")) #predict marginal effects from model for foundation spp. loss
ppHphyllogg<-as.data.frame(ppHphyllogg) #create dataframe 

ppHphyllogg<-ppHphyllogg %>% #output for values gives you an x for variable. rename variable to match
  rename(PhyllospadixLoss=x) #rename to join to rest of dataframe

ppHphyllogg<-left_join(ppHphyllogg,PhylloDayNightall) #rejoin with main dataframe for ggplot

phyllopH<-ggplot(ppHphyllogg, aes(x =PhyllospadixLoss, y=pH)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=PhyllospadixLoss, y=predicted), color="#006d2c",size =2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs(x ='Surfgrass % loss \n (Phyllospadix spp.)',y=expression("Average pH"[T]))
phyllopH

#3
#pH and nep because significant day/night interaction
PDNpHall<- lm(pH ~(NEP+PhyllospadixLoss+Volume+TideHeight)*Day_Night,data = PhylloDayNightall)
ppHnepgg<-ggpredict(PDNpHall, c("PhyllospadixLoss","Day_Night")) #predict marginal effects from model for foundation spp. loss
ppHnepgg<-as.data.frame(ppHnepgg) #create dataframe 

ppHnepgg<-ppHnepgg %>% #output for values gives you an x for variable. rename variable to match
  rename(PhyllospadixLoss=x,Day_Night=group) #rename to join to rest of dataframe

ppHnepgg<-left_join(ppHnepgg,PhylloDayNightall) #rejoin with main dataframe for ggplot

pneppH<-ggplot(ppHnepgg, aes(x =PhyllospadixLoss, y=pH,color=Day_Night)) +
  geom_point(size=8,aes(shape=Removal_Control,color=Day_Night),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  scale_colour_manual(values = c("#006d2c",'#bdbdbd'))+
  geom_line(aes(x=PhyllospadixLoss, y=predicted, color=Day_Night),size =2,linetype=2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs(x =expression("Average"~NEP~(mmol~C/m^"2"*hr)),y='')
pneppH

#4
#temp and phyllo
#significant day/night interaction 
PDNTempall<-lm(MaxTemp~ (PhyllospadixLoss+Volume+TideHeight)*Day_Night, data =PhylloDayNightall)
ptempgg<-ggpredict(PDNTempall, c("PhyllospadixLoss", "Day_Night")) #predict marginal effects from model for foundation spp. loss
ptempgg<-as.data.frame(ptempgg) #create dataframe 

ptempgg<-ptempgg%>% #output for values gives you an x for variable. rename variable to match
  rename(PhyllospadixLoss=x,Day_Night=group) #rename to join to rest of dataframe

ptempgg<-left_join(ptempgg,PhylloDayNightall) #rejoin with main dataframe for ggplot

phyllotemp<-ggplot(ptempgg, aes(x =PhyllospadixLoss, y=MaxTemp,color=Day_Night)) +
  geom_point(size=8,aes(shape=Removal_Control,color=Day_Night),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  scale_colour_manual(values = c("#006d2c",'#bdbdbd'))+
  geom_line(aes(x=PhyllospadixLoss, y=predicted, color=Day_Night),size =2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs(x ='Surfgrass % loss \n (Phyllospadix spp.)',y="Average maximum temperature (°C)")
phyllotemp

#5
#Temp and NEP
PDNNEPall<-lm(NEP ~MaxTemp +MicroMacroAlgaeCover+NtoPRatio+TideHeight, data =PhylloDayNightall) 

ptempnep<-ggpredict(PDNNEPall, c("MaxTemp[all]")) #predict marginal effects from model max temp

ptempnep<-as.data.frame(ptempnep) #create dataframe 

ptempnep<-ptempnep%>% #output for values gives you an x for variable. rename variable to match
  rename(MaxTempadj=x) #rename to join to rest of dataframe
PhylloDayNightall$MaxTempadj<-format(round(PhylloDayNightall$MaxTemp, 3), nsmall = 3)  # Apply format function
PhylloDayNightall$MaxTempadj<-as.numeric(PhylloDayNightall$MaxTempadj)

ptempnep<-left_join(ptempnep,PhylloDayNightall) #rejoin with main dataframe for ggplot

ptempnepplot<-ggplot(ptempnep, aes(x =MaxTempadj, y=NEP)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=MaxTempadj, y=predicted), color="#006d2c",size =2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs(x ="Average maximum Temperature (°C)",y=expression("Average"~NEP~(mmol~C/m^"2"*hr)))
ptempnepplot

#6
#N:P and NEP
pntopnep<-ggpredict(PDNNEPall, c("NtoPRatio[all]")) #predict marginal effects from model max temp

pntopnep<-as.data.frame(pntopnep) #create dataframe 

pntopnep<-pntopnep%>% #output for values gives you an x for variable. rename variable to match
  rename(NtoPRatioadj=x) #rename to join to rest of dataframe

PhylloDayNightall$NtoPRatioadj<-format(round(PhylloDayNightall$NtoPRatio, 2), nsmall = 2)  # Apply format function
PhylloDayNightall$NtoPRatioadj<-as.numeric(PhylloDayNightall$NtoPRatioadj) #make ratio smaller to match output of predict model

pntopnep<-left_join(pntopnep,PhylloDayNightall) #rejoin with main dataframe for ggplot

pntopnepplot<-ggplot(pntopnep, aes(x =NtoPRatioadj, y=NEP)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=NtoPRatioadj, y=predicted), color="#006d2c",size =2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs(x ="Average N to P Ratio",y='')
pntopnepplot

#7
#pH and NEC
PDNNECall<- lm(NEC~pH +MaxTemp+TideHeight,data =PhylloDayNightall) #removed surfgrass loss since very correlated with pH and maxTemp

pHnecgg<-ggpredict(PDNNECall, c("pH[all]")) #predict marginal effects from model for foundation spp. loss
pHnecgg<-as.data.frame(pHnecgg) #create dataframe 

pHnecgg<-pHnecgg%>% #output for values gives you an x for variable. rename variable to match
  rename(pHadj=x) #rename to join to rest of dataframe

PhylloDayNightall$pHadj<-format(round(PhylloDayNightall$pH, 3), nsmall = 3)  # Apply format function
PhylloDayNightall$pHadj<-as.numeric(PhylloDayNightall$pHadj) #make ratio smaller to match output of predict model
pHnecgg<-left_join(pHnecgg,PhylloDayNightall) #rejoin with main dataframe for ggplot

pHNEC<-ggplot(pHnecgg, aes(x =pHadj, y=NEC)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=pHadj, y=predicted), color="#006d2c",size =2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs(x =expression("Average pH"[T]),y=expression("Average"~NEC~(mmol~CaCO["3"]/m^"2"*hr)))
pHNEC

#8
#tide height and nec
THnecgg<-ggpredict(PDNNECall, c("TideHeight[all]")) #predict marginal effects from model for foundation spp. loss
THnecgg<-as.data.frame(THnecgg) #create dataframe 

THnecgg<-THnecgg%>% #output for values gives you an x for variable. rename variable to match
  rename(TideHeightadj=x) #rename to join to rest of dataframe

PhylloDayNightall$TideHeightadj<-format(round(PhylloDayNightall$TideHeight, 3), nsmall = 3)  # Apply format function
PhylloDayNightall$TideHeightadj<-as.numeric(PhylloDayNightall$TideHeightadj) #make ratio smaller to match output of predict model
THnecgg<-left_join(THnecgg,PhylloDayNightall) #rejoin with main dataframe for ggplot

THNEC<-ggplot(THnecgg, aes(x =TideHeightadj, y=NEC)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=TideHeightadj, y=predicted), color="#006d2c",size =2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs(x ='Tide height (m)',y='')
THNEC

#patchwork everything together
phyllosig<-(phyllotemp | phylloMMA)/
  (phyllopH | pneppH)/
  (ptempnepplot | pntopnepplot)/
  (pHNEC|THNEC) +
  plot_annotation(tag_levels = 'a') &         #label each individual plot with letters A-G
  theme(plot.tag = element_text(size = 50))   #edit the lettered text
phyllosig
#ggsave(filename = "Output/phyllosigsem.pdf", useDingbats =FALSE,dpi=600,device = "pdf", width = 40, height = 45)

####Mytilus mussel SEM#####
#1
#micro/macroalage cover and mytilus
MDNMMalgaeall<-lm(MicroMacroAlgaeCover ~ MytilusLoss + Volume+TideHeight, data = Mytilusdaynightall)

mytmmagg<-ggpredict(MDNMMalgaeall, c("MytilusLoss")) #predict marginal effects from model for foundation spp. loss
mytmmagg<-as.data.frame(mytmmagg) #create dataframe 

mytmmagg<-mytmmagg%>% #output for values gives you an x for variable. rename variable to match
  rename(MytilusLoss=x) #rename to join to rest of dataframe

mytmmagg<-left_join(mytmmagg,Mytilusdaynightall) #rejoin with main dataframe for ggplot

mytmma<-ggplot(mytmmagg, aes(x =MytilusLoss, y=MicroMacroAlgaeCover)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=MytilusLoss, y=predicted), color="#045a8d",size =2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs(x ='CA mussel % loss \n (Mytilus californianus)',y='Change in micro/macroalgae % cover')
mytmma

#2
#mma and tideheight
MDNMMalgaeall<-lm(MicroMacroAlgaeCover ~ MytilusLoss + Volume+TideHeight, data = Mytilusdaynightall)
thmmagg<-ggpredict(MDNMMalgaeall, c("TideHeight")) #predict marginal effects from model for foundation spp. loss
thmmagg<-as.data.frame(thmmagg) #create dataframe 

thmmagg<-thmmagg%>% #output for values gives you an x for variable. rename variable to match
  rename(TideHeight=x) #rename to join to rest of dataframe

thmmagg<-left_join(thmmagg,Mytilusdaynightall) #rejoin with main dataframe for ggplot
thmma<-ggplot(thmmagg, aes(x =TideHeight, y=MicroMacroAlgaeCover)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=TideHeight, y=predicted), color="#045a8d",size =2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs(x='Tide height (m)',y='')
thmma

#3
#mytilus and pH
MDNpHall<- lm(pH ~ MytilusLoss+NEP+Volume+TideHeight, data =Mytilusdaynightall)

mytpHgg<-ggpredict(MDNpHall, c("MytilusLoss")) #predict marginal effects from model for foundation spp. loss
mytpHgg<-as.data.frame(mytpHgg) #create dataframe 

mytpHgg<-mytpHgg%>% #output for values gives you an x for variable. rename variable to match
  rename(MytilusLoss=x) #rename to join to rest of dataframe

mytpHgg<-left_join(mytpHgg,Mytilusdaynightall) #rejoin with main dataframe for ggplot
mytpH<-ggplot(mytpHgg, aes(x =MytilusLoss, y=pH)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=MytilusLoss, y=predicted), color="#045a8d",size =2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs(x ='CA mussel % loss \n (Mytilus californianus)',y=expression("Average pH"[T]))
mytpH

#4
#nep and pH
mneppHgg<-ggpredict(MDNpHall, c("NEP[all]")) #predict marginal effects from model for foundation spp. loss
mneppHgg<-as.data.frame(mneppHgg) #create dataframe 

mneppHgg<-mneppHgg%>% #output for values gives you an x for variable. rename variable to match
  rename(NEPadj=x) #rename to join to rest of dataframe
Mytilusdaynightall$NEPadj<-format(round(Mytilusdaynightall$NEP, 2), nsmall = 2)  # Apply format function
Mytilusdaynightall$NEPadj<-as.numeric(Mytilusdaynightall$NEPadj) #make ratio smaller to match output of predict model
mneppHgg<-left_join(mneppHgg,Mytilusdaynightall) #rejoin with main dataframe for ggplot

mytpHnep<-ggplot(mneppHgg, aes(x =NEPadj, y=pH)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=NEPadj, y=predicted), color="#045a8d",size =2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs(x=expression("Average"~NEP~(mmol~C/m^"2"*hr)),y='')
mytpHnep

#5
#nep and mma
#signficant day/night interaction
MDNNEPallint<-lm(NEP ~(NtoPRatio+MicroMacroAlgaeCover+TideHeight)*Day_Night,data = Mytilusdaynightall) 

mmanepgg<-ggpredict(MDNNEPallint, c("MicroMacroAlgaeCover","Day_Night")) #predict marginal effects from model for foundation spp. loss
mmanepgg<-as.data.frame(mmanepgg) #create dataframe 

mmanepgg<-mmanepgg%>% #output for values gives you an x for variable. rename variable to match
  rename(MicroMacroAlgaeCover=x,Day_Night=group) #rename to join to rest of dataframe

mmanepgg<-left_join(mmanepgg,Mytilusdaynightall) #rejoin with main dataframe for ggplot
mmanep<-ggplot(mmanepgg, aes(y =NEP, x=MicroMacroAlgaeCover,color=Day_Night)) +
  geom_point(size=8,aes(shape=Removal_Control,color=Day_Night),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  scale_colour_manual(values = c("#045a8d",'#bdbdbd'))+
  geom_line(aes(x=MicroMacroAlgaeCover, y=predicted, color=Day_Night,linetype=Day_Night),size =2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs(y=expression("Average"~NEP~(mmol~C/m^"2"*hr)),x='Change in micro/macroalgae % cover')
mmanep

#6
#nep and n to p
MDNNEPall<-lm(NEP ~NtoPRatio+MicroMacroAlgaeCover+TideHeight,data = Mytilusdaynightall) 
ntpnepgg<-ggpredict(MDNNEPall, c("NtoPRatio[all]")) #predict marginal effects from model for foundation spp. loss
ntpnepgg<-as.data.frame(ntpnepgg) #create dataframe 

ntpnepgg<-ntpnepgg%>% #output for values gives you an x for variable. rename variable to match
  rename(NtoPRatioedited=x) #rename to join to rest of dataframe

Mytilusdaynightall$NtoPRatioedited<-format(round(Mytilusdaynightall$NtoPRatio, 2), nsmall = 2)  # Apply format function
Mytilusdaynightall$NtoPRatioedited<-as.numeric(Mytilusdaynightall$NtoPRatioedited) #make ratio smaller to match output of predict model
ntpnepgg<-left_join(ntpnepgg,Mytilusdaynightall) #rejoin with main dataframe for ggplot

ntpnep<-ggplot(ntpnepgg, aes(y =NEP, x=NtoPRatioedited)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=NtoPRatioedited, y=predicted), color="#045a8d",size =2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs(x='Average N to P Ratio',y='')
ntpnep

#7
#nec and max temp
MDNNECall<- lm(NEC~ pH+MaxTemp+TideHeight, data = Mytilusdaynightall) #removed mussel loss since v related to pH&maxtemp
mtempnecgg<-ggpredict(MDNNECall, c("MaxTemp[all]")) #predict marginal effects from model for foundation spp. loss
mtempnecgg<-as.data.frame(mtempnecgg) #create dataframe 

mtempnecgg<-mtempnecgg%>% #output for values gives you an x for variable. rename variable to match
  rename(MaxTempedited=x) #rename to join to rest of dataframe

Mytilusdaynightall$MaxTempedited<-format(round(Mytilusdaynightall$MaxTemp, 3), nsmall = 3)  # Apply format function
Mytilusdaynightall$MaxTempedited<-as.numeric(Mytilusdaynightall$MaxTempedited) #make ratio smaller to match output of predict model

mtempnecgg<-left_join(mtempnecgg,Mytilusdaynightall) #rejoin with main dataframe for ggplot
mtempnec<-ggplot(mtempnecgg, aes(x =MaxTempedited, y=NEC)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=MaxTempedited, y=predicted), color="#045a8d",size =2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs(x='Average maximum temperature (°C)',y=expression("Average"~NEC~(mmol~CaCO["3"]/m^"2"*hr)))
mtempnec

#patchwork mytilus figs
mytilussig<-(mytmma | thmma)/
  (mmanep | ntpnep) /
  (mytpH | mytpHnep |mtempnec)+
  plot_annotation(tag_levels = 'a') &         #label each individual plot with letters A-G
  theme(plot.tag = element_text(size = 50))   #edit the lettered text
mytilussig
ggsave(filename = "Output/mytilussigsem.pdf", useDingbats =FALSE,dpi=600,device = "pdf", width = 40, height = 45)

####Light and temp correlation#####

#supplemental figure
#Light and temp correlation in SEM
PhylloDay<-PhylloDayNightall%>%
  filter(Day_Night=='Day')
cor.test(PhylloDay$MaxTemp, PhylloDay$Light,  method = "pearson")
#t = 6.1464, df = 14, p-value = 2.534e-05
#95 percent confidence interval:
 #0.6217138 0.9483361
#sample estimates: cor 0.8541739 
Mytilusday<-Mytilusdaynightall%>%
  filter(Day_Night=='Day')
cor.test(Mytilusday$MaxTemp,Mytilusday$Light,  method = "pearson")
#t = 5.7644, df = 13, p-value = 6.551e-05
#95 percent confidence interval: 0.5931612 0.9482484
#sample estimates:cor 0.8478125 
#graphs
phyllolightandtemp<-ggplot(PhylloDay, aes(y =MaxTemp, x=Light)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_smooth(method = 'lm',color ="#006d2c",size =1.5,alpha =0.3)+
  theme_sleek() +
  theme(axis.title.y=element_text(color="black", size=45), 
        axis.title.x=element_text(color="black", size=40),
        axis.text.x =element_text(color="black", size=35),
        axis.text.y =element_text(color="black", size=35)) +
  theme(legend.position="none")+
  labs(y='Change in daily max temp (°C)',x='')

mytiluslightandtemp<-ggplot(Mytilusday, aes(y =MaxTemp, x=Light)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_smooth(method = 'lm',color ="#045a8d",size =1.5,alpha =0.3)+
  theme_sleek() +
  theme(axis.title.x=element_text(face="italic", color="black", size=45), 
        axis.title.y=element_blank(),
        axis.text.x =element_text(color="black",size=35),
        axis.text.y =element_text(color="black", size=35)) +
  theme(legend.position="none")+
  labs(x=expression('Average change in light' (PFD~µmol~photons~m^{-2}~s^{-1})), y='')

lighttempsem<-phyllolightandtemp+mytiluslightandtemp +      #patchwork to combine plots
  plot_annotation(tag_levels = 'a') &         #label each individual plot with letters A-G
  theme(plot.tag = element_text(size = 50))   #edit the lettered text

lighttempsem
#ggsave(filename = "Output/SEMsuppLightandTempgraphs.pdf", useDingbats =FALSE,dpi=300,device = "pdf", width = 30, height = 25)

