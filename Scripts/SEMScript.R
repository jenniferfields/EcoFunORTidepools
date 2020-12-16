## Data processing and Structural equation model for both surfgrass and mussel models
## By: Jenn Fields
## Last updated: 10.23.2020
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
  mutate(NtoP = (NH4_umol_L+NN_umol_L)/PO_umol_L) %>% #group by time point to get n:p ratio per tide pool per time point
  dplyr::group_by(PoolID, Foundation_spp, Removal_Control,Day_Night,Before_After)  %>%
  summarise(AvgpH =mean(pH_insitu),
            AvgNtoP =mean(NtoP),
            AvgNN=mean(NN_umol_L))


#summarise average max temp over the low tide sample (n=3-4 values per tide pool/time period)
SummaryTemp<-TempandLightSumall%>%
  dplyr::group_by(PoolID, Day_Night,Before_After)  %>%
  summarise(AvgTempmax = mean(Temp.max),
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

PP<-SEMall%>%
  dplyr::group_by(Foundation_spp,PoolID,Removal_Control,Day_Night) %>%
  summarise(SAVav = mean(SAtoV),THav = mean(TideHeight),SAav=mean(SurfaceArea),Vav=mean(Vol)) 


SEMcombined<-left_join(SummaryNEPNECcommpp,Sumbiogeochemtemp)
SEMcombined<-as.data.frame(SEMcombined)

#Takes the delta between after-before period
SEMallavg<-SEMcombined%>%
  dplyr::group_by(PoolID, Foundation_spp, Removal_Control,Day_Night)  %>%
  dplyr::summarise(NEPdelta = AvgNEP[Before_After == 'After'] - AvgNEP[Before_After == 'Before'],
                   NECdelta = AvgNEC[Before_After == 'After'] - AvgNEC[Before_After == 'Before'],
                   pHdelta= AvgpH[Before_After == 'After'] - AvgpH[Before_After == 'Before'],
                   NtoPdelta = AvgNtoP[Before_After == 'After'] - AvgNtoP[Before_After == 'Before'],
                   NNdelta = AvgNN[Before_After == 'After'] - AvgNN[Before_After == 'Before'],
                   Tempmaxdelta = AvgTempmax[Before_After == 'After'] - AvgTempmax[Before_After == 'Before'],
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

SEMallavg%>%
  filter(Foundation_spp=="Phyllospadix"&PoolID !='26') %>%
  ggplot(aes(x=Tempmaxdelta,y=NNdelta,color=Day_Night))+
    geom_point()+
    geom_smooth(method="lm")

SEMallavg%>%
  filter(Foundation_spp=="Mytilus") %>%
  ggplot(aes(x=Tempmaxdelta,y=NNdelta,color=Day_Night,shape=Removal_Control))+
  geom_point()+
  geom_smooth(method="lm")
####Surfgrass SEM####
PhylloDayNightall<-SEMallavg %>%
  filter(Foundation_spp == 'Phyllospadix') %>%
  dplyr::rename(NEP =NEPdelta, NEC=NECdelta, pH=pHdelta, NtoPRatio=NtoPdelta,
                MaxTemp = Tempmaxdelta,Light=Parmeandelta,MytilusLoss=Mytilusdelta,PhyllospadixLoss=Phyllodelta, MicroMacroAlgaeCover=micromacroalgaedelta,
                SAtoVRatio=SAVav,TideHeight=THav,SurfaceArea=SA,Volume=V, NN=NNdelta) #rename cols to match sem

###problem children round 1000XX
ggplot(noTP26, aes(x=Volume,y=MaxTemp, color=Day_Night))+
  geom_point()+
  geom_smooth(method="lm")

ggplot(noTP26, aes(x=MicroMacroAlgaeCover,y=NEP))+
  geom_point()+
  geom_smooth(method="lm")


noTP26<-PhylloDayNightall%>%
  filter(PoolID != 26) #removed tp 26 because outlier within N:P data (+4SD away) 
noTP26<-as.data.frame(noTP26)
ggpairs(noTP26[c(5:11,14:15,18:21)])
##max temp and sa:v
#surfgrass loss and sa to v ratio?

#surfgrass loss and sa to v ratio?
##max temp and sa:v
#max temp and pH
#sav and pH

#make day/night a factor for multigroup analysis
noTP26$Day_Night<-as.factor(noTP26$Day_Night)

#models based off hypotheses
PDNMMAlgaeall<-lm(MicroMacroAlgaeCover~ PhyllospadixLoss+Volume+TideHeight,  data = noTP26)
PDNTempall<-lm(MaxTemp~ PhyllospadixLoss+Volume+TideHeight, data = noTP26)
PDNNtoPall<-lm(NtoPRatio ~PhyllospadixLoss+Volume +TideHeight, data =noTP26)
PDNNECall<- lm(NEC~pH +MaxTemp+TideHeight,data =noTP26)
PDNpHall<- lm(pH ~NEP+PhyllospadixLoss+Volume+TideHeight,data = noTP26)
PDNNEPall<-lm(NEP ~MaxTemp +MicroMacroAlgaeCover+NtoPRatio+TideHeight, data =noTP26) 


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

Mytilusdaynightall$Day_Night<-as.factor(Mytilusdaynightall$Day_Night)

#linear model with tide height and mussel loss 
#get resids of tide height (so the effect of tide height above mussel loss since they are correlated)
MusselTideHeight<-lm(TideHeight~MytilusLoss,data=Mytilusdaynightall)
#summary(MusselTideHeight) #correlated
TH.resid<-resid(MusselTideHeight)
Mytilusdaynightall$TH.resid<-TH.resid
x<-lm(TH.resid~MytilusLoss,data=Mytilusdaynightall)
#summary(x) #not correlated

ggpairs(Mytilusdaynightall[c(5:11,13,15,18:21)])


#models based on hypotheses
MDNMMalgaeall<-lm(MicroMacroAlgaeCover ~ MytilusLoss + Volume+TH.resid, data = Mytilusdaynightall)
MDNTempall<-lm(MaxTemp~MytilusLoss +Volume+TH.resid , data = Mytilusdaynightall)
#Mytilusdaynightall$logLight<-sign(Mytilusdaynightall$Light)*log(abs(Mytilusdaynightall$Light+1))
#MDNLightall<-lm(Light ~ MytilusLoss +SAtoVRatio+TideHeight, data = Mytilusdaynightall)
MDNtoPall<-lm(NtoPRatio ~ MytilusLoss+Volume+TH.resid,  data = Mytilusdaynightall)
MDNNECall<- lm(NEC~ pH+MaxTemp+TH.resid, data = Mytilusdaynightall) #removed mussel loss since v related to pH
MDNpHall<- lm(pH ~ MytilusLoss+NEP+ Volume+TH.resid, data =Mytilusdaynightall)
MDNNEPall<-lm(NEP ~NtoPRatio+MicroMacroAlgaeCover+TH.resid,data = Mytilusdaynightall) 

Mdaynep<-Mytilusdaynightall%>%
  filter(Day_Night=="Day")
mytilusnep<-lm(NEP~Light+MicroMacroAlgaeCover+NtoPRatio+SAtoVRatio +TideHeight,data=Mdaynep)
summary(mytilusnep)


qqp(resid(MDNMMalgaeall),"norm") 
qqp(resid(MDNTempall),"norm") 
#qqp(resid(MDNLightall),"norm") #log light data for normality
#plot(MDNLightall) #one outlier with log light data
qqp(resid(MDNtoPall),"norm") 
qqp(resid(MDNNECall),"norm") 
qqp(resid(MDNpHall),"norm") 
qqp(resid(MDNNEPall),"norm") 


MytilusDNSEMall<-psem(MDNMMalgaeall,
                   #MDNLightall,
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
#Micro/macroalgae and phyllo

pmmagg<-ggpredict(PAvgMMAlgae, c("PhyllospadixLoss")) #predict marginal effects from model for foundation spp. loss
plot(pmmagg) #plot output 
pmmagg<-as.data.frame(pmmagg) #create dataframe 

pmmagg<-pmmagg %>% #output for values gives you an x for variable. rename variable to match
  rename(PhyllospadixLoss=x) #rename to join to rest of dataframe

pmmagg<-left_join(pmmagg,Phyllounscaled) #rejoin with main dataframe for ggplot

#display raw data but prediction line and confidence intervals are from ggpredict model
phylloMMA<-ggplot(pmmagg, aes(x =PhyllospadixLoss, y=MicroMacroAlgaeCover)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=PhyllospadixLoss, y=predicted), color="#006d2c",size =1.5)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs( x= '', y='Change micro/macroalgae cover') 
phylloMMA

#light and phyllo
plightgg<-ggpredict(PAvgLight, c("PhyllospadixLoss")) #predict marginal effects from model for foundation spp. loss
plightgg<-as.data.frame(plightgg) #create dataframe 

plightgg<-plightgg %>% #output for values gives you an x for variable. rename variable to match
  rename(PhyllospadixLoss=x) #rename to join to rest of dataframe

plightgg<-left_join(plightgg,Phyllounscaled) #rejoin with main dataframe for ggplot

phyllolight<-ggplot(plightgg, aes(x =PhyllospadixLoss, y=Light)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=PhyllospadixLoss, y=predicted), color="#006d2c",size =1.5)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs(x ='Surfgrass loss \n (Phyllospadix spp.)',y=expression('Average light' (PFD~µmol~photons~m^{-2}~s^{-1})))
phyllolight

#temp and phyllo
ptempgg<-ggpredict(PAvgTemp, c("PhyllospadixLoss")) #predict marginal effects from model for foundation spp. loss
ptempgg<-as.data.frame(ptempgg) #create dataframe 

ptempgg<-ptempgg%>% #output for values gives you an x for variable. rename variable to match
  rename(PhyllospadixLoss=x) #rename to join to rest of dataframe

ptempgg<-left_join(ptempgg,Phyllounscaled) #rejoin with main dataframe for ggplot

phyllotemp<-ggplot(ptempgg, aes(x =PhyllospadixLoss, y=MaxTemp)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=PhyllospadixLoss, y=predicted), color="#006d2c",size =1.5)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs(x ='Surfgrass loss \n (Phyllospadix spp.)',y="Maximum Temperature (°C)")
phyllotemp

#SA/V and pH
savphgg<-ggpredict(PAvgpH, c("SAtoVRatio")) #predict marginal effects from model for foundation spp. loss
savphgg<-as.data.frame(savphgg) #create dataframe 

savphgg<-savphgg %>% #output for values gives you an x for variable. rename variable to match
  rename(SAtoVRatio=x) #rename to join to rest of dataframe

savphgg<-left_join(savphgg,Phyllounscaled) #rejoin with main dataframe for ggplot

SAVpH<-ggplot(savphgg, aes(x =SAtoVRatio, y=pH)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=SAtoVRatio, y=predicted), color="#006d2c",size =1.5)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=45), 
        axis.title.y=element_text(color="black", size=45),
        axis.text.x =element_text(color="black", size=30),
        axis.text.y =element_text(color="black", size=30)) +
  theme(legend.position="none")+
  labs(x ='SA to V ratio',y='') 
SAVpH

#NEP and pH
phyllophgg<-ggpredict(PAvgpH, c("PhyllospadixLoss")) #predict marginal effects from model for foundation spp. loss
phyllophgg<-as.data.frame(phyllophgg) #create dataframe 

phyllophgg<-phyllophgg %>% #output for values gives you an x for variable. rename variable to match
  rename(PhyllospadixLoss=x) #rename to join to rest of dataframe

phyllophgg<-left_join(phyllophgg,Phyllounscaled) #rejoin with main dataframe for ggplot

phyllopH<-ggplot(phyllophgg, aes(x =PhyllospadixLoss, y=pH)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=PhyllospadixLoss, y=predicted), color="#006d2c",size =1.5)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=45), 
        axis.title.y=element_text(color="black", size=45),
        axis.text.x =element_text(color="black", size=30),
        axis.text.y =element_text(color="black", size=30)) +
  theme(legend.position="none")+
  labs(y ='pH',x='Surfgrass loss \n (Phyllospadix spp.)')
phyllopH

#SA to V and NEC
savnecgg<-ggpredict(PAvgNEC, c("SAtoVRatio")) #predict marginal effects from model for foundation spp. loss
savnecgg<-as.data.frame(savnecgg) #create dataframe 

savnecgg<-savnecgg%>% #output for values gives you an x for variable. rename variable to match
  rename(SAtoVRatio=x) #rename to join to rest of dataframe

savnecgg<-left_join(savnecgg,Phyllounscaled) #rejoin with main dataframe for ggplot
SAVNEC<-ggplot(savnecgg, aes(x =SAtoVRatio, y=NEC)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=SAtoVRatio, y=predicted), color="#006d2c",size =1.5)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40))+
  theme(legend.position="none")+
  labs(x ='SA to V ratio',y='')
SAVNEC

#temp and nec
tempnecgg<-ggpredict(PAvgNEC, c("MaxTemp")) #predict marginal effects from model for foundation spp. loss
tempnecgg<-as.data.frame(tempnecgg) #create dataframe 

tempnecgg<-tempnecgg%>% #output for values gives you an x for variable. rename variable to match
  rename(MaxTemp=x) #rename to join to rest of dataframe

tempnecgg<-left_join(tempnecgg,Phyllounscaled) #rejoin with main dataframe for ggplot
tempNEC<-ggplot(tempnecgg, aes(x =MaxTemp, y=NEC)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=MaxTemp, y=predicted), color="#006d2c",size =1.5)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40))+
  theme(legend.position="none")+
  labs(x ='Maximum temperature (°C)',y='')
tempNEC

#pH and NEC

pHnecgg<-ggpredict(PAvgNEC, c("pH")) #predict marginal effects from model for foundation spp. loss
pHnecgg<-as.data.frame(pHnecgg) #create dataframe 

pHnecgg<-pHnecgg%>% #output for values gives you an x for variable. rename variable to match
  rename(pH=x) #rename to join to rest of dataframe

pHnecgg<-left_join(pHnecgg,Phyllounscaled) #rejoin with main dataframe for ggplot
pHNEC<-ggplot(pHnecgg, aes(x =pH, y=NEC)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=pH, y=predicted), color="#006d2c",size =1.5)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs(x ='pH',y=expression(Average~NEC~(mmol~CaCO["3"]/m^"2"*hr)))
pHNEC

#tide height and nec
THnecgg<-ggpredict(PAvgNEC, c("TideHeight")) #predict marginal effects from model for foundation spp. loss
THnecgg<-as.data.frame(THnecgg) #create dataframe 

THnecgg<-THnecgg%>% #output for values gives you an x for variable. rename variable to match
  rename(TideHeight=x) #rename to join to rest of dataframe

THnecgg<-left_join(THnecgg,Phyllounscaled) #rejoin with main dataframe for ggplot
THNEC<-ggplot(THnecgg, aes(x =TideHeight, y=NEC)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=TideHeight, y=predicted), color="#006d2c",size =1.5)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs(x ='Tide height (m)',y='')
THNEC

#Light and NEP
lightnepgg<-ggpredict(PAvgNEP, c("Light")) #predict marginal effects from model for foundation spp. loss
lightnepgg<-as.data.frame(lightnepgg) #create dataframe 

lightnepgg<-lightnepgg%>% #output for values gives you an x for variable. rename variable to match
  rename(Light=x) #rename to join to rest of dataframe

lightnepgg<-left_join(lightnepgg,Phyllounscaled) #rejoin with main dataframe for ggplot
lightnep<-ggplot(lightnepgg, aes(x =Light, y=NEP)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=Light, y=predicted), color="#006d2c",size =1.5)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs(x =expression('Average light' (PFD~µmol~photons~m^{-2}~s^{-1})),y=expression(Average~NEP~(mmol~C/m^"2"*hr)))
lightnep

#patchwork everything together
phyllosig<-(phyllolight | phylloMMA | lightnep)/
  (phyllopH|SAVpH)/
  (pHNEC|tempNEC|THNEC|SAVNEC) +
  plot_annotation(tag_levels = 'a') &         #label each individual plot with letters A-G
  theme(plot.tag = element_text(size = 50, face = "bold"))   #edit the lettered text
phyllosig
ggsave(filename = "Output/phyllosigsem.pdf", useDingbats =FALSE,dpi=600,device = "pdf", width = 35, height = 45)

####Mytilus mussel SEM#####

#mma and mytilus
mytmmagg<-ggpredict(MAvgMMalgae, c("MytilusLoss")) #predict marginal effects from model for foundation spp. loss
mytmmagg<-as.data.frame(mytmmagg) #create dataframe 

mytmmagg<-mytmmagg%>% #output for values gives you an x for variable. rename variable to match
  rename(MytilusLoss=x) #rename to join to rest of dataframe

mytmmagg<-left_join(mytmmagg,Mytilusunscaled) #rejoin with main dataframe for ggplot
mytmma<-ggplot(mytmmagg, aes(x =MytilusLoss, y=MicroMacroAlgaeCover)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=MytilusLoss, y=predicted), color="#045a8d",size =1.5)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs(x ='CA mussel loss \n (Mytilus californianus)',y='Change in micro/macroalgae cover')
mytmma

#mytilus and pH

mytpHgg<-ggpredict(MAvgpH, c("MytilusLoss")) #predict marginal effects from model for foundation spp. loss
mytpHgg<-as.data.frame(mytpHgg) #create dataframe 

mytpHgg<-mytpHgg%>% #output for values gives you an x for variable. rename variable to match
  rename(MytilusLoss=x) #rename to join to rest of dataframe

mytpHgg<-left_join(mytpHgg,Mytilusunscaled) #rejoin with main dataframe for ggplot
mytpH<-ggplot(mytpHgg, aes(x =MytilusLoss, y=pH)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=MytilusLoss, y=predicted), color="#045a8d",size =1.5)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs(x ='',y='pH')

#nep and mma
mmanepgg<-ggpredict(MAvgNEP, c("MicroMacroAlgaeCover")) #predict marginal effects from model for foundation spp. loss
mmanepgg<-as.data.frame(mmanepgg) #create dataframe 

mmanepgg<-mmanepgg%>% #output for values gives you an x for variable. rename variable to match
  rename(MicroMacroAlgaeCover=x) #rename to join to rest of dataframe

mmanepgg<-left_join(mmanepgg,Mytilusunscaled) #rejoin with main dataframe for ggplot
mmanep<-ggplot(mmanepgg, aes(y =NEP, x=MicroMacroAlgaeCover)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=MicroMacroAlgaeCover, y=predicted), color="#045a8d",size =1.5)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs(y=expression(Average~NEP~(mmol~C/m^"2"*hr)),x='Change in micro/macroalgae cover')


#mma and tideheight
thmmagg<-ggpredict(MAvgMMalgae, c("TideHeight")) #predict marginal effects from model for foundation spp. loss
thmmagg<-as.data.frame(thmmagg) #create dataframe 

thmmagg<-thmmagg%>% #output for values gives you an x for variable. rename variable to match
  rename(TideHeight=x) #rename to join to rest of dataframe

thmmagg<-left_join(thmmagg,Mytilusunscaled) #rejoin with main dataframe for ggplot
thmma<-ggplot(thmmagg, aes(x =TideHeight, y=MicroMacroAlgaeCover)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=TideHeight, y=predicted), color="#045a8d",size =1.5)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs(x='Tide height (m)',y='Change in micro/macroalgae cover')

#nutrients and sa to v ratio
nutsavgg<-ggpredict(MAvgNtoP, c("SAtoVRatio")) #predict marginal effects from model for foundation spp. loss
nutsavgg<-as.data.frame(nutsavgg) #create dataframe 

nutsavgg<-nutsavgg%>% #output for values gives you an x for variable. rename variable to match
  rename(SAtoVRatio=x) #rename to join to rest of dataframe

nutsavgg<-left_join(nutsavgg,Mytilusunscaled) #rejoin with main dataframe for ggplot
nutsav<-ggplot(nutsavgg, aes(x =SAtoVRatio, y=NtoPRatio)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=SAtoVRatio, y=predicted), color="#045a8d",size =1.5)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs(x='SA to V ratio',y='N to P ratio')

#patchwork mytilus figs
mytilussig<-(mytmma|mytpH)/
  (mmanep|thmma|nutsav) +
  plot_annotation(tag_levels = 'a') &         #label each individual plot with letters A-G
  theme(plot.tag = element_text(size = 50, face = "bold"))   #edit the lettered text
mytilussig
ggsave(filename = "Output/mytilussigsem.pdf", useDingbats =FALSE,dpi=600,device = "pdf", width = 35, height = 45)

####Light and temp correlation#####
#supplemental figure
#Light and temp correlation in SEM
cor.test(Phyllounscaled$MaxTemp, Phyllounscaled$Light,  method = "pearson")
#t = 4.2129, df = 14, p-value = 0.0008685
#95 percent confidence interval:0.4003722 0.9071667\
#cor 0.7476867 
cor.test(Mytilusunscaled$MaxTemp, Mytilusunscaled$Light,  method = "pearson")
#t = 5.842, df = 13, p-value = 5.761e-05
#95 percent confidence interval:0.6004721 0.9493811
#cor 0.8509753 
#graphs
phyllolightandtemp<-ggplot(Phyllounscaled, aes(x =MaxTemp, y=Light)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_smooth(method = 'lm',color ="#006d2c",size =1.5,alpha =0.3)+
  theme_sleek() +
  theme(axis.title.x=element_text(face="italic", color="black", size=40), 
        axis.title.y=element_text(color="black", size=40),
        axis.text.x =element_text(color="black", size=35),
        axis.text.y =element_text(color="black", size=35)) +
  theme(legend.position="none")+
  labs( x= '', y='Change in daily max temp (°C)')

mytiluslightandtemp<-ggplot(Mytilusunscaled, aes(x =MaxTemp, y=Light)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_smooth(method = 'lm',color ="#045a8d",size =1.5,alpha =0.3)+
  theme_sleek() +
  theme(axis.title.x=element_text(face="italic", color="black", size=40), 
        axis.title.y=element_blank(),
        axis.text.x =element_text(color="black",size=35),
        axis.text.y =element_text(color="black", size=35)) +
  theme(legend.position="none")+
  labs(x=expression('Average change in light' (PFD~µmol~photons~m^{-2}~s^{-1})), y='Change in daily max temp (°C)')

lighttempsem<-phyllolightandtemp+mytiluslightandtemp +      #patchwork to combine plots
  plot_annotation(tag_levels = 'a') &         #label each individual plot with letters A-G
  theme(plot.tag = element_text(size = 50, face = "bold"))   #edit the lettered text

lighttempsem
ggsave(filename = "Output/SEMsuppLightandTempgraphs.pdf", useDingbats =FALSE,dpi=300,device = "pdf", width = 20, height = 20)

