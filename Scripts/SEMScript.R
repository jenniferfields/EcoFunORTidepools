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
devtools::install_github("seananderson/ggsidekick")
library(ggsidekick) #for theme sleek
library(tidyverse)
library(GGally)

#source code scripts
source("scripts/tidepoolphysicalparameters.R")
source("scripts/CleanCarbChem.R")
source("scripts/CommunityComp.R")
source("scripts/plot.psemfun.R") #plot draft SEM
#change poolid to character

SEMcommunitydata$PoolID<-as.character(SEMcommunitydata$PoolID)

#combine community comp with biogeochem
SEMdata<-left_join(DeltaSamples,SEMcommunitydata)

SEMdata<-as.data.frame(SEMdata)

#Calculate N:P ratio since collinearity between nutrient spp
SEMdata$NtoP<-(SEMdata$NNmean +SEMdata$NH4mean)/SEMdata$POmean

#only average day par values
SEMDataDL<-SEMdata%>%
  filter(Day_Night =="Day") %>%
  select(PoolID, Foundation_spp, Removal_Control, Before_After, Par.mean,Par.max) %>%
  group_by(PoolID, Foundation_spp, Removal_Control, Before_After) %>%
  summarise(AvgPar = mean(Par.mean), AvgParmax = mean(Par.max))

SEMdaynightAvg<- SEMdata %>%
  dplyr::group_by(PoolID, Foundation_spp, Removal_Control, Before_After)  %>%
  dplyr::summarise(AvgNEP = mean(NEP.mmol.m2.hr),
                   AvgNEC = mean(NEC.mmol.m2.hr),
                   AvgpH = mean(pHmean),
                   AvgNN = mean(NNmean),
                   AvgNH4 = mean(NH4mean),
                   AvgPO = mean(POmean),
                   AvgNtoP = mean(NtoP),
                   AvgTempmean = mean(Temp.mean),
                   AvgTempmax = mean(Temp.max),
                   AvgMytilus = mean(MusselCover),
                   AvgPhyllo = mean(SurfgrassCover),
                   AvgFleshyalgae = mean(macroalgae),
                   AvgCCA = mean(allCCA),
                   Avgproddom = mean(prodphyllodom),
                   AvgSAV = mean(SAtoV),
                   AvgTH = mean(TideHeight)) 

SEMdaynightAvg<-left_join(SEMdaynightAvg,SEMDataDL)

SEMdaynightAvg<-SEMdaynightAvg%>%
  dplyr::group_by(PoolID, Foundation_spp, Removal_Control)  %>%
  dplyr::summarise(NEPdeltaAvg = AvgNEP[Before_After == 'After'] - AvgNEP[Before_After == 'Before'],
                   NECdeltaAvg = AvgNEC[Before_After == 'After'] - AvgNEC[Before_After == 'Before'],
                   pHdeltaAvg = AvgpH[Before_After == 'After'] - AvgpH[Before_After == 'Before'],
                   NNdeltaAvg = AvgNN[Before_After == 'After'] - AvgNN[Before_After == 'Before'],
                   NH4deltaAvg = AvgNH4[Before_After == 'After'] - AvgNH4[Before_After == 'Before'],
                   POdeltaAvg = AvgPO[Before_After == 'After'] - AvgPO[Before_After == 'Before'],
                   NtoPdeltaAvg = AvgNtoP[Before_After == 'After'] - AvgNtoP[Before_After == 'Before'],
                   TempmeandeltaAvg = AvgTempmean[Before_After == 'After'] - AvgTempmean[Before_After == 'Before'],
                   Tempmaxdeltaavg = AvgTempmax[Before_After == 'After'] - AvgTempmax[Before_After == 'Before'],
                   Parmaxdeltaavg = AvgParmax[Before_After == 'After'] -AvgParmax[Before_After == 'Before'],
                   Percentmaxlightavg =100 * (AvgParmax[Before_After == 'After'] - AvgParmax[Before_After == 'Before']) /  AvgParmax[Before_After == 'Before'],
                   ParmeandeltaAvg = AvgPar[Before_After == 'After'] - AvgPar[Before_After == 'Before'],
                   Mytiluscover =AvgMytilus[Before_After == 'After'] - AvgMytilus[Before_After == 'Before'],
                   MytilusdeltaAvg = -1*(AvgMytilus[Before_After == 'After'] - AvgMytilus[Before_After == 'Before']),
                   PhyllodeltaAvg = -1*(AvgPhyllo[Before_After == 'After'] - AvgPhyllo[Before_After == 'Before']),
                   #multiply delta of mytilus and phyllospadix by -1 to change to Percent loss
                   micromacroalgaedeltaAvg = (AvgFleshyalgae[Before_After == 'After'] - AvgFleshyalgae[Before_After == 'Before']), 
                   #only micro/acroalgae cover
                   CCAdeltaAvg = (AvgCCA[Before_After == 'After'] - AvgCCA[Before_After == 'Before']), 
                   #crustose CCA and articulated CA
                   prodomdeltaAvg = (Avgproddom[Before_After == 'After'] -  Avgproddom[Before_After == 'Before']), 
                   #macroalgae cover - sessile consumers (including Mytilus)
                   SAVav = mean(AvgSAV), #ave between before and after since SA/V changed with fspp removal
                   THav = mean(AvgTH)) #tide height didn't change 

SEMdaynightAvg<-as.data.frame(SEMdaynightAvg)

#####Surfgrass SEM model#####
Phyllounscaled<-SEMdaynightAvg %>%
  filter(Foundation_spp == 'Phyllospadix') %>%
  dplyr::rename(NEP =NEPdeltaAvg, NEC=NECdeltaAvg, pH=pHdeltaAvg, NitriteNitrate =NNdeltaAvg, 
                Ammonium =NH4deltaAvg,Phosphate =POdeltaAvg,NtoPRatio=NtoPdeltaAvg,
                MaxLight =Parmaxdeltaavg, MaxTemp = Tempmaxdeltaavg, Light=ParmeandeltaAvg, 
                MytilusLoss=MytilusdeltaAvg,PhyllospadixLoss=PhyllodeltaAvg, MicroMacroAlgaeCover=micromacroalgaedeltaAvg,
                SAtoVRatio=SAVav,TideHeight=THav, RawTemp=TempmeandeltaAvg) #rename cols to match sem
#check collinearity
#ggpairs(Phyllounscaled[c(4:6,10,12,15,18:19,22:23)]) #good


PAvgMMAlgae<-lm(MicroMacroAlgaeCover~ PhyllospadixLoss+SAtoVRatio +TideHeight,  data = Phyllounscaled)
PAvgLight<-lm(Light ~PhyllospadixLoss+SAtoVRatio +TideHeight,data=Phyllounscaled)
PAvgTemp<-lm(MaxTemp~ PhyllospadixLoss+SAtoVRatio+TideHeight  ,data = Phyllounscaled)
PAvgNtoP<-lm(NtoPRatio ~PhyllospadixLoss+SAtoVRatio +TideHeight, data =Phyllounscaled)
PAvgNEC<- lm(NEC~pH + MaxTemp +PhyllospadixLoss +SAtoVRatio+TideHeight,data =Phyllounscaled)
PAvgpH<- lm(pH ~ PhyllospadixLoss+NEP + SAtoVRatio+TideHeight, data = Phyllounscaled)
PAvgNEP<-lm(NEP ~Light+ MicroMacroAlgaeCover
            + NtoPRatio+SAtoVRatio +TideHeight, data = Phyllounscaled) 

#normality and homogeneity of variance checks
plot(PAvgMMAlgae)
qqp(resid(PAvgMMAlgae),"norm")
#good

plot(PAvgTemp)
qqp(resid(PAvgTemp),"norm")
#good
plot(PAvgLight)
qqp(resid(PAvgLight),"norm")
#okay, one point out
plot(PAvgNtoP)
qqp(resid(PAvgNtoP),"norm")
#okay,two points barely out
plot(PAvgNEP)
qqp(resid(PAvgNEP),"norm")
#good
plot(PAvgpH)
qqp(resid(PAvgpH),"norm")
#okay one point out

plot(PAvgNEC)
qqp(resid(PAvgNEC),"norm")
#good
#piecewise SEM
PhylloSEM<-psem(PAvgMMAlgae,
                  PAvgLight,
                  PAvgTemp,
                  PAvgNtoP,
                  PAvgNEP, 
                  PAvgpH,
                  PAvgNEC,
                  MaxTemp%~~%Light)

PhylloSEMsum<-summary(PhylloSEM,standardize = 'scale',center = "TRUE") #scale data in sum
plot(PhylloSEM) #test plot to see what is significant

#output csvs
PhylloModels<-as.data.frame(PhylloSEMsum$call)
Pvalue<-as.data.frame(PhylloSEMsum$dTable)
Fishers<-as.data.frame(PhylloSEMsum$Cstat)
ICs<-as.data.frame(PhylloSEMsum$IC)
coeff<-as.data.frame(PhylloSEMsum$coefficients)
R2<-as.data.frame(PhylloSEMsum$R2)

write.csv(PhylloModels, 'Output/phyllosemmodelNEP.csv' )
write.csv(Pvalue, 'Output/phyllosempvalNEP.csv' )
write.csv(Fishers, 'Output/phyllosemfishersNEP.csv' )
write.csv(ICs, 'Output/phylloicsNEP.csv' )
write.csv(coeff, 'Output/phyllosemmodelcoeffNEP.csv' )
write.csv(R2, 'Output/phyllosemmodelr2NEP.csv' )

#####CA Mussel SEM#####
Mytilusunscaled<-SEMdaynightAvg %>%
  filter(Foundation_spp == 'Mytilus') %>%
  dplyr::rename(NEP =NEPdeltaAvg, NEC=NECdeltaAvg, pH=pHdeltaAvg, NitriteNitrate =NNdeltaAvg, 
                Ammonium =NH4deltaAvg,Phosphate =POdeltaAvg,NtoPRatio=NtoPdeltaAvg,
                MaxLight =Parmaxdeltaavg, MaxTemp = Tempmaxdeltaavg, Light=ParmeandeltaAvg, 
                MytilusLoss=MytilusdeltaAvg,PhyllospadixLoss=PhyllodeltaAvg, MicroMacroAlgaeCover=micromacroalgaedeltaAvg,
                SAtoVRatio=SAVav,TideHeight=THav, RawTemp=TempmeandeltaAvg) #rename cols to match sem
#check collinearity
#ggpairs(Mytilusunscaled[c(4:6,10,12,15,17,19,22:23)]) #good

#models based on hypotheses
MAvgMMalgae<-lm(MicroMacroAlgaeCover ~ MytilusLoss + SAtoVRatio+TideHeight, data = Mytilusunscaled)
MAvgTemp<-lm(MaxTemp~MytilusLoss +SAtoVRatio+TideHeight , data = Mytilusunscaled)
MAvgLight<-lm(Light ~ MytilusLoss +SAtoVRatio+TideHeight, data = Mytilusunscaled)
MAvgNtoP<-lm(NtoPRatio ~  MytilusLoss+SAtoVRatio+TideHeight,  data = Mytilusunscaled)
MAvgNEC<- lm(NEC~ MytilusLoss + pH +SAtoVRatio+MaxTemp+TideHeight, data = Mytilusunscaled)
MAvgpH<- lm(pH ~ MytilusLoss+NEP+ SAtoVRatio+TideHeight, data = Mytilusunscaled)
MAvgNEP<-lm(NEP ~Light+MicroMacroAlgaeCover +
              NtoPRatio +SAtoVRatio +TideHeight,data = Mytilusunscaled)

MytilusSEM<-psem(MAvgMMalgae,
                   MAvgLight,
                   MAvgTemp,
                   MAvgNtoP,
                   MAvgNEP,
                   MAvgpH,
                   MAvgNEC,
                   MaxTemp%~~%Light)
MytilusSEMSum<-summary(MytilusSEM,standardize = "scale",center = "TRUE")
plot(MytilusSEM) 

#output into dataframes
MytilusModelsNEP<-as.data.frame(MytilusSEMSum$call)
MytilusPvalueNEP<-as.data.frame(MytilusSEMSum$dTable)
MytilusFishersNEP<-as.data.frame(MytilusSEMSum$Cstat)
MytilusICsNEP<-as.data.frame(MytilusSEMSum$IC)
MytiluscoeffNEP<-as.data.frame(MytilusSEMSum$coefficients)
MytilusR2NEP<-as.data.frame(MytilusSEMSum$R2)
#save csvs
write.csv(MytilusModelsNEP, 'Output/mytilusSEMmodelNEP.csv' )
write.csv(MytilusPvalueNEP, 'Output/mytilusSEMpvalNEP.csv' )
write.csv(MytilusFishersNEP, 'Output/mytilussemSEMfishersNEP.csv' )
write.csv(MytilusICsNEP, 'Output/mytilusSEMicsNEP.csv' )
write.csv(MytiluscoeffNEP, 'Output/mytilusSEMcoeffNEP.csv' )
write.csv(MytilusR2NEP, 'Output/mytilusSEMr2NEP.csv' )

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

