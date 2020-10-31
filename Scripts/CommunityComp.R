##Community Composition data for Sessile and Mobiles from Surfgrass and Mussel
##Oregon Tide Pools
##By: Jenn Fields
##Last updated: 10.23.2020


rm(list=ls()) #Clears the environment

#load libraries
library(vegan)
library(broom)
library(cowplot)
library(ggpubr)
library(reshape2)
library(hrbrthemes)
library(plotrix)
library(patchwork)
library(car)
library(MASS)
library(lattice)
library(plotrix)
devtools::install_github("seananderson/ggsidekick")
library(ggsidekick) #for theme sleek
library(effects) #for ggeffect
library(ggeffects) #marginal effects on complex models
library(margins)
library(emmeans)
library(parameters)
library(tidyverse)


#load data from community comp surveys
#source scripts
source("scripts/tidepoolphysicalparameters.R")
#data
Sessiles <- read_csv("Data/CommunityComposition/SessilesAll.csv")
Mobiles <- read_csv("Data/CommunityComposition/Mobiles.csv")
SessilesGroupings<-read_csv("Data/CommunityComposition/SessilesFun.csv")
MobileGroupings<-read_csv("Data/CommunityComposition/MobilesFun.csv")


#replace NA with 0
Mobiles[is.na(Mobiles)]<-0


#starts_with("Diatom")

#######Sessiles#########
#convert characters to numeric in sessile sheet
Sessiles$Epiactis.prolifera<-as.numeric(Sessiles$Epiactis.prolifera)
Sessiles$Chaetomorpha.linum<-as.numeric(Sessiles$Chaetomorpha.linum)
Sessiles$Costaria.costata<-as.numeric(Sessiles$Costaria.costata)
Sessiles[is.na(Sessiles)]<-0 


# Make all the community data a relative percent
PercentSessile<-100*Sessiles[7:ncol(Sessiles)]/Sessiles$Squares #change to rock--end spp
#View(PercentSessile)
#normalize to the sum of the total cover (since it can be greater than 100%)
PercentSessile<- 100*PercentSessile/rowSums(PercentSessile)



Communitymetrics<-cbind(Sessiles$PoolID, Sessiles$Foundation_spp, Sessiles$Removal_Control, 
                        Sessiles$Before_After,PercentSessile)
colnames(Communitymetrics)[1:4]<- c("PoolID","Foundation_spp","Removal_Control","Before_After")

Communitymetrics<-Communitymetrics %>%
  dplyr::filter(Before_After != "Immediate") %>%
  dplyr::filter(PoolID != 30) %>% #remove pool 30 since no biogeochem values for this pool 
  dplyr:: mutate(allCCA = (Crustose.coralline + Bosiella.spp + Corallina.spp + Corallina.vancouveriensis +
                             Calliarthron.tuberculosum), #creates column for all CCA
                 macroalgae = (Diatoms + Algae.film + Turf.algae +Acrosiphonia.coalita + Analipus.japonicus + Chaetomorpha.linum +
                                 Centroceras.Ceramium + Cladophora.columbiana + Cryptosiphonia.wooddii + Costaria.costata +
                                 Cumagloia.andersonii + Erythrophyllum.delessrioides +	Endocladia.muricata	+
                                 Fucus.gardneri +	Halosaccion.glandiforme	+ Pyropia.spp	+ Polysiphonia.spp +
                                 Ptilota.spp +	Cryptopluera.spp + Odonthalia.floccosa +	Microcladia.borealis +	Mazzaella.splendens +
                                 Mazzaella.flaccida	+ Mazzaella.oregona +	Neorhodomela.larix + Odanthalia.washingtoniensis +Scytosiphon.lomentaria	
                               + Osmundea.spectabilis	+ Ulva.spp + Plocamium.pacificum + Smithura.naiadum + Mastocarpus	+ Farlowia.mollis	+
                                 Savoiea.robusta +	Palmaria.hecatensis	+ Melobesia.mediocris +	Leathesia.marina + Callithamnion.pikeanum	+
                                 Laminara.setchellii + Schizymenia.pacifica),
                 macrophytes = (Phyllospadix.spp + macroalgae), #includes phyllospadix & macroalgae
                 macroCCA = (macroalgae + allCCA), #includes macroalgae and CCA
                 consumers = (Chthamalus +	Semibalanus.cariosus +	Balanus.nibulis	+ Balanus.glandula +
                                Pollicipes.polymerus +	tube.worm	+ Ophlitaspongia.pennata + Halichondria +	
                                Haliclona.permollis	+ Anthropluera.elegantissima	+ Anthropluera.xanthogrammica	+ 
                                Urticina.coriacea	+ Epiactis.prolifera +Anthopleura.artemisia	+ Stylantheca.spp),
                 #includes all consumers (no mytilus)
                 allconsumers = (consumers + Mytilus.californianus), #consumers and mytilus
                 prodphyllodom = (macroalgae - (consumers + Mytilus.californianus)), #producer dominance for phyllo model (so subtract mytilus too)
                 allproddom = (macrophytes - allconsumers)) %>%#prod dominance with foundation spp
  dplyr::select(PoolID,Foundation_spp,Removal_Control,Before_After, Mytilus.californianus, Phyllospadix.spp,allCCA, macroalgae,macrophytes, macroCCA,
                consumers,allconsumers,  prodphyllodom, allproddom)

#create data sheet with only these columns for SEM
write_csv(Communitymetrics, path="Data/CommunityComposition/Communitymetrics.csv")


SessilesAll<-data.frame(Sessiles$PoolID, Sessiles$Foundation_spp, Sessiles$Removal_Control, Sessiles$Before_After,PercentSessile)
colnames(SessilesAll)[1:4]<- c("PoolID","Foundation_spp","Removal_Control","Before_After")

#Change in foundation species cover
Funsppcover<- SessilesAll%>%
  dplyr::group_by(PoolID,Foundation_spp, Removal_Control) %>%
  summarise(Mytilusdelta = -1*(Mytilus.californianus[Before_After == 'After'] - Mytilus.californianus[Before_After == 'Before']),
            Phyllodelta = -1*(Phyllospadix.spp[Before_After == 'After'] - Phyllospadix.spp[Before_After == 'Before']))

#####Sessile functional groups####
SessilesMytilusFun<-SessilesAll %>%
  dplyr::filter(Before_After != 'Immediate') %>%
  dplyr::group_by(PoolID, Foundation_spp, Before_After, Removal_Control) %>%
  tidyr::pivot_longer(
    cols = Phyllospadix.spp:Stylantheca.spp,
    names_to = "Species", #creates column with species in longformat
    values_to = "Cover", #adds column for % cover
    values_drop_na = TRUE
  ) 
SessilesphylloFun<-SessilesAll %>%
  dplyr::filter(Before_After != 'Immediate') %>%
  dplyr::group_by(PoolID, Foundation_spp, Before_After, Removal_Control) %>%
  tidyr::pivot_longer(
    cols = c(8,10:70), #excluding phyllo since taking delta cover
    names_to = "Species", #creates column with species in longformat
    values_to = "Cover", #adds column for % cover
    values_drop_na = TRUE
  ) 
SessilesMytilusStacked<-left_join(SessilesMytilusFun,SessilesGroupings) #combine with fun groups 
SessilesPhylloStacked<-left_join(SessilesphylloFun,SessilesGroupings) #combine with fun groups 

#before then take avg of tide height and sa to v ratio
PhysicalParameters<-TidePooldes[, c("PoolID","Removal_Control","MaxDepth","Perimeter","SurfaceArea", "Vol", "SAtoV", "TideHeight")] #pull out the necessary columns and treatment 

PhysicalParameters$PoolID<-as.factor(PhysicalParameters$PoolID)
Funsppcover$PoolID<-as.factor(Funsppcover$PoolID)
Funsppandpp<-left_join(Funsppcover,PhysicalParameters)

PP<-PhysicalParameters %>%
  dplyr::group_by(PoolID,Removal_Control) %>%
  dplyr::summarise(SAVav = mean(SAtoV), #ave between before and after since SA/V changed with fspp removal
                   THav = mean(TideHeight),SAav=mean(SurfaceArea)) #tide height didn't change 

PP$PoolID<-as.factor(PP$PoolID)
Funsppcover$PoolID<-as.factor(Funsppcover$PoolID)
Funsppandpp<-left_join(Funsppcover,PP)


deltaMV<-SessilesMytilusStacked %>%
  dplyr::group_by(PoolID,Removal_Control,Foundation_spp,Before_After,Functional_Group) %>%
  summarise(SumCover = sum(Cover)) %>%
  dplyr::group_by(PoolID,Removal_Control,Foundation_spp, Functional_Group) %>%
  summarise(Deltacover = SumCover[Before_After=="After"]-SumCover[Before_After =="Before"])

phyllodeltacover<-SessilesPhylloStacked%>%
  dplyr::group_by(PoolID,Removal_Control,Foundation_spp,Functional_Group, Before_After) %>%
  summarise(SumCover = sum(Cover)) %>%
  dplyr::group_by(PoolID,Removal_Control,Foundation_spp, Functional_Group) %>%
  summarise(Deltacover = SumCover[Before_After=="After"]-SumCover[Before_After =="Before"])

phyllodeltacover$PoolID<-as.factor(phyllodeltacover$PoolID)
deltaMV$PoolID<-as.factor(deltaMV$PoolID)
PdeltaMV<-left_join(phyllodeltacover,Funsppandpp)
MdeltaMV<-left_join(deltaMV,Funsppandpp)

MdeltaMV<-MdeltaMV %>%
  dplyr::filter(Foundation_spp =="Mytilus") %>%
  dplyr::select(PoolID,Removal_Control,Foundation_spp, Functional_Group,Deltacover,Mytilusdelta,SAVav,THav)
MdeltaMV<-as.data.frame(MdeltaMV)

PdeltaMV<-PdeltaMV %>%
  dplyr::filter(Foundation_spp =="Phyllospadix") %>%
  dplyr::select(PoolID,Removal_Control,Foundation_spp, Functional_Group,Deltacover,Phyllodelta,SAVav,THav)
PdeltaMV<-as.data.frame(PdeltaMV)

#colorblind friendly graph colors for all ggplots
Colors<-c(
  Anemone="#54278f",
  ArticulatedCorallines="#c51b8a",
  Crustose="#c51b7d",
  SuspensionFeeder="#636363",
  Microalgae = "#78c679",
  Filamentous = "#31a354", 
  FoiloseAlgae="#006837",
  CorticatedFoliose="#addd8e",
  CorticatedMacro="#31a354",
  LeatheryMacro="#993404")

MdeltaMV$QCover<-sign(MdeltaMV$Deltacover)*sqrt(sqrt(abs(MdeltaMV$Deltacover)))
#take sign so adding positive and negative back and take sqrt of positive data since can't have sqrt of negative

#model with Quad rooted cover data
#Run model to get scaled estimates for plot to determine magnitude of effects with different ranges
MdeltaMV$Functional_Group<-as.factor(MdeltaMV$Functional_Group)
Mytilussessmod<-lm(QCover~(Mytilusdelta +SAVav+THav)*Functional_Group, data =MdeltaMV)

#plot(Mytilussessmod)
#plot(Mytilussessmod) #good
qqp(resid(Mytilussessmod), "norm") #good

summary(Mytilussessmod)
MdeltaMV$Functional_Group<-as.factor(MdeltaMV$Functional_Group)

#Marginal Effects with emtrends
Mytilussessfuntrends<-emtrends(Mytilussessmod,pairwise ~Functional_Group, var="sqrt(sqrt(Mytilusdelta))") 
#since data is quad rooted want the response related to quadrooted data
MSFun<-summary(Mytilussessfuntrends, infer=c(TRUE, TRUE))
#gives you t-ratio and pvalue for emtrends as well as contrasts
# CI and p-values are adjusted with "tukey" using the studentized range distribution with the number of means (n=10)
MSFuntrends<-as.data.frame(MSFun$emtrends) #create data frame for graphs
#write.csv(MSFuntrends, file="Output/MSFuntrends.csv")

MSFuntrends$sig<-c("nonsig","nonsig","nonsig","sig","nonsig","nonsig","nonsig","nonsig","sig","sig") #add column for sig

MSFuntrends<-MSFuntrends %>%
  rename("Mytilusdeltatrend"="sqrt(sqrt(Mytilusdelta)).trend") #rename sqrt col output to something easier



Msessfunplot<-ggplot(MSFuntrends,aes(x=Functional_Group, y=Mytilusdeltatrend,alpha=sig, col = Functional_Group)) +
  geom_point(size = 7) +
  coord_flip()+
  geom_pointrange(aes(y=Mytilusdeltatrend, ymin=lower.CL, ymax=upper.CL), size=2) + #graph estimates and st err for all values
  geom_hline(yintercept = 0, lty = 2,size =1.5)  + #add vertical line to graph and make it dashed 
  scale_colour_manual(values = Colors) +
  guides(alpha= FALSE, col = FALSE) + #remove alpha legend
  scale_alpha_discrete(range = c(0.1, 5.0)) + #set range for sig vs non sig alphas
  theme_classic() +
  theme(axis.text.x=element_text(size = 35, color = "black"), 
        axis.text.y=element_text(size = 35, color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 40, color = "black"),
        plot.title = element_text(size =40, color ="black", hjust = 0.5)) +
  theme(legend.position="none")+
  scale_x_discrete(labels=c("Crustose"="Crustose", "ArticulatedCorallines"= "Articulated corallines", "Microalgae"="Microalgae", "Filamentous"="Filamentous", "FoiloseAlgae"="Foliose algae", 
                            "CorticatedFoliose"="Corticated foliose","CorticatedMacro"="Corticated macroalgae", "LeatheryMacro"="Leathery macrophytes","Anemone"="Anemone","SuspensionFeeder"="Suspension feeders")) +  #rename y axis tickmarks 
  labs(y="Effect size (cover)",x = "",title =expression(~italic(M.)~''~italic(californianus)~ 'percent loss')) 
Msessfunplot


SAVMsessfuntrends<-emtrends(Mytilussessmod,pairwise ~Functional_Group, var="sqrt(sqrt(SAVav))" )
SAVMsessfun<-summary(SAVMsessfuntrends, infer=c(TRUE, TRUE))
#gives you t-ratio and pvalue for emtrends as well as contrasts
# CI and p-values are adjusted with "tukey" using the studentized range distribution with the number of means (n=10)
SAVMSFuntrends<-as.data.frame(SAVMsessfun$emtrends) #create data frame for graphs
#write.csv(SAVMSFuntrends, file="Output/SAVMSFuntrends.csv")
SAVMSFuntrends$sig<-c("nonsig","nonsig","nonsig","nonsig","nonsig","nonsig","nonsig","nonsig","sig","nonsig") #add column 
SAVMSFuntrends<-SAVMSFuntrends%>%
  rename("SAVav.trend"="sqrt(sqrt(SAVav)).trend") #rename sqrt col output to something easier


SAVMsessplot<-ggplot(SAVMSFuntrends,aes(x=Functional_Group, y=SAVav.trend, col = Functional_Group,alpha=sig)) +
  geom_point(size = 7) +
  coord_flip()+
  geom_pointrange(aes(y=SAVav.trend, ymin=lower.CL, ymax=upper.CL), size=2) + #graph estimates and st err for all values
  geom_hline(yintercept = 0, lty = 2,size =1.5)  + #add vertical line to graph and make it dashed 
  scale_colour_manual(values = Colors) +
  guides(alpha= FALSE, col = FALSE) + #remove alpha legend
  scale_alpha_discrete(range = c(0.1, 5.0)) + #set range for sig vs non sig alphas
  theme_classic() +
  theme(axis.text.x=element_text(size = 30, color = "black"), 
        axis.text.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 40, color = "black"),
        plot.title = element_text(size =40, color ="black", hjust = 0.5)) +
  theme(legend.position="none")+
  scale_x_discrete(labels=c("Crustose"="Crustose", "ArticulatedCorallines"= "Articulated Corallines", 
                            "Microalgae"="Microalgae", "Filamentous"="Filamentous", "FoiloseAlgae"="Foliose Algae", 
                            "CorticatedFoliose"="Corticated Foliose","CorticatedMacro"="Corticated Macroalgae", "LeatheryMacro"="Leathery Macrophytes","Anemone"="Anemone","SuspensionFeeder"="Suspension Feeders")) +  #rename y axis tickmarks 
  labs(x="",y = "Effect size (percent cover)",title="Size of tide pool (SA:V)") 
SAVMsessplot
THMsessfun<-emtrends(Mytilussessmod,pairwise ~Functional_Group, var="sqrt(sqrt(THav))" )
THMsessfuntrends<-summary(THMsessfun,infer=c(TRUE, TRUE))
THMsfuntrends<-as.data.frame(THMsessfuntrends$emtrends) 
#write.csv(THMsfuntrends, file="Output/THMsfuntrends.csv")

THMsfuntrends$sig<-c("nonsig","nonsig","nonsig","nonsig","nonsig","nonsig","nonsig","nonsig","sig","nonsig") #add column 
THMsfuntrends<-THMsfuntrends%>%
  rename("THav.trend"="sqrt(sqrt(THav)).trend") #rename sqrt col output to something easier
THMsessplot<-ggplot(THMsfuntrends,aes(x=Functional_Group, y=THav.trend, col = Functional_Group, alpha=sig)) +
  geom_point(size = 7) +
  coord_flip()+
  geom_pointrange(aes(y=THav.trend, ymin=lower.CL, ymax=upper.CL), size=2) + #graph estimates and st err for all values
  geom_hline(yintercept = 0, lty = 2,size =1.5)  + #add vertical line to graph and make it dashed 
  scale_colour_manual(values = Colors) +
  guides(alpha= FALSE, col = FALSE) + #remove alpha legend
  scale_alpha_discrete(range = c(0.1, 5.0)) + #set range for sig vs non sig alphas
  theme_classic() +
  theme(axis.text.x=element_text(size = 26, color = "black"), 
        axis.text.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(face="bold",size = 30, color = "black"),
        plot.title = element_text(size= 40, color ="black", hjust = 0.5)) +
  theme(legend.position="none")+
  scale_x_discrete(labels=c("Crustose"="Crustose", "ArticulatedCorallines"= "Articulated Corallines", 
                            "Microalgae"="Microalgae", "Filamentous"="Filamentous", "FoiloseAlgae"="Foliose Algae", 
                            "CorticatedFoliose"="Corticated Foliose","CorticatedMacro"="Corticated Macroalgae", "LeatheryMacro"="Leathery Macrophytes","Anemone"="Anemone","SuspensionFeeder"="Suspension Feeders")) +  #rename y axis tickmarks 
  labs(y="",x = "",title="Tide height") 
THMsessplot
####Phyllo Sessile Analysis and ggplots####
PdeltaMV$QCover<-sign(PdeltaMV$Deltacover)*sqrt(sqrt(abs(PdeltaMV$Deltacover))) 
#don't need to do this below
PdeltaMV<-PdeltaMV %>%
  mutate_at(.vars= c('Phyllodelta','SAVav', 'THav'), 
            .funs = list(std = ~scale(.)))
PdeltaMV$Functional_Group<-as.factor(PdeltaMV$Functional_Group)
Phyllosessmod<-lm(QCover~(Phyllodelta +SAVav+THav) * Functional_Group, data =PdeltaMV)
#plot(Phyllosessmod) #good
qqp(resid(Phyllosessmod),"norm")
#okay few points out
summary(Phyllosessmod)

#all factors interaction with functional group
Phyllosessfun<-emtrends(Phyllosessmod,pairwise ~Functional_Group, var="sqrt(sqrt(Phyllodelta))")

PSFun<-summary(Phyllosessfun, infer=c(TRUE, TRUE))
#gives you t-ratio and pvalue for emtrends as well as contrasts
# CI and p-values are adjusted with "tukey" using the studentized range distribution with the number of means (n=10)
PSFuntrends<-as.data.frame(PSFun$emtrends) #create data frame for graphs
#write.csv(PSFuntrends, file="Output/PSFuntrends.csv")

PSFuntrends$sig<-c("sig","nonsig","sig","nonsig","sig","nonsig","nonsig","nonsig","nonsig","nonsig") #add column for sig
PSFuntrends<-PSFuntrends%>%
  rename(Phyllodelta.trend="sqrt(sqrt(Phyllodelta)).trend")
#ggplots of Effect Sizes
Psessfunplot<-ggplot(PSFuntrends,aes(x=Functional_Group, y=Phyllodelta.trend, col = Functional_Group, alpha=sig)) +
  geom_point(size = 7) +
  coord_flip()+
  geom_pointrange(aes(y=Phyllodelta.trend, ymin=lower.CL, ymax=upper.CL), size=2) + #graph estimates and st err for all values
  geom_hline(yintercept = 0, lty = 2, size =1.5)  + #add vertical line to graph and make it dashed 
  scale_colour_manual(values = Colors) +
  guides(alpha= FALSE, col = FALSE) + #remove alpha legend
  scale_alpha_discrete(range = c(0.1, 5.0)) + #set range for sig vs non sig alphas
  theme_classic() +
  theme(axis.text.x=element_text(size = 35, color = "black"), 
        axis.text.y=element_text(size = 35, color = "black"),
        axis.title.y = element_text(size = 40, color = "black",face="bold"),
        axis.title.x = element_text(size = 40, color = "black"),
        plot.title = element_text(size =40, color ="black", hjust = 0.5)) +
  theme(legend.position="none")+
  scale_x_discrete(labels=c("Crustose"="Crustose", "ArticulatedCorallines"= "Articulated corallines", 
                            "Microalgae"="Microalgae", "Filamentous"="Filamentous", 
                            "FoiloseAlgae"="Foliose algae", 
                            "CorticatedFoliose"="Corticated foliose","CorticatedMacro"="Corticated macroalgae", "LeatheryMacro"="Leathery macrophytes","Anemone"="Anemone","SuspensionFeeder"="Suspension feeders")) +  #rename y axis tickmarks 
  labs(y="Effect size (cover)",x = "",title =expression(~italic(Phyllospadix)~"spp. percent loss" )) 


SAVPsessfun<-emtrends(Phyllosessmod,pairwise ~Functional_Group, var="sqrt(sqrt(SAVav))" )
SAVPsessfuntrends<-summary(SAVPsessfun,infer=c(TRUE, TRUE))
SAVPsfuntrends<-as.data.frame(SAVPsessfuntrends$emtrends) 
#write.csv(SAVPsfuntrends, file="Output/SAVPsfuntrendss.csv")

SAVPsfuntrends$sig<-c("nonsig","nonsig","nonsig","nonsig","nonsig","nonsig","nonsig","nonsig","nonsig","nonsig")
SAVPsfuntrends<-SAVPsfuntrends %>%
  rename(SAVav.trend="sqrt(sqrt(SAVav)).trend")
SAVPsessplot<-ggplot(SAVPsfuntrends,aes(x=Functional_Group, y=SAVav.trend, col = Functional_Group,alpha =sig)) +geom_point(size = 7, alpha=0.3) +
  geom_point(size=7)+
  coord_flip()+
  geom_pointrange(aes(y=SAVav.trend, ymin=lower.CL, ymax=upper.CL), size=2) + #graph estimates and st err for all values
  geom_hline(yintercept = 0, lty = 2, size =1.5)  + #add vertical line to graph and make it dashed 
  scale_colour_manual(values = Colors) +
  guides(alpha= FALSE, col = FALSE) + #remove alpha legend
  scale_alpha_discrete(range = c(0.1, 5.0)) + #set range for sig vs non sig alphas
  theme_classic() +
  theme(axis.text.x=element_text(size = 30, color = "black"), 
        axis.text.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 40, color = "black"),
        plot.title = element_text(size =40, color ="black", hjust = 0.5)) +
  labs(y = "Effect size (percent cover)",x="",title ="Size of tide pool (SA:V)") +
  theme(legend.position="none")+
  scale_x_discrete(labels=c("Crustose"="Crustose", "ArticulatedCorallines"= "Articulated Corallines", 
                            "Microalgae"="Microalgae", "Filamentous"="Filamentous", 
                            "FoiloseAlgae"="Foliose Algae", 
                            "CorticatedFoliose"="Corticated Foliose","CorticatedMacro"="Corticated Macroalgae", 
                            "LeatheryMacro"="Leathery Macrophytes","Anemone"="Anemone",
                            "SuspensionFeeder"="Suspension Feeders"))   #rename y axis tickmarks 


THPsessfun<-emtrends(Phyllosessmod,pairwise ~Functional_Group, var="sqrt(sqrt(THav))" )
THPsessfuntrends<-summary(THPsessfun,infer=c(TRUE, TRUE))
THPsfuntrends<-as.data.frame(THPsessfuntrends$emtrends) 
#write.csv(THPsfuntrends, file="Output/THPsfuntrends.csv")

THPsfuntrends$sig<-c("nonsig","nonsig","nonsig","nonsig","nonsig","nonsig","nonsig","nonsig","nonsig","nonsig")
THPsfuntrends<-THPsfuntrends%>%
  rename(THav.trend="sqrt(sqrt(THav)).trend")
THPsessplot<-ggplot(THPsfuntrends,aes(x=Functional_Group, y=THav.trend, col = Functional_Group, alpha =sig)) +
  geom_point(size = 7) +
  coord_flip()+
  geom_pointrange(aes(y=THav.trend, ymin=lower.CL, ymax=upper.CL), size=2) + #graph estimates and st err for all values
  geom_hline(yintercept = 0, lty = 2, size =1.5)  + #add vertical line to graph and make it dashed 
  scale_colour_manual(values = Colors) +
  guides(alpha= FALSE, col = FALSE) + #remove alpha legend
  scale_alpha_discrete(range = c(0.1, 5.0)) + #set range for sig vs non sig alphas
  theme_classic() +
  theme(axis.text.x=element_text(size = 30, color = "black"), 
        axis.text.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(face="bold",size = 30, color = "black"),
        plot.title = element_text(size =40, color ="black", hjust = 0.5)) +
  theme(legend.position="none")+
  scale_x_discrete(labels=c("Crustose"="Crustose", "ArticulatedCorallines"= "Articulated Corallines", 
                            "Microalgae"="Microalgae", "Filamentous"="Filamentous",
                            "FoiloseAlgae"="Foliose Algae", 
                            "CorticatedFoliose"="Corticated Foliose","CorticatedMacro"="Corticated Macroalgae", "LeatheryMacro"="Leathery Macrophytes","Anemone"="Anemone","SuspensionFeeder"="Suspension Feeders")) +  #rename y axis tickmarks 
  labs(y="",x = "", title ="Tide height")

#####Mobiles######

MobilesFun<-Mobiles %>%
  dplyr::filter(Before_After != 'Immediate') %>%
  dplyr::group_by(PoolID, Foundation_spp, Before_After, Removal_Control) %>%
  tidyr::pivot_longer(
    cols = Nuttalina.spp:Gunnel,
    names_to = "Species", #creates column with species in longformat
    values_to = "Count", #adds column for % cover
    values_drop_na = TRUE
  ) 
Mobstacked<-left_join(MobilesFun,MobileGroupings)
Mobstacked$PoolID<-as.factor(Mobstacked$PoolID)
Mobdata<-left_join(Mobstacked,Funsppandpp)
Mobdata$Std.Count<-Mobdata$Count/Mobdata$SAav #std counts by SA Count/m^2

Mobdata<- Mobdata%>%
  dplyr::group_by(PoolID,Removal_Control,Foundation_spp,Before_After,Functional_Group,Species) %>%
  summarise(SumDensity = sum(Std.Count)) %>%
  dplyr::group_by(PoolID,Removal_Control,Foundation_spp, Functional_Group,Species) %>%
  summarise(DeltaCount = SumDensity[Before_After=="After"]-SumDensity[Before_After =="Before"]) 

Mobdata<-left_join(Mobdata,Funsppandpp)

Mobdata<-as.data.frame(Mobdata)


Mobdata$logCount<-sign(Mobdata$DeltaCount)*log(abs(Mobdata$DeltaCount))+1

MytilusMobmod<-Mobdata%>%
  dplyr::filter(Foundation_spp == "Mytilus" & Functional_Group != "SuspensionFeeder") #no suspension feeders in mytilus pools

PhylloMobmod<-Mobdata%>%
  filter(Foundation_spp == "Phyllospadix"& Functional_Group != "SuspensionFeeder") 
#removed suspension feeder since only one in one tide pool in one time period 
#normality assumptions were met when suspension feeder (sea cucumber) was removed
Mytilusmobmod<-lm(logCount~(Mytilusdelta +SAVav+THav)*Functional_Group, data =MytilusMobmod)

#plot(Mytilusmobmod) #good
qqp(resid(Mytilusmobmod), "norm") #good
summary(Mytilusmobmod)

Phyllomobmod<-lm(logCount~(Phyllodelta +SAVav+THav)*Functional_Group, data =PhylloMobmod)
#plot(Phyllomobmod) #good

qqp(resid(Phyllomobmod), "norm") #okay some points out
summary(Phyllomobmod)

#marginal effects with emtrends

#colorblind friendly colors for all mobile plots
MobColors<-c(
  Carnivore = "#08519c",
  Herbivore = "#6baed6",
  Omnivores ="#3182bd")

#####Mussel mobile emtrend plots#####
Mytilusmobfun<-emtrends(Mytilusmobmod,pairwise ~Functional_Group, var="log(Mytilusdelta)+1")
Mytilusmobfuntrends<-summary(Mytilusmobfun,infer=c(TRUE, TRUE))
Mmobfuntrends<-as.data.frame(Mytilusmobfuntrends$emtrends) 
#write.csv(Mmobfuntrends,file="Output/Mmobfuntrends.csv")
Mmobfuntrends$sig<-c("nonsig","sig","nonsig")
Mmobfuntrends<-Mmobfuntrends%>%
  rename(Mytilusdelta.trend="log(Mytilusdelta)+1.trend")
Mmobfunplot<-ggplot(Mmobfuntrends,aes(x=Functional_Group, y=Mytilusdelta.trend, col = Functional_Group, alpha = sig)) +
  geom_point(size = 7) +
  coord_flip()+
  geom_pointrange(aes(y=Mytilusdelta.trend, ymin=lower.CL, ymax=upper.CL),size=2) + #graph estimates and st err for all values
  geom_hline(yintercept = 0, lty = 2, size = 1.5)  + #add vertical line to graph and make it dashed 
  guides(alpha= FALSE, col = FALSE) + #remove alpha legend
  scale_alpha_discrete(range = c(0.1, 5.0)) + #set range for sig vs non sig alphas
  scale_colour_manual(values = MobColors) +
  theme_classic() +
  theme(axis.text.x=element_text(size = 35, color = "black"), 
        axis.text.y=element_text(size = 35, color = "black"),
        axis.title.y = element_text(size = 40, color = "black",face="bold"),
        axis.title.x = element_text(size = 40, color = "black"))+
  scale_x_discrete(labels=c("Carnivore"="Carnivores", "Herbivore"= "Herbivores", 
                            "Omnivores"="Omnivores")) +
  theme(legend.position="none")+
  labs(y="Effect size (density)", x = "")
Mmobfunplot

SAVMmobfun<-emtrends(Mytilusmobmod,pairwise ~Functional_Group, var="log(SAVav)+1" )
SAVMmobfuntrends<-summary(SAVMmobfun,infer=c(TRUE, TRUE))
SAVMmobfun<-as.data.frame(SAVMmobfuntrends$emtrends) 
#write.csv(SAVMmobfun,file="Output/SAVMmobfun.csv")

SAVMmobfun$sig<-c("nonsig","sig","nonsig")
SAVMmobfun<-SAVMmobfun%>%
  rename(SAVav.trend="log(SAVav)+1.trend")
SAVMmobplot<-ggplot(SAVMmobfun,aes(x=Functional_Group, y=SAVav.trend, col = Functional_Group,alpha = sig)) +
  geom_point(size = 7) +
  coord_flip()+
  geom_pointrange(aes(y=SAVav.trend, ymin=lower.CL, ymax=upper.CL), size=2) + #graph estimates and st err for all values
  geom_hline(yintercept = 0, lty = 2,size =1.5)  + #add vertical line to graph and make it dashed 
  scale_colour_manual(values = MobColors) +
  guides(alpha= FALSE, col = FALSE) + #remove alpha legend
  scale_alpha_discrete(range = c(0.1, 5.0)) + #set range for sig vs non sig alphas
  theme_classic() +
  theme(axis.text.x=element_text(size = 30, color = "black"), 
        axis.text.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(face="bold",size = 40, color = "black"))+
  scale_x_discrete(labels=c("Carnivore"="Carnivores", "Herbivore"= "Herbivores", 
                            "Omnivores"="Omnivores")) +
  theme(legend.position="none")+
  labs(y="Effect size" ~(count/m^2), x = "")

THMmobfun<-emtrends(Mytilusmobmod,pairwise ~Functional_Group, var="log(THav)+1" )
THMmobfuntrends<-summary(THMmobfun,infer=c(TRUE, TRUE))
THMmobfunt<-as.data.frame(THMmobfuntrends$emtrends) 
#write.csv(THMmobfunt,file="Output/THMmobfunt.csv")
THMmobfunt$sig<-c("nonsig","sig","nonsig")

THMmobfunt<-THMmobfunt %>%
  rename(THav.trend="log(THav)+1.trend")
THMmobplot<-ggplot(THMmobfunt,aes(x=Functional_Group, y=THav.trend, col = Functional_Group,alpha = sig)) +
  geom_point(size = 7) +
  coord_flip()+
  geom_pointrange(aes(y=THav.trend, ymin=lower.CL, ymax=upper.CL),size=2) + #graph estimates and st err for all values
  geom_hline(yintercept = 0, lty = 2, size =1.5)  + #add vertical line to graph and make it dashed 
  scale_colour_manual(values = MobColors) +
  guides(alpha= FALSE, col = FALSE) + #remove alpha legend
  scale_alpha_discrete(range = c(0.1, 5.0)) + #set range for sig vs non sig alphas
  theme_classic() +
  theme(axis.text.x=element_text(size =30, color = "black"), 
        axis.text.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(face="bold",size = 40, color = "black"))+
  scale_x_discrete(labels=c("Carnivore"="Carnivores", "Herbivore"= "Herbivores", 
                            "Omnivores"="Omnivores")) +
  theme(legend.position="none")+
  labs(y="", x = "")
THMmobplot
#######Surfgrass mobile emtrend plots######
Phyllomobfun<-emtrends(Phyllomobmod,pairwise ~Functional_Group, var="log(Phyllodelta)+1")
Phyllomobfuntrends<-summary(Phyllomobfun,infer=c(TRUE, TRUE))
Phyllomobfunnt<-as.data.frame(Phyllomobfuntrends$emtrends) 
#write.csv(Phyllomobfunnt,file="Output/Phyllomobfunnt.csv")

Phyllomobfunnt$sig<-c("nonsig","nonsig","nonsig")
Phyllomobfunnt<-Phyllomobfunnt %>%
  rename(Phyllodelta.trend="log(Phyllodelta)+1.trend")
Pmobfunplot<-ggplot(Phyllomobfunnt,aes(x=Functional_Group, y=Phyllodelta.trend, col = Functional_Group, alpha = sig)) +
  geom_point(size = 7) +
  coord_flip()+
  geom_pointrange(aes(y=Phyllodelta.trend, ymin=lower.CL, ymax=upper.CL),size=2) + #graph estimates and st err for all values
  geom_hline(yintercept = 0, lty = 2, size = 1.5)  + #add vertical line to graph and make it dashed 
  guides(alpha= FALSE, col = FALSE) + #remove alpha legend
  scale_alpha_discrete(range = c(0.1, 5.0)) + #set range for sig vs non sig alphas
  scale_colour_manual(values = MobColors) +
  theme_classic() +
  theme(axis.text.x=element_text(size = 35, color = "black"), 
        axis.text.y=element_text(size =35, color = "black"),
        axis.title.y = element_text(size = 40, color = "black"),
        axis.title.x = element_text(size = 40, color = "black"))+
  scale_x_discrete(labels=c("Carnivore"="Carnivores", "Herbivore"= "Herbivores", 
                            "Omnivores"="Omnivores")) +
  theme(legend.position="none")+
  labs(y="Effect size (density)", x = "")

SAVPmobfun<-emtrends(Phyllomobmod,pairwise ~Functional_Group, var="log(SAVav)+1" )
SAVPmobfuntrends<-summary(SAVPmobfun,infer=c(TRUE, TRUE))
SAVPmobfunt<-as.data.frame(SAVPmobfuntrends$emtrends) 
#write.csv(SAVPmobfunt,file="Output/SAVPmobfunt.csv")

SAVPmobfunt$sig<-c("nonsig","nonsig","nonsig")
SAVPmobfunt<-SAVPmobfunt%>%
  rename(SAVav.trend="log(SAVav)+1.trend")
SAVPmobplot<-ggplot(SAVPmobfunt,aes(x=Functional_Group, y=SAVav.trend, col = Functional_Group,alpha = sig)) +
  geom_point(size = 7) +
  coord_flip()+
  geom_pointrange(aes(y=SAVav.trend, ymin=lower.CL, ymax=upper.CL), size=2) + #graph estimates and st err for all values
  geom_hline(yintercept = 0, lty = 2, size =1.5)  + #add vertical line to graph and make it dashed 
  scale_colour_manual(values = MobColors) +
  guides(alpha= FALSE, col = FALSE) + #remove alpha legend
  scale_alpha_discrete(range = c(0.1, 5.0)) + #set range for sig vs non sig alphas
  theme_classic() +
  theme(axis.text.x=element_text(size = 30, color = "black"), 
        axis.text.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 40, color = "black"))+
  scale_x_discrete(labels=c("Carnivore"="Carnivores", "Herbivore"= "Herbivores", 
                            "Omnivores"="Omnivores")) +
  theme(legend.position="none")+
  labs(y="Effect size" ~(count/m^2), x = "")

THPmobfun<-emtrends(Phyllomobmod,pairwise ~Functional_Group, var="log(THav)+1" )
THPmobfuntrends<-summary(THPmobfun,infer=c(TRUE, TRUE))
THPmobfunt<-as.data.frame(THPmobfuntrends$emtrends) 
#write.csv(SAVPmobfunt,file="Output/THPmobfunt.csv")
THPmobfunt$sig<-c("nonsig","nonsig","nonsig")
THPmobfunt<-THPmobfunt%>%
  rename(THav.trend="log(THav)+1.trend" )
THPmobplot<-ggplot(THPmobfunt,aes(x=Functional_Group, y=THav.trend, col = Functional_Group,alpha = sig)) +
  geom_point(size = 7) +
  coord_flip()+
  geom_pointrange(aes(y=THav.trend, ymin=lower.CL, ymax=upper.CL), size=2) + #graph estimates and st err for all values
  geom_hline(yintercept = 0, lty = 2, size =1.5)  + #add vertical line to graph and make it dashed 
  scale_colour_manual(values = MobColors) +
  guides(alpha= FALSE, col = FALSE) + #remove alpha legend
  scale_alpha_discrete(range = c(0.1, 5.0)) + #set range for sig vs non sig alphas
  theme_classic() +
  theme(axis.text.x=element_text(size = 30, color = "black"), 
        axis.text.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(face="bold",size = 30, color = "black"))+
  scale_x_discrete(labels=c("Carnivore"="Carnivores", "Herbivore"= "Herbivores", 
                            "Omnivores"="Omnivores")) +
  theme(legend.position="none")+
  labs(y="", x = "")


######patchwork for phyllo and mobile separated plots#####
Phyllosessmobplots<- Psessfunplot+SAVPsessplot+THPsessplot +Pmobfunplot+SAVPmobplot+THPmobplot+#patchwork to combine plots
  plot_annotation(tag_levels = 'A') &         #label each individual plot with letters A-G
  theme(plot.tag = element_text(size = 26, face = "bold"))   #edit the lettered text

Phyllosessmobplots
ggsave(filename = "Output/Phyllosessmobfunplots.pdf", useDingbats =FALSE,dpi=600,device = "pdf", width = 25, height = 20)

Mytilussessmobplots<-Msessfunplot+SAVMsessplot+THMsessplot +Mmobfunplot+SAVMmobplot+THMmobplot+     #patchwork to combine plots
  plot_annotation(tag_levels = 'A') &         #label each individual plot with letters A-G
  theme(plot.tag = element_text(size = 26, face = "bold"))   #edit the lettered text

Mytilussessmobplots
ggsave(filename = "Output/Mytilussessmobfunplots.pdf", useDingbats =FALSE,dpi=600,device = "pdf", width = 25, height = 20)


####Supplemental Figure####
Phyllodist<-Funsppcover %>%
  filter(Foundation_spp =='Phyllospadix') %>%
  ggplot(aes(x=Phyllodelta , fill = Removal_Control)) +
  geom_density(alpha=.5) +
  theme_sleek() +
  scale_fill_manual(values = c("#3182bd","#bdbdbd")) +
  theme(axis.text.x=element_text(size = 25, color = "black"), 
        axis.text.y=element_text(size = 25, color = "black"),
        axis.title.y = element_text(size = 30, color = "black"),
        axis.title.x = element_text(face="italic",size = 30, color = "black"))+
  theme(legend.position="none")+
  labs(y = 'Density of cover', x ='Surfgrass percent loss \n (Phyllospadix spp.)')

Mytilusdist<-Funsppcover %>%
  filter(Foundation_spp =='Mytilus') %>%
  ggplot(aes(x=Mytilusdelta , fill = Removal_Control)) +
  geom_density(alpha=.5) +
  theme_sleek() +
  scale_fill_manual(values = c("#3182bd","#bdbdbd")) +
  theme(axis.text.x=element_text(size = 25, color = "black"), 
        axis.text.y=element_text(size = 25, color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(face="italic",size = 30, color = "black"))+
  theme(legend.position="none")+
  labs(y = 'Density of cover', x ='CA mussel percent loss \n (M. californianus)',fill = 'Control or Removal Pools')


#distplots<-Phyllodist+ Mytilusdist+#patchwork to combine plots
# plot_annotation(tag_levels = 'a') &         #label each individual plot with letters A-G
#theme(plot.tag = element_text(size = 26, face = "bold"))   #edit the lettered text

#distplots
ggsave(filename = "Output/distplots.pdf", useDingbats =FALSE,dpi=600,device = "pdf", width = 25, height = 20)


