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
library(RColorBrewer)
devtools::install_github("seananderson/ggsidekick")
library(ggsidekick) #for theme sleek
library(effects) #for ggeffect
library(ggeffects) #marginal effects on complex models
library(margins)
library(emmeans)
library(ggrepel)
library(parameters)
library(tidyverse)


#load data from community comp surveys
#source scripts
source("scripts/tidepoolphysicalparameters.R")
#data
Sessiles <- read_csv("Data/CommunityComposition/SessilesAll.csv")
Mobiles <- read_csv("Data/CommunityComposition/Mobiles.csv")
SessilesGroupings<-read_csv("Data/CommunityComposition/SessilesFunwithcitations.csv")
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
StandardizedSessile<- 100*PercentSessile/rowSums(PercentSessile)


Communitymetrics<-cbind(Sessiles$PoolID, Sessiles$Foundation_spp, Sessiles$Removal_Control, 
                        Sessiles$Before_After,StandardizedSessile,PercentSessile$Mytilus.californianus,PercentSessile$Phyllospadix.spp)
Communitymetrics <- Communitymetrics %>%
  rename(PoolID = "Sessiles$PoolID", Foundation_spp = "Sessiles$Foundation_spp",Removal_Control ="Sessiles$Removal_Control",
         Before_After ="Sessiles$Before_After", AdjMusselCover=Mytilus.californianus, AdjSurfgrassCover=Phyllospadix.spp, MusselCover = "PercentSessile$Mytilus.californianus",SurfgrassCover = "PercentSessile$Phyllospadix.spp")
#rename joined columns


SEMcommunitydata<-Communitymetrics %>%
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
                 macrophytes = (AdjSurfgrassCover + macroalgae), #includes phyllospadix & macroalgae
                 macroCCA = (macroalgae + allCCA), #includes macroalgae and CCA
                 consumers = (Chthamalus +	Semibalanus.cariosus +	Balanus.nibulis	+ Balanus.glandula +
                                Pollicipes.polymerus +	tube.worm	+ Ophlitaspongia.pennata + Halichondria +	
                                Haliclona.permollis	+ Anthropluera.elegantissima	+ Anthropluera.xanthogrammica	+ 
                                Urticina.coriacea	+ Epiactis.prolifera +Anthopleura.artemisia	+ Stylantheca.spp),
                 #includes all consumers (no mytilus)
                 allconsumers = (consumers + AdjMusselCover), #consumers and mytilus
                 prodphyllodom = (macroalgae - (consumers + AdjMusselCover)), #producer dominance for phyllo model (so subtract mytilus too)
                 allproddom = (macrophytes - allconsumers)) %>%#prod dominance with foundation spp
  dplyr::select(PoolID,Foundation_spp,Removal_Control,Before_After,AdjMusselCover, AdjSurfgrassCover, MusselCover, SurfgrassCover,allCCA, macroalgae,macrophytes, macroCCA,
                consumers,allconsumers,  prodphyllodom, allproddom)


#Change in foundation species cover
Funsppcover<- Communitymetrics%>%
  filter(Before_After != 'Immediate') %>%
  dplyr::group_by(PoolID,Foundation_spp, Removal_Control) %>%
  summarise(Mytilusdelta = -1*(MusselCover[Before_After == 'After'] - MusselCover[Before_After == 'Before']),
            Phyllodelta = -1*(SurfgrassCover[Before_After == 'After'] - SurfgrassCover[Before_After == 'Before']))


PP<-TidePooldes %>%
  dplyr::group_by(PoolID,Removal_Control) %>%
  dplyr::summarise(SAVav = mean(SAtoV), #ave between before and after since SA/V changed with fspp removal
                   THav = mean(TideHeight),SAav=mean(SurfaceArea),Vav=mean(Vol),Depthav=mean(MaxDepth),
                   loggerdepth=mean(LoggerDepth)) #tide height didn't change 

PP$PoolID<-as.factor(PP$PoolID)
Funsppcover$PoolID<-as.factor(Funsppcover$PoolID)
Funsppandpp<-left_join(Funsppcover,PP)

#ggpairs(Funsppandpp[c(4:9)])
#mussel and surfgrass are not correlated with volume; t
#take residuals of mussel loss and tide height 


#####nMDS Surfgrass######
PhyllocommunitynMDS<-Communitymetrics%>%
  dplyr::filter(Before_After != "Immediate" & Foundation_spp == "Phyllospadix") 
MytiluscommunitynMDS<-Communitymetrics%>%
  dplyr::filter(Before_After != "Immediate" & Foundation_spp == "Mytilus") 
#create dataframe that has both before (baseline so low number in this case 1) and after
#surfgrass loss that can be used in plot later as color/size variable 
PFunsppcover<-Funsppcover%>%
  filter(Foundation_spp =="Phyllospadix") 
PFunsppafter<-PFunsppcover[c(1,5)]
PFunsppafter$Before_After<- "After"
PBeforespploss<-PhyllocommunitynMDS%>%
  select(PoolID,Before_After) %>%
  filter(Before_After =='Before')
PBeforespploss$Phyllodelta<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) #create column of baseline of phyllodelta
Phylloloss<-rbind(PBeforespploss,PFunsppafter) #combine before and after together
PhyllocommunitynMDS$PoolID<-as.character(PhyllocommunitynMDS$PoolID)
PhyllocommunitynMDS<-left_join(Phylloloss,PhyllocommunitynMDS) #combine with rest of dataframe by pool id

PhyllonMDS<-PhyllocommunitynMDS[-c(1:10,72:73)]  

PhyllonMDS<-PhyllonMDS%>%
  dplyr::select_if(colSums(.) != 0) #remove columns with 0s (spp found in only mussel pools)

set.seed(267)
PhylloSessiles<-metaMDS(sqrt(sqrt(PhyllonMDS)),k=2, distance='bray', trymax = 50, autotransform = FALSE) #add more iterations

PhylloSessiles$stress #0.17


#ordplot in ggplot
psdata.scores <- as.data.frame(scores(PhylloSessiles))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
psdata.scores$labels <- rownames(psdata.scores)  # create a column of site names, from the rownames of data.scores
psspecies.scores <- as.data.frame(scores(PhylloSessiles, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
psspecies.scores$Species <- rownames(psspecies.scores)  # create a column of species, from the rownames of species.scores

Mostabun<-PhyllonMDS%>%
  dplyr::select_if(colSums(.) > 50) #get top 8 spp (over 50% of all pools)

#colors for all sesssile spp graphs
Colors<-c(
  Anemone="#54278f",
  ArticulatedCorallines="#c51b8a",
  Crustose="#d4b9da",
  SuspensionFeeder="#636363",
  Microalgae = "#6baed6",
  Filamentous = "#c7e9c0", 
  FoiloseAlgae="#a1d99b",
  CorticatedFoliose="#31a354",
  CorticatedMacro="#a50f15",
  LeatheryMacro="#006d2c")

phyllospp<-left_join(psspecies.scores,SessilesGroupings) #combine with functional groups

#label for geom_label function
subset<-phyllospp%>%
  filter(Species =="Algae.film"|Species == "Diatoms" | Species == "Chaetomorpha.linum"| 
           Species == "Ptilota.spp"|Species == "Analipus.japonicus"|Species=="Crustose.coralline"|
           Species=="Noncoralline.crust"|Species=="Bosiella.spp"|Species=="Calliathron.tuberculosum"|
           Species=="Ulva.spp"|Species=="Smithura.naiadum") 
subset$Species<-c("Diatoms","CCA","Algae film","Noncoralline crust","Bossiella spp","A. japonicus",
                  "C. linum","Ptilota spp","Ulva spp","Smithora naiadum") #rename to easier labels

ordSesSurf<-ggplot(phyllospp)+
                   geom_point(aes(x=NMDS1,y=NMDS2,color=Functional_Group),size=8,shape =15) + 
  geom_label_repel(data=subset,aes(x=NMDS1,y=NMDS2,label=Species),
                   direction=c("both"),nudge_y=0.2,color="#006d2c",size = 12) +  # add the species labels
  #geom_text(data=psspecies.scores,aes(x=NMDS1,y=NMDS2,label=Species),color="#045a8d",size =8) +  # add the species labels
  scale_color_manual(values=Colors,guide = "legend",labels =c("Anemone","Articulated corallines","Corticated foliose","Corticated macroalgae",
                                                              "Crustose","Filamentous","Foliose","Leathery macrophytes","Microalgae","Suspension feeders"))+
  labs(color="Functional group")+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=6)),size=FALSE)+
  theme(axis.text.x = element_blank(),  # remove x-axis text
                    axis.text.y = element_blank(), # remove y-axis text
                    axis.ticks = element_blank(),  # remove axis ticks
                    axis.title.x = element_blank(), # remove x-axis labels
                    axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")
ordSesSurf       

phyllopp<-Funsppandpp%>%
  filter(Foundation_spp =="Phyllospadix")
phyllograph<-cbind(PhyllocommunitynMDS,phyllopp$Vav,phyllopp$THav)


phyllograph<-phyllograph%>%
  rename(Vol="phyllopp$Vav",TH="phyllopp$THav")

psessperm<-adonis(sqrt(sqrt(PhyllonMDS))~Phyllodelta+Vol+TH, phyllograph, permutations = 999, 
                               method="bray")
psessperm

####Phyllo Sessile nMDS graph####
PSessilesnMDSpts<-data.frame(PhylloSessiles$points)

PSessilesnMDSgraph<-cbind(phyllograph,
                     PSessilesnMDSpts)

PSessilesnMDSgraph$MDS1<-as.numeric(PSessilesnMDSgraph$MDS1)

PSessilesnMDSgraph$MDS2<-as.numeric(PSessilesnMDSgraph$MDS2)

#create dataframe for centroids with median from x and y axes
pscentroids <- aggregate(cbind(MDS1,MDS2)~Before_After*Removal_Control,
                           PSessilesnMDSgraph,median)

#for arrows fucnction:
x0 <- pscentroids %>%
  filter(Before_After == "Before") %>%
  select(MDS1)
x0<-as.matrix(x0)

y0 <- pscentroids %>%
  filter(Before_After == "Before") %>%
  select(MDS2)
y0<-as.matrix(y0)

x1<-pscentroids %>%
  filter(Before_After == "After") %>%
  select(MDS1)
x1<-as.matrix(x1)

y1<-pscentroids %>%
  filter(Before_After == "After") %>%
  select(MDS2)
y1<-as.matrix(y1)

pgroupings<-c("After" = 2,"Before" = 17)

Surfgrasssessilesplot<-ggplot(PSessilesnMDSgraph, aes(x = MDS1 , y= MDS2,shape = Before_After)) + #basic plot
  geom_point(aes(color =Phyllodelta, size =Phyllodelta, alpha=3,stroke=2), shape=16) +
  scale_color_distiller(palette = "Greens",guide = "legend")+
  scale_size(range = c(1,15)) +
  geom_point(data=pscentroids, size=10, stroke = 2.75) +
  theme_classic() +
  scale_shape_manual(values = c(pgroupings))+
  geom_segment(aes(x = x0[1], y = y0[1], xend = (x1[1]), yend = (y1[1])),size = 1,#segment with arrow for surfgrass before/after control
               colour = "#3182bd", arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  geom_segment(aes(x = x0[2], y = y0[2], xend = (x1[2]), yend = (y1[2])),linetype = 2,size = 1, #segment with arrow for surfgrass before/after removal
               colour = "#bdbdbd",arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  labs(x ='', y = 'Sessile community nMDS2', shape='', color='',size='', linetype ='Before or after') +
  annotate("text",  size=14, x=0.7, y=0.75, label= "2D stress = 0.12") +
  theme(axis.text = element_blank(), 
        axis.title.x = element_text(color="black", size=50), 
        axis.title.y = element_text(color="black", size=50), 
        axis.ticks = element_blank(),
        legend.title = element_text(color="black", size=40), 
        legend.text = element_text(color = "black", size = 40), 
        legend.position= "none",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  guides(shape = "none", alpha = "none")+
  guides(colour = guide_legend(nrow = 1))#makes legend only one row
Surfgrasssessilesplot

#ggsave("Output/surfgrassprocrustesswtich.pdf",useDingbats = FALSE, width=40, height=30,dpi=600, unit="cm")

####Mussel Sessile nMDS#####

#create dataframe that has both before (baseline so low number in this case 1) and after
#surfgrass loss that can be used in plot later as color/size variable 
MFunsppcover<-Funsppcover%>%
  filter(Foundation_spp =="Mytilus") 
MFunsppafter<-MFunsppcover[c(1,4)]
MFunsppafter$Before_After<- "After"
MBeforespploss<-MytiluscommunitynMDS%>%
  select(PoolID,Before_After) %>%
  filter(Before_After =='Before')
MBeforespploss$Mytilusdelta<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) #create column of baseline of phyllodelta
Musselloss<-rbind(MBeforespploss,MFunsppafter) #combine before and after together
MytiluscommunitynMDS$PoolID<-as.character(MytiluscommunitynMDS$PoolID)
MytiluscommunitynMDS<-left_join(Musselloss,MytiluscommunitynMDS) #combine with rest of dataframe by pool id

Myilussessspplist<-MytiluscommunitynMDS[-c(1:10,72:73)] 
Myilussessspplist<-Myilussessspplist%>%
  dplyr::select_if(colSums(.) != 0) #remove columns with 0s (spp found in surfgrass pools)

set.seed(267)
McombinednMDS<-metaMDS(sqrt(sqrt(Myilussessspplist)),k=2, distance='bray', trymax = 50, autotransform  = FALSE) #add more iterations
McombinednMDS$stress  #0.1736432

Mostabunm<-Myilussessspplist%>%
  dplyr::select_if(colSums(.) > 50) #get top 7 spp (over 50% of all pools)
#nMDS of surgrass
#ordiplot in ggplot
msdata.scores <- as.data.frame(scores(McombinednMDS))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
msdata.scores$site <- rownames(msdata.scores)  # create a column of site names, from the rownames of data.scores
msspecies.scores <- as.data.frame(scores(McombinednMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
msspecies.scores$Species <- rownames(msspecies.scores)  # create a column of species, from the rownames of species.scores

musselspp<-left_join(msspecies.scores,SessilesGroupings)

labels<-musselspp%>%
  filter(Species =="Diatoms"|Species =="Crustose.coralline"|Species =="Algae.film"|Species =="Noncoralline.crust"|Species =="Odonthalia.floccosa"|
           Species =="Chthamalus"|Species =="Anthropluera.xanthogrammica")

labels$Species<-c("Diatoms","CCA","Algae film","Noncoralline crust","Odonthalia floccosa","Chthamalus spp","Anthopluera xanthogrammica")
ordSesMussel<-ggplot(musselspp) + 
  geom_point(aes(x=NMDS1,y=NMDS2,color=Functional_Group),size=8, shape =15) + 
  #geom_text(data=msspecies.scores,aes(x=NMDS1,y=NMDS2,label=species),color="#045a8d",size =8) +  # add all species labels
  geom_label_repel(data=labels,aes(x=NMDS1,y=NMDS2,label=Species),
                   direction=c("both"),nudge_y=0.3,color="#045a8d",size = 12) +  # add the species labels
  scale_color_manual(values=Colors,guide = "legend",labels =c("Anemone","Articulated corallines","Corticated foliose","Corticated macroalgae",
                                                              "Crustose","Filamentous","Foliose","Leathery macrophytes","Microalgae","Suspension feeders"))+
  labs(color="Functional group")+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=6)),size=FALSE)+
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_blank(), # remove x-axis labels
        axis.title.y = element_blank(),
        legend.position = "none")
ordSesMussel

mpp<-Funsppandpp%>%
  filter(Foundation_spp =="Mytilus")
mgraph<-cbind(MytiluscommunitynMDS,mpp$Vav,mpp$THav)

mgraph<-mgraph%>%
  rename(Vol="mpp$Vav",TH="mpp$THav")


msessperm<-adonis(sqrt(sqrt(Myilussessspplist))~Mytilusdelta+Vol+TH, mgraph, permutations = 999, 
                  method="bray")
msessperm

####mussel sessile nMDS plot####

MSessilesnMDSpts<-data.frame(McombinednMDS$points)

MSessilesnMDSgraph<-cbind(mgraph,
                          MSessilesnMDSpts)

MSessilesnMDSgraph$MDS1<-as.numeric(MSessilesnMDSgraph$MDS1)

MSessilesnMDSgraph$MDS2<-as.numeric(MSessilesnMDSgraph$MDS2)

#create dataframe for centroids with median from x and y axes
mscentroids <- aggregate(cbind(MDS1,MDS2)~Before_After*Removal_Control,
                         MSessilesnMDSgraph,median)

#for arrows fucnction:
v0 <- mscentroids %>%
  filter(Before_After == "Before") %>%
  select(MDS1)
v0<-as.matrix(v0)

z0 <- mscentroids %>%
  filter(Before_After == "Before") %>%
  select(MDS2)
z0<-as.matrix(z0)

v1<-mscentroids %>%
  filter(Before_After == "After") %>%
  select(MDS1)
v1<-as.matrix(v1)

z1<-mscentroids %>%
  filter(Before_After == "After") %>%
  select(MDS2)
z1<-as.matrix(z1)

mgroupings<-c("After" = 2,"Before" = 17)

Musselsessilesplot<-ggplot(MSessilesnMDSgraph, aes(x = MDS1 , y= MDS2,shape = Before_After)) + #basic plot
  geom_point(aes(color =Mytilusdelta, size =Mytilusdelta, alpha=3,stroke=2), shape=16) +
  scale_color_distiller(palette = "Blues",guide = "legend")+
  scale_size(range = c(1,15)) +
  geom_point(data=mscentroids, size=10, stroke = 2.75) +
  theme_classic() +
  scale_shape_manual(values = c(mgroupings))+
  geom_segment(aes(x = v0[1], y = z0[1], xend = (v1[1]), yend = (z1[1])),size = 1,#segment with arrow for surfgrass before/after control
               colour = "#3182bd", arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  geom_segment(aes(x = v0[2], y = z0[2], xend = (v1[2]), yend = (z1[2])),linetype = 2,size = 1, #segment with arrow for surfgrass before/after removal
               colour = "#bdbdbd",arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  labs(x ='', y = 'Sessile community nMDS2', shape='', color='',size='', linetype ='Before or after') +
  annotate("text",  size=14, x=0.68, y=0.6, label= "2D stress = 0.17") +
  theme(axis.text = element_blank(), 
        axis.title.x = element_text(color="black", size=50), 
        axis.title.y = element_text(color="black", size=50), 
        legend.title = element_text(color="black", size=40), 
        axis.ticks = element_blank(),
        legend.text = element_text(color = "black", size = 40), 
        legend.position= "none",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  guides(shape = "none", alpha = "none")+
  guides(colour = guide_legend(nrow = 1))#makes legend only one row
Musselsessilesplot

######Mobile nMDS######
 
SurfgrassMobiles <-Mobiles %>%
  filter(Before_After !="Immediate" & Foundation_spp =="Phyllospadix")


MusselMobiles <-Mobiles %>%
  filter(Before_After !="Immediate" & Foundation_spp =="Mytilus") 
  

SurfgrassMobiles$PoolID<-as.character(SurfgrassMobiles$PoolID)
PhyllomobnMDS<-left_join(Phylloloss,SurfgrassMobiles) #combine with rest of dataframe by pool id

PhyllomoblistnMDS<-PhyllomobnMDS[-c(1:5)]  

Phyllomoblist<-PhyllomoblistnMDS%>%
  select_if(colSums(.) != 0) #remove columns with 0s (spp found in mussel pools)

set.seed(267)
Phyllomob<-metaMDS(sqrt(sqrt(Phyllomoblist)),k=2, distance='bray', trymax = 50, autotransform = FALSE) #add more iterations

Phyllomob$stress #0.2304031

Mostabunmobp<-Phyllomoblist%>%
  dplyr::select_if(colSums(.) > 100) #get top 9 spp 

#nMDS of surgrass
#ordplot in ggplot
pmdata.scores <- as.data.frame(scores(Phyllomob))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
pmdata.scores$site <- rownames(pmdata.scores)  # create a column of site names, from the rownames of data.scores
pmspecies.scores <- as.data.frame(scores(Phyllomob, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
pmspecies.scores$Species <- rownames(pmspecies.scores)  # create a column of species, from the rownames of species.scores

#mobile functional group colors
MobColors<-c(
  Carnivore = "#c6dbef",
  Herbivore = "#67a9cf",
  Omnivores ="#016c59")
phyllomobspp<-left_join(pmspecies.scores,MobileGroupings)
phyllomobspp<-phyllomobspp%>%
  filter(Functional_Group != "SuspensionFeeder")

plabels<-phyllomobspp%>%
  filter(Species =="Tonicella.spp"|Species =="Littorina.spp"|Species =="Tegula.funebralis"|Species =="Lottia.spp"|Species =="Pagarus.spp"|
           Species =="Pugettia.producta"|Species =="Heptacarpus.sitchensis"|Species =="Strongylocentrotus.purpuratus"|Species =="Sculpin")

plabels$Species<-c("Tonicella spp","Littorina spp","Tegula funebralis","Lottia spp","Pagurus spp","Pugettia producta","Heptacarpus sitchensis",
                   "Strongylocentrotus purpuratus","Sculpin")


ordMobSurf<-ggplot(phyllomobspp) + 
  #geom_text(data=pmspecies.scores,aes(x=NMDS1,y=NMDS2,label=species),color="#006d2c",size = 8) +  # add the species labels
  geom_point(aes(x=NMDS1,y=NMDS2,color=Functional_Group),size=8,shape=15) +
  geom_label_repel(data=plabels,aes(x=NMDS1,y=NMDS2,label=Species),
                   direction=c("both"),nudge_y=0.3,color="#006d2c",size = 12) +  # add the species labels
  scale_color_manual(values=MobColors,guide = "legend",labels =c("Carnivores","Herbivores","Omnivores"))+
  labs(x= "nMDS1",color="Functional group")+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=6)),size=FALSE)+
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(color="black", size=40),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")
ordMobSurf

phyllomobgraph<-cbind(PhyllomobnMDS,phyllopp$Vav,phyllopp$THav)

phyllomobgraph<-phyllomobgraph%>%
  rename(Vol="phyllopp$Vav",TH="phyllopp$THav")

pmobperm<-adonis(sqrt(sqrt(Phyllomoblist))~Phyllodelta+Vol+TH,phyllomobgraph, permutations = 999, 
                  method="bray")
pmobperm
####Phyllo mobile nMDS graph####
PmobnMDSpts<-data.frame(Phyllomob$points)

PmobnMDSgraph<-cbind(phyllomobgraph,
                     PmobnMDSpts)

PmobnMDSgraph$MDS1<-as.numeric(PmobnMDSgraph$MDS1)
PmobnMDSgraph$MDS2<-as.numeric(PmobnMDSgraph$MDS2)

#create dataframe for centroids with median from x and y axes
pmcentroids <- aggregate(cbind(MDS1,MDS2)~Before_After*Removal_Control,
                         PmobnMDSgraph,median)

#for arrows fucnction:
a0 <- pmcentroids %>%
  filter(Before_After == "Before") %>%
  select(MDS1)
a0<-as.matrix(a0)

b0 <- pmcentroids %>%
  filter(Before_After == "Before") %>%
  select(MDS2)
b0<-as.matrix(b0)

a1<-pmcentroids %>%
  filter(Before_After == "After") %>%
  select(MDS1)
a1<-as.matrix(a1)

b1<-pmcentroids %>%
  filter(Before_After == "After") %>%
  select(MDS2)
b1<-as.matrix(b1)


Surfgrassmobplot<-ggplot(PmobnMDSgraph, aes(x = MDS1 , y= MDS2,shape = Before_After)) + #basic plot
  geom_point(aes(color =Phyllodelta, size =Phyllodelta, alpha=3,stroke=2), shape=16) +
  scale_color_distiller(palette = "Greens",guide = "legend")+
  scale_size(range = c(1,15)) +
  geom_point(data=pmcentroids, size=10, stroke = 2.75) +
  theme_classic() +
  scale_shape_manual(values = c(pgroupings))+
  geom_segment(aes(x = a0[1], y = b0[1], xend = (a1[1]), yend = (b1[1])),size = 1,#segment with arrow for surfgrass before/after control
               colour = "#3182bd", arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  geom_segment(aes(x = a0[2], y = b0[2], xend = (a1[2]), yend = (b1[2])),linetype = 2,size = 1, #segment with arrow for surfgrass before/after removal
               colour = "#bdbdbd",arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  labs(x ='nMDS1', y = 'Mobile community nMDS2', shape='', color='Surfgrass loss',size='Surfgrass loss', linetype ='Before or after') +
  annotate("text",  size=14, x=0.45, y=0.75, label= "2D stress = 0.23") +
  theme(axis.text = element_blank(), 
        axis.title.x = element_text(color="black", size=50), 
        axis.title.y = element_text(color="black", size=50), 
        axis.ticks = element_blank(),  # remove axis ticks
        legend.title = element_text(color="black", size=40), 
        legend.text = element_text(color = "black", size = 40), 
        legend.position= "top",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  guides(shape = "none", alpha = "none")+
  guides(colour = guide_legend(nrow = 1))#makes legend only one row
Surfgrassmobplot
#####Mussel mobile nMDS#####
MusselMobiles$PoolID<-as.character(MusselMobiles$PoolID)
MytilusmobnMDS<-left_join(Musselloss,MusselMobiles) #combine with rest of dataframe by pool id

Mytilusmob<-MytilusmobnMDS[-c(1:5)] 
Mytilusmobspplist<-Mytilusmob%>%
  select_if(colSums(.) != 0) #remove columns with 0s (spp found in surfgrass pools)
set.seed(267)
MmobnMDS<-metaMDS(sqrt(sqrt(Mytilusmobspplist)),k=2, distance='bray', trymax = 50, autotransform  = FALSE) #add more iterations
MmobnMDS$stress  #.1962671

Mostabunmobm<-Mytilusmobspplist%>%
  dplyr::select_if(colSums(.) > 100) #get top 9 spp 

#ordiplot in ggplot
mmdata.scores <- as.data.frame(scores(MmobnMDS))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
mmdata.scores$site <- rownames(mmdata.scores)  # create a column of site names, from the rownames of data.scores
mmspecies.scores <- as.data.frame(scores(MmobnMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
mmspecies.scores$Species <- rownames(mmspecies.scores)  # create a column of species, from the rownames of species.scores

musselmobspp<-left_join(mmspecies.scores,MobileGroupings)

labelsm<-musselmobspp%>%
  filter(Species =="Nucella.ostrina"|Species =="Callistoma.canaliculata"|Species =="Lottorina.spp"|Species =="Tegula.funebralis"|Species =="Pagarus.spp"|
           Species =="Lottia.spp"|Species =="Strongylocentrotus.purpuratus"|Species =="Sculpin")

labelsm$Species<-c("Nucella ostrina","Callistoma canaliculata","Tegula funebralis","Lottia spp","Pagurus spp",
                   "Strongylocentrotus purpuratus","Sculpin")

ordMobMussel<-ggplot(musselmobspp) + 
  geom_point(aes(x=NMDS1,y=NMDS2,color=Functional_Group),size=8, shape=15) +
  scale_color_manual(values=MobColors,guide = "legend",labels =c("Carnivores","Herbivores","Omnivores"))+
  geom_label_repel(data=labelsm,aes(x=NMDS1,y=NMDS2,label=Species),
                   direction=c("both"),nudge_y=0.1,color="#045a8d",size = 12) +  # add the species labels
  labs(x= "nMDS1",color="Functional group")+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=6)),size=FALSE)+
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(color="black", size=40),
        axis.title.y = element_blank(),
        legend.title = element_blank(), 
        legend.position = "none")
ordMobMussel

mmobgraph<-cbind(MytilusmobnMDS,mpp$Vav,mpp$THav)

mmobgraph<-mmobgraph%>%
  rename(Vol="mpp$Vav",TH="mpp$THav")

set.seed(267)
mmobperm<-adonis(sqrt(sqrt(Mytilusmobspplist))~Mytilusdelta+Vol+TH,mmobgraph, permutations = 999, 
                 method="bray")
mmobperm
######Mussel mobile nMDS graph#####
MmobnMDSpts<-data.frame(MmobnMDS$points)

MmobnMDSgraph<-cbind(mmobgraph,
                     MmobnMDSpts)

MmobnMDSgraph$MDS1<-as.numeric(MmobnMDSgraph$MDS1)

MmobnMDSgraph$MDS2<-as.numeric(MmobnMDSgraph$MDS2)

#create dataframe for centroids with median from x and y axes
mmcentroids <- aggregate(cbind(MDS1,MDS2)~Before_After*Removal_Control,
                         MmobnMDSgraph,median)

#for arrows fucnction:
c0 <- mmcentroids %>%
  filter(Before_After == "Before") %>%
  select(MDS1)
c0<-as.matrix(c0)

d0 <- mmcentroids  %>%
  filter(Before_After == "Before") %>%
  select(MDS2)
d0<-as.matrix(d0)

c1<-mmcentroids  %>%
  filter(Before_After == "After") %>%
  select(MDS1)
c1<-as.matrix(c1)

d1<-mmcentroids  %>%
  filter(Before_After == "After") %>%
  select(MDS2)
d1<-as.matrix(d1)

Musselmobplot<-ggplot(MmobnMDSgraph, aes(x = MDS1 , y= MDS2,shape = Before_After)) + #basic plot
  geom_point(aes(color =Mytilusdelta, size =Mytilusdelta, alpha=3,stroke=2), shape=16) +
  scale_color_distiller(palette = "Blues",guide = "legend")+
  scale_size(range = c(1,15)) +
  geom_point(data=mmcentroids, size=10, stroke = 2.75) +
  theme_classic() +
  scale_shape_manual(values = c(mgroupings))+
  geom_segment(aes(x = c0[1], y = d0[1], xend = (c1[1]), yend = (d1[1])),size = 1,#segment with arrow for surfgrass before/after control
               colour = "#3182bd", arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  geom_segment(aes(x = c0[2], y = d0[2], xend = (c1[2]), yend = (d1[2])),linetype = 2,size = 1, #segment with arrow for surfgrass before/after removal
               colour = "#bdbdbd",arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  labs(x ='nMDS1', y = 'Mobile community nMDS2', shape='', color='Mussel loss',size='Mussel loss', linetype ='Before or after') +
  annotate("text",  size=14, x=0.45, y=0.6, label= "2D stress = 0.20") +
  theme(axis.text = element_blank(), 
        axis.title.x = element_text(color="black", size=50), 
        axis.title.y = element_text(color="black", size=50), 
        axis.ticks = element_blank(),  # remove axis ticks
        legend.title = element_text(color="black", size=40), 
        legend.text = element_text(color = "black", size = 40), 
        legend.position= "top",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  guides(shape = "none", alpha = "none")+
  guides(colour = guide_legend(nrow = 1))#makes legend only one row
Musselmobplot

#patchwork of ord plots and nMDS per foudnation spp

Mcommgraphs<-Musselsessilesplot+ordSesMussel+Musselmobplot+ordMobMussel+
  plot_annotation(tag_levels = 'a') &         #label each individual plot with letters A-G
  theme(plot.tag = element_text(size = 50, face = "bold"))   #edit the lettered text
Mcommgraphs
ggsave("Output/McommnMDS.pdf",useDingbats = FALSE, width=75, height=60,dpi=600, unit="cm")

Pcommgraphs<-Surfgrasssessilesplot+ordSesSurf+Surfgrassmobplot+ordMobSurf+
  plot_annotation(tag_levels = 'a') &         #label each individual plot with letters A-G
  theme(plot.tag = element_text(size = 50, face = "bold"))   #edit the lettered text
ggsave("Output/ScommnMDS.pdf",useDingbats = FALSE, width=75, height=60,dpi=600, unit="cm")

####Species Richness and diversity anovas#####

Sessilesall<-Communitymetrics%>%
  dplyr::filter(Before_After != "Immediate")
sessspplist<-Sessilesall[-c(1:9,71:72)] #remove grouping parameters and foundation spp.
sessspplist$Richness<-specnumber(sessspplist) #richness

sessspp<-cbind(Sessilesall$PoolID,Sessilesall$Foundation_spp,Sessilesall$Removal_Control,
               Sessilesall$Before_After,sessspplist$Richness)
sessspp<-as.data.frame(sessspp)
sessspp<-sessspp%>%
  rename(PoolID = "V1",Foundation_spp ="V2", Removal_Control ="V3", Before_After ="V4",
         Richness ="V5")
sessspp$Richness<-as.numeric(sessspp$Richness)

deltasessspp<-sessspp%>%
  dplyr::group_by(PoolID, Foundation_spp,Removal_Control)%>%
  summarise(DeltaRich = Richness[Before_After =="After"]-Richness[Before_After =="Before"])

deltasessrichfunpp<-left_join(deltasessspp,Funsppandpp)

phyllosessrich<-deltasessrichfunpp %>%
  filter(Foundation_spp =="Phyllospadix")

mytilusessrich<-deltasessrichfunpp %>%
  filter(Foundation_spp =="Mytilus")
#ggpairs(mytilusessrich[c(5:11)])

#ggpairs(phyllosessrich[c(6:11)])
#richness and diversity phyllo
phyllosessrichmod<-lm(DeltaRich~Phyllodelta +Vav+THav, data=phyllosessrich)
#plot(phyllosessrichmod) #good
qqp(resid(phyllosessrichmod),"norm") #good

summary(phyllosessrichmod)

library(ggeffects)
phyllospp<-ggpredict(phyllosessrichmod, c("Phyllodelta")) #predict marginal effects from model for foundation spp. loss
plot(phyllospp) #plot output 
phyllospp<-as.data.frame(phyllospp) #create dataframe 

phyllospp<-phyllospp %>% #output for values gives you an x for variable. rename variable to match
  rename(Phyllodelta=x) #rename to join to rest of dataframe

phyllospp<-left_join(phyllospp,phyllosessrich) #rejoin with main dataframe for ggplot

#display raw data but prediction line and confidence intervals are from ggpredict model
phyllospp<-ggplot(phyllospp, aes(x =Phyllodelta, y=DeltaRich)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  #geom_line(aes(x=Phyllodelta, y=predicted), color="#006d2c",size =2,linetype=2)+
  #geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs(x ='', y = 'Change in sessile species richness') 
phyllospp


#richness and diversity mytilus
mytilusessrich$logrichness<-sign(mytilusessrich$DeltaRich)*log(abs(mytilusessrich$DeltaRich))

mytilussessrichmod<-lm(logrichness~Mytilusdelta +Vav+THav, data=mytilusessrich)
#plot(mytilussessrichmod) #good
qqp(resid(mytilussessrichmod),"norm") #one point out log to meet normality

summary(mytilussessrichmod)
mytilussppr<-ggpredict(mytilussessrichmod, c("Mytilusdelta")) #predict marginal effects from model for foundation spp. loss
plot(mytilussppr) #plot output 
mytilussppr<-as.data.frame(mytilussppr) #create dataframe 

mytilussppr<-mytilussppr%>% #output for values gives you an x for variable. rename variable to match
  rename(Mytilusdelta=x) #rename to join to rest of dataframe

mytilussppr<-left_join(mytilussppr,mytilusessrich) #rejoin with main dataframe for ggplot
mytilussppr<-mytilussppr%>%
  mutate(transpredict=exp(predicted),rich=exp(logrichness),md=exp(Mytilusdelta), trancl= exp(conf.low),tranch=exp(conf.high))
#library(confidence)
#confidence::backtransform(-0.02305,type=c("log"))
#display raw data but prediction line and confidence intervals are from ggpredict model
mspr<-ggplot(mytilussppr, aes(x =Mytilusdelta, y=logrichness)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=Mytilusdelta, y=predicted), color="#045a8d",size =2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs(x ='', y = '') 
mspr

######Mobiles Species richness anovas######
MobileBA<-Mobiles %>%
  dplyr::filter(Before_After != 'Immediate') 
Mobilelist<-MobileBA[-c(1:4)]
Mobilelist$Richness<-specnumber(Mobilelist) #richness
Mobilespprich<-cbind(MobileBA$PoolID,MobileBA$Foundation_spp,MobileBA$Removal_Control,
                     MobileBA$Before_After,Mobilelist$Richness,Mobilelist$H)
Mobilespprich<-as.data.frame(Mobilespprich)
Mobilespprich<-Mobilespprich%>%
  rename(PoolID = "V1",Foundation_spp ="V2", Removal_Control ="V3", Before_After ="V4",
         Richness ="V5")
Mobilespprich$Richness<-as.numeric(Mobilespprich$Richness)
deltamobrich<-Mobilespprich %>%
  dplyr::group_by(PoolID, Foundation_spp,Removal_Control)%>%
  summarise(DeltaRich = Richness[Before_After =="After"]-Richness[Before_After =="Before"])

deltamobrichfunpp<-left_join(deltamobrich,Funsppandpp)

phyllomobrich<-deltamobrichfunpp %>%
  filter(Foundation_spp =="Phyllospadix")

mytilusmobrich<-deltamobrichfunpp %>%
  filter(Foundation_spp =="Mytilus")
#richness and diversity phyllo
phyllomobrichmod<-lm(DeltaRich~Phyllodelta +Vav+THav, data=phyllomobrich)
#plot(phyllomobrichmod) #good
qqp(resid(phyllomobrichmod),"norm") #good
summary(phyllomobrichmod)

phyllosppmobr<-ggpredict(phyllomobrichmod, c("Phyllodelta")) #predict marginal effects from model for foundation spp. loss
#plot(phyllosppmobr) #plot output 
phyllosppmobr<-as.data.frame(phyllosppmobr) #create dataframe 

phyllosppmobr<-phyllosppmobr%>% #output for values gives you an x for variable. rename variable to match
  rename(Phyllodelta=x) #rename to join to rest of dataframe

phyllosppmobr<-left_join(phyllosppmobr,phyllomobrich) #rejoin with main dataframe for ggplot

#display raw data but prediction line and confidence intervals are from ggpredict model
phyllospprm<-ggplot(phyllosppmobr, aes(x =Phyllodelta, y=DeltaRich)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  #geom_line(aes(x=Phyllodelta, y=predicted), color="#006d2c",size =2,linetype=2)+
 # geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40)) +
  theme(legend.position="none")+
  labs(x ='Surfgrass loss \n (Phyllospadix spp.)', y = 'Change in moblie species richness') 
phyllospprm


#richness and diversity mytilus
mytilusmobrichmod<-lm(DeltaRich~Mytilusdelta +Vav+THav, data=mytilusmobrich)
#plot(mytilusmobrichmod) #good
qqp(resid(mytilusmobrichmod),"norm") #good
summary(mytilusmobrichmod)
mytilussppmr<-ggpredict(mytilusmobrichmod, c("Mytilusdelta")) #predict marginal effects from model for foundation spp. loss
plot(mytilussppmr) #plot output 
mytilussppmr<-as.data.frame(mytilussppmr) #create dataframe 

mytilussppmr<-mytilussppmr%>% #output for values gives you an x for variable. rename variable to match
  rename(Mytilusdelta=x) #rename to join to rest of dataframe

mytilussppr<-left_join(mytilussppmr,mytilusmobrich) #rejoin with main dataframe for ggplot

#display raw data but prediction line and confidence intervals are from ggpredict model
mspmr<-ggplot(mytilussppr, aes(x =Mytilusdelta, y=DeltaRich)) +
  geom_point(size=8,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  #geom_line(aes(x=Mytilusdelta, y=predicted), color="#045a8d",size =2, linetype=2)+
  #geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=50), 
        axis.title.y=element_text(color="black", size=50),
        axis.text.x =element_text(color="black", size=40),
        axis.text.y =element_text(color="black", size=40))+
  theme(legend.position="none")+
  labs(x ='CA mussel loss \n (Mytilus californianus)', y = '') 
mspmr



#patchwork of richness plots
richpm<-phyllospp+mspr+phyllospprm+mspmr+
  plot_annotation(tag_levels = 'a') &         #label each individual plot with letters A-D
  theme(plot.tag = element_text(size = 40, face = "bold"))   #edit the lettered text
richpm
ggsave(filename = "Output/richnessplots.pdf", useDingbats =FALSE,dpi=600,device = "pdf",width = 28, height = 30)


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


