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
         Before_After ="Sessiles$Before_After", MusselCover = "PercentSessile$Mytilus.californianus",SurfgrassCover = "PercentSessile$Phyllospadix.spp")
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
  dplyr::select(PoolID,Foundation_spp,Removal_Control,Before_After,MusselCover, SurfgrassCover, Mytilus.californianus, Phyllospadix.spp,allCCA, macroalgae,macrophytes, macroCCA,
                consumers,allconsumers,  prodphyllodom, allproddom)


#Change in foundation species cover
Funsppcover<- Communitymetrics%>%
  filter(Before_After != 'Immediate') %>%
  dplyr::group_by(PoolID,Foundation_spp, Removal_Control) %>%
  summarise(Mytilusdelta = -1*(MusselCover[Before_After == 'After'] - MusselCover[Before_After == 'Before']),
            Phyllodelta = -1*(SurfgrassCover[Before_After == 'After'] - SurfgrassCover[Before_After == 'Before']))

#before then take avg of tide height and sa to v ratio
PhysicalParameters<-TidePooldes[, c("PoolID","Removal_Control","MaxDepth","Perimeter","SurfaceArea", "Vol", "SAtoV", "TideHeight")] #pull out the necessary columns and treatment 

PhysicalParameters$PoolID<-as.factor(PhysicalParameters$PoolID)
Funsppcover$PoolID<-as.factor(Funsppcover$PoolID)
Funsppandpp<-left_join(Funsppcover,PhysicalParameters)

PP<-PhysicalParameters %>%
  dplyr::group_by(PoolID,Removal_Control) %>%
  dplyr::summarise(SAVav = mean(SAtoV), #ave between before and after since SA/V changed with fspp removal
                   THav = mean(TideHeight),SAav=mean(SurfaceArea),Depthav=mean(MaxDepth)) #tide height didn't change 

PP$PoolID<-as.factor(PP$PoolID)
Funsppcover$PoolID<-as.factor(Funsppcover$PoolID)
Funsppandpp<-left_join(Funsppcover,PP)


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

set.seed(267)
PhylloSessiles<-metaMDS(sqrt(sqrt(PhyllonMDS)),k=2, distance='bray', trymax = 50, autotransform = FALSE) #add more iterations

PhylloSessiles$stress #0.1696881
phyllopp<-Funsppandpp%>%
  filter(Foundation_spp =="Phyllospadix")
phyllograph<-cbind(PhyllocommunitynMDS,phyllopp$SAVav,phyllopp$THav)

phyllograph<-phyllograph%>%
  rename(SAV="phyllopp$SAVav",TH="phyllopp$THav")

psessperm<-adonis(log(PhyllonMDS+1)~Phyllodelta+SAV+TH, phyllograph, permutations = 999, 
                               method="bray")
psessperm

#nMDS of surgrass
ordiplot(PhylloSessiles, type = 'text')

PSessilesnMDSpts<-data.frame(PhylloSessiles$points)
Pphypar<-TidePooldes%>%
  filter(Foundation_spp =="Phyllospadix")
PSessilesnMDS<-cbind(PhyllocommunitynMDS$PoolID,
                     PhyllocommunitynMDS$Removal_Control,
                     PhyllocommunitynMDS$Before_After,
                     PhyllocommunitynMDS$Phyllodelta,
                     Pphypar$TideHeight,
                     Pphypar$SAtoV,
                     PSessilesnMDSpts)


colnames(PSessilesnMDS)[1:6]<- c("PoolID","Removal_Control","Before_After","SurfgrassLoss","TideHeight","SAV")

set.seed(267)
permanovaSessilemodel<-adonis(PsessnMDS~SurfgrassLoss +SAV +TideHeight, PSessilesnMDS, permutations = 999, 
                              method="bray")
permanovaSessilemodel
#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#SurfgrassLoss  1    1.2622 1.26223  7.5809 0.18411  0.001 ***
# SAV            1    0.4005 0.40051  2.4054 0.05842  0.031 *  
# TideHeight     1    0.5311 0.53107  3.1896 0.07746  0.007 ** 
# Residuals     28    4.6620 0.16650         0.68001           
#Total         31    6.8558                 1.00000 
PSessilesnMDS$MDS1<-as.numeric(PSessilesnMDS$MDS1)
PSessilesnMDS$MDS2<-as.numeric(PSessilesnMDS$MDS2)


## add a column combining after before and foundation species
SessilesnMDS$AB_F<-factor(paste(SessilesnMDS$Before_After, SessilesnMDS$Foundation_spp))

#create dataframe for centroids with median from x and y axes
centroids <- aggregate(cbind(MDS1,MDS2)~Foundation_spp*Before_After*Removal_Control*AB_F,SessilesnMDS,median)


#separate out diatoms points to make different color
diatomsyes<-SessilesnMDS %>%
  filter(DiatomPresent =="Yes" & Before_After =="After" & Removal_Control =="Removal")

diatomsno<-SessilesnMDS %>%
  filter(DiatomPresent =="No" & Before_After =="After" & Removal_Control =="Removal")

#for arrows fucnction:
#no diatoms
x0 <- centroids %>%
  filter(Before_After == "Before") %>%
  select(MDS1)
x0<-as.matrix(x0)

y0 <- centroids %>%
  filter(Before_After == "Before") %>%
  select(MDS2)
y0<-as.matrix(y0)

x1<-centroids %>%
  filter(Before_After == "After") %>%
  select(MDS1)
x1<-as.matrix(x1)

y1<-centroids %>%
  filter(Before_After == "After") %>%
  select(MDS2)
y1<-as.matrix(y1)

#create dataframe for centroids with median from x and y axes
pcentroids <- aggregate(cbind(NMDS1,NMDS2)~Before_After*Removal_Control,phyllograph,median)
#use median since nonnormal data

#create groupings for shape labels
phyllocommgraph$RC_BA<-factor(paste(phyllocommgraph$Before_After, phyllocommgraph$Removal_Control))

pgroupings<-c("After" = 2,"Before" = 17)

#for arrows fucnction:
x0 <- pcomcentroids %>%
  filter(Before_After=="Before") %>%
  select(NMDS1)
x0<-as.matrix(x0)
y0 <- pcomcentroids %>%
  filter(Before_After == "Before") %>%
  select(NMDS2)
y0<-as.matrix(y0)
x1<-pcomcentroids %>%
  filter(Before_After == "After") %>%
  select(NMDS1)
x1<-as.matrix(x1)
y1<-pcomcentroids %>%
  filter(Before_After == "After") %>%
  select(NMDS2)
y1<-as.matrix(y1)


Surfgrassplot<-ggplot(phyllocommgraph, aes(x = NMDS1 , y= NMDS2,shape = Before_After)) + #basic plot
  geom_point(aes(color =SurfgrassLoss, size =SurfgrassLoss, alpha=3,stroke=2), shape=16) +
  scale_color_distiller(palette = "Greens",guide = "legend")+
  scale_size(range = c(1,15)) +
  geom_point(data=pcomcentroids, size=10, stroke = 2.75) +
  #stat_ellipse(aes(linetype=RC_BA))+
  #ordiellipse(PhylloSessilesafter, phyllocommgraph$RC_BA, display = "sites", 
  #kind = "se", conf = 0.95, label = T)+
  theme_classic() +
  scale_shape_manual(values = c(pgroupings))+
  geom_segment(aes(x = x0[1], y = y0[1], xend = (x1[1]), yend = (y1[1])),size = 1,#segment with arrow for surfgrass before/after control
               colour = "#3182bd", arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  geom_segment(aes(x = x0[2], y = y0[2], xend = (x1[2]), yend = (y1[2])),linetype = 2,size = 1, #segment with arrow for surfgrass before/after removal
               colour = "#bdbdbd",arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  labs(x ='nMDS1', y = 'nMDS2', shape='', color='Surfgrass loss',size='Surfgrass loss', linetype ='Before or after') +
  theme(axis.text = element_text(color = "black", size = 35), 
        axis.title.x = element_text(color="black", size=40), 
        axis.title.y = element_text(color="black", size=40), 
        legend.title = element_text(color="black", size=40), 
        legend.text = element_text(color = "black", size = 35), 
        legend.position= "top",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  guides(shape = "none", alpha = "none")+
  guides(colour = guide_legend(nrow = 1))#makes legend only one row
Surfgrassplot
ggsave("Output/surfgrassprocrustesswtich.pdf",useDingbats = FALSE, width=40, height=30,dpi=600, unit="cm")



####Procrustes graphs IDk what is happening really####

Phyllobefore<-PhyllocommunitynMDS%>%
  filter(Before_After =='Before')
Phylloafter<-PhyllocommunitynMDS%>%
  filter(Before_After =='After')
PsessnMDSbefore<-Phyllobefore[-c(1:10,72:73)]  
PsessnMDSafter<-Phylloafter[-c(1:10,72:73)]  

#create dissimlarity matrices for procrustes analysis
#log transform community data to remove affect of rare spp
set.seed(267)
PhylloSessilesbefore<-metaMDS(log(PsessnMDSbefore+1),k=2, distance='bray', trymax = 50, autotransform = FALSE) #add more iterations

set.seed(267)
PhylloSessilesafter<-metaMDS(log(PsessnMDSafter+1),k=2, distance='bray', trymax = 50, autotransform = FALSE) #add more iterations
pro<-procrustes(PhylloSessilesbefore,PhylloSessilesafter) #procrustes rotation with before as target and after as rotation

Targetp<-as.data.frame(pro$X)
Yrotp<-as.data.frame(pro$Yrot)

Targetp$Before_After<-"After"
Yrotp$Before_After<-"Before"
Yrotp<-Yrotp%>%
  rename(NMDS1="V1", NMDS2= "V2") #rename columns for rbind
                         
combinedprophyllo<-rbind(Targetp,Yrotp)

phyllopp<-Funsppandpp%>%
  filter(Foundation_spp =="Phyllospadix")
phyllocommgraph<-cbind(combinedprophyllo,PhyllocommunitynMDS$PoolID,PhyllocommunitynMDS$Removal_Control,PhyllocommunitynMDS$Phyllodelta,
                 phyllopp$SAVav,phyllopp$THav)
                         
phyllocommgraph<-phyllocommgraph %>%
  rename(PoolID="PhyllocommunitynMDS$PoolID",Removal_Control="PhyllocommunitynMDS$Removal_Control",
         SurfgrassLoss="PhyllocommunitynMDS$Phyllodelta", SAV="phyllopp$SAVav",TideHeight="phyllopp$THav")

#unsure how to add procrustes analysis to permanova
set.seed(267)
permanovaSessilemodel<-adonis(formula=pro~SurfgrassLoss +SAV +TideHeight, phyllocommgraph, permutations = 999, 
                              method="bray")
permanovaSessilemodel




#create dataframe for centroids with median from x and y axes
pcomcentroids <- aggregate(cbind(NMDS1,NMDS2)~Before_After*Removal_Control,phyllocommgraph,median)
#use median since nonnormal data

#create groupings for shape labels
phyllocommgraph$RC_BA<-factor(paste(phyllocommgraph$Before_After, phyllocommgraph$Removal_Control))

pgroupings<-c("After" = 2,"Before" = 17)

#for arrows fucnction:
x0 <- pcomcentroids %>%
  filter(Before_After=="Before") %>%
  select(NMDS1)
x0<-as.matrix(x0)
y0 <- pcomcentroids %>%
  filter(Before_After == "Before") %>%
  select(NMDS2)
y0<-as.matrix(y0)
x1<-pcomcentroids %>%
  filter(Before_After == "After") %>%
  select(NMDS1)
x1<-as.matrix(x1)
y1<-pcomcentroids %>%
  filter(Before_After == "After") %>%
  select(NMDS2)
y1<-as.matrix(y1)


Surfgrassplot<-ggplot(phyllocommgraph, aes(x = NMDS1 , y= NMDS2,shape = Before_After)) + #basic plot
  geom_point(aes(color =SurfgrassLoss, size =SurfgrassLoss, alpha=3,stroke=2), shape=16) +
  scale_color_distiller(palette = "Greens",guide = "legend")+
  scale_size(range = c(1,15)) +
  geom_point(data=pcomcentroids, size=10, stroke = 2.75) +
  #stat_ellipse(aes(linetype=RC_BA))+
  #ordiellipse(PhylloSessilesafter, phyllocommgraph$RC_BA, display = "sites", 
              #kind = "se", conf = 0.95, label = T)+
  theme_classic() +
  scale_shape_manual(values = c(pgroupings))+
  geom_segment(aes(x = x0[1], y = y0[1], xend = (x1[1]), yend = (y1[1])),size = 1,#segment with arrow for surfgrass before/after control
               colour = "#3182bd", arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  geom_segment(aes(x = x0[2], y = y0[2], xend = (x1[2]), yend = (y1[2])),linetype = 2,size = 1, #segment with arrow for surfgrass before/after removal
             colour = "#bdbdbd",arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
 labs(x ='nMDS1', y = 'nMDS2', shape='', color='Surfgrass loss',size='Surfgrass loss', linetype ='Before or after') +
  theme(axis.text = element_text(color = "black", size = 35), 
        axis.title.x = element_text(color="black", size=40), 
        axis.title.y = element_text(color="black", size=40), 
        legend.title = element_text(color="black", size=40), 
        legend.text = element_text(color = "black", size = 35), 
        legend.position= "top",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  guides(shape = "none", alpha = "none")+
  guides(colour = guide_legend(nrow = 1))#makes legend only one row
Surfgrassplot
ggsave("Output/surfgrassprocrustesswtich.pdf",useDingbats = FALSE, width=40, height=30,dpi=600, unit="cm")




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

McombinednMDS<-metaMDS(log(Myilussessspplist+1),k=2, distance='bray', trymax = 50, autotransform = FALSE) #add more iterations
McombinednMDS$stress  #0.166

set.seed(267)
permanovaSessilemodel<-adonis(formula=pro~SurfgrassLoss +SAV +TideHeight, phyllocommgraph, permutations = 999, 
                              method="bray")
permanovaSessilemodel

###procrustes####

Mbefore<-MytiluscommunitynMDS%>%
  filter(Before_After =='Before')
Mafter<-MytiluscommunitynMDS%>%
  filter(Before_After =='After')
MsessnMDSbefore<-Mbefore[-c(1:10,72:73)]  
MsessnMDSafter<-Mafter[-c(1:10,72:73)]  

#create dissimlarity matrices for procrustes analysis
#log transform community data to remove affect of rare spp
set.seed(267)
MSessilesbefore<-metaMDS(log(MsessnMDSbefore+1),k=2, distance='bray', trymax = 50, autotransform = FALSE) #add more iterations

set.seed(267)
MSessileafter<-metaMDS(log(MsessnMDSafter+1),k=2, distance='bray', trymax = 50, autotransform = FALSE) #add more iterations

mpro<-procrustes(MSessilesbefore,MSessileafter) #procrustes rotation with before as target and after as
plot(mpro)
Targetm<-as.data.frame(mpro$X)
Yrotm<-as.data.frame(mpro$Yrot)

Targetm$Before_After<-"Before"
Yrotm$Before_After<-"After"
Yrotm<-Yrotm%>%
  rename(NMDS1="V1", NMDS2= "V2") #rename columns for rbind

combinedpromytilus<-rbind(Targetm,Yrotm)

mpp<-Funsppandpp%>%
  filter(Foundation_spp =="Mytilus")
mcommgraph<-cbind(combinedpromytilus,MytiluscommunitynMDS$PoolID,MytiluscommunitynMDS$Removal_Control,MytiluscommunitynMDS$Mytilusdelta,
                       mpp$SAVav,mpp$THav)

musselcommgraph<-mcommgraph %>%
  rename(PoolID="MytiluscommunitynMDS$PoolID",Removal_Control="MytiluscommunitynMDS$Removal_Control",
        MusselLoss="MytiluscommunitynMDS$Mytilusdelta", SAV="mpp$SAVav",TideHeight="mpp$THav")


#create dataframe for centroids with median from x and y axes
mcomcentroids <- aggregate(cbind(NMDS1,NMDS2)~Before_After*Removal_Control,musselcommgraph,median)
#use median since nonnormal data

#create groupings for shape labels

mgroupings<-c("After" = 2,"Before" = 17)

#for arrows fucnction:
x0 <- mcomcentroids %>%
  filter(Before_After=="Before") %>%
  select(NMDS1)
x0<-as.matrix(x0)
y0 <- mcomcentroids %>%
  filter(Before_After == "Before") %>%
  select(NMDS2)
y0<-as.matrix(y0)
x1<-mcomcentroids %>%
  filter(Before_After == "After") %>%
  select(NMDS1)
x1<-as.matrix(x1)
y1<-mcomcentroids %>%
  filter(Before_After == "After") %>%
  select(NMDS2)
y1<-as.matrix(y1)


musselplot<-ggplot(musselcommgraph, aes(x = NMDS1 , y= NMDS2,shape = Before_After)) + #basic plot
  geom_point(aes(color =MusselLoss, size =MusselLoss, alpha=3,stroke=2), shape=16) +
  scale_color_distiller(palette = "Blues",guide = "legend")+
  scale_size(range = c(1,15)) +
  geom_point(data=mcomcentroids, size=10, stroke = 2.75) +
  theme_classic() +
  scale_shape_manual(values = c(mgroupings))+
  geom_segment(aes(x = x0[1], y = y0[1], xend = (x1[1]), yend = (y1[1])),size = 1,#segment with arrow for surfgrass before/after control
               colour = "#3182bd", arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  geom_segment(aes(x = x0[2], y = y0[2], xend = (x1[2]), yend = (y1[2])),linetype = 2,size = 1, #segment with arrow for surfgrass before/after removal
               colour = "#bdbdbd",arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  # labs(x ='PC1 (32.88%)', y = 'PC2 (27.9%)', shape='', color='Surfgrass loss',size='Surfgrass loss', linetype ='Before or after') +
  #PC1  0.3288PC2 0.2790total 0.6079
  theme(axis.text = element_text(color = "black", size = 35), 
        axis.title.x = element_text(color="black", size=40), 
        axis.title.y = element_text(color="black", size=40), 
        legend.title = element_text(color="black", size=40), 
        legend.text = element_text(color = "black", size = 35), 
        legend.position= "top",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  guides(shape = "none", alpha = "none")+
  guides(colour = guide_legend(nrow = 1))#makes legend only one row
musselplot

#seems to be opposite of what it should be (before vs afer)
#Untransformed data
#Sessiles2D<-metaMDS(SessilesByPool[-c(1:4)],k=2, distance='bray', trymax = 50) #add more iterations

#let's look at the 2D stress. Is it < 0.2? 
PhylloSessiles$stress #0.171503


#stress plot for both transformed and untransformed data
stressplot(PhylloSessiles)

# basic plot
ordiplot(PhylloSessiles) # dots represent tide pools and 
#+ represents species
# add species names
ordiplot(PhylloSessiles, type = 'text')

PSessilesnMDSpts<-data.frame(PhylloSessiles$points)
Pphypar<-TidePooldes%>%
  filter(Foundation_spp =="Phyllospadix")
PSessilesnMDS<-cbind(PhyllocommunitynMDS$PoolID,
                    PhyllocommunitynMDS$Removal_Control,
                    PhyllocommunitynMDS$Before_After,
                    PhyllocommunitynMDS$Phyllodelta,
                    Pphypar$TideHeight,
                    Pphypar$SAtoV,
                    PSessilesnMDSpts)


colnames(PSessilesnMDS)[1:6]<- c("PoolID","Removal_Control","Before_After","SurfgrassLoss","TideHeight","SAV")

set.seed(267)
permanovaSessilemodel<-adonis(PsessnMDS~SurfgrassLoss +SAV +TideHeight, PSessilesnMDS, permutations = 999, 
                              method="bray")
permanovaSessilemodel
#Terms added sequentially (first to last)

              #Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#SurfgrassLoss  1    1.2622 1.26223  7.5809 0.18411  0.001 ***
 # SAV            1    0.4005 0.40051  2.4054 0.05842  0.031 *  
 # TideHeight     1    0.5311 0.53107  3.1896 0.07746  0.007 ** 
 # Residuals     28    4.6620 0.16650         0.68001           
#Total         31    6.8558                 1.00000 
PSessilesnMDS$MDS1<-as.numeric(PSessilesnMDS$MDS1)
PSessilesnMDS$MDS2<-as.numeric(PSessilesnMDS$MDS2)


## add a column combining after before and foundation species
SessilesnMDS$AB_F<-factor(paste(SessilesnMDS$Before_After, SessilesnMDS$Foundation_spp))

#create dataframe for centroids with median from x and y axes
centroids <- aggregate(cbind(MDS1,MDS2)~Foundation_spp*Before_After*Removal_Control*AB_F,SessilesnMDS,median)


#separate out diatoms points to make different color
diatomsyes<-SessilesnMDS %>%
  filter(DiatomPresent =="Yes" & Before_After =="After" & Removal_Control =="Removal")

diatomsno<-SessilesnMDS %>%
  filter(DiatomPresent =="No" & Before_After =="After" & Removal_Control =="Removal")

#for arrows fucnction:
#no diatoms
x0 <- centroids %>%
  filter(Before_After == "Before") %>%
  select(MDS1)
x0<-as.matrix(x0)

y0 <- centroids %>%
  filter(Before_After == "Before") %>%
  select(MDS2)
y0<-as.matrix(y0)

x1<-centroids %>%
  filter(Before_After == "After") %>%
  select(MDS1)
x1<-as.matrix(x1)

y1<-centroids %>%
  filter(Before_After == "After") %>%
  select(MDS2)
y1<-as.matrix(y1)


#create groupings for shape labels
groupings<-c("Before Phyllospadix" = 17,"After Phyllospadix" = 2,"Before Mytilus" = 16,"After Mytilus" = 1)
colors<-c("Removal" = "#bdbdbd", "Control" = "#3182bd")

SessilesnMDSplot<- ggplot(PSessilesnMDS, aes(x = MDS1 , y= MDS2, shape=Before_After,color = SurfgrassLoss, size=SurfgrassLoss)) + #basic plot
  geom_point(stroke = 2) + #geom_point(data=centroids,size=8, stroke = 2.75) + 
  #geom_point(data=diatomsyes, size =3,alpha = 0.6,stroke = 2, color ="#de2d26") +
  #geom_point(data=diatomsno, size =3,alpha = 0.6, stroke =2, color = "#756bb1") +
  theme_classic() +
  #xlim(-0.75, 1) + 
  labs(x ='nMDS1', y = 'nMDS2', color ='SurfgrassLoss') +
 # geom_segment(aes(x = x0[1], y = y0[1], xend = x1[1]+0.01, yend = y1[1]-0.04),size = 1,#segment with arrow for Mussels before/after control
               #colour = "#3182bd", arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  #geom_segment(aes(x = x0[2], y = y0[2], xend = x1[2]+0.01, yend = y1[2]-0.04),linetype = 2,size = 1, #segment with arrow for Mussels before/after removal
               #colour = "#bdbdbd",arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  #geom_segment(aes(x = x0[3], y = y0[3], xend = x1[3], yend = y1[3]-0.04),size = 1, #segment with arrow for phyllospadix before/after control
               #colour ="#3182bd", arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  #geom_segment(aes(x = x0[4], y = y0[4], xend = x1[4], yend = y1[4]-0.04),linetype = 2,size = 1, #segment with arrow for phyllospadix before/after removal
               #colour = "#bdbdbd",arrow = arrow(length = unit(0.3, "cm"), type = "closed")) +
  # theme(legend.text = element_text(size=22, face ="italic"),
  #legend.title = element_text(size = 22),
  #legend.spacing = unit(1, unit = "cm")) +
  theme(legend.position="none") +
  #scale_shape_manual(values=c(groupings)) +
  #scale_color_manual(values = c("#3182bd","#bdbdbd")) +
  #annotate(geom="text", x=0.9, y=-0.2, label="Before Removal",
  #color="black", size = 10, fontface = 2) +
  #annotate(geom="text", x=0.9, y=0.25, label="After Removal",
  #color="black", size = 10, fontface = 2) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(color = "black", size = 18), 
        axis.title.x = element_text(color="black", size=24, face="bold"), 
        axis.title.y = element_text(color="black", size=24, face="bold"), 
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) 
  #ggsave("Output/CommonsppSessileplot.pdf",useDingbats = FALSE, width=40, height=30,dpi=300, unit="cm")
SessilesnMDSplot
######nMDS Mussel######
####Species Richness and diversity anovas#####

Sessilesall<-Communitymetrics%>%
  dplyr::filter(Before_After != "Immediate")
sessspplist<-Sessilesall[-c(1:9,71:72)] #remove grouping parameters and foundation spp.
sessspplist$Richness<-specnumber(sessspplist) #richness
sessspplist$Diversity<-diversity(sessspplist, index = "shannon", MARGIN = 1, base = exp(1)) #diversity

sessspp<-cbind(Sessilesall$PoolID,Sessilesall$Foundation_spp,Sessilesall$Removal_Control,
               Sessilesall$Before_After,sessspplist$Richness,sessspplist$Diversity)
sessspp<-as.data.frame(sessspp)
sessspp<-sessspp%>%
  rename(PoolID = "V1",Foundation_spp ="V2", Removal_Control ="V3", Before_After ="V4",
         Richness ="V5",Diversity="V6")
sessspp$Richness<-as.numeric(sessspp$Richness)
sessspp$Diversity<-as.numeric(sessspp$Diversity)
deltasessspp<-sessspp%>%
  dplyr::group_by(PoolID, Foundation_spp,Removal_Control)%>%
  summarise(DeltaRich = Richness[Before_After =="After"]-Richness[Before_After =="Before"],
            DeltaDiversity=Diversity[Before_After =="After"]-Diversity[Before_After =="Before"])

deltasessrichfunpp<-left_join(deltasessspp,Funsppandpp)

phyllosessrich<-deltasessrichfunpp %>%
  filter(Foundation_spp =="Phyllospadix")

mytilusessrich<-deltasessrichfunpp %>%
  filter(Foundation_spp =="Mytilus")
#richness and diversity phyllo
phyllosessrichmod<-lm(DeltaRich~Phyllodelta +SAav+THav, data=phyllosessrich)
#plot(phyllosessrichmod) #good
qqp(resid(phyllosessrichmod),"norm") #good

anova(phyllosessrichmod)

library(ggeffects)
phyllospp<-ggpredict(phyllosessrichmod, c("Phyllodelta")) #predict marginal effects from model for foundation spp. loss
plot(phyllospp) #plot output 
phyllospp<-as.data.frame(phyllospp) #create dataframe 

phyllospp<-phyllospp %>% #output for values gives you an x for variable. rename variable to match
  rename(Phyllodelta=x) #rename to join to rest of dataframe

phyllospp<-left_join(phyllospp,phyllosessrich) #rejoin with main dataframe for ggplot

#display raw data but prediction line and confidence intervals are from ggpredict model
phyllospp<-ggplot(phyllospp, aes(x =Phyllodelta, y=DeltaRich)) +
  geom_point(size=5,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=Phyllodelta, y=predicted), color="#006d2c",size =2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=40), 
        axis.title.y=element_text(color="black", size=40),
        axis.text.x =element_text(color="black", size=30),
        axis.text.y =element_text(color="black", size=30)) +
  theme(legend.position="none")+
  labs(x ='Surfgrass loss \n (Phyllospadix spp.)', y = 'Delta Species Richness') 
phyllospp


phyllosessrichdiv<-lm(DeltaDiversity~Phyllodelta +SAav+THav, data=phyllosessrich)
#plot(phyllosessrichdiv) #good
qqp(resid(phyllosessrichdiv),"norm") #good

anova(phyllosessrichdiv)
phyllosppd<-ggpredict(phyllosessrichdiv, c("Phyllodelta")) #predict marginal effects from model for foundation spp. loss
plot(phyllosppd) #plot output 
phyllosppd<-as.data.frame(phyllosppd) #create dataframe 

phyllosppd<-phyllosppd %>% #output for values gives you an x for variable. rename variable to match
  rename(Phyllodelta=x) #rename to join to rest of dataframe

phyllosppd<-left_join(phyllosppd,phyllosessrich) #rejoin with main dataframe for ggplot

#display raw data but prediction line and confidence intervals are from ggpredict model
phyllosppd<-ggplot(phyllosppd, aes(x =Phyllodelta, y=DeltaDiversity)) +
  geom_point(size=5,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=Phyllodelta, y=predicted), color="#006d2c",size =2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=40), 
        axis.title.y=element_text(color="black", size=40),
        axis.text.x =element_text(color="black", size=30),
        axis.text.y =element_text(color="black", size=30)) +
  theme(legend.position="none")+
  labs(x ='Surfgrass loss \n (Phyllospadix spp.)', y = 'Delta Species diversity') 
phyllosppd
#richness and diversity mytilus
mytilusessrich$logrichness<-sign(mytilusessrich$DeltaRich)*log(abs(mytilusessrich$DeltaRich))

mytilussessrichmod<-lm(logrichness~Mytilusdelta +SAav+THav, data=mytilusessrich)
#plot(mytilussessrichmod) #good
qqp(resid(mytilussessrichmod),"norm") #one point out log to meet normality
anova(mytilussessrichmod)
mytilussppr<-ggpredict(mytilussessrichmod, c("Mytilusdelta")) #predict marginal effects from model for foundation spp. loss
plot(mytilussppr) #plot output 
mytilussppr<-as.data.frame(mytilussppr) #create dataframe 

mytilussppr<-mytilussppr%>% #output for values gives you an x for variable. rename variable to match
  rename(Mytilusdelta=x) #rename to join to rest of dataframe

mytilussppr<-left_join(mytilussppr,mytilusessrich) #rejoin with main dataframe for ggplot

#display raw data but prediction line and confidence intervals are from ggpredict model
mspr<-ggplot(mytilussppr, aes(x =Mytilusdelta, y=logrichness)) +
  geom_point(size=5,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=Mytilusdelta, y=predicted), color="#006d2c",size =2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=40), 
        axis.title.y=element_text(color="black", size=40),
        axis.text.x =element_text(color="black", size=30),
        axis.text.y =element_text(color="black", size=30)) +
  theme(legend.position="none")+
  labs(x ='Mussel Loss', y = 'Delta Species richness') 
mspr

mytilussessdivmod<-lm(DeltaDiversity~Mytilusdelta +SAav+THav, data=mytilusessrich)
#plot(mytilussessdivmod) #good
qqp(resid(mytilussessdivmod),"norm") #logged to meet normality
mytilussppd<-ggpredict(mytilussessdivmod, c("Mytilusdelta")) #predict marginal effects from model for foundation spp. loss
plot(mytilussppd) #plot output 
mytilussppd<-as.data.frame(mytilussppd) #create dataframe 

mytilussppd<-mytilussppd%>% #output for values gives you an x for variable. rename variable to match
  rename(Mytilusdelta=x) #rename to join to rest of dataframe

mytilussppd<-left_join(mytilussppd,mytilusessrich) #rejoin with main dataframe for ggplot

#display raw data but prediction line and confidence intervals are from ggpredict model
mspd<-ggplot(mytilussppd, aes(x =Mytilusdelta, y=DeltaDiversity)) +
  geom_point(size=5,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=Mytilusdelta, y=predicted), color="#006d2c",size =2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=40), 
        axis.title.y=element_text(color="black", size=40),
        axis.text.x =element_text(color="black", size=30),
        axis.text.y =element_text(color="black", size=30)) +
  theme(legend.position="none")+
  labs(x ='Mussel Loss', y = 'Delta Species diversity') 
mspd
anova(mytilussessdivmod)


#####Sessile functional groups####
SessilesMytilusFun<-Communitymetrics %>%
  dplyr::filter(Before_After != 'Immediate') %>%
  dplyr::group_by(PoolID, Foundation_spp, Before_After, Removal_Control) %>%
  tidyr::pivot_longer(
    cols = Phyllospadix.spp:Stylantheca.spp,
    names_to = "Species", #creates column with species in longformat
    values_to = "Cover", #adds column for % cover
    values_drop_na = TRUE
  ) 

SessilesphylloFun<-Communitymetrics %>%
  dplyr::filter(Before_After != 'Immediate') %>%
  dplyr::group_by(PoolID, Foundation_spp, Before_After, Removal_Control) %>%
  tidyr::pivot_longer(
    cols = c(8,10:70), #excluding phyllo since taking delta cover
    names_to = "Species", #creates column with species in longformat
    values_to = "Cover", #adds column for % cover
    values_drop_na = TRUE
  ) 
####Functional group log ratio graphs#####
SessilesMytilusStacked<-left_join(SessilesMytilusFun,SessilesGroupings) #combine with fun groups 
SessilesPhylloStacked<-left_join(SessilesphylloFun,SessilesGroupings) #combine with fun groups 

deltaMV<-SessilesMytilusStacked %>%
  dplyr::group_by(PoolID,Removal_Control,Foundation_spp,Before_After,Functional_Group) %>%
  summarise(SumCover = sum(Cover)) %>%
  dplyr::group_by(PoolID,Removal_Control,Foundation_spp, Functional_Group) %>%
  mutate(pluoneSumcover= SumCover+1) %>% #get rid of 0s
  summarise(Ratiocover = pluoneSumcover[Before_After=="After"]/pluoneSumcover[Before_After =="Before"]) %>%
  mutate(logRatiocover = log(Ratiocover)) 
  
phyllodeltacover<-SessilesPhylloStacked%>%
  dplyr::group_by(PoolID,Removal_Control,Foundation_spp,Functional_Group, Before_After) %>%
  summarise(SumCover = sum(Cover)) %>%
  dplyr::group_by(PoolID,Removal_Control,Foundation_spp, Functional_Group) %>%
  mutate(pluoneSumcover= SumCover+1) %>% #get rid of 0s
  summarise(Ratiocover = pluoneSumcover[Before_After=="After"]/pluoneSumcover[Before_After =="Before"]) %>%
  mutate(logRatiocover = log(Ratiocover)) 
  

phyllodeltacover$PoolID<-as.factor(phyllodeltacover$PoolID)
deltaMV$PoolID<-as.factor(deltaMV$PoolID)
PdeltaMV<-left_join(phyllodeltacover,Funsppandpp)
MdeltaMV<-left_join(deltaMV,Funsppandpp)

MdeltaMV<-MdeltaMV %>%
  dplyr::filter(Foundation_spp =="Mytilus") %>%
  dplyr::select(PoolID,Removal_Control,Foundation_spp, Functional_Group,Ratiocover,logRatiocover,Mytilusdelta,SAVav,THav,Depthav)
MdeltaMV<-as.data.frame(MdeltaMV)

PdeltaMV<-PdeltaMV %>%
  dplyr::filter(Foundation_spp =="Phyllospadix") %>% 
  dplyr::select(PoolID,Removal_Control,Foundation_spp, Functional_Group,Ratiocover,logRatiocover,Phyllodelta,SAVav,THav,Depthav)
PdeltaMV<-as.data.frame(PdeltaMV)


###log ratio plots phyllo#####
Plogratioplot<-ggplot(PdeltaMV,aes(y=Functional_Group, x=logRatiocover)) +
  geom_point(aes(color=Phyllodelta, size=Phyllodelta)) +
  scale_color_distiller(palette = "Greens",guide = "legend")+
  scale_size(range = c(3,30)) +
  geom_vline(xintercept = 0, lty = 2,size =1.5)  + #add vertical line to graph and make it dashed 
  theme_classic() +
  theme(axis.text.x=element_text(size = 35, color = "black"), 
        axis.text.y=element_text(size = 35, color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 40, color = "black"),
        legend.text = element_text(size = 35, color = "black"),
        legend.title = element_text(size = 35, color = "black"),
        legend.position="bottom",
        plot.title = element_text(size =40, color ="black", hjust = 0.5)) +
  guides(colour = guide_legend(nrow = 1))+ #makes legend only one row
  scale_y_discrete(labels=c("Crustose"="Crustose", "ArticulatedCorallines"= "Articulated corallines", "Microalgae"="Microalgae", "Filamentous"="Filamentous", "FoiloseAlgae"="Foliose algae", 
                            "CorticatedFoliose"="Corticated foliose","CorticatedMacro"="Corticated macroalgae", "LeatheryMacro"="Leathery macrophytes","Anemone"="Anemone","SuspensionFeeder"="Suspension feeders")) +  #rename y axis tickmarks 
  labs(x="Log ratio (percent cover)",y = "",color= "Surfgrass Loss",size ="Surfgrass Loss",title =expression(~italic(Phyllospadix)~'spp. percent loss')) 
Plogratioplot

PSAVlogratioplot<-ggplot(PdeltaMV,aes(y=Functional_Group, x=logRatiocover)) +
  geom_point(aes(color=SAVav, size=SAVav)) +
  scale_color_distiller(palette = "Purples",guide = "legend")+
  scale_size(range = c(3,30)) +
  geom_vline(xintercept = 0, lty = 2,size =1.5)  + #add vertical line to graph and make it dashed 
  theme_classic() +
  theme(axis.text.x=element_text(size = 35, color = "black"), 
        axis.text.y=element_text(size = 35, color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 40, color = "black"),
        legend.text = element_text(size = 35, color = "black"),
        legend.title = element_text(size = 35, color = "black"),
        legend.position="bottom",
        plot.title = element_text(size =40, color ="black", hjust = 0.5)) +
  scale_y_discrete(labels=c("Crustose"="Crustose", "ArticulatedCorallines"= "Articulated corallines", "Microalgae"="Microalgae", "Filamentous"="Filamentous", "FoiloseAlgae"="Foliose algae", 
                            "CorticatedFoliose"="Corticated foliose","CorticatedMacro"="Corticated macroalgae", "LeatheryMacro"="Leathery macrophytes","Anemone"="Anemone","SuspensionFeeder"="Suspension feeders")) +  #rename y axis tickmarks 
  labs(x="Log ratio (percent cover)",y = "",color= "SA:V",size ="SA:V",title ="Surface area to volume (SA:V)")
PSAVlogratioplot

PTHlogratioplot<-ggplot(PdeltaMV,aes(y=Functional_Group, x=logRatiocover)) +
  geom_point(aes(color=THav, size=THav)) +
  scale_color_distiller(palette = "Oranges",guide = "legend")+
  scale_size(range = c(3,30)) +
  geom_vline(xintercept = 0, lty = 2,size =1.5)  + #add vertical line to graph and make it dashed 
  theme_classic() +
  theme(axis.text.x=element_text(size = 35, color = "black"), 
        axis.text.y=element_text(size = 35, color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 40, color = "black"),
        legend.text = element_text(size = 35, color = "black"),
        legend.title = element_text(size = 35, color = "black"),
        legend.position="bottom",
        plot.title = element_text(size =40, color ="black", hjust = 0.5)) +
  scale_y_discrete(labels=c("Crustose"="Crustose", "ArticulatedCorallines"= "Articulated corallines", "Microalgae"="Microalgae", "Filamentous"="Filamentous", "FoiloseAlgae"="Foliose algae", 
                            "CorticatedFoliose"="Corticated foliose","CorticatedMacro"="Corticated macroalgae", "LeatheryMacro"="Leathery macrophytes","Anemone"="Anemone","SuspensionFeeder"="Suspension feeders")) +  #rename y axis tickmarks 
  labs(x="Log ratio (percent cover)",y = "",color= "Tide Height (m)",size ="Tide Height (m)",title ="Tide Height (m)")
PTHlogratioplot

###mussel log ratio plots###
Mlogratioplot<-ggplot(MdeltaMV,aes(y=Functional_Group, x=logRatiocover)) +
  geom_point(aes(color=Mytilusdelta, size=Mytilusdelta)) +
  scale_color_distiller(palette = "Blues",guide = "legend")+
  scale_size(range = c(3,30)) +
  geom_vline(xintercept = 0, lty = 2,size =1.5)  + #add vertical line to graph and make it dashed 
  theme_classic() +
  theme(axis.text.x=element_text(size = 35, color = "black"), 
        axis.text.y=element_text(size = 35, color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 40, color = "black"),
        legend.text = element_text(size = 35, color = "black"),
        legend.title = element_text(size = 35, color = "black"),
        legend.position="bottom",
        plot.title = element_text(size =40, color ="black", hjust = 0.5)) +
  scale_y_discrete(labels=c("Crustose"="Crustose", "ArticulatedCorallines"= "Articulated corallines", "Microalgae"="Microalgae", "Filamentous"="Filamentous", "FoiloseAlgae"="Foliose algae", 
                            "CorticatedFoliose"="Corticated foliose","CorticatedMacro"="Corticated macroalgae", "LeatheryMacro"="Leathery macrophytes","Anemone"="Anemone","SuspensionFeeder"="Suspension feeders")) +  #rename y axis tickmarks 
  labs(x="Log ratio (percent cover)",y = "",color= "Mussel Loss",size ="Mussel Loss",title =expression(~italic(M.)~''~italic(californianus)~ 'percent loss')) 
Mlogratioplot

MSAVlogratioplot<-ggplot(MdeltaMV,aes(y=Functional_Group, x=logRatiocover)) +
  geom_point(aes(color=SAVav, size=SAVav)) +
  scale_color_distiller(palette = "Purples",guide = "legend")+
  scale_size(range = c(3,30)) +
  geom_vline(xintercept = 0, lty = 2,size =1.5)  + #add vertical line to graph and make it dashed 
  theme_classic() +
  theme(axis.text.x=element_text(size = 35, color = "black"), 
        axis.text.y=element_text(size = 35, color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 40, color = "black"),
        legend.text = element_text(size = 35, color = "black"),
        legend.title = element_text(size = 35, color = "black"),
        legend.position="bottom",
        plot.title = element_text(size =40, color ="black", hjust = 0.5)) +
  scale_y_discrete(labels=c("Crustose"="Crustose", "ArticulatedCorallines"= "Articulated corallines", "Microalgae"="Microalgae", "Filamentous"="Filamentous", "FoiloseAlgae"="Foliose algae", 
                            "CorticatedFoliose"="Corticated foliose","CorticatedMacro"="Corticated macroalgae", "LeatheryMacro"="Leathery macrophytes","Anemone"="Anemone","SuspensionFeeder"="Suspension feeders")) +  #rename y axis tickmarks 
  labs(x="Log ratio (percent cover)",y = "",color= "SA:V",size ="SA:V",title ="Surface area to volume (SA:V)")
MSAVlogratioplot

MTHlogratioplot<-ggplot(MdeltaMV,aes(y=Functional_Group, x=logRatiocover)) +
  geom_point(aes(color=THav, size=THav)) +
  scale_color_distiller(palette = "Oranges",guide = "legend")+
  scale_size(range = c(3,30)) +
  geom_vline(xintercept = 0, lty = 2,size =1.5)  + #add vertical line to graph and make it dashed 
  theme_classic() +
  theme(axis.text.x=element_text(size = 35, color = "black"), 
        axis.text.y=element_text(size = 35, color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 40, color = "black"),
        legend.text = element_text(size = 35, color = "black"),
        legend.title = element_text(size = 35, color = "black"),
        legend.position="bottom",
        plot.title = element_text(size =40, color ="black", hjust = 0.5)) +
  scale_y_discrete(labels=c("Crustose"="Crustose", "ArticulatedCorallines"= "Articulated corallines", "Microalgae"="Microalgae", "Filamentous"="Filamentous", "FoiloseAlgae"="Foliose algae", 
                            "CorticatedFoliose"="Corticated foliose","CorticatedMacro"="Corticated macroalgae", "LeatheryMacro"="Leathery macrophytes","Anemone"="Anemone","SuspensionFeeder"="Suspension feeders")) +  #rename y axis tickmarks 
  labs(x="Log ratio (percent cover)",y = "",color= "Tide Height (m)",size ="Tide Height (m)",title ="Tide Height (m)")
MTHlogratioplot

#####Mobiles######

####Species richness & diversity anovas####
MobileBA<-Mobiles %>%
  dplyr::filter(Before_After != 'Immediate') 
Mobilelist<-MobileBA[-c(1:4)]
Mobilelist$Richness<-specnumber(Mobilelist) #richness
Mobilelist$H<-diversity(Mobilelist, index = "shannon", MARGIN = 1, base = exp(1)) #diversity
Mobilespprich<-cbind(MobileBA$PoolID,MobileBA$Foundation_spp,MobileBA$Removal_Control,
                     MobileBA$Before_After,Mobilelist$Richness,Mobilelist$H)
Mobilespprich<-as.data.frame(Mobilespprich)
Mobilespprich<-Mobilespprich%>%
  rename(PoolID = "V1",Foundation_spp ="V2", Removal_Control ="V3", Before_After ="V4",
         Richness ="V5",Diversity="V6")
Mobilespprich$Richness<-as.numeric(Mobilespprich$Richness)
Mobilespprich$Diversity<-as.numeric(Mobilespprich$Diversity)#change richness to numeric
deltamobrich<-Mobilespprich %>%
  dplyr::group_by(PoolID, Foundation_spp,Removal_Control)%>%
  summarise(DeltaRich = Richness[Before_After =="After"]-Richness[Before_After =="Before"],
            DeltaDiversity=Diversity[Before_After =="After"]-Diversity[Before_After =="Before"])

deltamobrichfunpp<-left_join(deltamobrich,Funsppandpp)

phyllomobrich<-deltamobrichfunpp %>%
  filter(Foundation_spp =="Phyllospadix")

mytilusmobrich<-deltamobrichfunpp %>%
  filter(Foundation_spp =="Mytilus")
#richness and diversity phyllo
phyllomobrichmod<-lm(DeltaRich~Phyllodelta +SAav+THav, data=phyllomobrich)
#plot(phyllomobrichmod) #good
qqp(resid(phyllomobrichmod),"norm") #good

phyllosppmobr<-ggpredict(phyllomobrichmod, c("Phyllodelta")) #predict marginal effects from model for foundation spp. loss
#plot(phyllosppmobr) #plot output 
phyllosppmobr<-as.data.frame(phyllosppmobr) #create dataframe 

phyllosppmobr<-phyllosppmobr%>% #output for values gives you an x for variable. rename variable to match
  rename(Phyllodelta=x) #rename to join to rest of dataframe

phyllosppmobr<-left_join(phyllosppmobr,phyllomobrich) #rejoin with main dataframe for ggplot

#display raw data but prediction line and confidence intervals are from ggpredict model
phyllospprm<-ggplot(phyllosppmobr, aes(x =Phyllodelta, y=DeltaRich)) +
  geom_point(size=5,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=Phyllodelta, y=predicted), color="#006d2c",size =2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=40), 
        axis.title.y=element_text(color="black", size=40),
        axis.text.x =element_text(color="black", size=30),
        axis.text.y =element_text(color="black", size=30)) +
  theme(legend.position="none")+
  labs(x ='Surfgrass loss \n (Phyllospadix spp.)', y = 'Delta Species richness') 
phyllospprm
anova(phyllomobrichmod)

phyllomobrichdiv<-lm(DeltaDiversity~Phyllodelta +SAav+THav, data=phyllomobrich)
#plot(phyllomobrichdiv) #good
qqp(resid(phyllomobrichdiv),"norm") #good

phyllosppmobd<-ggpredict(phyllomobrichdiv, c("Phyllodelta")) #predict marginal effects from model for foundation spp. loss
plot(phyllosppmobd) #plot output 
phyllosppmobd<-as.data.frame(phyllosppmobd) #create dataframe 

phyllosppmobd<-phyllosppmobd%>% #output for values gives you an x for variable. rename variable to match
  rename(Phyllodelta=x) #rename to join to rest of dataframe

phyllosppmobd<-left_join(phyllosppmobd,phyllomobrich) #rejoin with main dataframe for ggplot

#display raw data but prediction line and confidence intervals are from ggpredict model
phyllosppdm<-ggplot(phyllosppmobd, aes(x =Phyllodelta, y=DeltaDiversity)) +
  geom_point(size=5,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=Phyllodelta, y=predicted), color="#006d2c",size =2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=40), 
        axis.title.y=element_text(color="black", size=40),
        axis.text.x =element_text(color="black", size=30),
        axis.text.y =element_text(color="black", size=30)) +
  theme(legend.position="none")+
  labs(x ='Surfgrass loss \n (Phyllospadix spp.)', y = 'Delta Species diversity') 
phyllosppdm
anova(phyllomobrichdiv)

#richness and diversity mytilus
mytilusmobrichmod<-lm(DeltaRich~Mytilusdelta +SAav+THav, data=mytilusmobrich)
#plot(mytilusmobrichmod) #good
qqp(resid(mytilusmobrichmod),"norm") #good
mytilussppmr<-ggpredict(mytilusmobrichmod, c("Mytilusdelta")) #predict marginal effects from model for foundation spp. loss
plot(mytilussppmr) #plot output 
mytilussppmr<-as.data.frame(mytilussppmr) #create dataframe 

mytilussppmr<-mytilussppmr%>% #output for values gives you an x for variable. rename variable to match
  rename(Mytilusdelta=x) #rename to join to rest of dataframe

mytilussppr<-left_join(mytilussppmr,mytilusmobrich) #rejoin with main dataframe for ggplot

#display raw data but prediction line and confidence intervals are from ggpredict model
mspmr<-ggplot(mytilussppr, aes(x =Mytilusdelta, y=DeltaRich)) +
  geom_point(size=5,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=Mytilusdelta, y=predicted), color="#006d2c",size =2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=40), 
        axis.title.y=element_text(color="black", size=40),
        axis.text.x =element_text(color="black", size=30),
        axis.text.y =element_text(color="black", size=30)) +
  theme(legend.position="none")+
  labs(x ='Mussel Loss', y = 'Delta Species richness') 
mspmr
anova(mytilusmobrichmod)

mytilusmobrich$logdiversity<-sign(mytilusmobrich$DeltaDiversity)*log(abs(mytilusmobrich$DeltaDiversity))
mytilusmobdivmod<-lm(logdiversity~Mytilusdelta +SAav+THav, data=mytilusmobrich)
#plot(mytilusmobdivmod) #good
qqp(resid(mytilusmobdivmod),"norm") #logged to meet normality
mytilussppmd<-ggpredict(mytilusmobdivmod, c("Mytilusdelta")) #predict marginal effects from model for foundation spp. loss
plot(mytilussppmd) #plot output 
mytilussppmd<-as.data.frame(mytilussppmd) #create dataframe 

mytilussppmd<-mytilussppmd%>% #output for values gives you an x for variable. rename variable to match
  rename(Mytilusdelta=x) #rename to join to rest of dataframe

mytilussppd<-left_join(mytilussppmd,mytilusmobrich) #rejoin with main dataframe for ggplot

#display raw data but prediction line and confidence intervals are from ggpredict model
mspmd<-ggplot(mytilussppd, aes(x =Mytilusdelta, y=logdiversity)) +
  geom_point(size=5,aes(shape=Removal_Control),stroke=2) +
  scale_shape_manual(values = c(19,1)) +
  geom_line(aes(x=Mytilusdelta, y=predicted), color="#006d2c",size =2)+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.2) +
  theme_classic()+
  theme(axis.title.x=element_text(color="black", size=40), 
        axis.title.y=element_text(color="black", size=40),
        axis.text.x =element_text(color="black", size=30),
        axis.text.y =element_text(color="black", size=30)) +
  theme(legend.position="none")+
  labs(x ='Mussel Loss', y = 'log Delta Species diversity') 
mspmd

anova(mytilusmobdivmod)
####log ratio data####
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


LRMobdata<- Mobdata%>%
  dplyr::group_by(PoolID,Removal_Control,Foundation_spp,Before_After,Functional_Group) %>%
  summarise(SumDensity = sum(Std.Count)) %>%
  dplyr::group_by(PoolID,Removal_Control,Foundation_spp, Functional_Group) %>%
  mutate(plusonesumdens= SumDensity+1) %>% #get rid of 0s
  summarise(Ratiocount = plusonesumdens[Before_After=="After"]/plusonesumdens[Before_After =="Before"]) %>%
  mutate(logRatiocount =log(Ratiocount))

LRMobdata<-left_join(LRMobdata,Funsppandpp)

LRMobdata<-as.data.frame(LRMobdata)

MytilusMobmod<-LRMobdata%>%
  dplyr::filter(Foundation_spp == "Mytilus" & Functional_Group != "SuspensionFeeder") #no suspension feeders in mytilus pools

PhylloMobmod<-LRMobdata%>%
  filter(Foundation_spp == "Phyllospadix"& Functional_Group != "SuspensionFeeder") 
#removed suspension feeder since only one in one tide pool in one time period 

####log ratio mobile phyllo####
Pmobfunplot<-ggplot(PhylloMobmod,aes(y=Functional_Group, x=logRatiocount)) +
  geom_point(aes(color=Phyllodelta, size=Phyllodelta)) +
  scale_color_distiller(palette = "Greens",guide = "legend")+
  scale_size(range = c(3,30)) +
  geom_vline(xintercept = 0, lty = 2, size = 1.5)  + #add vertical line to graph and make it dashed 
  theme_classic() +
  theme(axis.text.x=element_text(size = 35, color = "black"), 
        axis.text.y=element_text(size =35, color = "black"),
        axis.title.y = element_text(size = 40, color = "black"),
        axis.title.x = element_text(size = 40, color = "black"),
        legend.position='none')+
  scale_y_discrete(labels=c("Carnivore"="Carnivores", "Herbivore"= "Herbivores", 
                            "Omnivores"="Omnivores")) +
  labs(x="Log ratio (density" ~count/m^2~")", y= "",color = "Surfgrass Loss", size = "Surfgrass Loss")
Pmobfunplot

PSAVmobfunplot<-ggplot(PhylloMobmod,aes(y=Functional_Group, x=logRatiocount)) +
  geom_point(aes(color=SAVav, size=SAVav)) +
  scale_color_distiller(palette = "Purples",guide = "legend")+
  scale_size(range = c(3,30)) +
  geom_vline(xintercept = 0, lty = 2, size = 1.5)  + #add vertical line to graph and make it dashed 
  theme_classic() +
  theme(axis.text.x=element_text(size = 35, color = "black"), 
        axis.text.y=element_text(size =35, color = "black"),
        axis.title.y = element_text(size = 40, color = "black"),
        axis.title.x = element_text(size = 40, color = "black"),
        legend.position='none')+
  scale_y_discrete(labels=c("Carnivore"="Carnivores", "Herbivore"= "Herbivores", 
                            "Omnivores"="Omnivores")) +
  labs(x="Log ratio (density" ~count/m^2~")", y= "",color = "SA:V", size = "SA:V")
PSAVmobfunplot

PTHmobfunplot<-ggplot(PhylloMobmod,aes(y=Functional_Group, x=logRatiocount)) +
  geom_point(aes(color=THav, size=THav)) +
  scale_color_distiller(palette = "Oranges",guide = "legend")+
  scale_size(range = c(3,30)) +
  geom_vline(xintercept = 0, lty = 2, size = 1.5)  + #add vertical line to graph and make it dashed 
  theme_classic() +
  theme(axis.text.x=element_text(size = 35, color = "black"), 
        axis.text.y=element_text(size =35, color = "black"),
        axis.title.y = element_text(size = 40, color = "black"),
        axis.title.x = element_text(size = 40, color = "black"),
        legend.position='none')+
  scale_y_discrete(labels=c("Carnivore"="Carnivores", "Herbivore"= "Herbivores", 
                            "Omnivores"="Omnivores")) +
  labs(x="Log ratio (density" ~count/m^2~")", y= "",color = "Tide Height (m)", size = "Tide Height (m)")
PTHmobfunplot
###log ratio mobile mussels####
Mmobfunplot<-ggplot(MytilusMobmod,aes(y=Functional_Group, x=logRatiocount)) +
  geom_point(aes(color=Mytilusdelta, size=Mytilusdelta)) +
  scale_color_distiller(palette = "Blues",guide = "legend")+
  scale_size(range = c(3,30)) +
  geom_vline(xintercept = 0, lty = 2, size = 1.5)  + #add vertical line to graph and make it dashed 
  theme_classic() +
  theme(axis.text.x=element_text(size = 35, color = "black"), 
        axis.text.y=element_text(size =35, color = "black"),
        axis.title.y = element_text(size = 40, color = "black"),
        axis.title.x = element_text(size = 40, color = "black"),
        legend.position='none')+
  scale_y_discrete(labels=c("Carnivore"="Carnivores", "Herbivore"= "Herbivores", 
                            "Omnivores"="Omnivores")) +
  labs(x="Log ratio (density" ~count/m^2~")", y= "",color = "Mussel Loss", size = "Mussel Loss")
Mmobfunplot

MSAVmobfunplot<-ggplot(MytilusMobmod,aes(y=Functional_Group, x=logRatiocount)) +
  geom_point(aes(color=SAVav, size=SAVav)) +
  scale_color_distiller(palette = "Purples",guide = "legend")+
  scale_size(range = c(3,30)) +
  geom_vline(xintercept = 0, lty = 2, size = 1.5)  + #add vertical line to graph and make it dashed 
  theme_classic() +
  theme(axis.text.x=element_text(size = 35, color = "black"), 
        axis.text.y=element_text(size =35, color = "black"),
        axis.title.y = element_text(size = 40, color = "black"),
        axis.title.x = element_text(size = 40, color = "black"),
        legend.position='none')+
  scale_y_discrete(labels=c("Carnivore"="Carnivores", "Herbivore"= "Herbivores", 
                            "Omnivores"="Omnivores")) +
  labs(x="Log ratio (density" ~count/m^2~")", y= "",color = "SA:V", size = "SA:V")
MSAVmobfunplot

MTHmobfunplot<-ggplot(MytilusMobmod,aes(y=Functional_Group, x=logRatiocount)) +
  geom_point(aes(color=THav, size=THav)) +
  scale_color_distiller(palette = "Oranges",guide = "legend")+
  scale_size(range = c(3,30)) +
  geom_vline(xintercept = 0, lty = 2, size = 1.5)  + #add vertical line to graph and make it dashed 
  theme_classic() +
  theme(axis.text.x=element_text(size = 35, color = "black"), 
        axis.text.y=element_text(size =35, color = "black"),
        axis.title.y = element_text(size = 40, color = "black"),
        axis.title.x = element_text(size = 40, color = "black"),
        legend.position='none')+
  scale_y_discrete(labels=c("Carnivore"="Carnivores", "Herbivore"= "Herbivores", 
                            "Omnivores"="Omnivores")) +
  labs(x="Log ratio (density" ~count/m^2~")", y= "",color = "Tide Height (m)", size = "Tide Height (m)")
MTHmobfunplot

######patchwork for phyllo and mobile separated plots#####
Phylloplots<-Plogratioplot+PSAVlogratioplot+PTHlogratioplot+Pmobfunplot+PSAVmobfunplot+PTHmobfunplot+#patchwork to combine plots
  plot_annotation(tag_levels = 'A') &         #label each individual plot with letters A-G
  theme(plot.tag = element_text(size = 26, face = "bold"))   #edit the lettered text

Phylloplots
ggsave(filename = "Output/Phyllologratioplots.pdf", useDingbats =FALSE,dpi=600,device = "pdf", width = 45, height = 35)

Mytilusplots<-Mlogratioplot+MSAVlogratioplot+MTHlogratioplot +Mmobfunplot+MSAVmobfunplot+MTHmobfunplot+     #patchwork to combine plots
  plot_annotation(tag_levels = 'A') &         #label each individual plot with letters A-G
  theme(plot.tag = element_text(size = 26, face = "bold"))   #edit the lettered text
Mytilusplots
ggsave(filename = "Output/Mytiluslogratioplots.pdf", useDingbats =FALSE,dpi=600,device = "pdf", width = 45, height = 35)

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


