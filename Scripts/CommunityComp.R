##Community Composition data for Sessile and Mobiles from Surfgrass and Mussel
##Oregon Tide Pools
##By: Jenn Fields
##Last updated: 10.3.2021

###########################################
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
library(GGally)
library(ggrepel)
library(ecodist) #for pcoa
library(ape)
library(janitor)
library(parameters)
library(tidyverse)


#source scripts
source("scripts/tidepoolphysicalparameters.R")
#load data from community comp surveys
Sessiles <- read_csv("Data/CommunityComposition/SessilesAll.csv")
Mobiles <- read_csv("Data/CommunityComposition/Mobiles.csv")
SessilesGroupings<-read_csv("Data/CommunityComposition/SessilesFunwithcitations.csv")
MobileGroupings<-read_csv("Data/CommunityComposition/MobilesFun.csv")


#replace NA with 0
Mobiles[is.na(Mobiles)]<-0

#######Sessiles#########
#convert characters to numeric in sessile sheet
Sessiles$Epiactis.prolifera<-as.numeric(Sessiles$Epiactis.prolifera)
Sessiles$Chaetomorpha.linum<-as.numeric(Sessiles$Chaetomorpha.linum)
Sessiles$Costaria.costata<-as.numeric(Sessiles$Costaria.costata)
#replace NA with 0
Sessiles[is.na(Sessiles)]<-0 


# Make all the community data a relative percent
PercentSessile<-100*Sessiles[7:ncol(Sessiles)]/Sessiles$Squares #change to rock--end spp

#normalize to the sum of the total cover (since it can be greater than 100%)
StandardizedSessile<- 100*PercentSessile/rowSums(PercentSessile)


Communitymetrics<-cbind(Sessiles$PoolID, Sessiles$Foundation_spp, Sessiles$Removal_Control, 
                        Sessiles$Before_After,StandardizedSessile,PercentSessile$Mytilus.californianus,PercentSessile$Phyllospadix.spp)
Communitymetrics <- Communitymetrics %>%
  rename(PoolID = "Sessiles$PoolID", Foundation_spp = "Sessiles$Foundation_spp",Removal_Control ="Sessiles$Removal_Control",
         Before_After ="Sessiles$Before_After", AdjMusselCover=Mytilus.californianus, AdjSurfgrassCover=Phyllospadix.spp, MusselCover = "PercentSessile$Mytilus.californianus",SurfgrassCover = "PercentSessile$Phyllospadix.spp")
#rename joined columns

#combine data for community metrics for Strcutural equation models
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

#physical parameters averaged over time periods
PP<-TidePooldes %>%
  dplyr::group_by(PoolID,Removal_Control) %>%
  dplyr::summarise(SAVav = mean(SAtoV), #ave between before and after since SA/V changed with fspp removal
                   THav = mean(TideHeight),SAav=mean(SurfaceArea),Vav=mean(Vol),Depthav=mean(MaxDepth)) #tide height didn't change 
PP$PoolID<-as.factor(PP$PoolID)
Funsppcover$PoolID<-as.factor(Funsppcover$PoolID)
Funsppandpp<-left_join(Funsppcover,PP)
#ggpairs(Funsppandpp[4:9])
#mussel and surfgrass are not correlated with volume

#combine with functional groups
SessilesMytilusFun<-Communitymetrics%>%
  dplyr::filter(Before_After != 'Immediate') %>%
  dplyr::group_by(PoolID, Foundation_spp, Before_After, Removal_Control) %>%
  tidyr::pivot_longer(
    cols = AdjSurfgrassCover:Stylantheca.spp,
    names_to = "Species", #creates column with species in longformat
    values_to = "Cover", #adds column for % cover
    values_drop_na = TRUE
  ) 
SessilesphylloFun<-Communitymetrics%>%
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

#find change in functional groups between after and before period in each tide pool
DeltaSessMussels<-SessilesMytilusStacked %>%
  filter(Foundation_spp == 'Mytilus') %>%
  dplyr::group_by(PoolID,Removal_Control,Foundation_spp,Before_After,Functional_Group) %>%
  summarise(SumCover = sum(Cover)) %>%
  dplyr::group_by(PoolID,Removal_Control,Foundation_spp, Functional_Group) %>%
  dplyr::summarise(Deltacover = SumCover[Before_After=="After"]-SumCover[Before_After =="Before"])

DeltaSessSurf<-SessilesPhylloStacked%>%
  filter(Foundation_spp == 'Phyllospadix') %>%
  dplyr::group_by(PoolID,Removal_Control,Foundation_spp,Functional_Group, Before_After) %>%
  summarise(SumCover = sum(Cover)) %>%
  dplyr::group_by(PoolID,Removal_Control,Foundation_spp, Functional_Group) %>%
  dplyr::summarise(Deltacover = SumCover[Before_After=="After"]-SumCover[Before_After =="Before"])

DeltaSessSurf$PoolID<-as.factor(DeltaSessSurf$PoolID)
DeltaSessMussels$PoolID<-as.factor(DeltaSessMussels$PoolID)
DeltaSessSurf<-left_join(DeltaSessSurf,Funsppandpp)
DeltaSessMussels<-left_join(DeltaSessMussels,Funsppandpp)

#####Mobile data filtering####
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
count<-Mobdata%>%
  select(!Species) %>%  # remove the species names
  group_by(Foundation_spp, Removal_Control, Fun_group) %>% # summarize by group
  summarize(avg = mean(Std.Count, na.rm = TRUE))

Mobdata<- Mobdata%>%
  dplyr::group_by(PoolID,Removal_Control,Foundation_spp,Before_After,Fun_group,Species) %>%
  summarise(SumDensity = sum(Std.Count)) %>%
  dplyr::group_by(PoolID,Removal_Control,Foundation_spp, Fun_group,Species) %>%
  summarise(DeltaCount = SumDensity[Before_After=="After"]-SumDensity[Before_After =="Before"]) 

Mobdata<-left_join(Mobdata,Funsppandpp)
#remove suspension feeder since only 1 present in the functional group
Musselmob<-Mobdata%>%
  dplyr::filter(Foundation_spp=='Mytilus'& Fun_group != 'SuspensionFeeder')

Surfgrassmob<-Mobdata%>%
  dplyr::filter(Foundation_spp=='Phyllospadix'& Fun_group != 'SuspensionFeeder')

#### Functional groups PCoA-------------------------

# surfgrass sessile #############################

fungroups_wideSS <-pivot_wider(DeltaSessSurf,names_from = Functional_Group, values_from = Deltacover)

Enviro2SS<-fungroups_wideSS %>%
  ungroup() %>% #ungroup from factor categories
  select(PoolID:Depthav)

fungroup2SS<-fungroups_wideSS %>%
  ungroup() %>% #ungroup from factor categories
  select(Anemone:SuspensionFeeder)

# try PCOA

distSS <- vegdist(fungroup2SS,  method = "euclidean")
# Some distance measures may result in negative eigenvalues. In that case, add a correction:
PCOA_SS <- pcoa(distSS, correction = "cailliez")

biplot<-biplot.pcoa(PCOA_SS, fungroup2SS)
PCOAaxesSS <- PCOA_SS$vectors[,c(1,2)]
allPCOSS<-bind_cols(Enviro2SS,data.frame(PCOAaxesSS)) 


compute_arrows = function(given_pcoa, trait_df) {
  
  # Keep only quantitative or ordinal variables
  # /!\ Change this line for different dataset
  #     or select only quantitative/ordinal var. /!\
  
  n <- nrow(trait_df)
  points.stand <- scale(given_pcoa$vectors)
  
  # Compute covariance of variables with all axes
  S <- cov(trait_df, points.stand)
  
  # Select only positive eigenvalues
  pos_eigen = given_pcoa$values$Eigenvalues[seq(ncol(S))]
  
  # Standardize value of covariance (see Legendre & Legendre 1998)
  U <- S %*% diag((pos_eigen/(n - 1))^(-0.5))
  colnames(U) <- colnames(given_pcoa$vectors)
  
  # Add values of covariances inside object
  given_pcoa$U <- U
  
  return(given_pcoa)
}

trait_pcoa_arrows <- compute_arrows(PCOA_SS, fungroup2SS)

# Basic plot with individuals

# Now let's add the arrows
# Each arrow begins at the origin of the plot (x = 0, y = 0) and ends at the
# values of covariances of each variable

arrows_df <- as.data.frame(trait_pcoa_arrows$U*50)
arrows_df$variable <- rownames(arrows_df)

# Make pretty names
arrows_df$variable_pretty<-c("Anemone"," Articulated \nCorallines", 
                             "Corticated \nFoliose","Corticated \nMacroalgae",
                             "Crustose", "Filamentous", "Foliose", "Leathery \nMacroalgae",
                             "Microalgae", "Suspension \nFeeders") 

plotPCO_SS<-ggplot(allPCOSS, aes(Axis.1, Axis.2))+
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(aes(fill=Phyllodelta, size =Phyllodelta,stroke=1), shape=21,colour = "black")+
  geom_segment(data = arrows_df,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               size = 1.2,
               # Add arrow head
               arrow = arrow(length = unit(3, "mm"))) +
  ggrepel::geom_text_repel(data = arrows_df, aes(label = variable_pretty), size =9, max.overlaps =25)+
  theme_bw()+ # Add axes
  # Keep coord equal for each axis
  # coord_equal() +
  scale_fill_distiller(expression("Surfgrass percent loss"), palette = 'Greens',
                      direction = -1,
                      breaks=seq(-25, 100, by=25),guide="legend")+
  scale_size_continuous(expression("Surfgrass percent loss"),range = c(1,10),
                        breaks=seq(-25, 100, by=25))+
  guides(fill= guide_legend(nrow =1)) +
  theme(axis.text = element_text(color = "black", size = 35), 
        axis.title.x = element_text(color="black", size=40), 
        axis.title.y = element_text(color="black", size=40), 
        legend.title = element_text(color="black", size=35), 
        legend.text = element_text(color = "black", size = 35), 
        legend.position= "top",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  # Labs
  labs(#subtitle = "Arrows scale arbitrarily",
    # Add Explained Variance per axis
    x = paste0("Axis 1 (", round(trait_pcoa_arrows$values$Relative_eig[1] * 100, 2), "%)"),
    y = paste0("Axis 2 (", round(trait_pcoa_arrows$values$Relative_eig[2] * 100, 2), "%)")) 
# Theme change
#labs(subtitle = "Arrows scale arbitrarily") # make sure to add in the legend that the arrows are arbitrarily scaled
plotPCO_SS

# permanova
perm_phyllo_SS<-adonis(distSS~Enviro2SS$Phyllodelta+Enviro2SS$Vav+Enviro2SS$THav)
perm_phyllo_SS

#### Surfrass mobile ####################

fungroups_SM<-Surfgrassmob %>%
  select(!Species) %>%  # remove the species names
  group_by(PoolID,Foundation_spp, Removal_Control, Fun_group) %>% # summarize by group
  summarize(DeltaCount = sum(DeltaCount, na.rm = TRUE),
            Phyllodelta = mean(Phyllodelta),
            THav = mean(THav),
            Vav = mean(Vav)) %>%
  ungroup() %>%
  filter(PoolID != 18) %>% #removed tide pool 18 due to outlier with micrograzers (tide height influence)
  mutate(DeltaCount = sign(DeltaCount)*sqrt(abs(DeltaCount))) 
#  select(-DeltaCount)


fungroups_wideSM <-pivot_wider(fungroups_SM,names_from = Fun_group, values_from = DeltaCount)

Enviro2SM<-fungroups_wideSM %>%
  select(PoolID:Vav)

fungroup2SM<-fungroups_wideSM %>%
  select(Carnivore:Omnivore)

# try PCOA

distSM <- vegdist(fungroup2SM,  method = "euclidean")
# Some distance measures may result in negative eigenvalues. In that case, add a correction:
PCOASM <- pcoa(distSM, correction = "cailliez")

biplot<-biplot.pcoa(PCOASM, fungroup2SM)
PCOAaxesSM <- PCOASM$vectors[,c(1,2)]
allPCOSM<-bind_cols(Enviro2SM,data.frame(PCOAaxesSM)) 

trait_pcoa_arrows <- compute_arrows(PCOASM, fungroup2SM)

# Now let's add the arrows
# Each arrow begins at the origin of the plot (x = 0, y = 0) and ends at the
# values of covariances of each variable

arrows_df <- as.data.frame(trait_pcoa_arrows$U*5)
arrows_df$variable <- rownames(arrows_df)

plotPCO_SM<-ggplot(allPCOSM, aes(Axis.1, Axis.2))+
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(aes(fill=Phyllodelta, size =Phyllodelta,stroke=1), shape=21,colour = "black")+
  geom_segment(data = arrows_df,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               size = 1.2,
               # Add arrow head
               arrow = arrow(length = unit(3, "mm"))) +
  ggrepel::geom_text_repel(data = arrows_df, aes(label = variable), size =9,max.overlaps = 20)+
  theme_bw()+ # Add axes
  guides(fill= guide_legend(nrow =1)) +
  theme(axis.text = element_text(color = "black", size = 35), 
        axis.title.x = element_text(color="black", size=40), 
        axis.title.y = element_text(color="black", size=40), 
        legend.title = element_text(color="black", size=35), 
        legend.text = element_text(color = "black", size = 35), 
        legend.position= "none",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  # Keep coord equal for each axis
  # coord_equal(ratio = 1) +
  scale_fill_distiller(expression("Surfgrass percent loss"), palette = 'Greens',
                       direction = -1,
                       breaks=seq(-25, 100, by=25),guide="legend")+
  scale_size_continuous(expression("Surfgrass percent loss"),range = c(1,10),
                        breaks=seq(-25, 100, by=25))+
  # Labs
  labs(#subtitle = "Arrows scale arbitrarily",
    # Add Explained Variance per axis
    x = paste0("Axis 1 (", round(trait_pcoa_arrows$values$Relative_eig[1] * 100, 2), "%)"),
    y = paste0("Axis 2 (", round(trait_pcoa_arrows$values$Relative_eig[2] * 100, 2), "%)")) 
# Theme change
plotPCO_SM

# permanova
perm_phylloSM<-adonis(distSM~Enviro2SM$Phyllodelta+Enviro2SM$Vav+Enviro2SM$THav)
perm_phylloSM

#####Surfgrass Pcoa with tide pool 18####
fungroups_SM18<-Surfgrassmob %>%
  select(!Species) %>%  # remove the species names
  group_by(PoolID,Foundation_spp, Removal_Control, Fun_group) %>% # summarize by group
  summarize(DeltaCount = sum(DeltaCount, na.rm = TRUE),
            Phyllodelta = mean(Phyllodelta),
            THav = mean(THav),
            Vav = mean(Vav)) %>%
  ungroup() %>%
  mutate(DeltaCount = sign(DeltaCount)*sqrt(abs(DeltaCount))) 
#  select(-DeltaCount)


fungroups_wideSM18 <-pivot_wider(fungroups_SM18,names_from = Fun_group, values_from = DeltaCount)

Enviro2SM18<-fungroups_wideSM18 %>%
  select(PoolID:Vav)

fungroup2SM18<-fungroups_wideSM18 %>%
  select(Carnivore:Omnivore)

# try PCOA

distSM18 <- vegdist(fungroup2SM18,  method = "euclidean")
# Some distance measures may result in negative eigenvalues. In that case, add a correction:
PCOASM18 <- pcoa(distSM18, correction = "cailliez")

biplot<-biplot.pcoa(PCOASM18, fungroup2SM18)
PCOAaxesSM18 <- PCOASM18$vectors[,c(1,2)]
allPCO18<-bind_cols(Enviro2SM18,data.frame(PCOAaxesSM18)) 

trait_pcoa_arrows <- compute_arrows(PCOASM18, fungroup2SM18)

# Now let's add the arrows
# Each arrow begins at the origin of the plot (x = 0, y = 0) and ends at the
# values of covariances of each variable

arrows_df <- as.data.frame(trait_pcoa_arrows$U*5)
arrows_df$variable <- rownames(arrows_df)

plotPCO_SM18<-ggplot(allPCO18, aes(Axis.1, Axis.2))+
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(aes(fill=Phyllodelta, size =Phyllodelta,stroke=1), shape=21,colour = "black")+
  geom_segment(data = arrows_df,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               size = 1.2,
               # Add arrow head
               arrow = arrow(length = unit(3, "mm"))) +
  ggrepel::geom_text_repel(data = arrows_df, aes(label = variable),size =9, max.overlaps = 20)+
  theme_bw()+ # Add axes
  guides(fill= guide_legend(nrow =1)) +
  theme(axis.text = element_text(color = "black", size = 35), 
        axis.title.x = element_text(color="black", size=40), 
        axis.title.y = element_text(color="black", size=40), 
        legend.title = element_text(color="black", size=35), 
        legend.text = element_text(color = "black", size = 35), 
        legend.position= "top",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  # Keep coord equal for each axis
  # coord_equal(ratio = 1) +
  scale_fill_distiller(expression("Surfgrass percent loss"), palette = 'Greens',
                       direction = -1,
                       breaks=seq(-25, 100, by=25),guide="legend")+
  scale_size_continuous(expression("Surfgrass percent loss"),range = c(1,10),
                        breaks=seq(-25, 100, by=25))+
  # Labs
  labs(#subtitle = "Arrows scale arbitrarily",
    # Add Explained Variance per axis
    x = paste0("Axis 1 (", round(trait_pcoa_arrows$values$Relative_eig[1] * 100, 2), "%)"),
    y = paste0("Axis 2 (", round(trait_pcoa_arrows$values$Relative_eig[2] * 100, 2), "%)")) 
# Theme change
plotPCO_SM18
ggsave(filename = "Output/PCOA_pool18.pdf", useDingbats =FALSE,dpi=600,device = "pdf", width = 15, height = 12)

# permanova with tide pool 18
perm_phylloSM18<-adonis(distSM18~Enviro2SM18$Phyllodelta+Enviro2SM18$Vav+Enviro2SM18$THav)
perm_phylloSM18
###--------------- Mussels sessile #############
fungroups_MS<-DeltaSessMussels %>%
  mutate(Deltacover = sign(Deltacover)*sqrt(abs(Deltacover))) %>%
  select(-Phyllodelta)

fungroups_wideMS <-pivot_wider(fungroups_MS,names_from = Functional_Group, values_from = Deltacover  )

Enviro2MS<-fungroups_wideMS %>%
  ungroup()%>%
  select(PoolID:Depthav)

fungroup2MS<-fungroups_wideMS %>%
  ungroup()%>%
  select(Anemone:SuspensionFeeder)

# try PCOA

distMS <- vegdist(fungroup2MS,  method = "euclidean")
# Some distance measures may result in negative eigenvalues. In that case, add a correction:
PCOAMS <- pcoa(distMS, correction = "cailliez")

biplot<-biplot.pcoa(PCOAMS, fungroup2MS)
PCOAaxesMS <- PCOAMS$vectors[,c(1,2)]
allPCOMS<-bind_cols(Enviro2MS,data.frame(PCOAaxesMS)) 



trait_pcoa_arrows <- compute_arrows(PCOAMS, fungroup2MS)

# Basic plot with individuals

# Now let's add the arrows
# Each arrow begins at the origin of the plot (x = 0, y = 0) and ends at the
# values of covariances of each variable

arrows_df <- as.data.frame(trait_pcoa_arrows$U*5)
arrows_df$variable <- rownames(arrows_df)

# Make pretty names
arrows_df$variable_pretty<-c("Anemone"," Articulated \nCorallines", 
                             "Corticated \nFoliose","Corticated \nMacroalgae",
                             "Crustose", "Filamentous", "Foliose", "Leathery \nMacroalgae",
                             "Microalgae", "Suspension \nFeeders") 

plotPCO_MS<-ggplot(allPCOMS, aes(Axis.1, Axis.2))+
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(aes(fill=Mytilusdelta, size =Mytilusdelta,stroke=1), shape=21,colour = "black")+
  geom_segment(data = arrows_df,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               size = 1.2,
               # Add arrow head
               arrow = arrow(length = unit(3, "mm"))) +
  ggrepel::geom_text_repel(data = arrows_df, aes(label = variable_pretty),size =9, max.overlaps = 20)+
  theme_bw()+ # Add axes
  # Keep coord equal for each axis
  # coord_equal(1.5) +
  theme(axis.text = element_text(color = "black", size = 35), 
        axis.title.x = element_text(color="black", size=40), 
        axis.title.y = element_text(color="black", size=40), 
        legend.title = element_text(color="black", size=35), 
        legend.text = element_text(color = "black", size = 35), 
        legend.position= "top",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_fill_distiller(expression("Mussel percent loss"), palette = 'Blues',
                                  direction= -1,limits=c(-25,100),
                      breaks=seq(-25, 100, by=25),guide="legend")+
  scale_size_continuous(expression("Mussel percent loss"),range = c(1,10),limits=c(-25,100),
                        breaks=seq(-25, 100, by=25))+
  guides(fill= guide_legend(nrow =1)) +
  # Labs
  labs(#subtitle = "Arrows scale arbitrarily",
    # Add Explained Variance per axis
    x = paste0("Axis 1 (", round(trait_pcoa_arrows$values$Relative_eig[1] * 100, 2), "%)"),
    y = paste0("Axis 2 (", round(trait_pcoa_arrows$values$Relative_eig[2] * 100, 2), "%)")) 
# Theme change
plotPCO_MS

# permanova
perm_mussel_S<-adonis(distMS~Enviro2MS$Mytilusdelta+Enviro2MS$Vav+Enviro2MS$THav)
#labs(subtitle = "Arrows scale arbitrarily") # make sure to add in the legend that the arrows are arbitrarily scaled
perm_mussel_S

###---- Mussel Mobile ---###########
fungroups_MM<-Musselmob %>%
  select(!Species) %>%  # remove the species names
  group_by(PoolID,Foundation_spp, Removal_Control, Fun_group) %>% # summarize by group
  summarize(DeltaCount = sum(DeltaCount, na.rm = TRUE),
            Mytilusdelta = mean(Mytilusdelta),
            THav = mean(THav),
            Vav = mean(Vav)) %>%
  ungroup() %>%
  mutate(DeltaCount = sign(DeltaCount)*sqrt(abs(DeltaCount))) 


fungroups_wideMM <-pivot_wider(fungroups_MM,names_from = Fun_group, values_from = DeltaCount)

Enviro2MM<-fungroups_wideMM %>%
  select(PoolID:Vav)

fungroup2MM<-fungroups_wideMM %>%
  select(Carnivore:Omnivore)

# try PCOA

distMM <- vegdist(fungroup2MM,  method = "euclidean")
# Some distance measures may result in negative eigenvalues. In that case, add a correction:
PCOAMM <- pcoa(distMM, correction = "cailliez")

biplot<-biplot.pcoa(PCOAMM, fungroup2MM)
PCOAaxesMM <- PCOAMM$vectors[,c(1,2)]
allPCOMM<-bind_cols(Enviro2MM,data.frame(PCOAaxesMM)) 

trait_pcoa_arrows <- compute_arrows(PCOAMM, fungroup2MM)

# Now let's add the arrows
# Each arrow begins at the origin of the plot (x = 0, y = 0) and ends at the
# values of covariances of each variable

arrows_df <- as.data.frame(trait_pcoa_arrows$U*10)
arrows_df$variable <- rownames(arrows_df)

plotPCO_MM<-ggplot(allPCOMM, aes(Axis.1, Axis.2))+
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(aes(fill=Mytilusdelta, size =Mytilusdelta,stroke=1), shape=21,colour = "black")+
  geom_segment(data = arrows_df,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               size = 1.2,
               # Add arrow head
               arrow = arrow(length = unit(3, "mm"))) +
  ggrepel::geom_text_repel(data = arrows_df, aes(label = variable),size =9, max.overlaps = 20)+
  theme_bw()+ # Add axes
  # Keep coord equal for each axis
  coord_equal(ratio = 2) +
  theme(axis.text = element_text(color = "black", size = 35), 
        axis.title.x = element_text(color="black", size=40), 
        axis.title.y = element_text(color="black", size=40), 
        legend.title = element_text(color="black", size=35), 
        legend.text = element_text(color = "black", size = 35), 
        legend.position= "none",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_fill_distiller(expression("Mussel percent loss"), palette = 'Blues',
                       direction= -1,
                       breaks=seq(-25, 100, by=25),guide="legend")+
  scale_size_continuous(expression("Mussel percent loss"),range = c(1,10),
                        breaks=seq(-25, 100, by=25)) +
  # Labs
  labs(#subtitle = "Arrows scale arbitrarily",
    # Add Explained Variance per axis
    x = paste0("Axis 1 (", round(trait_pcoa_arrows$values$Relative_eig[1] * 100, 2), "%)"),
    y = paste0("Axis 2 (", round(trait_pcoa_arrows$values$Relative_eig[2] * 100, 2), "%)")) 
# Theme change
plotPCO_MM
# permanova
perm_musselMM<-adonis(distMM~Enviro2MM$Mytilusdelta+Enviro2MM$Vav+Enviro2MM$THav)
perm_musselMM

## bring together

PCOAall<-
  (plotPCO_SS+ plotPCO_SM)/
  (plotPCO_MS+plotPCO_MM) +
  plot_annotation(tag_levels = 'A') &         #label each individual plot with letters A-G
  theme(plot.tag = element_text(size =40))
  #plot_layout(guides = "collect")#edit the lettered text
PCOAall
ggsave(filename = "Output/PCOA_all.pdf", useDingbats =FALSE,dpi=600,device = "pdf", width = 25, height = 20)


#####PCoA surfgrass####
SurfgrassSessFunGroup <- read_csv("Data/SurfgrassSessFunGroup.csv")
DeltaSessSurf<-as.data.frame(DeltaSessSurf)
Surffungroups_wide <-pivot_wider(DeltaSessSurf,
                                 names_from = Functional_Group, 
                                 values_from = Deltacover) #change to wide format for analysis

Enviro2<-Surffungroups_wide %>%
  select(PoolID:Depthav)

fungroup2<-Surffungroups_wide %>%
  ungroup() %>% #ungroup from factor categories
  dplyr::select(Anemone:SuspensionFeeder)
# try PCOA
dist <- vegdist(fungroup2,  method = "euclidean")
# Some distance measures may result in negative eigenvalues. In that case, add a correction:
PCOA <- pcoa(dist, correction = "cailliez")
biplot<-biplot.pcoa(PCOA, fungroup2)
PCOAaxes <- PCOA$vectors[,c(1,2)]
allPCO<-bind_cols(Enviro2,data.frame(PCOAaxes))
ggplot(allPCO, aes(Axis.1, Axis.2))+
  geom_point(aes(size = Phyllodelta))
compute_arrows = function(given_pcoa, trait_df) {
  # Keep only quantitative or ordinal variables
  # /!\ Change this line for different dataset
  #     or select only quantitative/ordinal var. /!\
  n <- nrow(trait_df)
  points.stand <- scale(given_pcoa$vectors)
  # Compute covariance of variables with all axes
  S <- cov(trait_df, points.stand)
  # Select only positive eigenvalues
  pos_eigen = given_pcoa$values$Eigenvalues[seq(ncol(S))]
  # Standardize value of covariance (see Legendre & Legendre 1998)
  U <- S %*% diag((pos_eigen/(n - 1))^(-0.5))
  colnames(U) <- colnames(given_pcoa$vectors)
  # Add values of covariances inside object
  given_pcoa$U <- U
  return(given_pcoa)
}
trait_pcoa_arrows <- compute_arrows(PCOA, fungroup2)
# Basic plot with individuals
plot_pcoa <-  ggplot(allPCO, aes(Axis.1, Axis.2))+
  geom_point(aes(size = Phyllodelta))
plot_pcoa
# Now let's add the arrows
# Each arrow begins at the origin of the plot (x = 0, y = 0) and ends at the
# values of covariances of each variable
ggplot(allPCO, aes(Axis.1, Axis.2))+
  geom_point(aes(size = Phyllodelta))
arrows_df = as.data.frame(trait_pcoa_arrows$U*50)
arrows_df$variable = rownames(arrows_df)
plot_pcoa +
  geom_segment(data = arrows_df,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               # Add arrow head
               arrow = arrow(length = unit(3, "mm"))) +
  ggrepel::geom_label_repel(data = arrows_df, aes(label = variable))
#labs(subtitle = "Arrows scale arbitrarily") # make sure to add in the legend that the arrows are arbitrarily scaled
plot_pcoa +
  geom_segment(data = as.data.frame(trait_pcoa_arrows$U*50),
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               # Add arrow head
               arrow = arrow(length = unit(3, "mm"))) +
  # Add Axes Labels
  ggrepel::geom_text_repel(data = arrows_df, aes(label = variable), min.segment.length = .01,max.overlaps = 20 )  +
  # Add axes
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  # Keep coord equal for each axis
  coord_equal() +
  # Labs
  labs(subtitle = "Arrows scale arbitrarily",
       # Add Explained Variance per axis
       x = paste0("Axis 1 (", round(trait_pcoa_arrows$values$Relative_eig[1] * 100, 2), "%)"),
       y = paste0("Axis 2 (", round(trait_pcoa_arrows$values$Relative_eig[2] * 100, 2), "%)")) +
  # Theme change
  theme_bw()
# permanova
set.seed(456)
perm_phyllo<-adonis(dist~Enviro2$Phyllodelta+Enviro2$Vav+Enviro2$THav)
perm_phyllo

####PCoA Mussels####
Musselfungroups_wide <-pivot_wider(DeltaSessMussels,
                                 names_from = Functional_Group, 
                                 values_from = Deltacover) #change to wide format for analysis

Enviro2<-Musselfungroups_wide %>%
  select(PoolID:Depthav)

fungroup2<-Musselfungroups_wide %>%
  ungroup() %>% #ungroup from factor categories
  dplyr::select(Anemone:SuspensionFeeder)
# try PCOA
dist <- vegdist(fungroup2,  method = "euclidean")
# Some distance measures may result in negative eigenvalues. In that case, add a correction:
PCOA <- pcoa(dist, correction = "cailliez")
biplot<-biplot.pcoa(PCOA, fungroup2)
PCOAaxes <- PCOA$vectors[,c(1,2)]
allPCO<-bind_cols(Enviro2,data.frame(PCOAaxes))
ggplot(allPCO, aes(Axis.1, Axis.2))+
  geom_point(aes(size = Mytilusdelta))
compute_arrows = function(given_pcoa, trait_df) {
  # Keep only quantitative or ordinal variables
  # /!\ Change this line for different dataset
  #     or select only quantitative/ordinal var. /!\
  n <- nrow(trait_df)
  points.stand <- scale(given_pcoa$vectors)
  # Compute covariance of variables with all axes
  S <- cov(trait_df, points.stand)
  # Select only positive eigenvalues
  pos_eigen = given_pcoa$values$Eigenvalues[seq(ncol(S))]
  # Standardize value of covariance (see Legendre & Legendre 1998)
  U <- S %*% diag((pos_eigen/(n - 1))^(-0.5))
  colnames(U) <- colnames(given_pcoa$vectors)
  # Add values of covariances inside object
  given_pcoa$U <- U
  return(given_pcoa)
}
trait_pcoa_arrows <- compute_arrows(PCOA, fungroup2)
# Basic plot with individuals
plot_pcoa <-  ggplot(allPCO, aes(Axis.1, Axis.2))+
  geom_point(aes(size = Mytilusdelta))
plot_pcoa
# Now let's add the arrows
# Each arrow begins at the origin of the plot (x = 0, y = 0) and ends at the
# values of covariances of each variable
ggplot(allPCO, aes(Axis.1, Axis.2))+
  geom_point(aes(size = Mytilusdelta))
arrows_df = as.data.frame(trait_pcoa_arrows$U*50)
arrows_df$variable = rownames(arrows_df)
plot_pcoa +
  geom_segment(data = arrows_df,
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               # Add arrow head
               arrow = arrow(length = unit(3, "mm"))) +
  ggrepel::geom_label_repel(data = arrows_df, aes(label = variable))
#labs(subtitle = "Arrows scale arbitrarily") # make sure to add in the legend that the arrows are arbitrarily scaled
plot_pcoa +
  geom_segment(data = as.data.frame(trait_pcoa_arrows$U*50),
               x = 0, y = 0, alpha = 0.7,
               mapping = aes(xend = Axis.1, yend = Axis.2),
               # Add arrow head
               arrow = arrow(length = unit(3, "mm"))) +
  # Add Axes Labels
  ggrepel::geom_text_repel(data = arrows_df, aes(label = variable), min.segment.length = .01,max.overlaps = 20 )  +
  # Add axes
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  # Keep coord equal for each axis
  coord_equal() +
  # Labs
  labs(subtitle = "Arrows scale arbitrarily",
       # Add Explained Variance per axis
       x = paste0("Axis 1 (", round(trait_pcoa_arrows$values$Relative_eig[1] * 100, 2), "%)"),
       y = paste0("Axis 2 (", round(trait_pcoa_arrows$values$Relative_eig[2] * 100, 2), "%)")) +
  # Theme change
  theme_bw()
# permanova
set.seed(456)
perm_mytilus<-adonis(dist~Enviro2$Mytilusdelta+Enviro2$Vav+Enviro2$THav)
perm_mytilus



Musselmobfungroups_wide <-pivot_wider(Musselmob,
                                   names_from = Fun_group, 
                                   values_from = DeltaCount) #change to wide format for analysis


Musselmobfungroups_wide[is.na(Musselmobfungroups_wide)] <- 0


# permanova
#set.seed(456)
#perm_mytilus<-adonis(dist~Enviro2$Mytilusdelta+Enviro2$Vav+Enviro2$THav)
#perm_mytilus

#####nMDS Sessile Surfgrass######
PhyllocommunitynMDS<-Communitymetrics%>%
  dplyr::filter(Before_After != "Immediate" & Foundation_spp == "Phyllospadix") 

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
  filter(Species =="Algae.film"|Species == "Diatoms" |
           Species=="Crustose.coralline"|Species=="Noncoralline.crust"|Species=="Bosiella.spp"|Species=="Calliathron.tuberculosum"|
           Species=="Ulva.spp"|Species=="Smithura.naiadum") 
subset$Species<-c("Diatoms","CCA","Algae film","Noncoralline crust","Bossiella spp",
                  "Ulva spp","Smithora naiadum") #rename to easier labels

ordSesSurf<-ggplot(phyllospp)+
                   geom_point(aes(x=NMDS1,y=NMDS2,color=Functional_Group),size=8,shape =15) + 
  geom_label_repel(data=subset,aes(x=NMDS1,y=NMDS2,label=Species),
                   direction=c("both"),nudge_y=0.2,color="black",size = 12) +  # add the species labels
  #geom_text(data=psspecies.scores,aes(x=NMDS1,y=NMDS2,label=Species),color="#045a8d",size =8) +  # add the species labels
  scale_color_manual(values=Colors)+
  labs(x='',color="Functional group",y='Surfgrass species scores nMDS2',title='Sessile community')+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=6)),size=FALSE)+
  theme(axis.text.x = element_blank(),  # remove x-axis text
                    axis.text.y =element_blank(),
                    axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(color="black", size=50),
        axis.title.y = element_text(color="black", size=50),
        plot.title = element_text(color="black", size=50,hjust = 0.5,face="bold"),#centered and bolded title
        legend.title = element_blank(),
        legend.position = "none")
ordSesSurf       

phyllopp<-Funsppandpp%>%
  filter(Foundation_spp =="Phyllospadix")
phyllograph<-cbind(PhyllocommunitynMDS,phyllopp$Vav,phyllopp$THav)
phyllograph<-phyllograph%>%
  rename(Vol="phyllopp$Vav",TH="phyllopp$THav")
phyllographNyssa<-phyllograph[-c(6:8,10,72:73)] 
write.csv(phyllographNyssa, file="output/SurfgrassSessileComm.csv")
set.seed(267)
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
  geom_point(aes(fill=Phyllodelta, size =Phyllodelta,stroke=1), shape=21,colour = "black") + #black line surrounding pts
  scale_fill_distiller(palette = "Greys",guide = "legend",direction =1)+
  scale_size(range = c(1,15)) +
  geom_point(data=pscentroids, size=10, stroke = 2.75) +
  theme_classic() +
  scale_shape_manual(values = c(pgroupings))+
  geom_segment(aes(x = x0[1], y = y0[1], xend = (x1[1]), yend = (y1[1])),size = 1,#segment with arrow for surfgrass before/after control
               colour = "#3182bd", arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  geom_segment(aes(x = x0[2], y = y0[2], xend = (x1[2]), yend = (y1[2])),linetype = 2,size = 1, #segment with arrow for surfgrass before/after removal
               colour = "#252525",arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  labs(x ='', y = 'Surfgrass community nMDS2', shape='', title='Sessile community',fill='Foundation species percent loss',size='Foundation species percent loss', linetype ='Before or after') +
  theme(axis.text = element_blank(), 
        axis.title.x = element_text(color="black", size=50), 
        axis.title.y = element_text(color="black", size=50), 
        axis.ticks = element_blank(),
        legend.title = element_text(color="black", size=40), 
        legend.text = element_text(color = "black", size = 40), 
        legend.position = 'top',
        plot.title = element_text(color="black", size=50,hjust = 0.5,face="bold"),#centered and bolded title
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  guides(shape = "none", alpha = "none")+
  guides(fill= guide_legend(nrow =1))#makes legend only one row
Surfgrasssessilesplot


####Mussel Sessile nMDS#####
MytiluscommunitynMDS<-Communitymetrics%>%
  dplyr::filter(Before_After != "Immediate" & Foundation_spp == "Mytilus") 
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

labels$Species<-c("Diatoms","CCA","Algae film","Noncoralline crust","Odonthalia floccosa","Chthamalus spp","Anthopleura xanthogrammica")
ordSesMussel<-ggplot(musselspp) + 
  geom_point(aes(x=NMDS1,y=NMDS2,color=Functional_Group),size=8, shape =15) + 
  #geom_text(data=msspecies.scores,aes(x=NMDS1,y=NMDS2,label=species),color="#045a8d",size =8) +  # add all species labels
  geom_label_repel(data=labels,aes(x=NMDS1,y=NMDS2,label=Species),
                   direction=c("both"),nudge_y=0.3,color="black",size = 12) +  # add the species labels
  scale_color_manual(values=Colors,guide = "legend",labels =c("Anemone","Articulated corallines","Corticated foliose","Corticated macroalgae",
                                                              "Crustose","Filamentous","Foliose","Leathery macrophytes","Microalgae","Suspension feeders"))+
  labs(x='nMDS1',color="Functional group",y='Mussel species scores nMDS2')+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=6)),size=FALSE)+
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(color="black", size=50),
        axis.title.y = element_text(color="black", size=50),
        legend.position = "none")
ordSesMussel

mpp<-Funsppandpp%>%
  filter(Foundation_spp =="Mytilus")
mgraph<-cbind(MytiluscommunitynMDS,mpp$Vav,mpp$THav)

mgraph<-mgraph%>%
  rename(Vol="mpp$Vav",TH="mpp$THav")

set.seed(267)
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
  geom_point(aes(fill=Mytilusdelta, size =Mytilusdelta,stroke=1), shape=21, color='black') +
  scale_fill_distiller(palette = "Greys",guide = "legend",direction=1,limits=c(-25,100))+
  scale_size(range = c(1,15),limits=c(-25,100)) +
  geom_point(data=mscentroids, size=10, stroke = 2.75) +
  theme_classic() +
  scale_shape_manual(values = c(mgroupings))+
  geom_segment(aes(x = v0[1], y = z0[1], xend = (v1[1]), yend = (z1[1])),size = 1,#segment with arrow for surfgrass before/after control
               colour = "#3182bd", arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  geom_segment(aes(x = v0[2], y = z0[2], xend = (v1[2]), yend = (z1[2])),linetype = 2,size = 1, #segment with arrow for surfgrass before/after removal
               colour = "#252525",arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  labs(x ='nMDS1', y = 'Mussel community nMDS2', shape='', fill='Foundation species percent loss',size='Foundation species percent loss', linetype ='Before or after') +
  theme(axis.text = element_blank(), 
        axis.title.x = element_text(color="black", size=50), 
        axis.title.y = element_text(color="black", size=50), 
        legend.title = element_text(color="black", size=40), 
        axis.ticks = element_blank(),
        legend.text = element_text(color = "black", size = 40), 
        legend.position = 'none',
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  guides(shape = "none")+
  guides(fill= guide_legend(nrow =1))#makes legend only one row
Musselsessilesplot

######Mobile nMDS######
 
SurfgrassMobiles <-Mobiles %>%
  filter(Before_After !="Immediate" & Foundation_spp =="Phyllospadix")


MusselMobiles <-Mobiles %>%
  filter(Before_After !="Immediate" & Foundation_spp =="Mytilus") 
  
####Surfgrass mobile nMDS#####
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
                   direction=c("both"),nudge_y=0.3,color="black",size = 12) +  # add the species labels
  scale_color_manual(values=MobColors,guide = "legend",labels =c("Carnivores","Herbivores","Omnivores"))+
  labs(x= "",y='',title='Mobile community',color="Functional group")+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=6)),size=FALSE)+
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(color="black", size=50),
        axis.title.y = element_text(color="black", size=50),
        plot.title = element_text(color="black", size=50,hjust = 0.5,face="bold"),#centered and bolded title
        legend.title = element_blank(),
        legend.position = "none")
ordMobSurf

phyllomobgraph<-cbind(PhyllomobnMDS,phyllopp$Vav,phyllopp$THav)
write.csv(phyllomobgraph, file="output/SurfgrassMobileComm.csv")
phyllomobgraph<-phyllomobgraph%>%
  rename(Vol="phyllopp$Vav",TH="phyllopp$THav")
set.seed(267)
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
  geom_point(aes(fill =Phyllodelta, size =Phyllodelta, stroke=1), shape=21,color='black') +
  scale_fill_distiller(palette = "Greys",guide = "legend",direction=1)+
  scale_size(range = c(1,15)) +
  geom_point(data=pmcentroids, size=10, stroke = 2.75) +
  theme_classic() +
  scale_shape_manual(values = c(pgroupings))+
  geom_segment(aes(x = a0[1], y = b0[1], xend = (a1[1]), yend = (b1[1])),size = 1,#segment with arrow for surfgrass before/after control
               colour = "#3182bd", arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  geom_segment(aes(x = a0[2], y = b0[2], xend = (a1[2]), yend = (b1[2])),linetype = 2,size = 1, #segment with arrow for surfgrass before/after removal
               colour = "#252525",arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  labs(x ='', y = '', title='Mobile community',shape='', color='',size='', linetype ='Before or after') +
  theme(axis.text = element_blank(), 
        axis.title.x = element_text(color="black", size=50), 
        axis.title.y = element_text(color="black", size=50), 
        plot.title = element_text(color="black", size=50,hjust = 0.5,face="bold"),#centered and bolded title
        axis.ticks = element_blank(),  # remove axis ticks
        legend.position= "none",
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
                   direction=c("both"),nudge_y=0.1,color="black",size = 12) +  # add the species labels
  labs(x= "nMDS1",y='',color="Functional group")+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=6)),size=FALSE)+
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(color="black", size=50),
        axis.title.y = element_text(color="black", size=50),
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
  geom_point(aes(fill =Mytilusdelta, size =Mytilusdelta,stroke=1), shape=21,color='black') +
  scale_fill_distiller(palette = "Greys",guide = "legend",direction=1,limits=c(-25,100))+
  scale_size(range = c(1,15),limits=c(-25,100)) +
  geom_point(data=mmcentroids, size=10, stroke = 2.75) +
  theme_classic() +
  scale_shape_manual(values = c(mgroupings))+
  geom_segment(aes(x = c0[1], y = d0[1], xend = (c1[1]), yend = (d1[1])),size = 1,#segment with arrow for surfgrass before/after control
               colour = "#3182bd", arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  geom_segment(aes(x = c0[2], y = d0[2], xend = (c1[2]), yend = (d1[2])),linetype = 2,size = 1, #segment with arrow for surfgrass before/after removal
               colour = "#252525",arrow = arrow(length = unit(0.3, "cm"),type = "closed")) +
  labs(x ='nMDS1', y = '', shape='', color='',size='', linetype ='Before or after') +
  theme(axis.text = element_blank(), 
        axis.title.x = element_text(color="black", size=50), 
        axis.title.y = element_text(color="black", size=50), 
        axis.ticks = element_blank(),  # remove axis ticks
        legend.position= "none",
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  guides(shape = "none", alpha = "none")
Musselmobplot

#patchwork of ord plots and nMDS per foudnation spp
Mcommgraphs<-Musselsessilesplot+Musselmobplot+ordSesMussel+ordMobMussel+
  plot_annotation(tag_levels = 'a') &         #label each individual plot with letters A-G
  theme(plot.tag = element_text(size = 50))   #edit the lettered text
Mcommgraphs
#ggsave("Output/McommnMDS.pdf",useDingbats = FALSE, width=75, height=60,dpi=600, unit="cm")

Pcommgraphs<-Surfgrasssessilesplot+Surfgrassmobplot+ordSesSurf+ordMobSurf+
  plot_annotation(tag_levels = 'a') &         #label each individual plot with letters A-G
  theme(plot.tag = element_text(size = 50))   #edit the lettered text
#ggsave("Output/ScommnMDS.pdf",useDingbats = FALSE, width=75, height=60,dpi=600, unit="cm")
#combined surfgrass and mussel community plots
ComnMDS<-Surfgrasssessilesplot+Surfgrassmobplot+Musselsessilesplot+Musselmobplot+
  plot_annotation(tag_levels = 'a') &         #label each individual plot with letters A-G
  theme(plot.tag = element_text(size = 50))   #edit the lettered text
ComnMDS
#ggsave("Output/CommunitycompnMDS.pdf",useDingbats = FALSE,  width = 35, height = 30,dpi=600)

#supplement spp scores
Sppscores<-ordSesSurf+ordMobSurf+ordSesMussel+ordMobMussel+
  plot_annotation(tag_levels = 'a') &         #label each individual plot with letters A-G
  theme(plot.tag = element_text(size = 50))   #edit the lettered text
Sppscores
#ggsave("Output/sppscoresupplement.pdf",useDingbats = FALSE, width = 35, height = 30,dpi=600)

####Species Richness Anovas#####
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
###Sessile species richness anova and plots####
#Phyllospadix pool sessile richness
phyllosessrichmod<-lm(DeltaRich~Phyllodelta +Vav+THav, data=phyllosessrich)
#plot(phyllosessrichmod) #good
qqp(resid(phyllosessrichmod),"norm") #good

summary(phyllosessrichmod)

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

#Mytilus pools sessile richness
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

###Mobiles Species richness anovas and plots####
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

#Phyllospadix pools mobile richness
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


#Mytilus pools mobile richness
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
#ggsave(filename = "Output/richnessplots.pdf", useDingbats =FALSE,dpi=600,device = "pdf",width = 28, height = 30)


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
#ggsave(filename = "Output/distplots.pdf", useDingbats =FALSE,dpi=600,device = "pdf", width = 25, height = 20)


