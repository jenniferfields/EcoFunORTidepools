## DAG for Structural Equation model for mussel and surfgrass models
##Code for SEM graph figures then further modified in Adobe Illustrator
## By: Jenn Fields
## Last updated: 10.23.2020
###################

rm(list=ls()) #Clears the environment

#load libraries
library(DiagrammeR)
library(rsvg) #to export to a pdf
library(DiagrammeRsvg) #to export to a pdf
#Visit http://rich-iannone.github.io/DiagrammeR/graphviz_and_mermaid.html for more options

#diagrams of SEM output for mussel and surfgrass day/night models
#kept labels in script as placeholders, but deleted in illustrator for final figures
#SEM dag for mussel SEM
#With temp in NEP model
MytilusDay <- grViz("
digraph boxes_and_circles { 
#this is the type of diagram you are using

# add node statements 
node [shape  = rectangle,
      penwidth =2.5,
      fontname = Helvetica,
      fontsize = 25]

MytilusLoss[label = 'CA Mussel loss \n(M. californianus)'];NEC; NEP; pH; 
MMC[label = 'Micro/macroalgae cover']; 
V[label = 'Volume']; 
MaxTemp[label = 'Temperature']; 
NtoPRatio[label = 'N:P']; 
TH[label = 'Tide height']

#End statements
# color red is significant negative; where as blue is significant positive; 
#grey is nonsignificant (p>0.05)
# size of path represents the standardized coefficients * 10 
# labels are significant standardized coefficients


MytilusLoss-> MMC[color = SteelBlue3, penwidth =7.9, label = '0.79',fontsize = 15,
fontname = Helvetica, minlen = 1.5]
MytilusLoss -> MaxTemp[color = grey, penwidth =3.4,fontsize = 15,
fontname = Helvetica, minlen = 1.5]
MytilusLoss-> NtoPRatio[color =  grey, penwidth = 3.6, minlen =1.5]
MytilusLoss->pH[color = SteelBlue3, penwidth =4.3, label = '0.43',fontsize = 15, fontname = Helvetica, minlen = 1.5]


V-> MMC[color =grey, penwidth = 0.3, minlen = 1.5]
TH-> MMC[color =OrangeRed2, penwidth = 3.9, label = '-0.39', fontsize = 15, fontname = Helvetica, minlen = 1.5]
V->MaxTemp[color = grey, penwidth = 0.41, minlen =1.5] 
TH->MaxTemp[color = grey, penwidth = 0.99, minlen =1.5]
V-> NtoPRatio[color = grey, penwidth = 4.0, minlen =1.5]
TH->NtoPRatio[color = grey, penwidth=4.0, minlen =1.5]
TH->NEP[color=grey,penwidth=2.4,minlen=1.5]
V-> pH[color = grey, penwidth =1.6, minlen = 1.5]
TH->pH[color = grey, penwidth =2.3, minlen = 1.5]
TH->NEC[color =grey, penwidth=2.7,minlen=1.5]

NtoPRatio-> NEP[color = OrangeRed2, penwidth = 2.0, label ='-0.20', fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]
MMC-> NEP[color = SteelBlue3, penwidth =8.9, label = '0.89', fontsize = 15, fontname = Helvetica, minlen =1.5]
NEP -> pH[color = SteelBlue3, penwidth =4.0,label='0.40' , fontsize = 15, fontname = Helvetica, minlen =1.5]
pH -> NEC[color = grey, penwidth =4.1, minlen =1.5] 
MaxTemp -> NEC[color = OrangeRed2, penwidth = 6.1, label ='-0.61', fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]

# add a graph statement
graph [nodesep = 0.1]
}")
MytilusDay

#export image to a PDF
MytilusDay %>%
  export_svg %>% 
  charToRaw %>%
  rsvg_pdf("Output/MytilusDayDAG.pdf", width =375, height =325)

#Dag sem for mussel night
MytilusNight <- grViz("
digraph boxes_and_circles { 
#this is the type of diagram you are using

# add node statements 
node [shape  = rectangle,
      penwidth =2.5,
      fontname = Helvetica,
      fontsize = 25]

MytilusLoss[label = 'CA mussel loss \n(M. californianus)'];NEC; NEP; pH; 
MMC[label = 'Micro/macroalgae cover']; 
V[label = 'Volume']; 
MaxTemp[label = 'Temperature']; 
NtoPRatio[label = 'N:P']; 
TH[label = 'Tide height']

#End statements
# color red is significant negative; where as blue is significant positive; 
#grey is nonsignificant (p>0.05)
# size of path represents the standardized coefficients * 10
# labels are significant standardized coefficients


MytilusLoss-> MMC[color = SteelBlue3, penwidth =7.9, label = '0.79',fontsize = 15,
fontname = Helvetica, minlen = 1.5]
MytilusLoss -> MaxTemp[color = grey penwidth =6.2, fontsize = 15,
fontname = Helvetica, minlen = 1.5]
MytilusLoss-> NtoPRatio[color =  grey, penwidth = 2.5, minlen =1.5]
MytilusLoss->pH[color = SteelBlue3, penwidth =7.5, label = '0.75',fontsize = 15, fontname = Helvetica, minlen = 1.5]

V-> MMC[color =grey, penwidth = 0.35, minlen = 1.5]
TH-> MMC[color =OrangeRed2, penwidth = 3.9, label = '-0.39', fontsize = 15, fontname = Helvetica, minlen = 1.5]
V->MaxTemp[color = grey, penwidth = 0.75, minlen =1.5] 
TH->MaxTemp[color = grey, penwidth = 1.8, minlen =1.5]
V-> NtoPRatio[color = grey, penwidth = 2.8, minlen =1.5]
TH->NtoPRatio[color = grey, penwidth=2.8, minlen =1.5]
TH->NEP[color=grey,penwidth=4.5,minlen=1.5]
V-> pH[color = grey, penwidth =2.9, minlen = 1.5]
TH->pH[color = grey, penwidth =4.0, minlen = 1.5]
TH->NEC[color =grey, penwidth=2.7,minlen=1.5]

MMC-> NEP[color = grey, penwidth =2.5, fontsize = 15, fontname = Helvetica, minlen =1.5]
NtoPRatio-> NEP[color = OrangeRed2, penwidth = 5.3, label ='-0.53', fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]
NEP -> pH[color = SteelBlue3, penwidth =3.9,label='0.39' , fontsize = 15, fontname = Helvetica, minlen =1.5]
MaxTemp -> NEC[color = OrangeRed2, penwidth = 3.4, label ='-0.34', fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]
pH -> NEC[color = grey, penwidth =2.4, minlen =1.5] 

# add a graph statement
graph [nodesep = 0.1]
}")
MytilusNight

#export image to a PDF
MytilusNight%>%
  export_svg %>% 
  charToRaw %>%
  rsvg_pdf("Output/MytilusNightDAG.pdf", width =375, height =325)


#Dag sem for surfgrass Day 
PhyllospadixDay<-
  grViz("
digraph boxes_and_circles {

# add node statements
node [shape  = rectangle,
      penwidth = 2.5,
      fontname = Helvetica,
      fontsize = 25]

PhyllospadixLoss[label = 'Surfgrass loss \n(Phyllospadix spp.)'];NEC; NEP; pH; 
MMC[label = 'Micro/macroalgae cover']; 
V[label = 'Volume'];
MaxTemp[label = 'Temperature']; 
NtoPRatio[label = 'N:P']; 
TH[label='Tide height']

#End statements
# color red is significant negative; where as black is significant positive; 
#grey is nonsignificant (p>0.05)
# size of path represents the standardized effect size * 9 labels are significant standardized ES
#** in label indicates significant day/night interaction
#check within funciton how to get labels on arrows
PhyllospadixLoss-> MMC[color = SteelBlue3, penwidth = 5.6, label = '0.55',fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]
PhyllospadixLoss-> MaxTemp[color = SteelBlue3, penwidth = 7.2, label = '0.72',fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]
PhyllospadixLoss-> NtoPRatio[color = grey, penwidth =2.3, minlen = 1.5]
PhyllospadixLoss->pH[color=SteelBlue3, penwidth = 3.1, label = '0.31',fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]

V-> MMC[color = grey, penwidth = 0.6, minlen = 1.5] 
TH->MMC[color=grey,penwidth=1.9,minlen=1.5]
V->MaxTemp[color = grey, penwidth = 0.3, minlen = 1.5]
TH->MaxTemp[color =grey, penwidth=2.1, minlen=1.5]
V-> NtoPRatio[color = grey, penwidth = 1.1, minlen = 1.5] 
TH->NtoPRatio[color=grey,penwidth=1.4,minlen=1.5]
TH->NEP[color=grey,penwidth=2.1,minlen=1.5]
V-> pH[color = grey, penwidth = 0.1,fontsize = 15, minlen = 1.5]
TH->pH[color=grey, penwidth=1.2,minlen=1.5]
TH->NEC[color = OrangeRed2, penwidth = 2.7, label ='-0.27', fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]


MMC-> NEP[color =  grey, penwidth =0.1, minlen = 1.5]
MaxTemp -> NEP[color = SteelBlue3, penwidth = 5.4, label ='0.54', fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]
NtoPRatio-> NEP[color = OrangeRed2, penwidth = 2.2, label ='0.22', fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]
NEP -> pH[color = grey, penwidth = 4.2, minlen = 1.5]
MaxTemp -> NEC[color = grey, penwidth = 2.6, minlen = 1.5]
pH -> NEC[color = SteelBlue3, penwidth = 6.4, label = '0.64',fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]

# add a graph statement
graph [nodesep = 0.1]
}")
PhyllospadixDay

#export image to a PDF
PhyllospadixDay %>%
  export_svg %>% 
  charToRaw %>% 
  rsvg_pdf("Output/PhyllospadixDay.pdf",width = 375, height=325)


#Dag sem for surfgrass night SEM
#** in label indicates significant day/night interaction
PhyllospadixNight<-
  grViz("
digraph boxes_and_circles {

# add node statements
node [shape  = rectangle,
      penwidth = 2.5,
      fontname = Helvetica,
      fontsize = 25]

PhyllospadixLoss[label = 'Surfgrass loss \n(Phyllospadix spp.)'];NEC; NEP; pH; 
MMC[label = 'Micro/macroalgae cover']; 
V[label = 'Volume'];
MaxTemp[label = 'Temperature']; 
NtoPRatio[label = 'N:P']; 
TH[label='Tide height']

#End statements
# color red is significant negative; where as black is significant positive; 
#grey is nonsignificant (p>0.05)
# size of path represents the standardized effect size * 9 labels are significant standardized ES
#check within funciton how to get labels on arrows
PhyllospadixLoss-> MMC[color = SteelBlue3, penwidth = 5.6, label = '0.56',fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]
PhyllospadixLoss-> MaxTemp[color = grey,penwidth = 2.4, minlen = 1.5]
PhyllospadixLoss-> NtoPRatio[color = grey, penwidth =2.5, minlen = 1.5]
PhyllospadixLoss->pH[color = SteelBlue3, penwidth = 9.9, label = '0.89',fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]

V-> MMC[color = grey, penwidth = 0.6, minlen = 1.5] 
TH->MMC[color=grey,penwidth=1.9,minlen=1.5]
V->MaxTemp[color = grey, penwidth = 0.8, minlen = 1.5]
TH->MaxTemp[color =grey, penwidth=5.0, minlen=1.5]
V-> NtoPRatio[color = grey, penwidth = 1.1, minlen = 1.5] 
TH->NtoPRatio[color=grey,penwidth=1.4,minlen=1.5]
TH->NEP[color=grey,penwidth=4.1,minlen=1.5]
V-> pH[color = grey, penwidth = 0.3,fontsize = 15, minlen = 1.5]
TH->pH[color=grey, penwidth=3.4,minlen=1.5]
TH->NEC[color = OrangeRed2, penwidth = 3.8, label ='-0.38', fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]

MMC-> NEP[color =  grey, penwidth =0.2, minlen = 1.5]
MaxTemp -> NEP[color = SteelBlue3, penwidth = 4.3, label ='0.43', fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]
NtoPRatio-> NEP[color = OrangeRed2, penwidth = 4.1, label ='-0.41', fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]
NEP->pH[color = grey, penwidth = 3.3, minlen = 1.5]
MaxTemp -> NEC[color = grey, penwidth = 1.5,fontsize = 15,minlen = 1.5]
pH -> NEC[color = SteelBlue3, penwidth = 3.1, label = '0.31',fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]

# add a graph statement
graph [nodesep = 0.1]
}")
PhyllospadixNight

#export image to a PDF
PhyllospadixNight %>%
  export_svg %>% 
  charToRaw %>% 
  rsvg_pdf("Output/PhyllospadixNight.pdf",width = 375, height=325)

