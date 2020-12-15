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

#SEM dag for mussel SEM
#With temp in NEP model
MytilusDay <- grViz("
digraph boxes_and_circles { 
#this is the type of diagram you are using

# add node statements 
node [shape  = rectangle,
      penwidth =2.5,
      fontname = Helvetica,
      fontsize = 20]

MytilusLoss[label = 'CA Mussel Loss \n(M. californianus)'];NEC; NEP; pH; 
MMC[label = 'Micro/Macroalgae Cover']; 
SAtoVRatio[label = 'SA:V']; 
MaxTemp[label = 'Temperature']; 
NtoPRatio[label = 'N:P']; 
TH[label = 'Tide Height']

#End statements
# color red is significant negative; where as blue is significant positive; 
#grey is nonsignificant (p>0.05)
# size of path represents the standardized coefficients * 9 
# labels are significant standardized coefficients


MytilusLoss-> MMC[color = SteelBlue3, penwidth =7.2, label = '0.80',fontsize = 15,
fontname = Helvetica, minlen = 1.5]
MytilusLoss -> MaxTemp[color = SteelBlue3, penwidth =3.13, label = '0.35',fontsize = 15,
fontname = Helvetica, minlen = 1.5]
MytilusLoss-> NtoPRatio[color =  grey, penwidth = 3.8, minlen =1.5]
MytilusLoss->pH[color = SteelBlue3, penwidth =4.27, label = '0.47',fontsize = 15, fontname = Helvetica, minlen = 1.5]
MytilusLoss -> NEC[color = grey, penwidth =1.3, minlen =1.5] 

SAtoVRatio-> MMC[color =OrangeRed2, penwidth = 3.04, label = '-0.34', fontsize = 15, fontname = Helvetica, minlen = 1.5]
TH-> MMC[color =OrangeRed2, penwidth = 4.46, label = '-0.50', fontsize = 15, fontname = Helvetica, minlen = 1.5]
SAtoVRatio->MaxTemp[color = grey, penwidth = 1.06, minlen =1.5] 
TH->MaxTemp[color = grey, penwidth = 0.76, minlen =1.5]
SAtoVRatio-> NtoPRatio[color = OrangeRed2, penwidth = 5.62, label = '-0.62',fontsize = 15,fontname = Helvetica, minlen =1.5]
TH->NtoPRatio[color = grey, penwidth=3.9, minlen =1.5]
SAtoVRatio-> NEP[color = grey, penwidth = 0.87, minlen =1.5] 
TH->NEP[color=grey,penwidth=2.62,minlen=1.5]
SAtoVRatio-> pH[color = SteelBlue3, penwidth =2.13, label = '0.24',fontsize = 15, fontname = Helvetica, minlen = 1.5]
TH->pH[color = SteelBlue3, penwidth =2.0, label = '0.22',fontsize = 15, fontname = Helvetica, minlen = 1.5]
SAtoVRatio-> NEC[color = grey, penwidth = 2.45, minlen =1.5] 
TH->NEC[color =grey, penwidth=3.44,minlen=1.5]

MaxTemp -> NEP[color = grey, penwidth =1.2, label = '**', fontsize = 15, fontname = Helvetica, minlen =1.5]
MMC-> NEP[color = SteelBlue3, penwidth =8.63, label = '0.96**', fontsize = 15, fontname = Helvetica, minlen =1.5]
NtoPRatio-> NEP[color = OrangeRed2, penwidth = 2.24, label ='0.25', fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]
NEP -> pH[color = SteelBlue3, penwidth =3.61,label='0.40' , fontsize = 15, fontname = Helvetica, minlen =1.5]
MaxTemp -> NEC[color = OrangeRed2, penwidth = 5.47, label ='-0.61', fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]
pH -> NEC[color = grey, penwidth =3.28, minlen =1.5] 

# add a graph statement
graph [nodesep = 0.1]
}")
MytilusDay

#export image to a PDF
MytilusDay %>%
  export_svg %>% 
  charToRaw %>%
  rsvg_pdf("Output/MytilusDayDAG.pdf", width =450, height =350)

#Dag sem for mussel night
MytilusNight <- grViz("
digraph boxes_and_circles { 
#this is the type of diagram you are using

# add node statements 
node [shape  = rectangle,
      penwidth =2.5,
      fontname = Helvetica,
      fontsize = 20]

MytilusLoss[label = 'CA Mussel Loss \n(M. californianus)'];NEC; NEP; pH; 
MMC[label = 'Micro/Macroalgae Cover']; 
SAtoVRatio[label = 'SA:V']; 
MaxTemp[label = 'Temperature']; 
NtoPRatio[label = 'N:P']; 
TH[label = 'Tide Height']

#End statements
# color red is significant negative; where as blue is significant positive; 
#grey is nonsignificant (p>0.05)
# size of path represents the standardized coefficients * 9 
# labels are significant standardized coefficients


MytilusLoss-> MMC[color = SteelBlue3, penwidth =7.2, label = '0.80',fontsize = 15,
fontname = Helvetica, minlen = 1.5]
MytilusLoss -> MaxTemp[color = SteelBlue3, penwidth =5.7, label = '0.64',fontsize = 15,
fontname = Helvetica, minlen = 1.5]
MytilusLoss-> NtoPRatio[color =  grey, penwidth = 2.2, minlen =1.5]
MytilusLoss->pH[color = SteelBlue3, penwidth =7.5, label = '0.84',fontsize = 15, fontname = Helvetica, minlen = 1.5]
MytilusLoss -> NEC[color = grey, penwidth =1.3, minlen =1.5] 

SAtoVRatio-> MMC[color =OrangeRed2, penwidth = 3.04, label = '-0.34', fontsize = 15, fontname = Helvetica, minlen = 1.5]
TH-> MMC[color =OrangeRed2, penwidth = 4.46, label = '-0.50', fontsize = 15, fontname = Helvetica, minlen = 1.5]
SAtoVRatio->MaxTemp[color = grey, penwidth = 1.95, minlen =1.5] 
TH->MaxTemp[color = grey, penwidth = 1.4, minlen =1.5]
SAtoVRatio-> NtoPRatio[color = OrangeRed2, penwidth = 3.31, label = '-0.37',fontsize = 15,fontname = Helvetica, minlen =1.5]
TH->NtoPRatio[color = grey, penwidth=2.27, minlen =1.5]
SAtoVRatio-> NEP[color = grey, penwidth = 1.6, minlen =1.5] 
TH->NEP[color=grey,penwidth=4.8,minlen=1.5]
SAtoVRatio-> pH[color = SteelBlue3, penwidth =3.8, label = '0.42',fontsize = 15, fontname = Helvetica, minlen = 1.5]
TH->pH[color = SteelBlue3, penwidth =3.5, label = '0.39',fontsize = 15, fontname = Helvetica, minlen = 1.5]
SAtoVRatio-> NEC[color = grey, penwidth = 2.45, minlen =1.5] 
TH->NEC[color =grey, penwidth=3.49,minlen=1.5]

MaxTemp -> NEP[color = grey, penwidth =6.23, label = '0.69**', fontsize = 15, fontname = Helvetica, minlen =1.5]
MMC-> NEP[color = grey, penwidth =2.4, label = '**', fontsize = 15, fontname = Helvetica, minlen =1.5]
NtoPRatio-> NEP[color = OrangeRed2, penwidth = 6.92, label ='0.77', fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]
NEP -> pH[color = SteelBlue3, penwidth =3.48,label='0.39' , fontsize = 15, fontname = Helvetica, minlen =1.5]
MaxTemp -> NEC[color = OrangeRed2, penwidth = 3.03, label ='-0.34', fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]
pH -> NEC[color = grey, penwidth =1.89, minlen =1.5] 

# add a graph statement
graph [nodesep = 0.1]
}")
MytilusNight

#export image to a PDF
MytilusNight%>%
  export_svg %>% 
  charToRaw %>%
  rsvg_pdf("Output/MytilusNightDAG.pdf", width =450, height =350)


#Dag sem for surfgrass Day 
PhyllospadixDay<-
  grViz("
digraph boxes_and_circles {

# add node statements
node [shape  = rectangle,
      penwidth = 2.5,
      fontname = Helvetica,
      fontsize = 20]

PhyllospadixLoss[label = 'Surfgrass Loss \n(Phyllospadix spp.)'];NEC; NEP; pH; 
MMC[label = 'Micro/Macroalgae Cover']; 
SAtoVRatio[label = 'SA:V'];
MaxTemp[label = 'Temperature']; 
NtoPRatio[label = 'N:P']; 
TH[label='Tide Height']

#End statements
# color red is significant negative; where as black is significant positive; 
#grey is nonsignificant (p>0.05)
# size of path represents the standardized effect size * 9 labels are significant standardized ES
#** in label indicates significant day/night interaction
#check within funciton how to get labels on arrows
PhyllospadixLoss-> MMC[color = SteelBlue3, penwidth = 5.9, label = '0.65',fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]
PhyllospadixLoss-> MaxTemp[color = SteelBlue3, penwidth = 6.7, label = '0.75**',fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]
PhyllospadixLoss-> NtoPRatio[color = grey, penwidth =3.1, minlen = 1.5]
PhyllospadixLoss->pH[color=grey, penwidth = 1.6, minlen = 1.5]
PhyllospadixLoss->NEC[color=grey, penwidth =0.05, minlen =1.5]

SAtoVRatio-> MMC[color = grey, penwidth = 1.5, minlen = 1.5] 
SAtoVRatio->MaxTemp[color = grey, penwidth = 0.19, minlen = 1.5]
SAtoVRatio-> NtoPRatio[color = grey, penwidth = 0.45, minlen = 1.5] 
SAtoVRatio-> NEP[color = grey, penwidth = 0.11, minlen = 1]
SAtoVRatio-> pH[color = SteelBlue3, penwidth = 3.8, label = '0.42**',fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]
SAtoVRatio-> NEC[color = OrangeRed2, penwidth = 2.68, label = '-0.30',fontsize = 15, fontface = bold,
fontname = Helvetica, minlen = 1.5]

TH->MMC[color=grey,penwidth=1.2,minlen=1.5]
TH->MaxTemp[color =grey, penwidth=1.9, minlen=1.5]
TH->NtoPRatio[color=grey,penwidth=1.6,minlen=1.5]
TH->NEP[color=grey,penwidth=1.8,minlen=1.5]
TH->pH[color=grey, penwidth=0.63,minlen=1.5]
TH->NEC[color = grey, penwidth = 1.8,  minlen = 1.5]

MMC-> NEP[color =  grey, penwidth =0.27, minlen = 1.5]
MaxTemp -> NEP[color = SteelBlue3, penwidth = 4.6, label ='0.51', fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]
NtoPRatio-> NEP[color = OrangeRed2, penwidth = 1.83, label ='0.20', fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]
NEP -> pH[color = SteelBlue3, penwidth = 5.22, label ='0.58**', fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]
MaxTemp -> NEC[color = grey, penwidth = 0.33,label ='**', minlen = 1.5]
pH -> NEC[color = SteelBlue3, penwidth = 7.12, label = '0.79**',fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]

# add a graph statement
graph [nodesep = 0.1]
}")
PhyllospadixDay

#export image to a PDF
PhyllospadixDay %>%
  export_svg %>% 
  charToRaw %>% 
  rsvg_pdf("Output/PhyllospadixDay.pdf",width = 450, height=350)


#Dag sem for surfgrass night SEM
#** in label indicates significant day/night interaction
PhyllospadixNight<-
  grViz("
digraph boxes_and_circles {

# add node statements
node [shape  = rectangle,
      penwidth = 2.5,
      fontname = Helvetica,
      fontsize = 20]

PhyllospadixLoss[label = 'Surfgrass Loss \n(Phyllospadix spp.)'];NEC; NEP; pH; 
MMC[label = 'Micro/Macroalgae Cover']; 
SAtoVRatio[label = 'SA:V'];
MaxTemp[label = 'Temperature']; 
NtoPRatio[label = 'N:P']; 
TH[label='Tide Height']

#End statements
# color red is significant negative; where as black is significant positive; 
#grey is nonsignificant (p>0.05)
# size of path represents the standardized effect size * 9 labels are significant standardized ES
#check within funciton how to get labels on arrows
PhyllospadixLoss-> MMC[color = SteelBlue3, penwidth = 5.9, label = '0.65',fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]
PhyllospadixLoss-> MaxTemp[color = grey, label ='**',fontsize = 15,penwidth = 1.0, minlen = 1.5]
PhyllospadixLoss-> NtoPRatio[color = grey, penwidth =2.2, minlen = 1.5]
PhyllospadixLoss->pH[color=grey, penwidth = 5.0, minlen = 1.5]
PhyllospadixLoss->NEC[color=grey, penwidth =0.06, minlen =1.5]

SAtoVRatio-> MMC[color = grey, penwidth = 1.5, minlen = 1.5] 
TH->MMC[color=grey,penwidth=1.2,minlen=1.5]
SAtoVRatio->MaxTemp[color = grey, penwidth = 0.46, minlen = 1.5]
TH->MaxTemp[color =grey, penwidth=4.6, minlen=1.5]
SAtoVRatio-> NtoPRatio[color = grey, penwidth = 0.32, minlen = 1.5] 
TH->NtoPRatio[color=grey,penwidth=1.6,minlen=1.5]
SAtoVRatio-> NEP[color = grey, penwidth = 0.23, minlen = 1]
TH->NEP[color=grey,penwidth=3.5,minlen=1.5]
SAtoVRatio-> pH[color = grey, penwidth = 0.07,label ='**',fontsize = 15, minlen = 1.5]
TH->pH[color=grey, penwidth=1.99,minlen=1.5]
SAtoVRatio-> NEC[color = OrangeRed2, penwidth = 3.95, label = '-0.44',fontsize = 15, fontface = bold,
fontname = Helvetica, minlen = 1.5]
TH->NEC[color = grey, penwidth = 2.7,  minlen = 1.5]

MMC-> NEP[color =  grey, penwidth =0.51, minlen = 1.5]
MaxTemp -> NEP[color = SteelBlue3, penwidth = 3.6, label ='0.40', fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]
NtoPRatio-> NEP[color = OrangeRed2, penwidth = 4.76, label ='-0.53', fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]
NEP -> pH[color = grey, penwidth =3.0,label ='**',fontsize = 15, minlen = 1.5]
MaxTemp -> NEC[color = grey, penwidth = 6.92,label ='**',fontsize = 15,minlen = 1.5]
pH -> NEC[color = SteelBlue3, penwidth = 3.37, label = '0.37',fontsize = 15, fontface = bold,fontname = Helvetica, minlen = 1.5]

# add a graph statement
graph [nodesep = 0.1]
}")
PhyllospadixNight

#export image to a PDF
PhyllospadixNight %>%
  export_svg %>% 
  charToRaw %>% 
  rsvg_pdf("Output/PhyllospadixNight.pdf",width = 450, height=350)
