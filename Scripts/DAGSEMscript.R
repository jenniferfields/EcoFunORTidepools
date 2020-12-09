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
Mytilus <- grViz("
digraph boxes_and_circles { 
#this is the type of diagram you are using

# add node statements 
node [shape  = rectangle,
      penwidth =2.5,
      fontname = Helvetica,
      fontsize = 18]

MytilusLoss[label = 'CA Mussel ( M. californianus ) Loss'];NEC; NEP; pH; 
MMC[label = 'Micro/Macroalgae Cover']; 
SAtoVRatio[label = 'SA:V']; Light;
MaxTemp[label = 'Temperature']; 
NtoPRatio[label = 'N:P']; 
TH[label = 'Tide Height']

#End statements
# color red is significant negative; where as blue is significant positive; 
#grey is nonsignificant (p>0.05)
# size of path represents the standardized coefficients * 9 
# labels are significant standardized coefficients


MytilusLoss-> MMC[color = SteelBlue3, penwidth =7.2, label = '0.80',fontsize = 13,
fontname = Helvetica, minlen = 1.5]
MytilusLoss -> Light[color = grey, penwidth = 4.5, minlen = 1.5]
MytilusLoss -> MaxTemp[color = grey, penwidth = 5.0, minlen = 1.5]
MytilusLoss-> NtoPRatio[color =  grey, penwidth = 3.5, minlen =1.5]
MytilusLoss->pH[color = SteelBlue3, penwidth =6.83, label = '0.76',fontsize = 13, fontname = Helvetica, minlen = 1.5]
MytilusLoss -> NEC[color = grey, penwidth =5.0, minlen =1.5] 
SAtoVRatio-> MMC[color = grey, penwidth = 3.0, minlen =1.5] 
SAtoVRatio->Light[color = grey, penwidth = 1.7, minlen =1.5]
SAtoVRatio->MaxTemp[color = grey, penwidth = 2.4, minlen =1.5] 
SAtoVRatio-> NtoPRatio[color = OrangeRed2, penwidth = 6.21, label = '-0.69',fontsize = 13,fontname = Helvetica, minlen =1.5]
SAtoVRatio-> NEP[color = grey, penwidth = 1.2, minlen =1.5] 
SAtoVRatio-> pH[color = grey, penwidth = 2.6, minlen =1.5] 
SAtoVRatio-> NEC[color = grey, penwidth = 3.6, minlen =1.5] 
TH-> MMC[color =OrangeRed2, penwidth = 4.45, label = '-0.50', fontsize = 13, fontname = Helvetica, minlen = 1.5]
TH-> Light[color =grey, penwidth = 1.1, minlen=1.5]
TH->MaxTemp[color = grey, penwidth = 1.1, minlen =1.5]
TH->NtoPRatio[color = grey, penwidth=3.1, minlen =1.5]
TH->NEP[color=grey,penwidth=3.7,minlen=1.5]
TH->pH[color=grey,penwidth=1.9,minlen=1.5]
TH->NEC[color =grey, penwidth=4.3,minlen=1.5]
MMC-> NEP[color = SteelBlue3, penwidth =6.6, label = '0.74', fontsize = 13, fontname = Helvetica, minlen =1.5]
Light -> NEP[color = grey, penwidth = 1.1, minlen =1.5] 
NtoPRatio-> NEP[color = grey, penwidth = 1.9, minlen =1.5] 
NEP -> pH[color = grey, penwidth = 1.6, minlen =1.5] 
MaxTemp -> NEC[color = grey, penwidth = 2.4, minlen =1.5]
pH -> NEC[color = grey, penwidth =1.6, minlen =1.5] 
Light->MaxTemp[color = black, penwidth =1, label = 'CE 0.83', fontsize = 13, fontname = Helvetica,minlen = 1.5]

# add a graph statement
graph [nodesep = 0.1]
}")
Mytilus

#export image to a PDF
M<-Mytilus %>%
  export_svg %>% 
  charToRaw %>%
  rsvg_pdf("Output/MytilusDAG.pdf", width =450, height =450)


#Dag sem for surfgrass SEM
#With light in NEP model
Phyllospadix<-
  grViz("
digraph boxes_and_circles {

# add node statements
node [shape  = rectangle,
      penwidth = 2.5,
      fontname = Helvetica,
      fontsize = 18]

PhyllospadixLoss[label = 'Surfgrass (Phyllospadix spp.) Loss'];NEC; NEP; pH; 
MMC[label = 'Micro/Macroalgae Cover']; 
SAtoVRatio[label = 'SA:V']; Light;
MaxTemp[label = 'Temperature']; 
NtoPRatio[label = 'N:P']; 
TH[label='Tide Height']

#End statements
# color red is significant negative; where as black is significant positive; 
#grey is nonsignificant (p>0.05)
# size of path represents the standardized effect size * 8 labels are significant standardized ES
#check within funciton how to get labels on arrows
PhyllospadixLoss-> MMC[color = SteelBlue3, penwidth = 5.9, label = '0.65',fontsize = 13, fontface = bold,fontname = Helvetica, minlen = 1.5]
PhyllospadixLoss -> Light[color = SteelBlue3, penwidth = 7.9, label = '0.87',fontsize = 13, fontface = bold,fontname = Helvetica, minlen = 1.5]
PhyllospadixLoss-> MaxTemp[color = SteelBlue3, penwidth = 5.9, label = '0.66',fontsize = 13, fontface = bold,fontname = Helvetica, minlen = 1.5]
PhyllospadixLoss-> NtoPRatio[color = grey, penwidth =0.4, minlen = 1.5]
PhyllospadixLoss->pH[color=SteelBlue3, penwidth = 3.6, label = '0.40',fontsize = 13, fontface = bold,fontname = Helvetica, minlen = 1.5]
PhyllospadixLoss->NEC[color=grey, penwidth =3.9, minlen =1.5]
MMC-> NEP[color =  grey, penwidth =0.11, minlen = 1.5]
Light -> NEP[color = SteelBlue3, penwidth = 5.1, label =0.57, fontsize = 13, fontface = bold,fontname = Helvetica, minlen = 1.5]
NtoPRatio-> NEP[color = grey, penwidth = 0.92, minlen = 1.5]
NEP -> pH[color = grey, penwidth = 2.7,  minlen = 1.5]
MaxTemp -> NEC[color = SteelBlue3, penwidth = 5.7, label=0.63,fontsize = 13, fontface = bold,fontname = Helvetica,minlen = 1.5] #0.37*5 =1.8
pH -> NEC[color = SteelBlue3, penwidth = 10.3, label = 1.1,fontsize = 13, fontface = bold,fontname = Helvetica, minlen = 1.5]
SAtoVRatio-> MMC[color = grey, penwidth = 1.5, minlen = 1.5] 
SAtoVRatio->Light[color = grey, penwidth = 2.0, minlen = 1.5] 
SAtoVRatio->MaxTemp[color = grey, penwidth = 0.21, minlen = 1.5]
SAtoVRatio-> NtoPRatio[color = grey, penwidth = 2.3, minlen = 1.5] 
SAtoVRatio-> NEP[color = grey, penwidth = 0.69, minlen = 1]
SAtoVRatio-> pH[color = SteelBlue3, penwidth = 3.8, label = '0.42',fontsize = 13, fontface = bold,fontname = Helvetica, minlen = 1.5]
SAtoVRatio-> NEC[color = OrangeRed2, penwidth = 5.9, label = '-0.65',fontsize = 13, fontface = bold,
fontname = Helvetica, minlen = 1.5]
TH->Light[color =grey, penwidth=0.54, minlen =1.5]
TH->MMC[color=grey,penwidth=1.2,minlen=1.5]
TH->MaxTemp[color =grey, penwidth=2.8, minlen=1.5]
TH->NtoPRatio[color=grey,penwidth=2.5,minlen=1.5]
TH->NEP[color=grey,penwidth=1.1,minlen=1.5]
TH->pH[color=grey, penwidth=1.0,minlen=1.5]
TH->NEC[color = OrangeRed2, penwidth = 4.1, label = '-0.46',fontsize = 13, fontface = bold,
fontname = Helvetica, minlen = 1.5]
Light ->MaxTemp[color = black, penwidth = 1, label = 'CE 0.66',fontsize = 13, fontface = bold,
fontname = Helvetica, minlen = 1.5]
# add a graph statement
graph [nodesep = 0.1]
}")
Phyllospadix

#export image to a PDF
Phyllospadix %>%
  export_svg %>% 
  charToRaw %>% 
  rsvg_pdf("Output/PhyllospadixDAG.pdf",width = 450, height=450)
