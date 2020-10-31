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
# size of path represents the standardized coefficients * 8 
# labels are significant standardized coefficients


MytilusLoss-> MMC[color = SteelBlue3, penwidth = 5.8, label = 0.72,fontsize = 13,
fontname = Helvetica, minlen = 1.5]
MytilusLoss -> Light[color = SteelBlue3, penwidth = 4.7, label = 0.58,fontsize = 13,fontname = Helvetica, minlen = 1.5]
MytilusLoss -> MaxTemp[color = SteelBlue3, penwidth = 5.2, label = 0.64,fontsize = 13,fontname = Helvetica, minlen = 1.5]
MytilusLoss-> NtoPRatio[color =  grey, penwidth = 3, fontsize = 13,fontname = Helvetica, minlen =1.5]
MytilusLoss -> NEC[color = grey, penwidth =4.1, minlen =1.5] #0.2068*5=1.1
SAtoVRatio-> MMC[color = grey, penwidth = 3.2, minlen =1.5] #-0.2808*5=1.4
SAtoVRatio->Light[color = grey, penwidth = 1.1, minlen =1.5] #0.1117*5=0.56
SAtoVRatio->MaxTemp[color = grey, penwidth = 1.7, minlen =1.5] #0.2534*5=1.3
SAtoVRatio-> NtoPRatio[color = OrangeRed2, penwidth = 5.3, label = -0.66,fontsize = 13,fontname = Helvetica, minlen =1.5]
SAtoVRatio-> NEP[color = grey, penwidth = 1.1, minlen =1.5] #0.1958*5=0.979
SAtoVRatio-> pH[color = grey, penwidth = 2.0, minlen =1.5] #0.1526*5=0.763
SAtoVRatio-> NEC[color = grey, penwidth = 3.6, minlen =1.5] #0.1319 *5 =0.66
TH-> MMC[color =OrangeRed2, penwidth = 3.7, label = 0.47, fontsize = 13, fontname = Helvetica, minlen = 1.5]
TH-> Light[color =grey, penwidth = 0.8, minlen=1.5]
TH->MaxTemp[color = grey, penwidth = 1.2, minlen =1.5]
TH->NtoPRatio[color = grey, penwidth=2.8, minlen =1.5]
TH->NEP[color=grey,penwidth=3,minlen=1.5]
TH->pH[color=grey,penwidth=2.3,minlen=1.5]
TH->NEC[color =grey, penwidth=4,minlen=1.5]
MMC-> NEP[color = SteelBlue3, penwidth =5.8, label = 0.73, fontsize = 13, fontname = Helvetica, minlen =1.5]
MaxTemp -> NEP[color = grey, penwidth = 0.7, minlen =1.5] #0.0725*5=0.36
NtoPRatio-> NEP[color = grey, penwidth = 1.6, minlen =1.5] #.108*5=0.54
NEP -> pH[color = grey, penwidth = 2.0, minlen =1.5] #0.0381*5=0.1905
MaxTemp -> NEC[color = grey, penwidth = 2.6, minlen =1.5] #.122*5=0.61
pH -> NEC[color = grey, penwidth =0.75, minlen =1.5] #-.505*5=-2.5
MytilusLoss->pH[color = SteelBlue3, penwidth =6.88, label = '0.86',fontsize = 13, fontname = Helvetica, minlen = 1.5]
Light->MaxTemp[color = black, penwidth =6.5, label = 'CE 0.81', fontsize = 13, fontname = Helvetica,minlen = 1.5]

# add a graph statement
graph [nodesep = 0.1]
}")
Mytilus

#export image to a PDF
M<-Mytilus %>%
  export_svg %>% 
  charToRaw %>%
  rsvg_pdf("Output/MytilusDAG.pdf", width =700, height =700)


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
PhyllospadixLoss-> MMC[color = SteelBlue3, penwidth = 5.4, label = 0.67,fontsize = 13, fontface = bold,fontname = Helvetica, minlen = 1.5]
PhyllospadixLoss -> Light[color = SteelBlue3, penwidth = 5.3, label = 0.67,fontsize = 13, fontface = bold,fontname = Helvetica, minlen = 1.5]
PhyllospadixLoss-> MaxTemp[color = grey, penwidth =3.6, minlen = 1.5]#0.3397*5=1.3885
PhyllospadixLoss-> NtoPRatio[color = grey, penwidth = 1.0, minlen = 1.5] #0.0308*5=0.154
PhyllospadixLoss->pH[color =grey, penwidth = 2.6, minlen =1.5] #0.2898*5 =1.449
PhyllospadixLoss->NEC[color=grey, penwidth =3.3, minlen =1.5] #0.1425*5=0.712
MMC-> NEP[color =  grey, penwidth =0.19, minlen = 1.5]
Light -> NEP[color = SteelBlue3, penwidth = 4.6, label =0.57, fontsize = 13, fontface = bold,fontname = Helvetica, minlen = 1.5]
NtoPRatio-> NEP[color = grey, penwidth = 0.8, minlen = 1.5]
NEP -> pH[color = SteelBlue3, penwidth = 3.0, label = 0.37,fontsize = 13, fontface = bold,fontname = Helvetica, minlen = 1.5]
MaxTemp -> NEC[color = SteelBlue3, penwidth = 4.5, label=0.56,fontsize = 13, fontface = bold,fontname = Helvetica,minlen = 1.5] #0.37*5 =1.8
pH -> NEC[color = SteelBlue3, penwidth = 9.0, label = 1.1,fontsize = 13, fontface = bold,fontname = Helvetica, minlen = 1.5]
SAtoVRatio-> MMC[color = grey, penwidth = 1.1, minlen = 1.5] #0.0562*5=0.281
SAtoVRatio->Light[color = grey, penwidth = 0.5, minlen = 1.5] #0.1602*5=0.801
SAtoVRatio->MaxTemp[color = grey, penwidth = 1.0, minlen = 1.5] #0.2767*5=1.38
SAtoVRatio-> NtoPRatio[color = grey, penwidth = 2.8, minlen = 1.5] #0.2235*5= 1.1175
SAtoVRatio-> NEP[color = grey, penwidth = 0.6, minlen = 1]#0.0316*5 =0.158
SAtoVRatio-> pH[color = SteelBlue3, penwidth = 3.7, label = 0.46,fontsize = 13, fontface = bold,fontname = Helvetica, minlen = 1.5]
SAtoVRatio-> NEC[color = OrangeRed2, penwidth = 5.2, label = -0.65,fontsize = 13, fontface = bold,
fontname = Helvetica, minlen = 1.5]
TH->Light[color =grey, penwidth=0.07, minlen =1.5]
TH->MMC[color=grey,penwidth=0.3,minlen=1.5]
TH->MaxTemp[color =grey, penwidth=2.8, minlen=1.5]
TH->NtoPRatio[color=grey,penwidth=2.5,minlen=1.5]
TH->NEP[color=grey,penwidth=1,minlen=1.5]
TH->pH[color=grey, penwidth=1.3,minlen=1.5]
TH->NEC[color = OrangeRed2, penwidth = 3.9, label = -0.49,fontsize = 13, fontface = bold,
fontname = Helvetica, minlen = 1.5]
Light ->MaxTemp[color = black, penwidth = 6.1, label = 'CE 0.76',fontsize = 13, fontface = bold,
fontname = Helvetica, minlen = 1.5]
# add a graph statement
graph [nodesep = 0.1]
}")
Phyllospadix

#export image to a PDF
Phyllospadix %>%
  export_svg %>% 
  charToRaw %>% 
  rsvg_pdf("Output/PhyllospadixDAG.pdf",width = 700, height=700)