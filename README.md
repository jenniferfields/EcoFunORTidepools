# Foundation species loss alters ecosystem functioning within temperate tidepool communities 

This repository includes data, analysis, and methods for the following manuscript

**Authors:** Jennifer B. Fields*, Nyssa J. Silbiger

**Journal:** Marine Ecology Progress Series

**Link:** https://doi.org/10.3354/meps13978

[![DOI](https://zenodo.org/badge/307201313.svg)](https://zenodo.org/badge/latestdoi/307201313)


**Funding and permit support:** Funding for this project was provided by California State University, Northridge (CSUN) Research and Graduate Studies, CSUN Department of Biology, CSU Council on Ocean Affairs, Science & Technology Graduate Research Grant, Association of Retired Faculty, Sigma Xi Grants-in-Aid of Research, Dr. Julie Gorchynski research grant, and NSF BIO-OCE #2044837 to NJS. This project was conducted under ODFW Scientific Taking Permit #22930 and Oregon Parks and Recreation Department Scientific Take Permit #017-19 as well as California DFW permit #13047 and CA State Parks permit #18-820-50 for preliminary testing. This is CSUN contribution #367.


**Contents:** There are four folders, R project, and a README.md file.

**Folders:**

**[Data/Clean Data](https://github.com/jenniferfields/EcoFunORTidepools/tree/master/Data)**

Contains all clean and raw data from project. Within this folder there are four folders for data from 
1. Tide pool physical parameters
2. Community composition
3. Light and temperature loggers
4. Biogeochemistry sampling

**[Coding scripts](https://github.com/jenniferfields/EcoFunORTidepools/tree/master/Scripts)**

Contains all code for project. Within this folder there are six files:

*[Tide pool physical parameters script](https://github.com/jenniferfields/EcoFunORTidepools/blob/master/Scripts/tidepoolphysicalparameters.R)*

Contains calculations for tide pools volumes of n=16 surfgrass and n=16 mussel tide pools at Otter Rock, OR utilizing methods and code from *[Silbiger Lab Volume Dye Method](https://github.com/SilbigerLab/Protocols/tree/master/Environmental_Parameter_Protocols/Protocols/Volume_Dye_Method)*

*[Community compositon script](https://github.com/jenniferfields/EcoFunORTidepools/blob/master/Scripts/CommunityComp.R)*

Code to investigate how tide pool sessile and mobile communities were altered by surfgrass and mussel removal. Contains PCoA plot code and PERMANOVA analysis of community structure and community composition supplemental material plots

*[Temperature and light script](https://github.com/jenniferfields/EcoFunORTidepools/blob/master/Scripts/TemperatureandLight.R)*

Code to test how surfgrass and mussel loss affected maximum temperature and percent light values of tide pools. Contains temperature and light logger data cleaning and month-long temperature and light analyses. Includes nearshore ocean temperature over study period from ODFW Marine Reserves mooring sensor at 1 meter depth within Otter Rock Marine Reserve

*[Map script](https://github.com/jenniferfields/EcoFunORTidepools/blob/master/Scripts/OregonMaps.R)*

Code for both West Coast and Otter Rock Marine Reserve and Marine Garden maps along the Oregon Coast for the project location in Figure 1 of the manuscript.

*[Carbonate and biogeochemistry script](https://github.com/jenniferfields/EcoFunORTidepools/blob/master/Scripts/CleanCarbChem.R)* 

Code to characterize differences between time points (before and after foundation species removal) of local tide pool and adjacent ocean biogeochemistry and temperature. Contains carbonate chemistry and ecosystem metabolism (Net ecosystem calcification and production) calculations. Includes Principal Component Analyses and plots of maximum, average, and variance values of biogeochemistry (dissolved oxygen, pH, nutrients) and temperature tide pools and adjacent ocean

*[Structural Equation Model script](https://github.com/jenniferfields/EcoFunORTidepools/blob/master/Scripts/SEMScript.R)*

Code to test the direct and indirect effects of mussel and surfgrass foundation species removal on ecosystem metabolism mediated by changes in temperature, nutrients, algal cover, and pH. Contains structural equation model (SEM) input data, PiecewiseSEM model framework, and significant marginal effect plots. Note: This script ultilizes output from both carbonate and biogeochemistry and community composition scripts

*[DAG of Structural Equation Model script](https://github.com/jenniferfields/EcoFunORTidepools/blob/master/Scripts/DAGSEMscript.R)* 

Contains code using DiagrammeR package to create SEM figures

**[Output](https://github.com/jenniferfields/EcoFunORTidepools/tree/master/Output)**

Contains all figures of manuscript

**[Protocols](https://github.com/jenniferfields/EcoFunORTidepools/tree/master/Protocols)**

Contains protocols from field experiment including:

[Tide pool dye method](https://github.com/jenniferfields/EcoFunORTidepools/blob/master/Protocols/Dye_Method_Protocol.md)

[General biogeochemical sampling protocol](https://github.com/jenniferfields/EcoFunORTidepools/blob/master/Protocols/TidePoolSampling_SOP.md)

[Group biogeochemical sampling protocol](https://github.com/jenniferfields/EcoFunORTidepools/blob/master/Protocols/GroupWaterSampling_SOP.md)

[Community composition surveys](https://github.com/jenniferfields/EcoFunORTidepools/blob/master/Protocols/CommunityComposition_SOP.md)

[Foundation species removal](https://github.com/jenniferfields/EcoFunORTidepools/blob/master/Protocols/FoundationSpeciesRemoval_SOP.md)



