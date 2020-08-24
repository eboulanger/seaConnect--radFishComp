**ce README n'est pas mis à jour. Pour l'instant ce dépot git sers à vous partager les résultats/figures des analyses que je fait pour SEACONNECT au fur et à mesure.**

# seaConnect--radFishComp
Comparison of fish RAD sequence data: compare different species and different marker types for population/fisheries genomic

## Data
RADseq SNP calling was performed by [P. Guerrin using a dDocent-based pipeline](https://github.com/Grelot/seaConnect--dDocent)    
SNPs were then filtered for missing data, maf, ... and outliers were detected to retain 
only adaptive or neutral loci for both species. 
Finally, datasets are converted to the correct formats for different analyses.  
The scripts can be found in the project [seaConnect--dataPrep](https://github.com/eboulanger/seaConnect--dataPrep)      

## Dependencies
You will need to install the following software:  
- [Genodive](http://bentleydrummer.nl/software/software/GenoDive.html) 
- [NeEstimator](http://www.molecularfisherieslaboratory.com.au/neestimator-software/) 

And you will need to have the following R libraries:  
ggplot2, reshape2, pophelper, maps, mapdata, dplyr, plotrix, adegenet, ...

## 01-Kinship

## 02-DAPC

Infer population structure with DAPC from [adegenet](http://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf).

Interactive R-scripts `dip_dapc.R` and `mul_dapc.R
This script runs two types of DAPC:  
- with the Mediterranean ecoregions set as prior groups
- with the detection of best-fit number of populations  

For both a scatterplot is created. For the latter, a pie map with the posterior membership 
probabilities is created comparable to the ancestry pie maps created in `01-ADMIXTURE`,
as well as a barplot of the posterior membership probabilities with individuals ordered by longitude.

## 03-RDA

Constrained ordination to infer if genetic structure is driven by certain 
environmental variables & look for signals of local adaptation.

Environmental variables considered:  
- Sea Surface Salinity  
- Sea Surface Temperature  
- habitat: substrate PC axes  
- geographic distance: dbMEMs (18 for D. sargus, 20 for M. surmuletus)  
- larval connectivity: AEMs   (97 for D. sargus and M. surmuletus)  

We have too many potential explanatory variables so we'll go through several variable selection steps.  
- `compile_explanatory_vars.R` : run separate RDAs with dbMEMs or AEMs as expl variables.
Use ordi2step to perform (forward?) variable selection.
- `run_rda_neutral|adaptive.R` : run RDA with previously selected dbMEMs and AEMs, as well as SST, SSS and habitat PCs. 









