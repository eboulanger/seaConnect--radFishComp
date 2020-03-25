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

## 02-DAPC

Infer population structure with DAPC from [adegenet](http://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf).

Interactive R-scripts `dip_dapc.R` and `mul_dapc.R
This script runs two types of DAPC:  
- with the Mediterranean ecoregions set as prior groups
- with the detection of best-fit number of populations  

For both a scatterplot is created. For the latter, a pie map with the posterior membership 
probabilities is created comparable to the ancestry pie maps created in `01-ADMIXTURE`,
as well as a barplot of the posterior membership probabilities with individuals ordered by longitude.

## 05-kinship

## 07-Ne







