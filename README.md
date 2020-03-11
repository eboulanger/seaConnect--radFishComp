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
- [Admixture](http://software.genetics.ucla.edu/admixture/download.html)  
- [Genodive] (http://bentleydrummer.nl/software/software/GenoDive.html) 

And you will need to have the following R libraries:  
ggplot2, reshape2, pophelper, maps, mapdata, dplyr, plotrix, adegenet, ...


## 01-ADMIXTURE

Script to run admixture analysis on adaptive and neutral SNPs to detect population structure.  
After running the model for different values of K, it also returns a cross-validation plot made in R.

Before starting, create subdirectories for each species and type of marker.  
For example: `01-ADMIXTURE/01-Diplodus/c-adaptive/` and run your admixture analyses within these directories.
(Do the same for `02-Mullus` and `b-neutral`.)  

ADMIXTURE needs a dataset in plink format, all in the same directory: `.bed`, `.bim` and `.fam`  
However, ADMIXTURE does not accept chromosome names that are not human chromosomes.
To solve this, we change all chromosome names to zero in the `.bim` file:
- copy the plink files from the project directory `seaConnect--dataPrep/04-finalPrep`.
- rename the first column of the `.bim`file

```
cd 01-ADMIXTURE/
cp ../seaConnect--dataPrep/04-finalPrep/02-Mullus/mul_*.bed 02-Mullus/admx_mul_*.bed
cp ../seaConnect--dataPrep/04-finalPrep/02-Mullus/mul_*.bim 02-Mullus/admx_mul_*.bim
cp ../seaConnect--dataPrep/04-finalPrep/02-Mullus/mul_*.fam 02-Mullus/admx_mul_*.fam

awk '{$1=0;print $0}' admx_mul_*.bim > admx_mul_*.bim.tmp
mv admx_mul_*.bim.tmp admx_mul_*.bim
```
To run the bash admixture script: 

copy admixture.sh in the directory you want to work in (e.g. 01-Diplodus/c-adaptive/)
run the script with the following arguments: 

  $1 = inputFolder, species-specific folder within `seaConnect--dataPrep/04-finalPrep`  
  $2 = species code  
  $3 = marker type (neutral or adaptive)  

```
bash admixture.sh 01-Diplodus dip adaptive
bash admixture.sh 01-Diplodus dip neutral

bash admixture.sh 02-Mullus mul adaptive
bash admixture.sh 02-Mullus mul neutral
```

When the most fitting number of populations is detected, individual-based and sampling site-based
results (ancestry matrix Q) can be visualised with ancestry barplots and ancestry pie maps respectively.  
The code for this can be found in the interactive R-script called `ancestry_barplots.R`.  

Once you have your first results, you can go further and subsample individuals to do a second iteration,
for example all individuals belonging for more than 75% to cluster 1, to look for further sub-structuring.  

To do this, create a new sub-directory within your species directory, e.g. `01-Diplodus/c-adaptive-iter2` 
and run the following interactive bash script:

```
admixture_subset.sh
```

This script will 
- open R to extract the list of cluster 1 individuals
- subset your vcf file with the just-created individual list
- convert the subsetted vcf file to the correct formats
- run admixture for K = 1-5 and 2000 bootstraps for cross-validation  

Results can again be visualised with the R-script `ancestry_barplots.R`

## 02-DAPC

Infer population structure with DAPC from [adegenet](http://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf).

Interactive R-script `run_dapc.R`  
This script runs two types of DAPC:  
- with the Mediterranean ecoregions set as prior groups
- with the detection of best-fit number of populations  

For both a scatterplot is created. For the latter, a pie map with the posterior membership 
probabilities is created comparable to the ancestry pie maps created in `01-ADMIXTURE`,
as well as a barplot of the posterior membership probabilities with individuals ordered by longitude.







