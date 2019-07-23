# seaConnect--radFishComp
Comparison of fish RAD sequence data: compare different species and different marker types for population/fisheries genomic

## Data
RADseq SNP calling was performed by [P. Guerrin using a dDocent-based pipeline](https://github.com/Grelot/seaConnect--dDocent)  
SNPs were then filtered for missing data, maf, ... and outliers were detected to retain 
only adaptive or neutral loci for both species. By [Emilie Boulanger](https://github.com/eboulanger/seaConnect--dataPrep)
Finally, datasets are converted to the correct formats for different analyses.  
Move the necessary files to your radFishComp data folder
```
cp ../seaConnect--dataPrep/04-finalPrep/01-Diplodus/*.bed 00-Data/01-Diplodus/
cp ../seaConnect--dataPrep/04-finalPrep/01-Diplodus/*.bim 00-Data/01-Diplodus/
```

## Dependencies
You will need to install the following software:  
- [Admixture](http://software.genetics.ucla.edu/admixture/download.html)
And you will need to install the following R libraries:  
- ggplot2
- reshape2
- pophelper
- maps
- mapdata
- dplyr
- plotrix


## 01-ADMIXTURE

Script to run admixture analysis on adaptive and neutral SNPs to detect population structure.  
After running the model for different values of K, it also returns a cross-validation plot made in R.

Define arguments:  
  $1 = inputFolder
  $2 = spCode
  $3 = locusType

```
bash admixture.sh 01-Diplodus dip adaptive
bash admixture.sh 01-Diplodus dip neutral

bash admixture.sh 02-Mullus mul adaptive
bash admixture.sh 02-Mullus mul neutral
```

