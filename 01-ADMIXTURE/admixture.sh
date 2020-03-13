### script to run ADMIXTURE for SNP datasets

# define arguments
spCode=$1
locusType=$2

# e.g. 
# inputFolder=02-Mullus
# spCode=mul
# locusType=adaptive

inputFile=../admx_"$spCode"_"$locusType"_*.bed
numLoci=`basename $inputFile |  awk -F'[_.]' '{print $(NF-1)}'` # retrieves the file name, sets _ and . as separators and retrieves the information in the second last column, in this case the number of SNPs in our dataset


### launch admixture for K varying from 1 to 5, 
### -B speciefies the amount of bootstrapping iterations, here 2000 
### -j speciefies the amount of processors used (here 4)
### -cv= specifies the amount of cross-validation procedures performed (here 10)
### tee copies the cross-validation log/output to a .out file 
for K in 1 2 3 4 5; do admixture --cv=10 -B2000 -j4 $inputFile $K | tee log${K}.out; done

# for K in 1 2 3 4 5; do admixture --cv=5 -B200 $inputFile $K | tee log${K}.out; done # test if works


### Merge the cross validation information obtained from the log files
grep -h CV log*.out>cross_validation.txt

### Take the right order for individual id using the tfam file
cut -f 1 ../admx_"$spCode"_"$locusType"_"$numLoci".fam > id_admixture.txt

### plot cross-validation results in R
mkdir figures

R --vanilla
library(ggplot2)
library(reshape2)
crossval <- read.table("cross_validation.txt")
K <- c(1: nrow(crossval))
crossval_gg <- cbind(crossval, K) 
pdf("figures/cross_validation.pdf", width = 5, height = 5)
ggplot(data = crossval_gg, aes(x = K, y = V4)) +
  geom_line() +
  geom_point() +
  ylab("Cross-validation error")
dev.off()  
q()

### plot ancestry matrix Q for the optimal K
  # see interactive R script ancestry_barplots.R


