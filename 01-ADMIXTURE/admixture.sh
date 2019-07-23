
# define arguments
inputFolder=$1
spCode=$2
locusType=$3

# e.g. 
# inputFolder=02-Mullus
# spCode=mul
# locusType=adaptive

inputFile=../../../../seaConnect--dataPrep/04-finalPrep/"$inputFolder"/"$spCode"_"$locusType"_*.bed
numLoci=`basename $inputFile |  awk -F'[_.]' '{print $(NF-1)}'` # retrieves the file name, sets _ and . as separators and retrieves the information in the second last column


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
cut -f 1 ../../../../seaConnect--dataPrep/04-finalPrep/"$inputFolder"/"$spCode"_"$locusType"_"$numLoci".tfam > id_admixture.txt

### plot cross-validation results in R
R --vanilla
library(ggplot2)
library(reshape2)
crossval <- read.table("cross_validation.txt")
K <- c(1: nrow(crossval))
crossval_gg <- cbind(crossval, K) 
pdf("cross_validation.pdf", width = 5, height = 5)
ggplot(data = crossval_gg, aes(x = K, y = V4)) +
  geom_line() +
  geom_point() +
  ylab("Cross-validation error")
dev.off()  
q()

### plot ancestry matrix Q for the optimal K
  # see interactive R script ancestry_barplots.R


##### re-run ADMIXTURE on subset after first iteration #####

mkdir ../c-adaptive-iter2

### extract individual/sample list by ancestry fraction
### extract individuals that belong for more than 75% to cluster 1
R --vanilla
ancestry <- read.table("dip_adaptive_494.2.Q", col.names = c("Cluster1","Cluster2"))
id <- read.table("id_admixture.txt", col.names = "ID",stringsAsFactors = F)
ancestry <- cbind(id, ancestry)
c1id <- ancestry$ID[ancestry$Cluster1 > 0.75]
write.table(c1id, file = "../c-adaptive-iter2/cl1_ind.txt", quote = F, row.names = F, col.names = F)
q()

### tab-index vcf file for use with bcftools
cd ../c-adaptive-iter2
cp ../../../../seaConnect--dataPrep/04-finalPrep/01-Diplodus/dip_adaptive_*.vcf 
bgzip dip_adaptive_*.vcf # zips your file to sar_adaptive.recode.vcf.gz
tabix -p vcf dip_adaptive_*.vcf.gz

### subset the adaptive vcf file using the cluster 1 individuals list
bcftools view -S cl1_ind.txt -o dip_adaptive_494_iter2.vcf dip_adaptive_494.vcf.gz

#### now re-convert iter2 .vcf
# .vcf to .tped
vcftools --vcf dip_adaptive_494_iter2.vcf --plink-tped --out dip_adaptive_494_iter2
# .tped to .bed
plink --tped dip_adaptive_494_iter2.tped --tfam dip_adaptive_494_iter2.tfam --make-bed --out dip_adaptive_494_iter2 --noweb

### run ADMIXTURE
for K in 1 2 3 4 5; do admixture --cv=10 -B2000 -j4 dip_adaptive_494_iter2.bed $K | tee log${K}.out; done
