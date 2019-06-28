### move ADMIXTURE executable to bin so you don't have to copy it to the working directory of each analysis
   # cp ~/programmes/admixture_macosx-1.3.0/admixture /usr/local/bin/
   # note: only do this once 

# move to working directory
# this is where you want your analysis output to be sotred
cd ~/Documents/project_SEACONNECT/seaConnect--radFishComp/01-ADMIXTURE/01-Diplodus/c-adaptive/

### Add the path of the location of your .bed file in your terminal by typing:
### (so you don't have to copy - paste all files in working directory)
file="../../../00-Data/01-Diplodus/c-adaptive/sar_adaptive.bed"

### launch admixture for K varying from 1 to 5, 
### -B speciefies the amount of bootstrapping iterations, here 2000 
### -j speciefies the amount of processors used (here 4)
### -cv= specifies the amount of cross-validation procedures performed (here 10)
### tee copies the cross-validation log/output to a .out file 
for K in 1 2 3 4 5; do admixture --cv=10 -B2000 -j4 $file $K | tee log${K}.out; done

for K in 1 2 3 4 5; do admixture --cv=5 -B200 $file $K | tee log${K}.out; done


### Merge the cross validation information obtained from the log files
grep -h CV log*.out>cross_validation.txt

### Take the right order for individual id using the tfam file
cut -f 1 ../../../00-Data/01-Diplodus/c-adaptive/sar_adaptive.tfam > id_admixture.txt

### plot cross-validation results in R
R
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
  # see R script ancestry_barplots.R


##### re-run ADMIXTURE on subset after first iteration #####

### extract individual/sample list by ancestry fraction
R
ancestry <- read.table("sar_adaptive.2.Q", col.names = c("Cluster1","Cluster2"))
id <- read.table("id_admixture.txt", col.names = "ID",stringsAsFactors = F)
ancestry <- cbind(id, ancestry)
c1id <- ancestry$ID[ancestry$Cluster1 > 0.75]
write.table(c1id, file = "../c-adaptive-iter2/cl1_ind.txt", quote = F, row.names = F, col.names = F)
q()

### tab-index vcf file for use with bcftools
bgzip ../../../00-Data/01-Diplodus/c-adaptive/sar_adaptive.recode.vcf # zips your file to sar_adaptive.recode.vcf.gz
tabix -p vcf ../../../00-Data/01-Diplodus/c-adaptive/sar_adaptive.recode.vcf.gz

bcftools view -S ../c-adaptive-iter2/cl1_ind.txt -o ../c-adaptive-iter2/dip_adaptive_iter2.vcf ../../../00-Data/01-Diplodus/c-adaptive/sar_adaptive.recode.vcf.gz

# now re-convert iter2 .vcf and re-run admixture