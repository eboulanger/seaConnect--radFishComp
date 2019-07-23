### rename individual ID of vcf file to avoid downstream complications
### (names should never begin with a number)

### go to directory
cd ~/Documents/project_SEACONNECT/seaConnect--radFishComp/00-Data/new-Mullus

### get vcf file 
cp ../../../seaConnect--dataPrep/04-finalPrep/02-Mullus/mul_adaptive_*.vcf .

### extract sample IDs from vcf file
bcftools query -l mul_adaptive_*.vcf > mullus_id.txt

### in R, add "C" at beginning of each id
R --vanilla
id <- read.table("mullus_id.txt", stringsAsFactors=F)
id[1,1]
# currently our Mullus sample code is e.g.  "100_10E2"
# in sar_ddocent it is formatted as  "C100_D100i10E2"
# try to format Mullus id's in the same way to make admixture work
id1 <- gsub("_.*","", id$V1)
id2 <- gsub(".*_","", id$V1)
id3 <- gsub("E.*","", id2)
Cid <- paste("C", id1, "_", "D", id1, "i", id2, sep ="")
Cid[1]
newid <- paste0("C",id1,"i",id3)
write.table(newid, file = "new_mullus_id.txt", quote = F, row.names = F, col.names = F)
q()

### verify new names
paste mullus_id.txt new_mullus_id.txt 

### replace old ID with new ID and output to .vcf
bcftools reheader --samples new_mullus_id.txt -o mul_adaptive_newid.vcf mul_adaptive_2680.vcf

# re-extract sample IDs to verify
bcftools query -l mul_adaptive_newid.vcf 

### subset vcf to test if bootstrapping works with less individuals and less markers

# extract list of 10 individuals 
sed -n 1,10p new_mullus_id.txt > subset_mullus_id.txt 
# tab-index vcf file to use with bcftools
bgzip mul_adaptive_newid.vcf  # zips your file to .gz
tabix -p vcf mul_adaptive_newid.vcf.gz
# subset individuals from vcf
bcftools view -S subset_mullus_id.txt -o mul_adaptive_subset.vcf mul_adaptive_newid.vcf.gz

# extract list of 10 loci
cp ../../../seaConnect--dataPrep/04-finalPrep/02-Mullus/mul_outl_pos_2680.txt .
sed -n 1,10p mul_outl_pos_2680.txt > subset_mullus_loci.txt
# subset loci from vcf
vcftools --vcf mul_adaptive_subset.vcf --positions subset_mullus_loci.txt --recode --out mul_adaptive_subset
mv mul_adaptive_subset.recode.vcf mul_adaptive_subset.vcf

### convert for admixture formats
# .vcf to .tped
vcftools --vcf mul_adaptive_subset.vcf --plink-tped --out mul_adaptive_subset
# .tped to .bed
plink --tped mul_adaptive_subset.tped --tfam mul_adaptive_subset.tfam --make-bed --out mul_adaptive_subset --noweb

### run ADMIXTURE
for K in 1 2 3 4 5; do admixture --cv=5 -B20 -j4 mul_adaptive_subset.bed $K | tee log${K}.out; done
