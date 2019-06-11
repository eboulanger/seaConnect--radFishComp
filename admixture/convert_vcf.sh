
# convert vcf file to plink .bed : two steps
# .vcf to .tped
vcftools --vcf mullus.vcf --plink-tped --out mullus

# .tped to .bed
plink --tped mullus.tped --tfam mullus.tfam --make-bed --out mullus --noweb

# merci Pierre-Edouard
