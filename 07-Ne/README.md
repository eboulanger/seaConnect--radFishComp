Calculate effective populations size for D. sargus and M. surmuletus.

Use the populations detected with neutral markers by DAPC. 

Mullus: one main population in the whole med sea and a second population comprising of 
        Alboran sea individuals and scattered individuals

Diplodus: one main population in the whole med sea and a second population comprising of 
          Alboran sea and south-western individuals
          
Calculate Ne for both populations for each species. 

# step 1: 
create genepop file with DAPC populations using PGDSpider. input: neutral .vcf and 
population map obtained through DAPC. output: .gen.txt file. add info on first line and 
remove .txt extension so NeEstimator can read it. 

# step 2:  
calculate Ne using Ne estimator.


Diplodus: 
```
input file name: Documents/project_SEACONNECT/seaConnect--radFishComp/07-Ne/01_create_genepop/dip_neutral_7570_DapcPopK2.gen
output file name: Documents/project_SEACONNECT/seaConnect--radFishComp/07-Ne/02_estimate_Ne/dip_neutral_7570_DapcPopK2_NeLD.txt
```

Mullus:
```
input file name: Documents/project_SEACONNECT/seaConnect--radFishComp/07-Ne/01_create_genepop/mul_neutral_13548_DapcPopK2.gen
output file name: Documents/project_SEACONNECT/seaConnect--radFishComp/07-Ne/02_estimate_Ne/mul_neutral_13548_NeLD.txt
```
