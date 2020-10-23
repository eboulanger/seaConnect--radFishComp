Objective:

- calculate Kinship between individuals
- explore kinship patterns and compare them between species

individual-pairs with high kinship (>0.1) :
  what are the geographical distances? different patterns between species?   			

## ANALYSES	

# kinship

calculate with genodive. input format needed: genepop
use PGDspider to convert .vcf to genepop (.gen.txt) using a population map
population map = df with column 1 = ind identifiers, column 2 = "population" identifier 
                                                                 (e.g. sampling cell)
In .gen.txt file: add some text on the first line, otherwise genodive won't be able to read it.
                                                                 
 # kinship patterns
 
 training_kinship.R
 
 # kinship SG
 
 training_kinSG.R
                                                                 