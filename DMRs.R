#
# find DMRs

#####################################
## Prepubescent (cases)
#####################################
rm(list=ls())

library(minfi)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")

## remove samples with a a mutant form of IDH
GRset=GRset[,-grep('IDH|MUT|Mutation', pData(GRset)$IDH_status)]

## find DMRs
GRset_cases=GRset[,which(pData(GRset)$Dx_simplified=='HGG')]
indexes=c(which(pData(GRset_cases)$age < 9.5 & pData(GRset_cases)$sex=='Female'), 
	which(pData(GRset_cases)$age < 10.5 & pData(GRset_cases)$sex=='Male'))
GRset_prepubescent_cases=GRset_cases[,indexes]

input_GRset=GRset_prepubescent_cases

DMRs_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/DMRs_prepubescent_cases_1000bootstraps.rda"
B=1000 # number permutations

DMRs_v2(input_GRset=input_GRset, DMRs_filename=DMRs_filename, B=B)

#####################################
## Postpubescent (cases)
#####################################
rm(list=ls())

library(minfi)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")

## remove samples with a a mutant form of IDH
GRset=GRset[,-grep('IDH|MUT|Mutation', pData(GRset)$IDH_status)]

## find DMRs
GRset_cases=GRset[,which(pData(GRset)$Dx_simplified=='HGG')]
indexes=c(which(pData(GRset_cases)$age >= 9.5 & pData(GRset_cases)$age <= 21 & pData(GRset_cases)$sex=='Female'), 
	which(pData(GRset_cases)$age >= 10.5 & pData(GRset_cases)$age <= 21 & pData(GRset_cases)$sex=='Male'))
GRset_postpubescent_cases=GRset_cases[,indexes]

input_GRset=GRset_postpubescent_cases

DMRs_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/DMRs_postpubescent_cases_1000bootstraps.rda"
B=1000 # number permutations

DMRs_v2(input_GRset=input_GRset, DMRs_filename=DMRs_filename, B=B)

#####################################
## >21 (cases)
#####################################
rm(list=ls())

library(minfi)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")

## remove samples with a a mutant form of IDH
GRset=GRset[,-grep('IDH|MUT|Mutation', pData(GRset)$IDH_status)]

## find DMRs
GRset_cases=GRset[,which(pData(GRset)$Dx_simplified=='HGG')]
GRset_adult_cases=GRset_cases[,which(pData(GRset_cases)$age >21)]

input_GRset=GRset_adult_cases

DMRs_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/DMRs_adult_cases_1000bootstraps.rda"
B=1000 # number permutations

DMRs_v2(input_GRset=input_GRset, DMRs_filename=DMRs_filename, B=B)

#####################################
## 0-21 (cases)
#####################################
rm(list=ls())

library(minfi)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")

## remove samples with a a mutant form of IDH
GRset=GRset[,-grep('IDH|MUT|Mutation', pData(GRset)$IDH_status)]

## find DMRs
GRset_cases=GRset[,which(pData(GRset)$Dx_simplified=='HGG')]
GRset_allPeds_cases=GRset_cases[,which(pData(GRset_cases)$age <=21)]

input_GRset=GRset_allPeds_cases

DMRs_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/DMRs_allPeds_cases_1000bootstraps.rda"
B=1000 # number permutations

DMRs_v2(input_GRset=input_GRset, DMRs_filename=DMRs_filename, B=B)

#####################################
## Prepubescent (controls)
#####################################
rm(list=ls())

library(minfi)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")

## remove samples with a a mutant form of IDH
GRset=GRset[,-grep('IDH|MUT|Mutation', pData(GRset)$IDH_status)]

## find DMRs
GRset_controls=GRset[,which(pData(GRset)$Dx_simplified=='Control')]
indexes=c(which(pData(GRset_controls)$age < 9.5 & pData(GRset_controls)$sex=='Female'), 
	which(pData(GRset_controls)$age < 10.5 & pData(GRset_controls)$sex=='Male'))
GRset_prepubescent_controls=GRset_controls[,indexes]

input_GRset=GRset_prepubescent_controls

DMRs_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/DMRs_prepubescent_controls_1000bootstraps.rda"
B=1000 # number permutations

DMRs_v2(input_GRset=input_GRset, DMRs_filename=DMRs_filename, B=B)

#####################################
## Postpubescent (controls)
#####################################
rm(list=ls())

library(minfi)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")

## remove samples with a a mutant form of IDH
GRset=GRset[,-grep('IDH|MUT|Mutation', pData(GRset)$IDH_status)]

## find DMRs
GRset_controls=GRset[,which(pData(GRset)$Dx_simplified=='Control')]
indexes=c(which(pData(GRset_controls)$age >= 9.5 & pData(GRset_controls)$age <= 21 & pData(GRset_controls)$sex=='Female'), 
	which(pData(GRset_controls)$age >= 10.5 & pData(GRset_controls)$age <= 21 & pData(GRset_controls)$sex=='Male'))
GRset_postpubescent_controls=GRset_controls[,indexes]

input_GRset=GRset_postpubescent_controls

DMRs_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/DMRs_postpubescent_controls_1000bootstraps.rda"
B=1000 # number permutations

DMRs_v2(input_GRset=input_GRset, DMRs_filename=DMRs_filename, B=B)

#####################################
## >21 (controls)
#####################################
rm(list=ls())

library(minfi)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")

## remove samples with a a mutant form of IDH
GRset=GRset[,-grep('IDH|MUT|Mutation', pData(GRset)$IDH_status)]

## find DMRs
GRset_controls=GRset[,which(pData(GRset)$Dx_simplified=='Control')]
GRset_adult_controls=GRset_controls[,which(pData(GRset_controls)$age >21)]

input_GRset=GRset_adult_controls

DMRs_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/DMRs_adult_controls_1000bootstraps.rda"
B=1000 # number permutations

DMRs_v2(input_GRset=input_GRset, DMRs_filename=DMRs_filename, B=B)

#####################################
## 0-21 (controls)
#####################################
rm(list=ls())

library(minfi)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")

## remove samples with a a mutant form of IDH
GRset=GRset[,-grep('IDH|MUT|Mutation', pData(GRset)$IDH_status)]

## find DMRs
GRset_controls=GRset[,which(pData(GRset)$Dx_simplified=='Control')]
GRset_allPeds_controls=GRset_controls[,which(pData(GRset_controls)$age <=21)]

input_GRset=GRset_allPeds_controls

DMRs_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/DMRs_allPeds_controls_1000bootstraps.rda"
B=1000 # number permutations

DMRs_v2(input_GRset=input_GRset, DMRs_filename=DMRs_filename, B=B)

























