#
# find blocks

#####################################
## <= 6 yo (cases)
#####################################
library(minfi)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")

## remove samples with a a mutant form of IDH

GRset=GRset[,-grep('IDH|MUT|Mutation', pData(GRset)$IDH_status)] # 36/707 samples have an IDH mutation

## find blocks
GRset_prepubescent_cases=GRset[,which(pData(GRset)$age <=6 & pData(GRset)$Dx_simplified=='HGG')]
input_GRset=GRset_prepubescent_cases
blocks_filename="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/blocks_prepubescent_cases_500bootstraps.rda"
B=500 # number permutations

blocks_v2(input_GRset=input_GRset, blocks_filename=blocks_filename, B=B)

#####################################
## 6-21 (cases)
#####################################
library(minfi)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")

## remove samples with a a mutant form of IDH

GRset=GRset[,-grep('IDH|MUT|Mutation', pData(GRset)$IDH_status)] # 36/707 samples have an IDH mutation

## find blocks
GRset_6to21_cases=GRset[,which(pData(GRset)$age >6 & pData(GRset)$age<=21 & pData(GRset)$Dx_simplified=='HGG')]
input_GRset=GRset_6to21_cases
blocks_filename="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/blocks_6to21_cases_500bootstraps.rda"
B=500 # number permutations

blocks_v2(input_GRset=input_GRset, blocks_filename=blocks_filename, B=B)

#####################################
## >21 (cases)
#####################################
library(minfi)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")

## remove samples with a a mutant form of IDH

GRset=GRset[,-grep('IDH|MUT|Mutation', pData(GRset)$IDH_status)] # 36/707 samples have an IDH mutation

## find blocks
GRset_adult_cases=GRset[,which(pData(GRset)$age >21 & pData(GRset)$Dx_simplified=='HGG')]
input_GRset=GRset_adult_cases
blocks_filename="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/blocks_adult_cases_500bootstraps.rda"
B=500 # number permutations

blocks_v2(input_GRset=input_GRset, blocks_filename=blocks_filename, B=B)

#####################################
## 0-21 (cases)
#####################################
rm(list=ls())

library(minfi)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")

## remove samples with a a mutant form of IDH

GRset=GRset[,-grep('IDH|MUT|Mutation', pData(GRset)$IDH_status)] # 36/707 samples have an IDH mutation

## find blocks
GRset_allPeds_cases=GRset[,which(pData(GRset)$age <=21 & pData(GRset)$Dx_simplified=='HGG')]
input_GRset=GRset_allPeds_cases
blocks_filename="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/blocks_0to21_cases_500bootstraps.rda"
B=500 # number permutations

blocks_v2(input_GRset=input_GRset, blocks_filename=blocks_filename, B=B)

#####################################
## <= 6 yo (controls)
#####################################
library(minfi)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")

## remove samples with a a mutant form of IDH

GRset=GRset[,-grep('IDH|MUT|Mutation', pData(GRset)$IDH_status)] # 36/707 samples have an IDH mutation

## find blocks
GRset_prepubescent_controls=GRset[,which(pData(GRset)$age <=6 & pData(GRset)$Dx_simplified=='Control')]
input_GRset=GRset_prepubescent_controls
blocks_filename="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/blocks_prepubescent_controls_500bootstraps.rda"
B=500 # number permutations

blocks_v2(input_GRset=input_GRset, blocks_filename=blocks_filename, B=B)

#####################################
## 6-21 (controls)
#####################################
library(minfi)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")

## remove samples with a a mutant form of IDH

GRset=GRset[,-grep('IDH|MUT|Mutation', pData(GRset)$IDH_status)] # 36/707 samples have an IDH mutation

## find blocks
GRset_6to21_controls=GRset[,which(pData(GRset)$age >6 & pData(GRset)$age<=21 & pData(GRset)$Dx_simplified=='Control')]
input_GRset=GRset_6to21_controls
blocks_filename="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/blocks_6to21_controls_500bootstraps.rda"
B=500 # number permutations

blocks_v2(input_GRset=input_GRset, blocks_filename=blocks_filename, B=B)

#####################################
## >21 (controls)
#####################################
library(minfi)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")

## remove samples with a a mutant form of IDH

GRset=GRset[,-grep('IDH|MUT|Mutation', pData(GRset)$IDH_status)] # 36/707 samples have an IDH mutation

## find blocks
GRset_adult_controls=GRset[,which(pData(GRset)$age >21 & pData(GRset)$Dx_simplified=='Control')]
input_GRset=GRset_adult_controls
blocks_filename="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/blocks_adult_controls_500bootstraps.rda"
B=500 # number permutations

blocks_v2(input_GRset=input_GRset, blocks_filename=blocks_filename, B=B)

#####################################
## 0-21 (controls)
#####################################
rm(list=ls())

library(minfi)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")

## remove samples with a a mutant form of IDH

GRset=GRset[,-grep('IDH|MUT|Mutation', pData(GRset)$IDH_status)] # 36/707 samples have an IDH mutation

## find blocks
GRset_allPeds_controls=GRset[,which(pData(GRset)$age <=21 & pData(GRset)$Dx_simplified=='Control')]
input_GRset=GRset_allPeds_controls
blocks_filename="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/blocks_0to21_controls_500bootstraps.rda"
B=500 # number permutations

blocks_v2(input_GRset=input_GRset, blocks_filename=blocks_filename, B=B)

























