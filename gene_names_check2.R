#

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/gene_universe_450k.rda") # universe.450k

load('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/cases_0to21_DEseq2_results.rda')
sig_results=as.data.frame(rr[which(rr$padj<=0.1),])
grrs=sig_results$Ensembl

mm=match(grrs, universe.450k)
length(which(is.na(mm))) # 45

##

rm(list=ls())

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/gene_universe_450k.rda") # universe.450k

load('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/controls_hippo_0to21_DEseq2_results.rda')
sig_results=as.data.frame(rr[which(rr$padj<=0.1),])
grrs=sig_results$Ensembl

mm=match(grrs, universe.450k)
length(which(is.na(mm))) # 62

##

rm(list=ls())

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/gene_universe_450k.rda") # universe.450k

load('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/controls_dlpfc_0to21_DEseq2_results.rda')
sig_results=as.data.frame(rr[which(rr$padj<=0.1),])
grrs=sig_results$Ensembl

mm=match(grrs, universe.450k)
length(which(is.na(mm))) # 63

