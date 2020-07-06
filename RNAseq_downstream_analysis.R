#

###################################################################################################
## Downstream anlaysis
###################################################################################################

load("/athena/masonlab/scratch/users/nai2008/genes.rda")
genes$Gene.Synonym=as.vector(genes$Gene.Synonym)
genes$Gene.Name=as.vector(genes$Gene.Name)

ig=read.csv('/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/imprinted_human_genes.csv')
ig$Gene=as.vector(ig$Gene)
ig$Gene[which(ig$Gene=='TNDM')]='DMTN'
ig$Gene[which(ig$Gene=='INPP5F V2')]='INPP5F'

mm=match(ig$Gene, genes$Gene.Name)
ig_pt1=ig[which(!is.na(mm)),]
ig_pt2=ig[which(is.na(mm)),]

mm=match(ig_pt2$Gene, genes$Gene.Synonym)
ig_pt2A=ig_pt2[which(!is.na(mm)),]
ig_pt2B=ig_pt2[which(is.na(mm)),]

for(i in 1:nrow(ig_pt2A)){
	m=match(ig_pt2A$Gene[i], genes$Gene.Synonym)
	ig_pt2A$Gene[i]=genes$Gene.Name[m]
}

ig=rbind(ig_pt1,ig_pt2A)

#####################################
## is there overlap b/w genes differentially expressed by sex in cases vs controls
#####################################

load('/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/rdas/controls_0to21_DE_autosomal_genes.rda')
controls_DE=out_autosomal

load('/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/rdas/HGG_0to21_DE_autosomal_genes.rda')
HGG_DE=out_autosomal

mm=match(HGG_DE$gene, controls_DE$gene) 

length(which(!is.na(mm))) # 4
HGG_DE$gene[which(!is.na(mm))] # ZNF883 FAR2P2 ITIH5 DEFA3

#####################################
## IGSF3 exloration for Dr. Greenfield
#####################################

load('/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/rdas/cases_0to21_DEseq2_results.rda')
rr_cases=rr

print=rr_cases[grep('IGSF',rownames(rr_cases)),]
write.csv(print[,c(2,6,7)],file="/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/data_tables/IGSF_cases_exploration.csv")

load('/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/rdas/controls_0to21_DEseq2_results.rda')
rr_controls=rr

print=rr_controls[grep('IGSF',rownames(rr_controls)),]
write.csv(print[,c(2,6,7)],file="/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/data_tables/IGSF_controls_exploration.csv")

#####################################
## Exploration of DE genes in the CONTROL dataset
#####################################

## How many DE genes have |log2FC| >2 in the CONTROL dataset?

load('/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/rdas/controls_0to21_DE_autosomal_genes.rda')
controls_DE=out_autosomal

length(which(abs(controls_DE$log2FoldChange)>2)) # 1
controls_DE[which(abs(controls_DE$log2FoldChange)>2),]
#       gene chr log2FoldChange          FDR
# AC022730.4   8       2.248818 7.826913e-06     (a lncRNA)

## How many DE genes are imprinted genes:
mm=match(controls_DE$gene,ig$Gene)
length(which(!is.na(mm)))

#####################################
## Exploration of DE genes in the pHGG dataset
#####################################

## How many DE genes have |log2FC| >2 in the pHGG dataset?
load('/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/rdas/HGG_0to21_DE_autosomal_genes.rda')
HGG_DE=out_autosomal

length(which(abs(HGG_DE$log2FoldChange)>2)) # 29
HGG_DE$gene[which(abs(HGG_DE$log2FoldChange)>2)]
#  [1] MTRNR2L8   TUFMP1     ZNF883     AC017053.1 RPL23AP47  PAX3
#  [7] OBP2A      CCDC140    TMEM190    DLK1       NWD2       HOXB5
# [13] FEZF1-AS1  PRSS33     MMP13      SIX3       SIX3-AS1   GATA6-AS1
# [19] KISS1R     IGF2-AS    ALOX15     MIR217HG   DEFA3      PITX2
# [25] HOXD11     SLC24A5    NKX2-1     BHLHE23    SCN10A

## How many DE genes are imprinted genes:
mm=match(HGG_DE$gene,ig$Gene)

HGG_DE[which(!is.na(mm)),]

# gene chr log2FoldChange          FDR
# DLK1  14      -2.176929 0.0005325927
# IGF2  11       1.977341 0.0556489705
# IGF2-AS  11       2.776508 0.0228121817



