#

###################################################################################################
## Downstream anlaysis
###################################################################################################

#####################################
## is there overlap b/w * autosomal * genes differentially expressed by sex in cases vs controls
#####################################

# HGG samples
load('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/HGG_0to21_DE_genes.rda') #out
nrow(out) # 83
length(unique(out$Ensembl)) # 83
HGG_DE_autosomal=out[! out$chr %in% c('X','Y'),]

# Controls hippo samples
load('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/controls_hippo_0to21_DE_genes.rda') #out
nrow(out) #106
length(unique(out$Ensembl)) #106
controls_hippo_DE_autosomal=out[! out$chr %in% c('X','Y'),]

# Controls dlpfc samples
load('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/controls_dlpfc_0to21_DE_genes.rda') #out
nrow(out) #98
length(unique(out$Ensembl)) #98
controls_dlpfc_DE_autosomal=out[! out$chr %in% c('X','Y'),]

#

mm=match(as.vector(HGG_DE_autosomal$Ensembl), as.vector(controls_hippo_DE_autosomal$Ensembl)) 
length(which(!is.na(mm))) # 0
# HGG_DE_autosomal$gene[which(!is.na(mm))]

mm=match(as.vector(HGG_DE_autosomal$Ensembl), as.vector(controls_dlpfc_DE_autosomal$Ensembl))
length(which(!is.na(mm))) # 0
# HGG_DE_autosomal$gene[which(!is.na(mm))]

#####################################
## is there overlap b/w * non-autosomal * genes differentially expressed by sex in cases vs controls
#####################################

# HGG samples
load('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/HGG_0to21_DE_genes.rda') #out
HGG_DE_xy=out[ out$chr %in% c('X','Y'),]

# Controls hippo samples
load('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/controls_hippo_0to21_DE_genes.rda') #out
controls_hippo_DE_xy=out[ out$chr %in% c('X','Y'),]

# Controls dlpfc samples
load('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/controls_dlpfc_0to21_DE_genes.rda') #out
controls_dlpfc_DE_xy=out[ out$chr %in% c('X','Y'),]

#

mm=match(HGG_DE_xy$Ensembl, controls_hippo_DE_xy$Ensembl) 
length(which(!is.na(mm))) # 32
length(which(!is.na(mm)))/length(HGG_DE_xy$Ensembl) # 0.8
# HGG_DE_xy$gene[which(!is.na(mm))]

mm=match(HGG_DE_xy$Ensembl, controls_dlpfc_DE_xy$Ensembl) 
length(which(!is.na(mm))) # 33
length(which(!is.na(mm)))/length(HGG_DE_xy$Ensembl) # 0.825
# HGG_DE_xy$gene[which(!is.na(mm))]

#####################################
## IGSF3 exloration for Dr. Greenfield
#####################################

#cases
load('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/cases_0to21_DEseq2_results.rda')
rr_cases=rr
print=rr_cases[grep('IGSF',rr_cases$gene, ignore.case = TRUE),]
#write.csv(print[,c(2,6,7,8,9)],file="/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/data_tables/IGSF_cases_exploration.csv")
length(which(print$padj<=0.1)) # 0

#controls (hippocampus)
load('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/controls_hippo_0to21_DEseq2_results.rda')
rr_controls_hippo=rr
print=rr_controls_hippo[grep('IGSF',rr_controls_hippo$gene, ignore.case = TRUE),]
#write.csv(print[,c(2,6,7,8,9)],file="/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/data_tables/IGSF_controls_hippo_exploration.csv")
length(which(print$padj<=0.1)) # 0

#controls (DLPFC)
load('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/controls_dlpfc_0to21_DEseq2_results.rda')
rr_controls_dlpfc=rr
print=rr_controls_dlpfc[grep('IGSF',rr_controls_dlpfc$gene, ignore.case = TRUE),]
#write.csv(print[,c(2,6,7,8,9)],file="/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/data_tables/IGSF_controls_dlpfc_exploration.csv")
length(which(print$padj<=0.1)) # 0

#####################################
## Exploration of autosomal DE genes in the CONTROL HIPPOCAMPUS dataset
#####################################

## How many DE genes have |log2FC| >2 in the CONTROL hippo dataset?
load('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/controls_hippo_0to21_DE_genes.rda') #out
nrow(out) # 106
length(unique(out$Ensembl)) # 106
controls_hippo_DE_autosomal=out[! out$chr %in% c('X','Y'),]

length(which(abs(controls_hippo_DE_autosomal$log2FoldChange)>2)) # 2
controls_hippo_DE_autosomal$gene[which(abs(controls_hippo_DE_autosomal$log2FoldChange)>2)] # CH507-528H12.1 RPS26P6

median(abs(controls_hippo_DE_autosomal$log2FoldChange))
mean(abs(controls_hippo_DE_autosomal$log2FoldChange))
range(abs(controls_hippo_DE_autosomal$log2FoldChange))
range(controls_hippo_DE_autosomal$log2FoldChange)
quantile(controls_hippo_DE_autosomal$log2FoldChange)

# > median(abs(controls_hippo_DE_autosomal$log2FoldChange))
# [1] 0.2709475
# > mean(abs(controls_hippo_DE_autosomal$log2FoldChange))
# [1] 0.478677
# > range(abs(controls_hippo_DE_autosomal$log2FoldChange))
# [1] 0.1280611 2.5067420
# > range(controls_hippo_DE_autosomal$log2FoldChange)
# [1] -0.7136051  2.5067420
# > quantile(controls_hippo_DE_autosomal$log2FoldChange)
#         0%        25%        50%        75%       100%
# -0.7136051 -0.1280611  0.2392128  0.4685019  2.5067420

quantile(abs(controls_hippo_DE_autosomal$log2FoldChange))
#        0%       25%       50%       75%      100%
# 0.1280611 0.1983379 0.2709475 0.5427018 2.5067420
# IQR=0.3443639

#####################################
## Exploration of autosomal DE genes in the CONTROL DLPFC dataset
#####################################

## How many DE genes have |log2FC| >2 in the CONTROL hippo dataset?
load('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/controls_dlpfc_0to21_DE_genes.rda')
nrow(out) # 98
length(unique(out$Ensembl)) # 98
controls_dlpfc_DE_autosomal=out[! out$chr %in% c('X','Y'),]

length(which(abs(controls_dlpfc_DE_autosomal$log2FoldChange)>2)) # 3
controls_dlpfc_DE_autosomal$gene[which(abs(controls_dlpfc_DE_autosomal$log2FoldChange)>2)]
# RP11-174O3.3  RP11-927P21.2 CH507-254M2.3

median(abs(controls_dlpfc_DE_autosomal$log2FoldChange))
mean(abs(controls_dlpfc_DE_autosomal$log2FoldChange))
range(abs(controls_dlpfc_DE_autosomal$log2FoldChange))
range(controls_dlpfc_DE_autosomal$log2FoldChange)
quantile(controls_dlpfc_DE_autosomal$log2FoldChange)

# > median(abs(controls_dlpfc_DE_autosomal$log2FoldChange))
# [1] 0.5870916
# > mean(abs(controls_dlpfc_DE_autosomal$log2FoldChange))
# [1] 0.8112447
# > range(abs(controls_dlpfc_DE_autosomal$log2FoldChange))
# [1] 0.09114172 2.79181485
# > range(controls_dlpfc_DE_autosomal$log2FoldChange)
# [1] -2.791815  2.528285
# > quantile(controls_dlpfc_DE_autosomal$log2FoldChange)
#         0%        25%        50%        75%       100%
# -2.7918149 -0.4467670  0.2638402  0.8674927  2.5282847

quantile(abs(controls_dlpfc_DE_autosomal$log2FoldChange))
#         0%        25%        50%        75%       100%
# 0.09114172 0.39262695 0.58709161 1.02839634 2.79181485
# IQR=0.6357694

#####################################
## Exploration of autosomal DE genes in the pHGG dataset
#####################################

## How many DE genes have |log2FC| >2 in the pHGG dataset?
load('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/HGG_0to21_DE_genes.rda') #out
nrow(out) #83
length(unique(out$Ensembl)) #83
HGG_DE_autosomal=out[! out$chr %in% c('X','Y'),]

length(which(abs(HGG_DE_autosomal$log2FoldChange)>2)) # 26
HGG_DE_autosomal$gene[which(abs(HGG_DE_autosomal$log2FoldChange)>2)]
#  [1] LDHBP3        MTRNR2L8      WIF1          TMEM233       SLC18A3
#  [6] BNC1          DCT           RP11-128M1.1  TP73          DHH
# [11] RP11-903H12.3 SELE          ADGRF2        AANAT         TFAP2A-AS1
# [16] CTD-2207P18.1 RP11-114H24.2 RP5-936J12.1  DDX43         SLC24A5
# [21] CBLN4         MUC19         PRL           HOXA9         TBC1D3H
# [26] TBC1D3G

median(abs(HGG_DE_autosomal$log2FoldChange))
mean(abs(HGG_DE_autosomal$log2FoldChange))
range(abs(HGG_DE_autosomal$log2FoldChange))
range(HGG_DE_autosomal$log2FoldChange)
quantile(HGG_DE_autosomal$log2FoldChange)

# > median(abs(HGG_DE_autosomal$log2FoldChange))
# [1] 2.146255
# > mean(abs(HGG_DE_autosomal$log2FoldChange))
# [1] 3.012383
# > range(abs(HGG_DE_autosomal$log2FoldChange))
# [1]  0.3581242 25.5012303
# > range(HGG_DE_autosomal$log2FoldChange)
# [1] -10.46916  25.50123
# > quantile(HGG_DE_autosomal$log2FoldChange)
#         0%        25%        50%        75%       100%
# -10.469159  -1.400537   1.113722   2.302499  25.501230

quantile(abs(HGG_DE_autosomal$log2FoldChange))
#         0%        25%        50%        75%       100%
#  0.3581242  1.1968844  2.1462547  2.7765965 25.5012303
# IQR=1.579712


