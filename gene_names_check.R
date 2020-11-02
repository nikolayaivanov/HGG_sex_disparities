load('/athena/masonlab/scratch/users/nai2008/gencode_v35_hg19_genes_and_genomicState.rda')
genes1=as.vector(genes$Geneid)

library(EnsDb.Hsapiens.v86)
edb = EnsDb.Hsapiens.v86
genes2=as.data.frame(genes(edb))$gene_id

mm=match(genes1,genes2)
length(which(is.na(mm))) # 6653

mm=match(genes2,genes1)
length(which(is.na(mm))) # 8228

##

rm(list=ls())

load('/athena/masonlab/scratch/users/nai2008/gencode_v35_hg19_genes_and_genomicState.rda')
genes=as.vector(genes$Geneid)

load('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/cases_0to21_DEseq2_results.rda')
sig_results=as.data.frame(rr[which(rr$padj<=0.1),])
grrs=sig_results$Ensembl

mm=match(grrs, genes)
length(which(is.na(mm))) # 3
sig_results[which(is.na(mm)),]

#                   baseMean log2FoldChange     lfcSE      stat       pvalue
# ENSG00000233864 645.361507       4.679592 0.2218268 21.095696 8.711428e-99
# ENSG00000260287   4.166453      25.501230 3.2537540  7.837479 4.596804e-15
# ENSG00000276399 127.239933       1.015777 0.2464106  4.122295 3.751167e-05
#                         padj chr         Ensembl     gene
# ENSG00000233864 3.145406e-95   Y ENSG00000233864   TTTY15
# ENSG00000260287 5.745298e-12  17 ENSG00000260287  TBC1D3G
# ENSG00000276399 2.257369e-02  17 ENSG00000276399 FLJ36000

##

rm(list=ls())

load('/athena/masonlab/scratch/users/nai2008/gencode_v35_hg19_genes_and_genomicState.rda')
genes=as.vector(genes$Geneid)

load('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/controls_hippo_0to21_DEseq2_results.rda')
sig_results=as.data.frame(rr[which(rr$padj<=0.1),])
grrs=sig_results$Ensembl

mm=match(grrs, genes)
length(which(is.na(mm))) # 7
sig_results[which(is.na(mm)),]

#                     baseMean log2FoldChange      lfcSE      stat        pvalue
# ENSG00000262621    77.612034     -0.2955301 0.07091361 -4.167466  3.080039e-05
# ENSG00000282826   330.170218      0.4305078 0.10748549  4.005265  6.194805e-05
# ENSG00000283443    25.991238      1.1940792 0.31054765  3.845076  1.205152e-04
# ENSG00000280441     4.617863      2.0498976 0.56662187  3.617752  2.971724e-04
# ENSG00000205664   119.875432     -0.6012117 0.10184217 -5.903367  3.561568e-09
# ENSG00000233864   388.311285      6.5346565 0.30331551 21.544089 6.015466e-103
# ENSG00000198712 59702.782256     -0.3708088 0.09814487 -3.778178  1.579800e-04
#                         padj chr         Ensembl           gene
# ENSG00000262621 1.235047e-02  16 ENSG00000262621  LA16c-306E5.2
# ENSG00000282826 2.408747e-02  20 ENSG00000282826         FRG1CP
# ENSG00000283443 4.069450e-02  20 ENSG00000283443   RP11-462H3.2
# ENSG00000280441 7.882315e-02  21 ENSG00000280441 CH507-528H12.1
# ENSG00000205664 1.904178e-06   X ENSG00000205664  RP11-706O15.1
# ENSG00000233864 1.102678e-99   Y ENSG00000233864         TTTY15
# ENSG00000198712 5.005236e-02  MT ENSG00000198712         MT-CO2

##

rm(list=ls())

load('/athena/masonlab/scratch/users/nai2008/gencode_v35_hg19_genes_and_genomicState.rda')
genes=as.vector(genes$Geneid)

load('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/controls_dlpfc_0to21_DEseq2_results.rda')
sig_results=as.data.frame(rr[which(rr$padj<=0.1),])
grrs=sig_results$Ensembl

mm=match(grrs, genes)
length(which(is.na(mm))) # 4
sig_results[which(is.na(mm)),]

#                   baseMean log2FoldChange      lfcSE      stat        pvalue
# ENSG00000283443  30.200086       1.402621 0.30095943  4.660498  3.154449e-06
# ENSG00000280164   3.158777       2.528285 0.69772205  3.623627  2.905001e-04
# ENSG00000205664 171.444692      -0.425392 0.09385985 -4.532204  5.837143e-06
# ENSG00000233864 544.050584       8.196597 0.23953726 34.218463 1.285002e-256
#                          padj chr         Ensembl          gene
# ENSG00000283443  1.495611e-03  20 ENSG00000283443  RP11-462H3.2
# ENSG00000280164  8.058906e-02  21 ENSG00000280164 CH507-254M2.3
# ENSG00000205664  2.642318e-03   X ENSG00000205664 RP11-706O15.1
# ENSG00000233864 4.787000e-253   Y ENSG00000233864        TTTY15


