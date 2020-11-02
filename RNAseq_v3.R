#

###################################################################################################
## Cavatica data (HGG cases 0-21 y.o.)
###################################################################################################
library(DESeq2)

# import data
library(tximport)

library(plyr)

kallistoFiles = list.files("/athena/masonlab/scratch/users/cem2009/home_backup/Projects/PediatricMeningioma/data/RNAseq_counts", "kallisto.abundance.tsv.gz", recursive=T, full.names=T)
metadata = read.csv("/athena/masonlab/scratch/users/cem2009/home_backup/Projects/PediatricMeningioma/data/metadata.csv")
metadata$sample_id = gsub("-", "_", metadata$sample_id)
metadata2 = metadata[ metadata$name %in% basename(kallistoFiles), ]

save(metadata2, file="/athena/masonlab/scratch/users/nai2008/RNAseq/cavatica/metadata2.rda")
#load("/athena/masonlab/scratch/users/nai2008/RNAseq/cavatica/metadata2.rda")

names(kallistoFiles) = mapvalues(basename(kallistoFiles), metadata2$name, metadata2$sample_id)

library(EnsDb.Hsapiens.v86)
edb = EnsDb.Hsapiens.v86
tx = transcripts(edb, columns = c("tx_name", "gene_id", "gene_name"), return.type="DataFrame")
tx2 = tx[, c(1,2)]

tx_kallisto = tximport(kallistoFiles, type = "kallisto", tx2gene = tx2, ignoreTxVersion=T)

save(tx_kallisto, file="/athena/masonlab/scratch/users/nai2008/RNAseq/cavatica/cavatica_kallistoSymbol_counts.rda")
#load("/athena/masonlab/scratch/users/nai2008/RNAseq/cavatica/cavatica_kallistoSymbol_counts.rda")

# make count matrix
counts = round(tx_kallisto$counts)
mode(counts) = "integer"

all(metadata2$sample_id==colnames(counts)) #TRUE

length(which(duplicated(rownames(counts))==TRUE)) # 0, all rownames of the count matrix are unique

# keep only 'high-grade glioma' samples
metadata2=metadata2[grep('High-grade glioma', metadata2$disease_type),]

# drop samples with unlabeled sex
drop=which(metadata2$gender =='Not Available')
metadata2=metadata2[-drop,]

# drop samples with unlabeled age
drop=which(is.na(metadata2$age_at_diagnosis))
metadata2=metadata2[-drop,]

# convert age from days to years
metadata2$age_at_Dx_yrs=signif(metadata2$age_at_diagnosis/365,3)

# remove samples >21 yo @ Dx
drop=which(metadata2$age_at_Dx_yrs > 21 )
metadata2=metadata2[-drop,]

# add info regading when samples was taken (eg, at diagnosis, recurrence, second malignancy)
mwst=read.csv('/athena/masonlab/scratch/users/nai2008/RNAseq/cavatica/metadata_pHGG_samples_2020-09-03_withSamplingTime.csv')
mm=match(metadata2$sample_id,mwst$sample_id)
metadata2$Diagnosis.Type=as.vector(mwst$Diagnosis.Type[mm])

# keep only samples in which the 'initial tumor' was sampled (remove 'proressive','recurrence', and 'second malignancy' diagnosis types)
metadata2=metadata2[which(metadata2$Diagnosis.Type %in% c('Initial CNS Tumor')),]
metadata2=metadata2[-which(duplicated(metadata2$Kids.First.Participant.ID)=='TRUE'),]

nrow(metadata2) #49
length(unique(metadata2$Kids.First.Participant.ID)) #49
table(as.vector(metadata2$gender))
# Female   Male
#     25     24

# make pheno table
phenoData=data.frame(age=as.vector(metadata2$age_at_Dx_yrs),
	sex=factor(as.vector(metadata2$gender)),
	Dx='HGG',
	Study_ID='Cavatica',
	my_unique_ID=paste0('Cavatica_', metadata2$sample_id), sample_id=metadata2$sample_id, case_id=metadata2$case_id, 
	primary_site=metadata2$primary_site)

mm=match(metadata2$sample_id,colnames(counts))
counts=counts[,mm]

all(metadata2$sample_id==colnames(counts)) #TRUE

colnames(counts)=phenoData$my_unique_ID

save(counts, phenoData, file="/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/cavatica_counts_cleaned.rda")

write.csv(phenoData, file="/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/data_tables/cavatica_phenoData.csv", row.names=FALSE)

## Differential expresson

library(DESeq2)

load("/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/cavatica_counts_cleaned.rda")

# How many samples
nrow(phenoData) # 49
table(phenoData$sex)
# Female   Male
#     25     24

range(phenoData$age) # 1.35 18.60

# make DESeq2 object
dds=DESeqDataSetFromMatrix(countData = counts, colData = phenoData, design = ~ sex)

# estimate the library size correction and save the normalized counts matrix
dds <- estimateSizeFactors(dds)
norm.cts <- counts(dds, normalized=TRUE)

# we want a normalized count of at least 10 in 4 or more samples
filter <- rowSums(norm.cts >= 10) >= 4
norm.cts <- norm.cts[filter,]
counts <- counts[filter,]

# run SVA
library(sva)

mm <- model.matrix(~ sex, colData(dds))
mm0 <- model.matrix(~ 1, colData(dds))

svaobj <- svaseq(norm.cts, mod=mm, mod0=mm0) 
# Number of significant surrogate variables is: 9
colnames(svaobj$sv)=paste0('SV_',1:ncol(svaobj$sv))

phenoData=cbind(phenoData, svaobj$sv)

## differetial expression analysis

# make DESeq2 object
dds=DESeqDataSetFromMatrix(countData = counts, colData = phenoData, 
	design = ~ SV_1 + SV_2 + SV_3 + SV_4 + SV_5 + SV_6 + SV_7 + SV_8 + SV_9 + sex)

# run DE analysis
dds=DESeq(dds)

# results
rr=results(dds, alpha=0.1)
summary(rr)
# out of 32496 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 60, 0.18%
# LFC < 0 (down)     : 23, 0.071%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 2)

pdf('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/pdfs/HGG_Cooks_distances.pdf')

par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2, cex.axis=.8)

dev.off()

# add gene symbols, chr, & Ensembl gene IDs
library(EnsDb.Hsapiens.v86)
edb = EnsDb.Hsapiens.v86
genes=as.data.frame(genes(edb))

mm=match(rownames(rr), genes$gene_id)
length(which(is.na(mm))) # 0
rr$chr=as.vector(genes$seqnames[mm])
rr$Ensembl=as.vector(rownames(rr))
rr$gene=as.vector(genes$gene_name[mm])

save(rr, file='/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/cases_0to21_DEseq2_results.rda')
# load('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/cases_0to21_DEseq2_results.rda')

# print results
sig_results=as.data.frame(rr[which(rr$padj<=0.1),])

# order results by LFC
oo=order(sig_results$log2FoldChange)
sig_results=sig_results[oo,]

# make output table
out=data.frame(gene=sig_results$gene, chr=sig_results$chr, Ensembl=sig_results$Ensembl, log2FoldChange=sig_results$log2FoldChange, FDR=sig_results$padj) 

nrow(out)
nrow(out[! out$chr %in% c('X','Y'),])
# 83 DE genes; 43 of them on autosomes

### Are any of the autosomal DE genes TFs or imprinted genes? Add this info to the output table

out$TF=FALSE
out$imprinted_gene=FALSE
out$expressed_allele=NA

# load imprinted genes
load('/athena/masonlab/scratch/users/nai2008/Imprinted_genes/imprinted_genes.rda') # ig

# drop imprinted genes without an Ensembl ID
ig=ig[-which(is.na(ig$Ensembl.ID)),]
ig$Ensembl.ID=as.vector(ig$Ensembl.ID)
nrow(ig) #102
length(unique(ig$Ensembl.ID)) #102

# load TFs
TFs=read.csv("/athena/masonlab/scratch/users/nai2008/Human_TFs_DatabaseExtract_v_1.01.csv")
TFs$Ensembl_ID=as.vector(TFs$Ensembl_ID)
nrow(TFs) # 2765
length(unique(TFs$Ensembl_ID)) # 2765
which(is.na(TFs$Ensembl_ID)) # 0

# How many autosomal DE genes are imprinted genes:
mm=match(out$Ensembl,ig$Ensembl.ID)
out$imprinted_gene[which(!is.na(mm))]=TRUE
length(which(out$imprinted_gene==TRUE)) # 1
out$gene[which(out$imprinted_gene==TRUE)] # TP73

index_ig=which(out$imprinted_gene==TRUE)

if (length(index_ig)>0){
	for (i in 1:length(index_ig)){
		m=match(out$Ensembl[index_ig[i]],ig$Ensembl.ID)
		out$expressed_allele[index_ig[i]]=as.vector(ig$ExpressedAllele[m])
	}
}

# How many DE genes are TFs:
mm=match(out$Ensembl,TFs$Ensembl_ID)
out$TF[which(!is.na(mm))]=TRUE
length(which(out$TF==TRUE)) # 11
out$gene[which(out$TF==TRUE)] # BNC1  TP73  GLIS1 ZFX   KDM5C HOXA9 HSFY2 KDM5D TBL1Y ZFY   SRY

save(out, file='/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/HGG_0to21_DE_genes.rda')
# load('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/HGG_0to21_DE_genes.rda')

# print table of autosomal genes only
write.csv(out[! out$chr %in% c('X','Y'),],file="/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/data_tables/DE_genes_autosomalOnly_HGG_samples_peds.csv", row.names=FALSE)

# print full table (autosomal and XY genes)
write.csv(out,file="/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/data_tables/DE_genes_autosomalandXY_HGG_samples_peds.csv", row.names=FALSE)


######### correlation b/w DNAm and gene expression

## match gene expression to DNAm [DMPs]
dmps_overlapping_genes=read.csv('/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/DMPs_allPeds_cases.csv')

mm=match(out$Ensembl, dmps_overlapping_genes$Ensembl)
length(out$gene[which(!is.na(mm))]) # 0 of the DE genes are overlapping ot are within 5kb of at least one DMP

if (length(out$gene[which(!is.na(mm))]) != 0) {

	gg=as.vector(out$Ensembl[which(!is.na(mm))])

	index=which(as.vector(dmps_overlapping_genes$Ensembl) == gg[1])
	meth_exprs_table=dmps_overlapping_genes[index,c(9,10,8,7,3,5,11)]

	for (i in 2:length(gg)){
		index=which(as.vector(dmps_overlapping_genes$Ensembl) == gg[i])
		aa_meth=dmps_overlapping_genes[index,c(9,10,8,7,3,5,11)]

		meth_exprs_table=rbind(meth_exprs_table,aa_meth)

	}


	mm=match(meth_exprs_table$Ensembl,out$Ensembl)

	if (all(as.vector(meth_exprs_table$Ensembl)==as.vector(out$Ensembl[mm]))) {	
		meth_exprs_table=cbind(meth_exprs_table, out[mm,c(5,6)])
		meth_exprs_table$gene=out$gene[mm]
		colnames(meth_exprs_table)[3]='DMP_pos'
		colnames(meth_exprs_table)[5]='DNAm_avg_change'
		colnames(meth_exprs_table)[6]='DNAm_FDR'
		colnames(meth_exprs_table)[8]='gene_expression_log2FC'
		colnames(meth_exprs_table)[9]='gene_expression_FDR'

	write.csv(meth_exprs_table,file="/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/data_tables/HGG_meth_exprs_table.csv", row.names=FALSE)

	}

}

## match gene expression to DNAm [DMRs] {0-21 year old cases}

# DMRs for 0-21 y.o. cases overlap genes: ENSG00000270123
dog='ENSG00000270123'
mm=match(dog, out$Ensembl)
length(which(!is.na(mm))) # 0

# volcano plot [autosomal genes only]
#library(EnhancedVolcano)
library(ggplot2)

rr_autosomal = rr[! rr$chr %in% c('X','Y'),]
which(is.na(rr_autosomal$padj)) # 0

pdf('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/pdfs/volcano_plot_HGG_0to21_autosomalGenes.pdf')

source('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/EnahancedVolcano_custom.R')

a='NS (FDR > 0.1)'
b=as.expression(bquote("|"~log[2]~"|"~"FC"~">"~"2"))
c=as.expression(bquote('FDR'<='0.1'))
d=as.expression(bquote("|"~log[2]~"|"~"FC"~">"~"2"~"&"~'FDR'<='0.1'))

EnhancedVolcano_custom(toptable=as.data.frame(rr_autosomal[,c(2,6)]),
	lab = rr_autosomal$gene,
	x = 'log2FoldChange',
	y = 'padj',
	ylab =  bquote(~-Log[10]~italic(FDR)),
	pCutoff = 0.1,
	FCcutoff=2,
	title='pHGG samples (0-21 y.o.)',
	subtitle='Autosomal genes differentially expressed by sex',
	caption=NULL,
	legendLabels=c(a,b,c,d),
	legendPosition="bottom",
	border='full',
	col=c('grey30', 'forestgreen', 'royalblue', 'red2'),
	colAlpha=1,
	xlim=c(-13,30),
	#ylim=c(0,12.25),
	shape=21,
	labhjust=0,
	legendLabSize=8,
	legendIconSize=15
	)

#plot the legend

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
a='Non-significant (FDR > 0.1)'
b=as.expression(bquote(log[2]~"FC"~">"~"|2|"))
c=as.expression(bquote('FDR'<='0.1'))
d=as.expression(bquote(log[2]~"FC"~">"~"|2|"~"&"~'FDR'<='0.1'))
title=expression(bold('Legend'))
legend("center", legend=c(a,b,c,d), pch=21, col='black', pt.bg=c("grey30", "forestgreen", "royalblue", "red2"),  cex=1.2, pt.cex=2, title=title)

dev.off()

# volcano plot [X & Y genes only]
#library(EnhancedVolcano)
library(ggplot2)

rr_xy = rr[ rr$chr %in% c('X','Y'),]
which(is.na(rr_xy$padj)) # 0

pdf('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/pdfs/volcano_plot_HGG_0to21_XYGenes.pdf')

source('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/EnahancedVolcano_custom.R')

a='NS (FDR > 0.1)'
b=as.expression(bquote("|"~log[2]~"|"~"FC"~">"~"2"))
c=as.expression(bquote('FDR'<='0.1'))
d=as.expression(bquote("|"~log[2]~"|"~"FC"~">"~"2"~"&"~'FDR'<='0.1'))

EnhancedVolcano_custom(toptable=as.data.frame(rr_xy[,c(2,6)]),
	lab = rr_xy$gene,
	x = 'log2FoldChange',
	y = 'padj',
	ylab =  bquote(~-Log[10]~italic(FDR)),
	pCutoff = 0.1,
	FCcutoff=2,
	title='pHGG samples (0-21 y.o.)',
	subtitle='X & Y genes differentially expressed by sex',
	caption=NULL,
	legendLabels=c(a,b,c,d),
	legendPosition="bottom",
	border='full',
	col=c('grey30', 'forestgreen', 'royalblue', 'red2'),
	colAlpha=1,
	xlim=c(-10,13),
	#ylim=c(X,X),
	shape=21,
	labhjust=0,
	legendLabSize=8,
	legendIconSize=15
	)

#plot the legend

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
a='Non-significant (FDR > 0.1)'
b=as.expression(bquote(log[2]~"FC"~">"~"|2|"))
c=as.expression(bquote('FDR'<='0.1'))
d=as.expression(bquote(log[2]~"FC"~">"~"|2|"~"&"~'FDR'<='0.1'))
title=expression(bold('Legend'))
legend("center", legend=c(a,b,c,d), pch=21, col='black', pt.bg=c("grey30", "forestgreen", "royalblue", "red2"),  cex=1.2, pt.cex=2, title=title)

dev.off()

###################################################################################################
## LIBD data: 0-21 y.o. controls {HIPPOCAMPUS}
###################################################################################################

# load data
library(SummarizedExperiment)

load("/athena/masonlab/scratch/users/nai2008/RNAseq/rse_gene_unfiltered.Rdata") #rse_gene

all(colData(rse_gene)$RNum==colnames(assays(rse_gene)$counts)) #TRUE

# make pheno data
pd=colData(rse_gene)

# load cell type proportions (computed by performing RNA deconvolution using an algorithm developed by Darmanis et al)
load("/athena/masonlab/scratch/users/nai2008/RNAseq/RNA_cell_proportions_brainSeq_phase2.rda") # propEsts

mm=match(pd$RNum,rownames(propEsts))

pd=cbind(pd,propEsts[mm,])

# make counts matrix
counts=assays(rse_gene)$counts
counts = round(counts)
mode(counts) = "integer"

# drop non-controls
pd=pd[which(pd$Dx=='Control'),]

# drop samples >21 yo
pd=pd[pd$Age<=21,]

# drop fetal samples
pd=pd[pd$Age>0,]

# keep only hippocampal samples
pd=pd[which(pd$Region=='HIPPO'),]

#
mm=match(pd$RNum,colnames(counts))
counts=counts[,mm]

sex=pd$Sex
	sex[which(sex=='M')]='Male'
	sex[which(sex=='F')]='Female'

nrow(pd) # 60
length(unique(pd$RNum)) # 60
length(unique(pd$BrNum)) # 60
table(sex)
# Female   Male
#     18     42

phenoData=data.frame(age=pd$Age,
sex=factor(as.vector(sex)),
Dx='Control',
Study_ID='LIBD',
my_unique_ID=paste0('LIBD_', pd$RNum), RNum=pd$RNum, case_id=pd$BrNum, site=pd$Region)

# add cell-type proportion estimates
phenoData=cbind(phenoData,pd[,55:62])

all(colnames(counts)==phenoData$RNum) #TRUE

colnames(counts)=as.vector(phenoData$my_unique_ID)

save(counts, phenoData, file="/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/LIBD_controls_hippo_counts_cleaned.rda")

write.csv(phenoData, file="/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/data_tables/LIBD_controls_hippo_phenoData.csv")

## Differential expresson

library(DESeq2)

load("/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/LIBD_controls_hippo_counts_cleaned.rda")

# How many samples
nrow(phenoData) # 60
length(unique(phenoData$RNum)) # 60
length(unique(phenoData$case_id)) # 60

nrow(phenoData) # 60
table(phenoData$sex)
# Female   Male
#     18     42

range(phenoData$age) # 0.06845 20.85000

range(phenoData$Astrocytes)
median(phenoData$Astrocytes)
range(phenoData$Neurons)
median(phenoData$Neurons)

# > range(phenoData$Astrocytes)
# [1] 0.2703236 0.5362156
# > median(phenoData$Astrocytes)
# [1] 0.3885961
# > range(phenoData$Neurons)
# [1] 0.1835502 0.8022841
# > median(phenoData$Neurons)
# [1] 0.5990485

# make DESeq2 object
dds=DESeqDataSetFromMatrix(countData = counts, colData = phenoData, design = ~ sex)

# estimate the library size correction and save the normalized counts matrix
dds <- estimateSizeFactors(dds)
norm.cts <- counts(dds, normalized=TRUE)

# we want a normalized count of at least 10 in 4 or more samples
filter <- rowSums(norm.cts >= 10) >= 4
norm.cts <- norm.cts[filter,]
counts <- counts[filter,]

# run SVA
library(sva)

mm <- model.matrix(~ sex, colData(dds))
mm0 <- model.matrix(~ 1, colData(dds))

svaobj <- svaseq(norm.cts, mod=mm, mod0=mm0, n.sv=11) # Number of significant surrogate variables varies b/w 11 and 12, so we set n.sv=11
colnames(svaobj$sv)=paste0('SV_',1:ncol(svaobj$sv))

phenoData=cbind(phenoData, svaobj$sv)

## differetial expression analysis

# make DESeq2 object
dds=DESeqDataSetFromMatrix(countData = counts, colData = phenoData, 
	design = ~ SV_1 + SV_2 + SV_3 + SV_4 + SV_5 + SV_6 + SV_7 + SV_8 + SV_9 + SV_10 + SV_11 + sex)

# run DE analysis
dds=DESeq(dds)

# results
rr=results(dds, alpha=0.1)
summary(rr)
# out of 25663 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 71, 0.28%
# LFC < 0 (down)     : 35, 0.14%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 1)

pdf('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/pdfs/Controls_hippocampus_Cooks_distances.pdf')

par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

dev.off()

source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')
tt=ss(rownames(rr),"\\.",1)
rownames(rr)=tt

# add gene symbols, chr, & Ensembl gene IDs
library(EnsDb.Hsapiens.v86)
edb = EnsDb.Hsapiens.v86
genes=as.data.frame(genes(edb))

mm=match(rownames(rr), genes$gene_id)
length(which(is.na(mm))) # 0
rr$chr=as.vector(genes$seqnames[mm])
rr$Ensembl=as.vector(rownames(rr))
rr$gene=as.vector(genes$gene_name[mm])

save(rr,file='/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/controls_hippo_0to21_DEseq2_results.rda')
# load('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/controls_hippo_0to21_DEseq2_results.rda')

# print results
sig_results=as.data.frame(rr[which(rr$padj<=0.1),])

# order results by LFC
oo=order(sig_results$log2FoldChange)
sig_results=sig_results[oo,]

# make output table
out=data.frame(gene=sig_results$gene, chr=sig_results$chr, Ensembl=sig_results$Ensembl, log2FoldChange=sig_results$log2FoldChange, FDR=sig_results$padj) 

nrow(out)
nrow(out[! out$chr %in% c('X','Y'),])
# 106 DE genes; 49 of them on autosomes

### Are any of the DE genes TFs or imprinted genes? Add this info to the output table

out$TF=FALSE
out$imprinted_gene=FALSE
out$expressed_allele=NA

# load imprinted genes
load('/athena/masonlab/scratch/users/nai2008/Imprinted_genes/imprinted_genes.rda') # ig

# drop imprinted genes without an Ensembl ID
length(which(is.na(ig$Ensembl.ID))) # 2
ig=ig[-which(is.na(ig$Ensembl.ID)),]
ig$Ensembl.ID=as.vector(ig$Ensembl.ID)
nrow(ig) #102
length(unique(ig$Ensembl.ID)) #102

# load TFs
TFs=read.csv("/athena/masonlab/scratch/users/nai2008/Human_TFs_DatabaseExtract_v_1.01.csv")
TFs$Ensembl_ID=as.vector(TFs$Ensembl_ID)
nrow(TFs) # 2765
length(unique(TFs$Ensembl_ID)) # 2765
which(is.na(TFs$Ensembl_ID)) # 0

# How many DE genes are imprinted genes:
mm=match(as.vector(out$Ensembl),ig$Ensembl.ID)
out$imprinted_gene[which(!is.na(mm))]=TRUE
length(which(out$imprinted_gene==TRUE)) # 1
out$gene[which(out$imprinted_gene==TRUE)] # NDN

index_ig=which(out$imprinted_gene==TRUE)

if (length(index_ig)>0){
	for (i in 1:length(index_ig)){
		m=match(out$Ensembl[index_ig[i]],ig$Ensembl.ID)
		out$expressed_allele[index_ig[i]]=as.vector(ig$ExpressedAllele[m])
	}
}

# How many DE genes are TFs:
mm=match(out$Ensembl,TFs$Ensembl_ID)
out$TF[which(!is.na(mm))]=TRUE
length(which(out$TF==TRUE)) # 14
out$gene[which(out$TF==TRUE)]
#  [1] ZFX     KDM5C   DDX3X   ZNF502  ZNF320  TRIP11  LARP7   NDN     HES6
# [10] ZMYND10 HSFY2   TBL1Y   ZFY     KDM5D

save(out, file='/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/controls_hippo_0to21_DE_genes.rda')
#load('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/controls_hippo_0to21_DE_genes.rda')

# print table of autosomal genes only
write.csv(out[! out$chr %in% c('X','Y'),],file="/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/data_tables/DE_genes_autosomalOnly_hippo_control_samples_peds.csv", row.names=FALSE)

# print full table (autosomal and XY genes)
write.csv(out,file="/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/data_tables/DE_genes_autosomalandXY_hippo_control_samples_peds.csv", row.names=FALSE)


######### correlation b/w DNAm and gene expression

## match gene expression to DNAm [DMPs]
dmps_overlapping_genes=read.csv('/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/DMPs_allPeds_controls.csv')

mm=match(out$Ensembl, dmps_overlapping_genes$Ensembl)
length(out$gene[which(!is.na(mm))]) # 6 of the DE genes are overlapping ot are within 5kb of at least one DMP

if (length(out$gene[which(!is.na(mm))]) != 0) {

	gg=as.vector(out$Ensembl[which(!is.na(mm))])

	index=which(as.vector(dmps_overlapping_genes$Ensembl) == gg[1])
	meth_exprs_table=dmps_overlapping_genes[index,c(9,10,8,7,3,5,11)]

	for (i in 2:length(gg)){
		index=which(as.vector(dmps_overlapping_genes$Ensembl) == gg[i])
		aa_meth=dmps_overlapping_genes[index,c(9,10,8,7,3,5,11)]

		meth_exprs_table=rbind(meth_exprs_table,aa_meth)

	}


	mm=match(meth_exprs_table$Ensembl,out$Ensembl)

	if (all(as.vector(meth_exprs_table$Ensembl)==as.vector(out$Ensembl[mm]))) {	
		meth_exprs_table=cbind(meth_exprs_table, out[mm,c(4,5)])
		meth_exprs_table$gene=out$gene[mm]
		colnames(meth_exprs_table)[3]='DMP_pos'
		colnames(meth_exprs_table)[5]='DNAm_avg_difference'
		colnames(meth_exprs_table)[6]='DNAm_FDR'
		colnames(meth_exprs_table)[8]='gene_expression_log2FC'
		colnames(meth_exprs_table)[9]='gene_expression_FDR'

	write.csv(meth_exprs_table,file="/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/data_tables/controls_hippocampus_meth_exprs_table.csv", row.names=FALSE)

	}

}

## match gene expression to DNAm [DMRs] {0-21 y.o. controls}

# DMRs for 0-21 y.o. controls overlap genes: ENSG00000271440; ENSG00000287279
dog=c('ENSG00000271440','ENSG00000287279')
mm=match(dog, out$Ensembl)
length(which(!is.na(mm))) # 0

# volcano plot [autosomal genes only]
#library(EnhancedVolcano)
library(ggplot2)

rr_autosomal = rr[! rr$chr %in% c('X','Y'),]
which(is.na(rr_autosomal$padj)) # 0

pdf('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/pdfs/volcano_plot_hippo_controls_0to21_autosomalGenes.pdf')

source('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/EnahancedVolcano_custom.R')

a='NS (FDR > 0.1)'
b=as.expression(bquote("|"~log[2]~"|"~"FC"~">"~"2"))
c=as.expression(bquote('FDR'<='0.1'))
d=as.expression(bquote("|"~log[2]~"|"~"FC"~">"~"2"~"&"~'FDR'<='0.1'))

EnhancedVolcano_custom(toptable=as.data.frame(rr_autosomal[,c(2,6)]),
	lab = rr_autosomal$gene,
	x = 'log2FoldChange',
	y = 'padj',
	ylab =  bquote(~-Log[10]~italic(FDR)),
	pCutoff = 0.1,
	FCcutoff=2,
	title='Controls (0-21 y.o.); hippocampal samples',
	subtitle='Autosomal genes differentially expressed by sex',
	caption=NULL,
	legendLabels=c(a,b,c,d),
	legendPosition="bottom",
	border='full',
	col=c('grey30', 'forestgreen', 'royalblue', 'red2'),
	colAlpha=1,
	#xlim=c(X,X),
	#ylim=c(X,X),
	shape=21,
	labhjust=0,
	legendLabSize=8,
	legendIconSize=15
	)

#plot the legend

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
a='Non-significant (FDR > 0.1)'
b=as.expression(bquote(log[2]~"FC"~">"~"|2|"))
c=as.expression(bquote('FDR'<='0.1'))
d=as.expression(bquote(log[2]~"FC"~">"~"|2|"~"&"~'FDR'<='0.1'))
title=expression(bold('Legend'))
legend("center", legend=c(a,b,c,d), pch=21, col='black', pt.bg=c("grey30", "forestgreen", "royalblue", "red2"),  cex=1.2, pt.cex=2, title=title)

dev.off()

# volcano plot [X & Y genes only]
#library(EnhancedVolcano)
library(ggplot2)

rr_xy = rr[ rr$chr %in% c('X','Y'),]
which(is.na(rr_xy$padj)) # 0

pdf('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/pdfs/volcano_plot_hippo_controls_0to21_XYGenes.pdf')

source('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/EnahancedVolcano_custom.R')

a='NS (FDR > 0.1)'
b=as.expression(bquote("|"~log[2]~"|"~"FC"~">"~"2"))
c=as.expression(bquote('FDR'<='0.1'))
d=as.expression(bquote("|"~log[2]~"|"~"FC"~">"~"2"~"&"~'FDR'<='0.1'))

EnhancedVolcano_custom(toptable=as.data.frame(rr_xy[,c(2,6)]),
	lab = rr_xy$gene,
	x = 'log2FoldChange',
	y = 'padj',
	ylab =  bquote(~-Log[10]~italic(FDR)),
	pCutoff = 0.1,
	FCcutoff=2,
	title='Controls (0-21 y.o.); hippocampal samples',
	subtitle='X & Y genes differentially expressed by sex',
	caption=NULL,
	legendLabels=c(a,b,c,d),
	legendPosition="bottom",
	border='full',
	col=c('grey30', 'forestgreen', 'royalblue', 'red2'),
	colAlpha=1,
	#xlim=c(X,X),
	#ylim=c(X,X),
	shape=21,
	labhjust=0,
	legendLabSize=8,
	legendIconSize=15
	)

#plot the legend

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
a='Non-significant (FDR > 0.1)'
b=as.expression(bquote(log[2]~"FC"~">"~"|2|"))
c=as.expression(bquote('FDR'<='0.1'))
d=as.expression(bquote(log[2]~"FC"~">"~"|2|"~"&"~'FDR'<='0.1'))
title=expression(bold('Legend'))
legend("center", legend=c(a,b,c,d), pch=21, col='black', pt.bg=c("grey30", "forestgreen", "royalblue", "red2"),  cex=1.2, pt.cex=2, title=title)

dev.off()

###################################################################################################
## LIBD data: 0-21 y.o. controls {DLPFC}
###################################################################################################

# load data
library(SummarizedExperiment)

load("/athena/masonlab/scratch/users/nai2008/RNAseq/rse_gene_unfiltered.Rdata") #rse_gene

all(colData(rse_gene)$RNum==colnames(assays(rse_gene)$counts)) #TRUE

# make pheno data
pd=colData(rse_gene)

# load cell type proportions (computed by performing RNA deconvolution using an algorithm developed by Darmanis et al)
load("/athena/masonlab/scratch/users/nai2008/RNAseq/RNA_cell_proportions_brainSeq_phase2.rda") # propEsts

mm=match(pd$RNum,rownames(propEsts))

pd=cbind(pd,propEsts[mm,])

# make counts matrix
counts=assays(rse_gene)$counts
counts = round(counts)
mode(counts) = "integer"

# drop non-controls
pd=pd[which(pd$Dx=='Control'),]

# drop samples >21 yo
pd=pd[pd$Age<=21,]

# drop fetal samples
pd=pd[pd$Age>0,]

# keep only DLPFC samples
pd=pd[which(pd$Region=='DLPFC'),] ######

#
mm=match(pd$RNum,colnames(counts))
counts=counts[,mm]

sex=pd$Sex
	sex[which(sex=='M')]='Male'
	sex[which(sex=='F')]='Female'

nrow(pd) # 65
length(unique(pd$RNum)) # 65
length(unique(pd$BrNum)) # 65
table(sex)
# Female   Male
#     24     41

phenoData=data.frame(age=pd$Age,
sex=factor(as.vector(sex)),
Dx='Control',
Study_ID='LIBD',
my_unique_ID=paste0('LIBD_', pd$RNum), RNum=pd$RNum, case_id=pd$BrNum, site=pd$Region)

# add cell-type proportion estimates
phenoData=cbind(phenoData,pd[,55:62])

all(colnames(counts)==phenoData$RNum) #TRUE

colnames(counts)=as.vector(phenoData$my_unique_ID)

save(counts, phenoData, file="/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/LIBD_controls_dlpfc_counts_cleaned.rda")

write.csv(phenoData, file="/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/data_tables/LIBD_controls_dlpfc_phenoData.csv")

## Differential expresson

library(DESeq2)

load("/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/LIBD_controls_dlpfc_counts_cleaned.rda")

# How many samples
nrow(phenoData) # 65
length(unique(phenoData$RNum)) # 65
length(unique(phenoData$case_id)) # 65

# How many samples
nrow(phenoData) # 65
table(phenoData$sex)
# Female   Male
#     24     41

range(phenoData$age)
# 0.005475 20.850000

range(phenoData$Astrocytes)
median(phenoData$Astrocytes)
range(phenoData$Neurons)
median(phenoData$Neurons)
# > range(phenoData$Astrocytes)
# [1] 0.1289110 0.4516653
# > median(phenoData$Astrocytes)
# [1] 0.360621
# > range(phenoData$Neurons)
# [1] 2.470641e-18 8.992461e-01
# > median(phenoData$Neurons)
# [1] 0.7086263

# make DESeq2 object
dds=DESeqDataSetFromMatrix(countData = counts, colData = phenoData, design = ~ sex)

# estimate the library size correction and save the normalized counts matrix
dds <- estimateSizeFactors(dds)
norm.cts <- counts(dds, normalized=TRUE)

# we want a normalized count of at least 10 in 4 or more samples
filter <- rowSums(norm.cts >= 10) >= 4
norm.cts <- norm.cts[filter,]
counts <- counts[filter,]

# run SVA
library(sva)

mm <- model.matrix(~ sex, colData(dds))
mm0 <- model.matrix(~ 1, colData(dds))

svaobj <- svaseq(norm.cts, mod=mm, mod0=mm0, n.sv=7) # Number of significant surrogate variables varies b/w 7 and 8, so we set n.sv=7
colnames(svaobj$sv)=paste0('SV_',1:ncol(svaobj$sv))

phenoData=cbind(phenoData, svaobj$sv)

## differetial expression analysis

# make DESeq2 object
dds=DESeqDataSetFromMatrix(countData = counts, colData = phenoData, 
	design = ~ SV_1 + SV_2 + SV_3 + SV_4 + SV_5 + SV_6 + SV_7 + sex)

# run DE analysis
dds=DESeq(dds)

# results
rr=results(dds, alpha=0.1)
summary(rr)
# out of 26077 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 62, 0.24%
# LFC < 0 (down)     : 36, 0.14%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 1)

pdf('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/pdfs/Controls_DLPFC_Cooks_distances.pdf')

par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

dev.off()

source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')
tt=ss(rownames(rr),"\\.",1)
rownames(rr)=tt

# add gene symbols, chr, & Ensembl gene IDs
library(EnsDb.Hsapiens.v86)
edb = EnsDb.Hsapiens.v86
genes=as.data.frame(genes(edb))

mm=match(rownames(rr), genes$gene_id)
length(which(is.na(mm))) # 0
rr$chr=as.vector(genes$seqnames[mm])
rr$Ensembl=as.vector(rownames(rr))
rr$gene=as.vector(genes$gene_name[mm])

save(rr,file='/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/controls_dlpfc_0to21_DEseq2_results.rda')
# load('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/controls_dlpfc_0to21_DEseq2_results.rda')

# print results
sig_results=as.data.frame(rr[which(rr$padj<=0.1),])

# order results by LFC
oo=order(sig_results$log2FoldChange)
sig_results=sig_results[oo,]

# make output table
out=data.frame(gene=sig_results$gene, chr=sig_results$chr, Ensembl=sig_results$Ensembl, log2FoldChange=sig_results$log2FoldChange, FDR=sig_results$padj) 

nrow(out)
nrow(out[! out$chr %in% c('X','Y'),])
# 98 DE genes; 42 of them on autosomes

### Are any of the DE genes TFs or imprinted genes? Add this info to the output table

out$TF=FALSE
out$imprinted_gene=FALSE
out$expressed_allele=NA

# load imprinted genes
load('/athena/masonlab/scratch/users/nai2008/Imprinted_genes/imprinted_genes.rda') # ig

# drop imprinted genes without an Ensembl ID
length(which(is.na(ig$Ensembl.ID))) # 2
ig=ig[-which(is.na(ig$Ensembl.ID)),]
ig$Ensembl.ID=as.vector(ig$Ensembl.ID)
nrow(ig) #102
length(unique(ig$Ensembl.ID)) #102

# load TFs
TFs=read.csv("/athena/masonlab/scratch/users/nai2008/Human_TFs_DatabaseExtract_v_1.01.csv")
TFs$Ensembl_ID=as.vector(TFs$Ensembl_ID)
nrow(TFs) # 2765
length(unique(TFs$Ensembl_ID)) # 2765
which(is.na(TFs$Ensembl_ID)) # 0

# How many DE genes are imprinted genes:
mm=match(as.vector(out$Ensembl),ig$Ensembl.ID)
out$imprinted_gene[which(!is.na(mm))]=TRUE
length(which(out$imprinted_gene==TRUE)) # 0
#out$gene[which(out$imprinted_gene==TRUE)]

index_ig=which(out$imprinted_gene==TRUE)

if (length(index_ig)>0){
	for (i in 1:length(index_ig)){
		m=match(out$Ensembl[index_ig[i]],ig$Ensembl.ID)
		out$expressed_allele[index_ig[i]]=as.vector(ig$ExpressedAllele[m])
	}
}

# How many DE genes are TFs:
mm=match(out$Ensembl,TFs$Ensembl_ID)
out$TF[which(!is.na(mm))]=TRUE
length(which(out$TF==TRUE)) # 9
out$gene[which(out$TF==TRUE)] # ZFX   KDM5C DDX3X ZRSR2 DAP   TBL1Y HSFY2 ZFY   KDM5D

save(out, file='/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/controls_dlpfc_0to21_DE_genes.rda')
#load('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/controls_dlpfc_0to21_DE_genes.rda')

# print table of autosomal genes only
write.csv(out[! out$chr %in% c('X','Y'),],file="/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/data_tables/DE_genes_autosomalOnly_dlpfc_control_samples_peds.csv", row.names=FALSE)

# print full table (autosomal and XY genes)
write.csv(out,file="/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/data_tables/DE_genes_autosomalandXY_dlpfc_control_samples_peds.csv", row.names=FALSE)

######### correlation b/w DNAm and gene expression

## match gene expression to DNAm [DMPs]
dmps_overlapping_genes=read.csv('/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/DMPs_allPeds_controls.csv')

mm=match(out$Ensembl, dmps_overlapping_genes$Ensembl)
length(out$gene[which(!is.na(mm))]) # 7 of the DE genes are overlapping ot are within 5kb of at least one DMP

if (length(out$gene[which(!is.na(mm))]) != 0) {

	gg=as.vector(out$Ensembl[which(!is.na(mm))])

	index=which(as.vector(dmps_overlapping_genes$Ensembl) == gg[1])
	meth_exprs_table=dmps_overlapping_genes[index,c(9,10,8,7,3,5,11)]

	for (i in 2:length(gg)){
		index=which(as.vector(dmps_overlapping_genes$Ensembl) == gg[i])
		aa_meth=dmps_overlapping_genes[index,c(9,10,8,7,3,5,11)]

		meth_exprs_table=rbind(meth_exprs_table,aa_meth)

	}


	mm=match(meth_exprs_table$Ensembl,out$Ensembl)

	if (all(as.vector(meth_exprs_table$Ensembl)==as.vector(out$Ensembl[mm]))) {	
		meth_exprs_table=cbind(meth_exprs_table, out[mm,c(4,5)])
		meth_exprs_table$gene=out$gene[mm]
		colnames(meth_exprs_table)[3]='DMP_pos'
		colnames(meth_exprs_table)[5]='DNAm_avg_difference'
		colnames(meth_exprs_table)[6]='DNAm_FDR'
		colnames(meth_exprs_table)[8]='gene_expression_log2FC'
		colnames(meth_exprs_table)[9]='gene_expression_FDR'

	write.csv(meth_exprs_table,file="/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/data_tables/controls_dlpfc_meth_exprs_table.csv", row.names=FALSE)

	}

}

## match gene expression to DNAm [DMRs] {0-21 y.o. controls}

# DMRs for 0-21 y.o. controls overlap genes: ENSG00000271440; ENSG00000287279
dog=c('ENSG00000271440','ENSG00000287279')
mm=match(dog, out$Ensembl)
length(which(!is.na(mm))) # 0

# volcano plot [autosomal genes only]
#library(EnhancedVolcano)
library(ggplot2)

rr_autosomal = rr[! rr$chr %in% c('X','Y'),]
which(is.na(rr_autosomal$padj)) # 0

pdf('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/pdfs/volcano_plot_dlpfc_controls_0to21_autosomalGenes.pdf')

source('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/EnahancedVolcano_custom.R')

a='NS (FDR > 0.1)'
b=as.expression(bquote("|"~log[2]~"|"~"FC"~">"~"2"))
c=as.expression(bquote('FDR'<='0.1'))
d=as.expression(bquote("|"~log[2]~"|"~"FC"~">"~"2"~"&"~'FDR'<='0.1'))

EnhancedVolcano_custom(toptable=as.data.frame(rr_autosomal[,c(2,6)]),
	lab = rr_autosomal$gene,
	x = 'log2FoldChange',
	y = 'padj',
	ylab =  bquote(~-Log[10]~italic(FDR)),
	pCutoff = 0.1,
	FCcutoff=2,
	title='Controls (0-21 y.o.); DLPFC samples',
	subtitle='Autosomal genes differentially expressed by sex',
	caption=NULL,
	legendLabels=c(a,b,c,d),
	legendPosition="bottom",
	border='full',
	col=c('grey30', 'forestgreen', 'royalblue', 'red2'),
	colAlpha=1,
	xlim=c(-3,3.5),
	#ylim=c(X,X),
	shape=21,
	labhjust=0,
	legendLabSize=8,
	legendIconSize=15
	)

#plot the legend

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
a='Non-significant (FDR > 0.1)'
b=as.expression(bquote(log[2]~"FC"~">"~"|2|"))
c=as.expression(bquote('FDR'<='0.1'))
d=as.expression(bquote(log[2]~"FC"~">"~"|2|"~"&"~'FDR'<='0.1'))
title=expression(bold('Legend'))
legend("center", legend=c(a,b,c,d), pch=21, col='black', pt.bg=c("grey30", "forestgreen", "royalblue", "red2"),  cex=1.2, pt.cex=2, title=title)

dev.off()

# volcano plot [X & Y genes only]
#library(EnhancedVolcano)
library(ggplot2)

rr_xy = rr[ rr$chr %in% c('X','Y'),]
which(is.na(rr_xy$padj)) # 0

pdf('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/pdfs/volcano_plot_dlpfc_controls_0to21_XYGenes.pdf')

source('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/EnahancedVolcano_custom.R')

a='NS (FDR > 0.1)'
b=as.expression(bquote("|"~log[2]~"|"~"FC"~">"~"2"))
c=as.expression(bquote('FDR'<='0.1'))
d=as.expression(bquote("|"~log[2]~"|"~"FC"~">"~"2"~"&"~'FDR'<='0.1'))

EnhancedVolcano_custom(toptable=as.data.frame(rr_xy[,c(2,6)]),
	lab = rr_xy$gene,
	x = 'log2FoldChange',
	y = 'padj',
	ylab =  bquote(~-Log[10]~italic(FDR)),
	pCutoff = 0.1,
	FCcutoff=2,
	title='Controls (0-21 y.o.); DLPFC samples',
	subtitle='X & Y genes differentially expressed by sex',
	caption=NULL,
	legendLabels=c(a,b,c,d),
	legendPosition="bottom",
	border='full',
	col=c('grey30', 'forestgreen', 'royalblue', 'red2'),
	colAlpha=1,
	xlim=c(-8,10),
	#ylim=c(X,X),
	shape=21,
	labhjust=0,
	legendLabSize=8,
	legendIconSize=15
	)

#plot the legend

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
a='Non-significant (FDR > 0.1)'
b=as.expression(bquote(log[2]~"FC"~">"~"|2|"))
c=as.expression(bquote('FDR'<='0.1'))
d=as.expression(bquote(log[2]~"FC"~">"~"|2|"~"&"~'FDR'<='0.1'))
title=expression(bold('Legend'))
legend("center", legend=c(a,b,c,d), pch=21, col='black', pt.bg=c("grey30", "forestgreen", "royalblue", "red2"),  cex=1.2, pt.cex=2, title=title)

dev.off()


#NA Ivanov




