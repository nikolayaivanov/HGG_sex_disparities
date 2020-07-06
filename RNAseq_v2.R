
###################################################################################################
## Cavatica data
###################################################################################################
library(DESeq2)

# import data
library(tximport)
​
library(plyr)

kallistoFiles = list.files("/athena/masonlab/scratch/users/cem2009/home_backup/Projects/PediatricMeningioma/data/RNAseq_counts", "kallisto.abundance.tsv.gz", recursive=T, full.names=T)
metadata = read.csv("/athena/masonlab/scratch/users/cem2009/home_backup/Projects/PediatricMeningioma/data/metadata.csv")
metadata$sample_id = gsub("-", "_", metadata$sample_id)
metadata2 = metadata[ metadata$name %in% basename(kallistoFiles), ]

save(metadata2, file="/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/cavatica/metadata2.rds")
load("/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/cavatica/metadata2.rds")

names(kallistoFiles) = mapvalues(basename(kallistoFiles), metadata2$name, metadata2$sample_id)

library(EnsDb.Hsapiens.v86)
edb = EnsDb.Hsapiens.v86
tx = transcripts(edb, columns = c("tx_name", "gene_id", "gene_name"), return.type="DataFrame")
tx2 = tx[, c(1,3)]

tx_kallisto = tximport(kallistoFiles, type = "kallisto", tx2gene = tx2, ignoreTxVersion=T)
tpm_kallisto = tximport(kallistoFiles, type = "kallisto", tx2gene = tx2, ignoreTxVersion=T, countsFromAbundance="lengthScaledTPM")
​
save(tx_kallisto, file="/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/cavatica/cavatica_kallistoSymbol_counts.rda")
save(tpm_kallisto, file="/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/cavatica/cavatica_kallistoSymbol_tpm.rda")

load("/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/cavatica/cavatica_kallistoSymbol_counts.rda")
load("/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/cavatica/cavatica_kallistoSymbol_tpm.rda")

# make count and tpm matrices
tpm = tpm_kallisto$counts

counts = round(tx_kallisto$counts)
mode(counts) = "integer"

all(metadata2$sample_id==colnames(counts)) #TRUE
all(colnames(counts)==colnames(tpm)) #TRUE

length(which(duplicated(rownames(tpm))==TRUE)) # 0, all rownames of the tpm and count matrices are unique

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

# make pheno table
phenoData=data.frame(age=as.vector(metadata2$age_at_Dx_yrs),
	sex=factor(as.vector(metadata2$gender)),
	Dx='HGG',
	Study_ID='Cavatica',
	my_unique_ID=paste0('Cavatica_', 1:nrow(metadata2)), sample_id=metadata2$sample_id)

mm=match(metadata2$sample_id,colnames(counts))
counts=counts[,mm]
tpm=tpm[,mm]

all(metadata2$sample_id==colnames(counts)) #TRUE
all(metadata2$sample_id==colnames(tpm)) # TRUE

phenoData=phenoData[,-6]

colnames(counts)=phenoData$my_unique_ID
colnames(tpm)=phenoData$my_unique_ID

all(rownames(counts)==rownames(tpm))
all(colnames(counts)==colnames(tpm))

load("/athena/masonlab/scratch/users/nai2008/genes.rda") # from Ensembl Biomart (http://useast.ensembl.org/biomart/martview/)
test=match(rownames(tpm),genes$Gene.Name)
index_na=which(is.na(test))

tpm1=tpm[-index_na,]
counts1=counts[-index_na,]

tpm2=tpm[index_na,]
counts2=counts[index_na,]

test=match(rownames(tpm2),genes$Gene.Synonym)
index_na=which(is.na(test))
tpm2=tpm2[-index_na,]
counts2=counts2[-index_na,]
mm=match(rownames(tpm2),genes$Gene.Synonym)
rownames(tpm2)=genes$Gene.Name[mm]
rownames(counts2)=genes$Gene.Name[mm]
tpm=rbind(tpm1,tpm2)
counts=rbind(counts1,counts2)

# add duplicates together (for TPM matrix)
length( which( duplicated(rownames(tpm)) == TRUE ) ) #56

index_of_duplictes=which(duplicated(rownames(tpm))==TRUE)
dup_genes=rownames(tpm)[index_of_duplictes]

for (i in 1:length(dup_genes)){

dup=which(rownames(tpm)==dup_genes[i])
if (length(dup)>1){
	colsum=as.vector(colSums(tpm[dup,]))
	cs=matrix(colsum, nrow=1, ncol=length(colsum))
	rownames(cs)=dup_genes[i]
	tpm=tpm[-dup,]
	tpm=rbind(tpm,cs)
	}
}

length(which(duplicated(rownames(tpm))==TRUE)) #0

# add duplicates together (for counts matrix)
length( which( duplicated(rownames(counts)) == TRUE ) ) #56

index_of_duplictes=which(duplicated(rownames(counts))==TRUE)
dup_genes=rownames(counts)[index_of_duplictes]

for (i in 1:length(dup_genes)){

dup=which(rownames(counts)==dup_genes[i])
if (length(dup)>1){
	colsum=as.vector(colSums(counts[dup,]))
	cs=matrix(colsum, nrow=1, ncol=length(colsum))
	rownames(cs)=dup_genes[i]
	counts=counts[-dup,]
	counts=rbind(counts,cs)
	}
}

length(which(duplicated(rownames(counts))==TRUE)) #0

## check that sex is labeled correctly

# find genes on chrY
load("/athena/masonlab/scratch/users/nai2008/genes.rda")
y_genes=genes$Gene.Name[which(genes$chr=='Y')]
y_genes=y_genes[!duplicated(y_genes)]

mm=match(y_genes,rownames(tpm))
y_genes=y_genes[-which(is.na(mm))]

mm=match(y_genes,rownames(tpm))
tpm_y=tpm[mm,]

csy=colSums(tpm_y)
sex=factor(as.vector(phenoData$sex), levels=c('Female','Male')) # Female = 0; Male=1;

color=sub('Female','pink',sex)
color=sub('Male','blue', color)

pdf('/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/pdfs/sex_verification_Cavatica_data.pdf')

barplot(csy, names.arg=names(csy), col=color, las=2, cex.names=.7, main='Sum of TPM values for genes on chrY')
legend(x='topright',legend=c('Female','Male'), col=c('pink','blue'), pch=15, cex=.7, title="Labeled sex")

dev.off()

save(counts, phenoData, file="/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/rdas/cavatica_counts_cleaned.rda")

########################################
## Differential expresson: cases 0-21
########################################
library(DESeq2)

load("/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/rdas/cavatica_counts_cleaned.rda")

# make DESeq2 object
dds=DESeqDataSetFromMatrix(countData = counts, colData = phenoData, design = ~ sex)

# estimate the library size correction and save the normalized counts matrix
dds <- estimateSizeFactors(dds)
norm.cts <- counts(dds, normalized=TRUE)

# we want a normalized count of at least 10 in 4 or more samples
dds <- estimateSizeFactors(dds)
nc <- counts(dds, normalized=TRUE)
filter <- rowSums(norm.cts >= 10) >= 4
norm.cts <- norm.cts[filter,]
counts <- counts[filter,]

# run SVA
library(sva)

mm <- model.matrix(~ sex, colData(dds))
mm0 <- model.matrix(~ 1, colData(dds))

svaobj <- svaseq(norm.cts, mod=mm, mod0=mm0)
colnames(svaobj$sv)=paste0('SV_',1:ncol(svaobj$sv))

phenoData=cbind(phenoData, svaobj$sv)

## differetial expression analysis

# make DESeq2 object
dds=DESeqDataSetFromMatrix(countData = counts, colData = phenoData, 
	design = ~ SV_1 + SV_2 + SV_3 + SV_4 + SV_5 + SV_6 + SV_7 + SV_8 + SV_9 + SV_10 + SV_11 + 
	SV_12 + SV_13 + SV_14 + SV_15 + SV_16 + sex)

# run DE analysis
#dds=DESeq(dds)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds, maxit=1000)

# 42 rows did not converge in beta
# omit rows that did not converge in beta (these are typically genes with very small counts and little power)
# see https://support.bioconductor.org/p/65091/
ddsClean <- dds[which(mcols(dds)$betaConv),]

# results
rr=results(ddsClean, alpha=0.1)
summary(rr)
# out of 24043 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 112, 0.47%
# LFC < 0 (down)     : 91, 0.38%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

load("/athena/masonlab/scratch/users/nai2008/genes.rda")
mm=match(rownames(rr), genes$Gene.Name)
rr$chr=genes$chr[mm]

save(rr, file='/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/rdas/cases_0to21_DEseq2_results.rda')

# print results
sig_results=as.data.frame(rr[which(rr$padj<=0.1),])

# order results by LFC
oo=order(sig_results$log2FoldChange)
sig_results=sig_results[oo,]

# add chr
mmg=match(rownames(sig_results),genes$Gene.Name)
chr=genes$chr[mmg]

# make output table
out=data.frame(gene=rownames(sig_results), chr=chr, log2FoldChange=sig_results$log2FoldChange, FDR=sig_results$padj) 

out_autosomal=out[! out$chr %in% c('X','Y'),]

nrow(out)
nrow(out_autosomal)
# 203 DE genes; 147 of them on autosomes

save(out_autosomal, file='/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/rdas/HGG_0to21_DE_autosomal_genes.rda')

# enrichr
write.csv(out_autosomal,file="/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/data_tables/DE_genes_HGG_samples_peds.csv", row.names=FALSE)

######### correlation b/w DNAm and gene expression

library(org.Hs.eg.db)
SYMBOL2EG=unlist(as.list(org.Hs.egSYMBOL2EG)) # Map between Entrez Gene Identifiers and Gene Symbols
SYMBOL2EG_df=data.frame(symbol=names(SYMBOL2EG), entrez=as.vector(SYMBOL2EG))

mm=match(out_autosomal$gene, SYMBOL2EG_df$symbol)
out_autosomal=out_autosomal[-which(is.na(mm)),]
mm=match(out_autosomal$gene, SYMBOL2EG_df$symbol)
out_autosomal$entrez.Gene.ID=SYMBOL2EG_df$entrez[mm]

## match gene expression to DNAm [DMPs]
load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_0to21_cases.rda")
DMPs_0to21_cases=all_dmps
dmps_overlapping_genes=DMPs_0to21_cases[! as.vector(DMPs_0to21_cases$region) %in% c('upstream','downstream'),]

mm=match(out_autosomal$entrez.Gene.ID, dmps_overlapping_genes$Geneid)
out_autosomal$gene[which(!is.na(mm))] # None of the DE genes are overlapping or within 2.5kb of a DMP

## match gene expression to DNAm [DMRs]

# there were no sex DMRs for 0-21 year old cases

# volcano plot
#library(EnhancedVolcano)

rr_autosomal = rr[! rr$chr %in% c('X','Y'),]

pdf('/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/pdfs/volcano_plot_HGG_0to21.pdf')

source('/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/EnahancedVolcano_custom.R')

a='NS (FDR > 0.1)'
b=as.expression(bquote(log[2]~"FC"~">"~"|2|"))
c=as.expression(bquote('FDR'<='0.1'))
d=as.expression(bquote(log[2]~"FC"~">"~"|2|"~"&"~'FDR'<='0.1'))

EnhancedVolcano_custom(toptable=as.data.frame(rr_autosomal[,c(2,6)]),
	lab = rownames(rr_autosomal),
	x = 'log2FoldChange',
	y = 'padj',
	ylab =  bquote(~-Log[10]~italic(FDR)),
	pCutoff = 0.1,
	FCcutoff=2,
	title='pHGG samples (0-21 y.o.)',
	subtitle='Genes differentially expressed by sex',
	caption=NULL,
	legendLabels=c(a,b,c,d),
	legendPosition="bottom",
	border='full',
	col=c('grey30', 'forestgreen', 'royalblue', 'red2'),
	colAlpha=1,
	xlim=c(-6.25,6.25),
	ylim=c(0,5.5),
	shape=21,
	labhjust=0,
	legendLabSize=10,
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
## LIBD control data
###################################################################################################

# load data
library(SummarizedExperiment)

load("/athena/masonlab/scratch/users/nai2008/RNAseq/rse_gene_unfiltered.Rdata") #rse_gene
#load("/athena/masonlab/scratch/users/nai2008/RNAseq/rse_tx_unfiltered.Rdata") #rse_tx

all(colData(rse_gene)$RNum==colnames(assays(rse_gene)$counts)) #TRUE
all(colData(rse_gene)$RNum==colnames(assays(rse_gene)$rpkm)) #TRUE

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

#make rpkm matrix
rpkm=assays(rse_gene)$rpkm

# drop non-controls
pd=pd[which(pd$Dx=='Control'),]

# drop samples >21 yo
pd=pd[pd$Age<=21,]

# drop fetal samples
pd=pd[pd$Age>0,]

#
mm=match(pd$RNum,colnames(counts))
counts=counts[,mm]
rpkm=rpkm[,mm]

sex=pd$Sex
	sex[which(sex=='M')]='Male'
	sex[which(sex=='F')]='Female'

phenoData=data.frame(age=pd$Age,
sex=factor(as.vector(sex)),
Dx='Control',
Study_ID='LIBD',
my_unique_ID=paste0('LIBD_', 1:nrow(pd)),
RNum=pd$RNum)

phenoData=cbind(phenoData,pd[,55:62])

all(colnames(counts)==phenoData$RNum) #TRUE
all(colnames(rpkm)==phenoData$RNum) #TRUE

colnames(counts)=as.vector(phenoData$my_unique_ID)
colnames(rpkm)=as.vector(phenoData$my_unique_ID)

all(rownames(rpkm)==rownames(counts)) # TRUE

# change name of genes from ensembl to symbol
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')
tt=ss(rownames(counts),"\\.",1)
rownames(rpkm)=tt
rownames(counts)=tt

load("/athena/masonlab/scratch/users/nai2008/genes.rda")
mm=match(rownames(counts), genes$Ensembl.ID)
counts=counts[-which(is.na(mm)),]
rpkm=rpkm[-which(is.na(mm)),]

mm=match(rownames(counts), genes$Ensembl.ID)
rownames(counts)=genes$Gene.Name[mm]
rownames(rpkm)=genes$Gene.Name[mm]

# add duplicates together (for rpkm matrix)
length( which( duplicated(rownames(rpkm)) == TRUE ) ) #61

index_of_duplictes=which(duplicated(rownames(rpkm))==TRUE)
dup_genes=rownames(rpkm)[index_of_duplictes]

for (i in 1:length(dup_genes)){

dup=which(rownames(rpkm)==dup_genes[i])
if (length(dup)>1){
	colsum=as.vector(colSums(rpkm[dup,]))
	cs=matrix(colsum, nrow=1, ncol=length(colsum))
	rownames(cs)=dup_genes[i]
	rpkm=rpkm[-dup,]
	rpkm=rbind(rpkm,cs)
	}
}

length(which(duplicated(rownames(rpkm))==TRUE)) #0

# add duplicates together (for counts matrix)
length( which( duplicated(rownames(counts)) == TRUE ) ) #61

index_of_duplictes=which(duplicated(rownames(counts))==TRUE)
dup_genes=rownames(counts)[index_of_duplictes]

for (i in 1:length(dup_genes)){

dup=which(rownames(counts)==dup_genes[i])
if (length(dup)>1){
	colsum=as.vector(colSums(counts[dup,]))
	cs=matrix(colsum, nrow=1, ncol=length(colsum))
	rownames(cs)=dup_genes[i]
	counts=counts[-dup,]
	counts=rbind(counts,cs)
	}
}

length(which(duplicated(rownames(counts))==TRUE)) #0

all(rownames(rpkm)==rownames(counts)) #TRUE

## check that sex is labeled correctly (by uzing RPKM)

# find genes on chrY
load("/athena/masonlab/scratch/users/nai2008/genes.rda")
y_genes=genes$Gene.Name[which(genes$chr=='Y')]
y_genes=y_genes[!duplicated(y_genes)]

mm=match(y_genes,rownames(rpkm))
y_genes=y_genes[-which(is.na(mm))]

mm=match(y_genes,rownames(rpkm))
rpkm_y=rpkm[mm,]

csy=colSums(rpkm_y)
sex=factor(as.vector(phenoData$sex), levels=c('Female','Male')) # Female = 0; Male=1;

color=sub('Female','pink',sex)
color=sub('Male','blue', color)

pdf('/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/pdfs/sex_verification_LIBD_data_RPKM.pdf')

barplot(csy, names.arg=names(csy), col=color, las=2, cex.names=.7, main='Sum of RPKM values for genes on chrY')
legend(x='topright',legend=c('Female','Male'), col=c('pink','blue'), pch=15, cex=.7, title="Labeled sex")

dev.off()

## check that sex is labeled correctly (by using TPM)

load("/athena/masonlab/scratch/users/nai2008/RNAseq/rse_tx_unfiltered.Rdata") #rse_tx
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')
mart_export=read.csv('/athena/masonlab/scratch/users/nai2008/mart_export_transcripts.csv')

all(colnames(assays(rse_gene)$rpkm)==colnames(assays(rse_tx)$tpm)) #TRUE

tpm=assays(rse_tx)$tpm

# mm=match(rownames(tpm), mart_export$Transcript.stable.ID.version)
# length(which(is.na(mm))) #22145

# mm=match(ss(rownames(tpm),'\\.',1),mart_export$Transcript.stable.ID)
# length(which(is.na(mm))) #1389

rownames(tpm)=ss(rownames(tpm),'\\.',1)

mm=match(rownames(tpm), mart_export$Transcript.stable.ID)
tpm=tpm[-which(is.na(mm)),]

mm=match(rownames(tpm), mart_export$Transcript.stable.ID)
chr=as.vector(mart_export$Chromosome.scaffold.name[mm])

tpm_y=tpm[which(chr=='Y'),]

mm=match(phenoData$RNum,colnames(tpm_y))
tpm_y=tpm_y[,mm]

all(phenoData$RNum==colnames(tpm_y)) #TRUE

colnames(tpm_y)=phenoData$my_unique_ID

csy=colSums(tpm_y)
sex=factor(as.vector(phenoData$sex), levels=c('Female','Male')) # Female = 0; Male=1;

color=sub('Female','pink',sex)
color=sub('Male','blue', color)

pdf('/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/pdfs/sex_verification_LIBD_data_TPM.pdf')

barplot(csy, names.arg=names(csy), col=color, las=2, cex.names=.7, main='Sum of TPM values for genes on chrY')
legend(x='topright',legend=c('Female','Male'), col=c('pink','blue'), pch=15, cex=.7, title="Labeled sex")

dev.off()


save(counts, phenoData, file="/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/rdas/LIBD_counts_cleaned.rda")

########################################
## Differential expresson: controls 0-21
########################################
library(DESeq2)

load("/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/rdas/LIBD_counts_cleaned.rda")
load("/athena/masonlab/scratch/users/nai2008/genes.rda")

# make DESeq2 object
dds=DESeqDataSetFromMatrix(countData = counts, colData = phenoData, design = ~ sex)

# estimate the library size correction and save the normalized counts matrix
dds <- estimateSizeFactors(dds)
norm.cts <- counts(dds, normalized=TRUE)

# we want a normalized count of at least 10 in 4 or more samples
dds <- estimateSizeFactors(dds)
nc <- counts(dds, normalized=TRUE)
filter <- rowSums(norm.cts >= 10) >= 4
norm.cts <- norm.cts[filter,]
counts <- counts[filter,]

# run SVA
library(sva)

mm <- model.matrix(~ sex, colData(dds))
mm0 <- model.matrix(~ 1, colData(dds))

svaobj <- svaseq(norm.cts, mod=mm, mod0=mm0)
colnames(svaobj$sv)=paste0('SV_',1:ncol(svaobj$sv))

phenoData=cbind(phenoData, svaobj$sv)

## differetial expression analysis

# make DESeq2 object
dds=DESeqDataSetFromMatrix(countData = counts, colData = phenoData, 
	design = ~ SV_1 + SV_2 + SV_3 + SV_4 + SV_5 + SV_6 + SV_7 + SV_8 + SV_9 + SV_10 + SV_11 + 
	SV_12 + SV_13 + SV_14 + SV_15 + sex)

# run DE analysis
#dds=DESeq(dds)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds, maxit=500)

# 2 rows did not converge in beta
# omit rows that did not converge in beta (these are typically genes with very small counts and little power)
# see https://support.bioconductor.org/p/65091/
ddsClean <- dds[which(mcols(dds)$betaConv),]

# results
rr=results(ddsClean, alpha=0.1)
summary(rr)
# out of 27265 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 276, 1%
# LFC < 0 (down)     : 249, 0.91%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

mm=match(rownames(rr), genes$Gene.Name)
rr$chr=genes$chr[mm]

save(rr,file='/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/rdas/controls_0to21_DEseq2_results.rda')

# print results
sig_results=as.data.frame(rr[which(rr$padj<=0.1),])

# order results by LFC
oo=order(sig_results$log2FoldChange)
sig_results=sig_results[oo,]

# make output table
out=data.frame(gene=rownames(sig_results), chr=sig_results$chr, log2FoldChange=sig_results$log2FoldChange, FDR=sig_results$padj) 
out_autosomal=out[! out$chr %in% c('X','Y'),]

nrow(out)
nrow(out_autosomal)
# 525 DE genes; 440 of them on autosomes

save(out_autosomal, file='/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/rdas/controls_0to21_DE_autosomal_genes.rda')
#load('/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/rdas/controls_0to21_DE_autosomal_genes.rda')

# enrichr
write.csv(out_autosomal,file="/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/data_tables/DE_genes_control_samples_peds.csv", row.names=FALSE)

######### correlation b/w DNAm and gene expression

library(org.Hs.eg.db)
SYMBOL2EG=unlist(as.list(org.Hs.egSYMBOL2EG)) # Map between Entrez Gene Identifiers and Gene Symbols
SYMBOL2EG_df=data.frame(symbol=names(SYMBOL2EG), entrez=as.vector(SYMBOL2EG))

mm=match(out_autosomal$gene, SYMBOL2EG_df$symbol)
out_autosomal=out_autosomal[-which(is.na(mm)),]

mm=match(out_autosomal$gene, SYMBOL2EG_df$symbol)
out_autosomal$entrez.Gene.ID=as.numeric(as.vector(SYMBOL2EG_df$entrez[mm]))

## match gene expression to DNAm [DMPs]
load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_0to21_controls.rda")
DMPs_0to21_controls=all_dmps
dmps_overlapping_genes=DMPs_0to21_controls[! as.vector(DMPs_0to21_controls$region) %in% c('upstream','downstream'),]
dmps_overlapping_genes$Geneid=as.numeric(as.vector(dmps_overlapping_genes$Geneid))

mm=match(out_autosomal$entrez.Gene.ID, dmps_overlapping_genes$Geneid)
length(out_autosomal$gene[which(!is.na(mm))]) # 13 of the DE genes are overlapping ot are within 2.5kb of at least one DMP

if (length(out_autosomal$gene[which(!is.na(mm))]) != 0) {

	gg=out_autosomal$entrez.Gene.ID[which(!is.na(mm))]
	n=which(!is.na(mm))

	index=which(dmps_overlapping_genes$Geneid == gg[1])
	meth_exprs_table=dmps_overlapping_genes[index,][,c(1,3,5,7:21)]

	for (i in 2:length(gg)){
		index=which(dmps_overlapping_genes$Geneid == gg[i])
		aa_meth=dmps_overlapping_genes[index,][,c(1,3,5,7:21)]

		meth_exprs_table=rbind(meth_exprs_table,aa_meth)

	}

}

mm=match(meth_exprs_table$Geneid,out_autosomal$entrez.Gene.ID)

if (all(as.vector(meth_exprs_table$name)==as.vector(out_autosomal$gene[mm]))) {	
	meth_exprs_table=cbind(meth_exprs_table, out_autosomal[mm,c(1,3,4)])
	meth_exprs_table=meth_exprs_table[,c(6,18,1,5,4,2,3,8,10,9,11,12,13,16,20,21)]
	colnames(meth_exprs_table)[6]='DNAm_avg_change'
	colnames(meth_exprs_table)[1]='gene'
	colnames(meth_exprs_table)[7]='DNAm_FDR'
	colnames(meth_exprs_table)[15]='gene_expression_log2FC'
	colnames(meth_exprs_table)[16]='gene_expression_FDR'
	colnames(meth_exprs_table)[2]='Entrez_GeneID'
}

write.csv(meth_exprs_table,file="/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/data_tables/controls_meth_exprs_table.csv", row.names=FALSE)


## match gene expression to DNAm [DMRs]

# the DMRs for 0-21 controls don't overlap any genes

# volcano plot
#library(EnhancedVolcano)

rr_autosomal = rr[! rr$chr %in% c('X','Y'),]

pdf('/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/pdfs/volcano_plot_controls_0to21.pdf')

source('/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/EnahancedVolcano_custom.R')

a='NS (FDR > 0.1)'
b=as.expression(bquote(log[2]~"FC"~">"~"|2|"))
c=as.expression(bquote('FDR'<='0.1'))
d=as.expression(bquote(log[2]~"FC"~">"~"|2|"~"&"~'FDR'<='0.1'))

EnhancedVolcano_custom(toptable=as.data.frame(rr_autosomal[,c(2,6)]),
	lab = rownames(rr_autosomal),
	x = 'log2FoldChange',
	y = 'padj',
	ylab =  bquote(~-Log[10]~italic(FDR)),
	pCutoff = 0.1,
	FCcutoff=2,
	title='Controls (0-21 y.o.)',
	subtitle='Genes differentially expressed by sex',
	caption=NULL,
	legendLabels=c(a,b,c,d),
	legendPosition="bottom",
	border='full',
	col=c('grey30', 'forestgreen', 'royalblue', 'red2'),
	colAlpha=1,
	xlim=c(-6.25,6.25),
	ylim=c(0,11),
	shape=21,
	labhjust=0,
	legendLabSize=10,
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
## Downstream anlaysis
###################################################################################################

## is there overlap b/w genes differentially expressed by sex in cases vs controls

load('/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/rdas/controls_0to21_DE_autosomal_genes.rda')
controls_DE=out_autosomal

load('/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/rdas/HGG_0to21_DE_autosomal_genes.rda')
HGG_DE=out_autosomal

mm=match(HGG_DE$gene, controls_DE$gene) 

length(which(!is.na(mm))) # 4
HGG_DE$gene[which(!is.na(mm))] # ZNF883 FAR2P2 ITIH5 DEFA3


## IGSF3 exloration for Dr. Greenfield

load('/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/rdas/cases_0to21_DEseq2_results.rda')
rr_cases=rr

print=rr_cases[grep('IGSF',rownames(rr_cases)),]
write.csv(print[,c(2,6,7)],file="/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/data_tables/IGSF_cases_exploration.csv")

load('/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/rdas/controls_0to21_DEseq2_results.rda')
rr_controls=rr

print=rr_controls[grep('IGSF',rownames(rr_controls)),]
write.csv(print[,c(2,6,7)],file="/athena/masonlab/scratch/users/nai2008/RNAseq/second_round_of_scripts/data_tables/IGSF_controls_exploration.csv")

## How many DE genes have |log2FC| >2 in the control dataset?



## How many DE genes have |log2FC| >2 in the pHGG dataset?


## How many imprinted genes are DE in the control dataset? Are they also DE in the pHGG dataset?


