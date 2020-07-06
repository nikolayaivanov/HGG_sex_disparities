##

# BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
# BiocManager::install("GenomicState")
# install.packages("devtools")
# install_github("ririzarr/rafalib")
# install.packages('RColorBrewer')

library(derfinder)
library(bumphunter)
library(BSgenome.Hsapiens.UCSC.hg19)
library(minfi)
library(devtools)
library(rafalib)
library(RColorBrewer)
library(GenomicState)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')
source('/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/dmrPlot_function.R')

# genes
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene, annotationPackage = 'org.Hs.eg.db', by='gene')


# GenomicState
# txdb=gencode_txdb(version = "31", genome = "hg19", chrs = paste0("chr", c(seq_len(22), "X", "Y", "M")))
# genomicState=gencode_genomic_state(txdb)

# save(txdb, genomicState, file='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/genomicState.rda')
load('/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/genomicState.rda')

#####################################
## >21 (cases)
#####################################

##
DMR_table_filename="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/data_tables/DMR_table_adult_cases.csv"
DMR_plots_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/pdfs/DMRs_adult_cases_500bootstraps.pdf'
load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/DMRs_adult_cases_500bootstraps.rda") #bumps
load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")

## stratify the dataset
GRset=GRset[,which(pData(GRset)$age >21 & pData(GRset)$Dx_simplified=='HGG')]
##

## remove samples with a a mutant form of IDH
GRset=GRset[,-grep('IDH|MUT|Mutation', pData(GRset)$IDH_status)]

annotation=getAnnotation(GRset)
annotation=as.data.frame(annotation)

length(which(bumps$table$fwer<=0.1)) # 2 sig DMRs

sigDMRs=bumps$table[bumps$table$fwer<=0.1,]

#length of the DMRs in genomic coordinates:
sigDMRs$end-sigDMRs$start+1 # 845 837

# make output table with DMR info, and all the genes they overlap

out_table=data.frame(chr=sigDMRs$chr, start=sigDMRs$start, end=sigDMRs$end, L=sigDMRs$L, p.value=sigDMRs$p.value, fwer=sigDMRs$fwer, overlapping_genes=NA, Entrez_GeneID=NA)

gr=GRanges(seqnames=sigDMRs$chr, 
	ranges=IRanges(sigDMRs$start, sigDMRs$end))

oo=findOverlaps(gr, genes, maxgap=0)
oo=as.data.frame(oo)

for (i in 1:nrow(sigDMRs)){

	ol=which(oo$queryHits==i)
	index=oo$subjectHits[ol]
	out_table$overlapping_genes[i]=paste(as.vector(genes$Gene[index,]), collapse="; ")
	out_table$Entrez_GeneID[i]=paste(as.vector(genes$Geneid[index,]), collapse="; ")

}

write.csv(out_table, DMR_table_filename, row.names=FALSE)

# match DMRs to genes (this is for the plot function)

gr=GRanges(seqnames=sigDMRs$chr, 
	ranges=IRanges(sigDMRs$start, sigDMRs$end))

match_DMRs_to_genes=matchGenes(gr, genes)

## plotting DMRs

pdf(DMR_plots_filename)

chr=as.vector(annotation$chr)
pos=as.vector(annotation$pos)
cluster=clusterMaker(chr,pos,maxGap=500)

dmrPlot(regions=bumps$tab[which(bumps$table$fwer <=.1),],
	p=getBeta(GRset),
	chr=as.vector(annotation$chr),
	pos=as.vector(annotation$pos),
	cluster=cluster,
	genes=match_DMRs_to_genes,
	coi=pData(GRset)$sex,
	build="hg19",
	number=length(which(bumps$table$fwer<=0.1)),
	gs=genomicState$fullGenome,
	species = "human",
	Jitter = TRUE,
	cols=c("red", "royalblue2"),
	lines = FALSE,
	linesSmooth = TRUE,
	title = TRUE,
	Legend = TRUE,
	colorRamp = FALSE, 
	meanSmooth=TRUE,
	plotCpG = TRUE,
	geneAnno = "gene" )

dev.off()


#####################################
## >21 (controls)
#####################################

##
DMR_table_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/data_tables/DMR_table_adult_controls.csv'
DMR_plots_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/pdfs/DMRs_adult_controls_500bootstraps.pdf'
load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/DMRs_adult_controls_500bootstraps.rda") #bumps
load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")

## stratify the dataset
GRset=GRset[,which(pData(GRset)$age >21 & pData(GRset)$Dx=='Control')]
##

annotation=getAnnotation(GRset)
annotation=as.data.frame(annotation)

length(which(bumps$table$fwer<=0.1)) # 2

sigDMRs=bumps$table[bumps$table$fwer<=0.1,]

#length of the DMRs in genomic coordinates:
sigDMRs$end-sigDMRs$start+1 # 1 1

# make output table with DMR info, and all the genes they overlap

out_table=data.frame(chr=sigDMRs$chr, start=sigDMRs$start, end=sigDMRs$end, L=sigDMRs$L, p.value=sigDMRs$p.value, fwer=sigDMRs$fwer, overlapping_genes=NA, Entrez_GeneID=NA)

gr=GRanges(seqnames=sigDMRs$chr, 
	ranges=IRanges(sigDMRs$start, sigDMRs$end))

oo=findOverlaps(gr, genes, maxgap=0)
oo=as.data.frame(oo)

for (i in 1:nrow(sigDMRs)){

	ol=which(oo$queryHits==i)
	index=oo$subjectHits[ol]

	if (length(index)==0){
		out_table$overlapping_genes[i]=NA
		out_table$Entrez_GeneID[i]=NA
	} else {
		out_table$overlapping_genes[i]=paste(as.vector(genes$Gene[index,]), collapse="; ")
		out_table$Entrez_GeneID[i]=paste(as.vector(genes$Geneid[index,]), collapse="; ")
	}

}

write.csv(out_table, DMR_table_filename, row.names=FALSE)

#match DMRs to genes (this is for the plot function)

gr=GRanges(seqnames=sigDMRs$chr, 
	ranges=IRanges(sigDMRs$start, sigDMRs$end))

match_DMRs_to_genes=matchGenes(gr, genes)

## plotting DMRs

pdf(DMR_plots_filename)

chr=as.vector(annotation$chr)
pos=as.vector(annotation$pos)
cluster=clusterMaker(chr,pos,maxGap=500)

dmrPlot(regions=bumps$tab[which(bumps$table$fwer <=.1),],
	p=getBeta(GRset),
	chr=as.vector(annotation$chr),
	pos=as.vector(annotation$pos),
	cluster=cluster,
	genes=match_DMRs_to_genes,
	coi=pData(GRset)$sex,
	build="hg19",
	number=length(which(bumps$table$fwer<=0.1)),
	gs=genomicState$fullGenome,
	species = "human",
	Jitter = TRUE,
	cols=c("red", "royalblue2"),
	lines = FALSE,
	linesSmooth = TRUE,
	title = TRUE,
	Legend = TRUE,
	colorRamp = FALSE, 
	meanSmooth=TRUE,
	plotCpG = TRUE,
	geneAnno = "gene" )

dev.off()

# 	regions=bumps$tab[which(bumps$table$fwer <=.1),]
# 	p=getBeta(GRset)
# 	chr=as.vector(annotation$chr)
# 	pos=as.vector(annotation$pos)
# 	cluster=cluster
# 	genes=match_DMRs_to_genes
# 	coi=pData(GRset)$sex
# 	build="hg19"
# 	number=length(which(bumps$table$fwer<=0.1))
# 	gs=genomicState$fullGenome
# 	species = "human"
# 	Jitter = TRUE
# 	cols=c("red", "royalblue2")
# 	lines = FALSE
# 	linesSmooth = TRUE
# 	title = TRUE
# 	Legend = TRUE
# 	colorRamp = FALSE
# 	meanSmooth=TRUE
# 	plotCpG = TRUE
# 	geneAnno = "gene"

#####################################
## 0-21 (cases)
#####################################

##
DMR_table_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/data_tables/DMR_table_0to21_cases.csv'
DMR_plots_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/pdfs/DMRs_0to21_cases_500bootstraps.pdf'
load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/DMRs_0to21_cases_500bootstraps.rda") #bumps
load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")

## stratify the dataset
GRset=GRset[,which(pData(GRset)$age <= 21 & pData(GRset)$Dx=='HGG')]
##

## remove samples with a a mutant form of IDH
GRset=GRset[,-grep('IDH|MUT|Mutation', pData(GRset)$IDH_status)]

annotation=getAnnotation(GRset)
annotation=as.data.frame(annotation)

length(which(bumps$table$fwer<=0.1)) # 0 (NO SIG DMRS) !!!!!!!

#####################################
## 0-21 (controls)
#####################################

##
DMR_table_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/data_tables/DMR_table_0to21_controls.csv'
DMR_plots_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/pdfs/DMRs_0to21_controls_500bootstraps.pdf'
load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/DMRs_0to21_controls_500bootstraps.rda") #bumps
load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")

## stratify the dataset
GRset=GRset[,which(pData(GRset)$age <= 21 & pData(GRset)$Dx=='Control')]
##

annotation=getAnnotation(GRset)
annotation=as.data.frame(annotation)

length(which(bumps$table$fwer<=0.1)) # 1

sigDMRs=bumps$table[bumps$table$fwer<=0.1,]

#length of the DMRs in genomic coordinates:
sigDMRs$end-sigDMRs$start+1 # 55

# make output table with DMR info, and all the genes they overlap

out_table=data.frame(chr=sigDMRs$chr, start=sigDMRs$start, end=sigDMRs$end, L=sigDMRs$L, p.value=sigDMRs$p.value, fwer=sigDMRs$fwer, overlapping_genes=NA, Entrez_GeneID=NA)

gr=GRanges(seqnames=sigDMRs$chr, 
	ranges=IRanges(sigDMRs$start, sigDMRs$end))

oo=findOverlaps(gr, genes, maxgap=0)
oo=as.data.frame(oo)

for (i in 1:nrow(sigDMRs)){

	ol=which(oo$queryHits==i)
	index=oo$subjectHits[ol]

	if (length(index)==0){
		out_table$overlapping_genes[i]=NA
		out_table$Entrez_GeneID[i]=NA
	} else {
		out_table$overlapping_genes[i]=paste(as.vector(genes$Gene[index,]), collapse="; ")
		out_table$Entrez_GeneID[i]=paste(as.vector(genes$Geneid[index,]), collapse="; ")
	}

}

write.csv(out_table, DMR_table_filename, row.names=FALSE)

#match DMRs to genes (this is for the plot function)

gr=GRanges(seqnames=sigDMRs$chr, 
	ranges=IRanges(sigDMRs$start, sigDMRs$end))

match_DMRs_to_genes=matchGenes(gr, genes)

## plotting DMRs

pdf(DMR_plots_filename)

chr=as.vector(annotation$chr)
pos=as.vector(annotation$pos)
cluster=clusterMaker(chr,pos,maxGap=500)

dmrPlot(regions=bumps$tab[which(bumps$table$fwer <=.1),],
	p=getBeta(GRset),
	chr=as.vector(annotation$chr),
	pos=as.vector(annotation$pos),
	cluster=cluster,
	genes=match_DMRs_to_genes,
	coi=pData(GRset)$sex,
	build="hg19",
	number=length(which(bumps$table$fwer<=0.1)),
	gs=genomicState$fullGenome,
	species = "human",
	Jitter = TRUE,
	cols=c("red", "royalblue2"),
	lines = FALSE,
	linesSmooth = TRUE,
	title = TRUE,
	Legend = TRUE,
	colorRamp = FALSE, 
	meanSmooth=TRUE,
	plotCpG = TRUE,
	geneAnno = "gene" )

dev.off()



