#
library(minfi)
library(bumphunter)
library(GenomicRanges)
library(rafalib)
source('/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/blockPlot_function.R')

# genes
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene, annotationPackage = 'org.Hs.eg.db', by='gene')

###################################
## 0-21 cases
###################################

blocks_table_filename="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/data_tables/blocks_table_0to21_cases.csv"
blocks_plots_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/pdfs/blocks_0to21_cases_500bootstraps.pdf'

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/blocks_0to21_cases_500bootstraps.rda") #blocks
load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")
GRset=GRset[,-grep('IDH|MUT|Mutation', pData(GRset)$IDH_status)] # remove samples with a a mutant form of IDH

# stratify
GRset=GRset[,which(pData(GRset)$age <=21 & pData(GRset)$Dx_simplified=='HGG')]
#

length(which(blocks$tab$fwer<=0.1))

sig_blocks=blocks$tab[which(blocks$tab$fwer<=0.1),]

#length of the blocks (in bp):
sig_blocks$end-sig_blocks$start+1 # 65085

# make output table with DMR info, and all the genes they overlap

out_table=data.frame(chr=sig_blocks$chr, start=sig_blocks$start, end=sig_blocks$end, L=sig_blocks$L, p.value=sig_blocks$p.value, fwer=sig_blocks$fwer, overlapping_genes=NA, Entrez_GeneID=NA)

gr=GRanges(seqnames=sig_blocks$chr, 
	ranges=IRanges(sig_blocks$start, sig_blocks$end))

oo=findOverlaps(gr, genes)
oo=as.data.frame(oo)

for (i in 1:nrow(sig_blocks)){

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

write.csv(out_table, blocks_table_filename, row.names=FALSE)


# plotting
collapsed_CpG_clusters=cpgCollapse(GRset, what='Beta')
cset=collapsed_CpG_clusters$object

blockPlot(cset=cset,
blocks450=blocks,
coi=pData(GRset)$sex,
N=length(which(blocks$tab$fwer<=0.1)),
blockname = "sex",
filename=blocks_plots_filename,
scale=10,
showMethPanel = TRUE,
showGenePanel=TRUE,
showDiffPanel=TRUE,
showCancerPanel = FALSE,
bty= "o" )

#####################################
## 0-21 (controls)
#####################################

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/blocks_0to21_controls_500bootstraps.rda")

length(which(blocks$tab$fwer<=0.1)) # 0 [NO SIG BLOCKS]

#####################################
## >21 (cases)
#####################################

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/blocks_adult_cases_500bootstraps.rda")

length(which(blocks$tab$fwer<=0.1)) # 0 [NO SIG BLOCKS]

#####################################
## >21 (controls)
#####################################

blocks_table_filename="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/data_tables/blocks_table_adults_controls.csv"
blocks_plots_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/pdfs/blocks_adults_controls_500bootstraps.pdf'

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/blocks_adult_controls_500bootstraps.rda")

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")
GRset=GRset[,-grep('IDH|MUT|Mutation', pData(GRset)$IDH_status)] # remove samples with a a mutant form of IDH

# stratify
GRset=GRset[,which(pData(GRset)$age >21 & pData(GRset)$Dx_simplified=='Control')]
#

length(which(blocks$tab$fwer<=0.1)) # 1

sig_blocks=blocks$tab[which(blocks$tab$fwer<=0.1),]

#length of the blocks (in bp):
as.integer(sig_blocks$end - sig_blocks$start) + 1 # 3397

# make output table with DMR info, and all the genes they overlap

out_table=data.frame(chr=sig_blocks$chr, start=sig_blocks$start, end=sig_blocks$end, L=sig_blocks$L, p.value=sig_blocks$p.value, fwer=sig_blocks$fwer, overlapping_genes=NA, Entrez_GeneID=NA)

gr=GRanges(seqnames=sig_blocks$chr, 
	ranges=IRanges(sig_blocks$start, sig_blocks$end))

oo=findOverlaps(gr, genes, maxgap=0)
oo=as.data.frame(oo)

for (i in 1:nrow(sig_blocks)){

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

write.csv(out_table, blocks_table_filename, row.names=FALSE)

# plotting
collapsed_CpG_clusters=cpgCollapse(GRset, what='Beta')
cset=collapsed_CpG_clusters$object

blockPlot(cset=cset,
blocks450=blocks,
coi=pData(GRset)$sex,
N=length(which(blocks$tab$fwer<=0.1)),
blockname = "sex",
filename=blocks_plots_filename,
scale=100,
showMethPanel = TRUE,
showGenePanel=TRUE,
showDiffPanel=TRUE,
showCancerPanel = FALSE,
bty= "o" )



















