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

# genes & GenomicState
# txdb=gencode_txdb(version = "35", genome = "hg19", chrs = paste0("chr", c(seq_len(22), "X", "Y", "M")))
# genomicState=gencode_genomic_state(txdb)
# genes=gencode_annotated_genes(txdb)

# save(txdb, genes, genomicState, file='/athena/masonlab/scratch/users/nai2008/gencode_v35_hg19_genes_and_genomicState.rda')
load('/athena/masonlab/scratch/users/nai2008/gencode_v35_hg19_genes_and_genomicState.rda')

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

#### load & stratify the dataset:
load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")

#remove samples with a a mutant form of IDH:
GRset=GRset[,-grep('IDH|MUT|Mutation', pData(GRset)$IDH_status)]

## CASES

GRset_cases=GRset[,which(pData(GRset)$Dx_simplified=='HGG')]

# prepubescent (cases) [girls <9.5; boys <10.5]
indexes=c(which(pData(GRset_cases)$age < 9.5 & pData(GRset_cases)$sex=='Female'), 
	which(pData(GRset_cases)$age < 10.5 & pData(GRset_cases)$sex=='Male'))
GRset_prepubescent_cases=GRset_cases[,indexes]

# postpubescent pediatric (cases) [girls >= 9.5 & <=21; boys >= 10.5 & <=21]
indexes=c(which(pData(GRset_cases)$age >= 9.5 & pData(GRset_cases)$age <= 21 & pData(GRset_cases)$sex=='Female'), 
	which(pData(GRset_cases)$age >= 10.5 & pData(GRset_cases)$age <= 21 & pData(GRset_cases)$sex=='Male'))
GRset_postpubescent_cases=GRset_cases[,indexes]

# >21 (cases)
GRset_adult_cases=GRset_cases[,which(pData(GRset_cases)$age >21)]

# 0-21 (cases)
GRset_allPeds_cases=GRset_cases[,which(pData(GRset_cases)$age <=21)]

## CONTROLS

GRset_controls=GRset[,which(pData(GRset)$Dx_simplified=='Control')]

# prepubescent (controls) [girls <9.5; boys <10.5]
indexes=c(which(pData(GRset_controls)$age < 9.5 & pData(GRset_controls)$sex=='Female'), 
	which(pData(GRset_controls)$age < 10.5 & pData(GRset_controls)$sex=='Male'))
GRset_prepubescent_controls=GRset_controls[,indexes]

# postpubescent pediatric (controls) [girls >= 9.5 & <=21; boys >= 10.5 & <=21]
indexes=c(which(pData(GRset_controls)$age >= 9.5 & pData(GRset_controls)$age <= 21 & pData(GRset_controls)$sex=='Female'), 
	which(pData(GRset_controls)$age >= 10.5 & pData(GRset_controls)$age <= 21 & pData(GRset_controls)$sex=='Male'))
GRset_postpubescent_controls=GRset_controls[,indexes]

# >21 (controls)
GRset_adult_controls=GRset_controls[,which(pData(GRset_controls)$age >21)]

# 0-21 (controls)
GRset_allPeds_controls=GRset_controls[,which(pData(GRset_controls)$age <=21)]

#####################################
## Function
#####################################

DMR_downstream_analysis = function(GRset, DMR_table_filename, DMR_plots_filename, bumps) {

annotation=getAnnotation(GRset)
annotation=as.data.frame(annotation)

sigDMRs=bumps$table[bumps$table$fwer<=0.1,]

#length of the DMRs in genomic coordinates:
gl=sigDMRs$end-sigDMRs$start+1

# make output table with DMR info, and all the genes they overlap
out_table=data.frame(chr=sigDMRs$chr, start=sigDMRs$start, end=sigDMRs$end, L=sigDMRs$L, genomic_length=gl, p.value=sigDMRs$p.value, fwer=sigDMRs$fwer, proximal_genes=NA, Ensembl_GeneID=NA, distance=NA)

gr=GRanges(seqnames=sigDMRs$chr, 
	ranges=IRanges(sigDMRs$start, sigDMRs$end))

oo=findOverlaps(gr, genes, maxgap=5000, ignore.strand=TRUE)
oo=as.data.frame(oo)

oo$dist=distance(gr[oo$queryHits],genes[oo$subjectHits],ignore.strand=TRUE)

for (i in 1:nrow(sigDMRs)){

	ol=which(oo$queryHits==i)

	if( length(ol) != 0 ){
		index=oo$subjectHits[ol]
		out_table$proximal_genes[i]=paste(as.vector(genes$Gene[index,]), collapse="; ")
		out_table$Ensembl_GeneID[i]=paste(as.vector(genes$Geneid[index,]), collapse="; ")
		out_table$distance[i]=paste(oo$dist[ol], collapse="; ")
	} else {
		out_table$proximal_genes[i]='No proximal genes'
	}

}

write.csv(out_table, DMR_table_filename, row.names=FALSE)

info=list()
olg=as.vector(genes$Geneid[oo$subjectHits])

# how many overlapping genes are imprinted genes?

mm_ig=match(olg,ig$Ensembl.ID)

if(length(which(!is.na(mm_ig))) !=0 ){

	if (length(which(is.na(mm_ig))) != 0) { mm_ig=mm_ig[-which(is.na(mm_ig))] }
	
	info[[1]]=as.vector(ig$Gene[mm_ig])

} else { info[[1]] = 'None of the DMRs overlap or are within 5kb of imprinted genes'}

names(info)[1]='Imprinted gene(s)'

# how many pverlapping genes are TFs?

mm_tf=match(olg,TFs$Ensembl_ID)

if(length(which(!is.na(mm_tf))) !=0 ){

	if ( length(which(is.na(mm_tf))) != 0 ) { mm_tf=mm_tf[-which(is.na(mm_tf))] }
	
	info[[2]]=as.vector(TFs$HGNC_symbol[mm_tf])

} else { info[[2]] = 'None of the DMRs overlap or are within 5kb of TFs'}

names(info)[2]='Transcription factors'

info[[3]]=out_table
names(info)[3]='Sig DMRs info'

# match DMRs to genes (this is for the plot function)

gr=GRanges(seqnames=sigDMRs$chr, 
	ranges=IRanges(sigDMRs$start, sigDMRs$end))

load('/athena/masonlab/scratch/users/nai2008/gencode_v35_hg19_genes_and_genomicState.rda')

match_DMRs_to_genes=matchGenes(gr, genes)

for (i in 1:nrow(match_DMRs_to_genes)){
	if(is.na(match_DMRs_to_genes$name[i])) { match_DMRs_to_genes$name[i]=as.vector(match_DMRs_to_genes$Geneid[i]) } 
}

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

return(info)

}

#####################################
## Prepubescent (cases)
#####################################

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/DMRs_prepubescent_cases_1000bootstraps.rda") #bumps

length(which(bumps$table$fwer<=0.1)) # 5 sig DMRs

DMR_downstream_analysis(GRset=GRset_prepubescent_cases, 
	DMR_table_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/DMR_table_prepubescent_cases.csv",
	DMR_plots_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/DMRs_prepubescent_cases.pdf",
	bumps=bumps)

# $`Imprinted gene(s)`
# [1] "None of the DMRs overlap or are within 5kb of imprinted genes"

# $`Transcription factors`
# [1] "None of the DMRs overlap or are within 5kb of TFs"

# $`Sig DMRs info`
#    chr     start       end  L genomic_length      p.value  fwer
# 1 chr1 151810485 151811364  8            880 1.130974e-05 0.011
# 2 chr6  30038859  30039600 35            742 5.045883e-05 0.051
# 3 chr6  30651825  30653799 27           1975 5.654869e-05 0.054
# 4 chr5   8457548   8458392  7            845 9.917770e-05 0.057
# 5 chr7   8481994   8483710 11           1717 1.374568e-04 0.099

#####################################
## Postpubescent (cases)
#####################################

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/DMRs_postpubescent_cases_1000bootstraps.rda")
length(which(bumps$table$fwer<=0.1)) # 0 sig DMRs

#####################################
## >21 (cases)
#####################################

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/DMRs_adult_cases_1000bootstraps.rda")
length(which(bumps$table$fwer<=0.1)) # 2 sig DMRs

DMR_downstream_analysis(GRset=GRset_adult_cases, 
	DMR_table_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/DMR_table_adult_cases.csv",
	DMR_plots_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/DMRs_adult_cases.pdf",
	bumps=bumps)

# $`Imprinted gene(s)`
# [1] "None of the DMRs overlap or are within 5kb of imprinted genes"

# $`Transcription factors`
# [1] "ZSCAN1"

# $`Sig DMRs info`
#     chr    start      end  L genomic_length     p.value  fwer proximal_genes
# 1 chr19 58545001 58545837 11            837 0.003413123 0.065         ZSCAN1
# 2 chr17 73286267 73286282  2             16 0.004490951 0.080       SLC25A19

#####################################
## 0-21 (cases)
#####################################

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/DMRs_allPeds_cases_1000bootstraps.rda")
length(which(bumps$table$fwer<=0.1)) # 1 sig DMRs

DMR_downstream_analysis(GRset=GRset_allPeds_cases, 
	DMR_table_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/DMR_table_allPeds_cases.csv",
	DMR_plots_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/DMRs_allPeds_cases.pdf",
	bumps=bumps)

# $`Imprinted gene(s)`
# [1] "VTRNA2-1"

# $`Transcription factors`
# [1] "None of the DMRs overlap or are within 5kb of TFs"

# $`Sig DMRs info`
#    chr     start       end  L genomic_length      p.value  fwer proximal_genes
# 1 chr5 135415693 135416613 16            921 0.0003906962 0.034             NA
#    Ensembl_GeneID
# 1 ENSG00000270123

#####################################
## Prepubescent (controls)
#####################################

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/DMRs_prepubescent_controls_1000bootstraps.rda")
length(which(bumps$table$fwer<=0.1)) # 0 sig DMRs

#####################################
## Postpubescent (controls)
#####################################

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/DMRs_postpubescent_controls_1000bootstraps.rda")
length(which(bumps$table$fwer<=0.1)) # 1 sig DMRs

DMR_downstream_analysis(GRset=GRset_postpubescent_controls, 
	DMR_table_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/DMR_table_postpubescent_controls.csv",
	DMR_plots_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/DMRs_postpubescent_controls.pdf",
	bumps=bumps)

# $`Imprinted gene(s)`
# [1] "None of the DMRs overlap or are within 5kb of imprinted genes"

# $`Transcription factors`
# [1] "None of the DMRs overlap or are within 5kb of TFs"

# $`Sig DMRs info`
#     chr    start      end L genomic_length     p.value  fwer proximal_genes
# 1 chr15 85203259 85203259 1              1 0.008151338 0.051            NMB
#    Ensembl_GeneID
# 1 ENSG00000197696

#####################################
## >21 (controls)
#####################################

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/DMRs_adult_controls_1000bootstraps.rda")
length(which(bumps$table$fwer<=0.1)) # 1 sig DMRs

DMR_downstream_analysis(GRset=GRset_adult_controls, 
	DMR_table_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/DMR_table_adult_controls.csv",
	DMR_plots_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/DMRs_adult_controls.pdf",
	bumps=bumps)

# $`Imprinted gene(s)`
# [1] "None of the DMRs overlap or are within 5kb of imprinted genes"

# $`Transcription factors`
# [1] "None of the DMRs overlap or are within 5kb of TFs"

# $`Sig DMRs info`
#     chr    start      end L genomic_length   p.value  fwer
# 1 chr17 45401833 45401833 1              1 0.7058824 0.047
#          proximal_genes                                    Ensembl_GeneID
# 1 EFCAB13; NA; THCAT158 ENSG00000178852; ENSG00000259753; ENSG00000263293

#####################################
## 0-21 (controls)
#####################################

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/DMRs_allPeds_controls_1000bootstraps.rda")
length(which(bumps$table$fwer<=0.1)) # 1 sig DMRs

DMR_downstream_analysis(GRset=GRset_allPeds_controls, 
	DMR_table_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/DMR_table_allPeds_controls.csv",
	DMR_plots_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/DMRs_allPeds_controls.pdf",
	bumps=bumps)

# $`Imprinted gene(s)`
# [1] "None of the DMRs overlap or are within 5kb of imprinted genes"

# $`Transcription factors`
# [1] "None of the DMRs overlap or are within 5kb of TFs"

# $`Sig DMRs info`
#    chr    start      end L genomic_length    p.value  fwer proximal_genes
# 1 chr6 28601365 28601443 7             79 0.05300353 0.084         NA; NA
#                     Ensembl_GeneID
# 1 ENSG00000271440; ENSG00000287279

################################################################################################################

## Are there overlaps in DMRs b/w cases & controls?

# Prepubescent
No

# Postpubescent
No

# 0-21
load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/DMRs_allPeds_cases_1000bootstraps.rda")
sigDMRs_cases=bumps$table[bumps$table$fwer<=0.1,]
gr_cases=GRanges(seqnames=sigDMRs_cases$chr, ranges=IRanges(sigDMRs_cases$start, sigDMRs_cases$end))

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/DMRs_allPeds_controls_1000bootstraps.rda")
sigDMRs_controls=bumps$table[bumps$table$fwer<=0.1,]
gr_controls=GRanges(seqnames=sigDMRs_controls$chr, ranges=IRanges(sigDMRs_controls$start, sigDMRs_controls$end))

oo=findOverlaps(gr_cases, gr_controls, maxgap=0, ignore.strand=TRUE)
# The 2 combined objects have no sequence levels in common.

# > 21 (adults)
load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/DMRs_adult_cases_1000bootstraps.rda")
sigDMRs_cases=bumps$table[bumps$table$fwer<=0.1,]
gr_cases=GRanges(seqnames=sigDMRs_cases$chr, ranges=IRanges(sigDMRs_cases$start, sigDMRs_cases$end))

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/DMRs_adult_controls_1000bootstraps.rda")
sigDMRs_controls=bumps$table[bumps$table$fwer<=0.1,]
gr_controls=GRanges(seqnames=sigDMRs_controls$chr, ranges=IRanges(sigDMRs_controls$start, sigDMRs_controls$end))

oo=findOverlaps(gr_cases, gr_controls, maxgap=0, ignore.strand=TRUE)
# No overlaps










