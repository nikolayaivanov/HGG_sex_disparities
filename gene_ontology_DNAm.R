## Determine the 'universe' - all genes that either overlap or are within 5kbs of the Illumina 450k array probes
library(minfi)
load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")
anno=getAnnotation(GRset)

anno.GRanges=GRanges(seqnames=anno$chr, ranges=IRanges(anno$pos, anno$pos))

load('/athena/masonlab/scratch/users/nai2008/gencode_v35_hg19_genes_and_genomicState.rda')

oo=as.data.frame(findOverlaps(anno.GRanges, genes, maxgap=5000, select='all', ignore.strand=TRUE))
hits=genes[oo$subjectHits,]
universe.450k=unique(as.vector(hits$Geneid))
length(universe.450k) # 43749

save(universe.450k, file="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/gene_universe_450k.rda")
# load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/gene_universe_450k.rda") # universe.450k

# how many probes is each gene associated with
mm=match(universe.450k, as.vector(genes$Geneid))
genes_u=genes[mm,]

oo=as.data.frame(findOverlaps(genes_u, anno.GRanges, maxgap=5000, select='all', ignore.strand=TRUE))

aa=universe.450k[oo$queryHits]

num_probes_per_gene=data.frame(universe.450k.genes=universe.450k, probes.per.gene=NA)

for (i in 1:length(universe.450k)) {
	num_probes_per_gene[i,2]=length(which(aa==universe.450k[i]))
}

length(which(is.na(num_probes_per_gene[,2]))) # 0

save(num_probes_per_gene, file="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/num_probes_per_gene.rda")
# load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/num_probes_per_gene.rda")

# KEGG
library("KEGG.db")
KEGG.PATHID2NAME=unlist(as.list(KEGGPATHID2NAME))
KEGG.PATHID2NAME=data.frame(kegg.path.ID=names(KEGG.PATHID2NAME), kegg.path.name=as.vector(KEGG.PATHID2NAME))

############ Drafts:

# From https://www.biostars.org/p/52101/

library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
#gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
gene.data <- getBM(attributes=c('ensembl_gene_id', 'go_id'), mart = ensembl, filters = 'go_parent_term', values = 'GO:0007507')



gene.data <- getBM(attributes=c('ensembl_gene_id', 'go_id'), mart = ensembl, filters = 'go', values = 'GO:0007507')



gene.data <- getBM(attributes=c('go_id','hgnc_symbol', 'ensembl_transcript_id'), filters = 'go', values = 'GO:0007507', mart = ensembl)

gene.data <- getBM(attributes=c('go_id'), 
	filters = 'go_id', values = 'GO:0007507', mart = ensembl)


tta=listAttributes(ensembl)

grep('GO term accession',tta$description, ignore.case=TRUE)

ttf=listFilters(ensembl)

grep('GO ID',ttf$description, ignore.case=TRUE)

grep('kegg',ttf$description, ignore.case=TRUE)

library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
gene.data <- getBM(attributes=c('hgnc_symbol', 
    'ensembl_transcript_id', 'go_id'), filters = 'go', 
    values = 'GO:0072599', mart = ensembl)











