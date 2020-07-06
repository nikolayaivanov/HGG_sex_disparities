library(org.Hs.eg.db)

#
ENSEMBL2EG=unlist(as.list(org.Hs.egENSEMBL2EG)) # Map Ensembl gene accession numbers with Entrez Gene identifiers
ENSEMBL2EG_df=data.frame(ensembl=names(ENSEMBL2EG), entrez=as.vector(ENSEMBL2EG))

#
SYMBOL2EG=unlist(as.list(org.Hs.egSYMBOL2EG)) # Map between Entrez Gene Identifiers and Gene Symbols
SYMBOL2EG_df=data.frame(symbol=names(SYMBOL2EG), entrez=as.vector(SYMBOL2EG))

#
EG2REFSEQ=unlist(as.list(org.Hs.egREFSEQ)) #Map between Entrez Gene Identifiers and RefSeq Identifiers
EG2REFSEQ_df=data.frame(entrez=names(EG2REFSEQ), refseq=as.vector(EG2REFSEQ))

#
EG2CHR=unlist(as.list(org.Hs.egCHR)) # Map Entrez Gene IDs to Chromosomes
EG2CHR_df=data.frame(chr=as.vector(EG2CHR), entrez=names(EG2CHR))

## make master df that has corresponding ensembl, entrez, refseq, gene symbol
master_gene_match_df=ENSEMBL2EG_df

m=match(master_gene_match_df$entrez,SYMBOL2EG_df$entrez)
SYMBOL2EG_df=SYMBOL2EG_df[m,]
master_gene_match_df$symbol=SYMBOL2EG_df$symbol

m=match(master_gene_match_df$entrez, EG2REFSEQ_df$entrez)
EG2REFSEQ_df=EG2REFSEQ_df[m,]
master_gene_match_df$refseq=EG2REFSEQ_df$refseq

m=match(master_gene_match_df$entrez, EG2CHR_df$entrez)
EG2CHR_df=EG2CHR_df[m,]
master_gene_match_df$chr=paste0('chr',EG2CHR_df$chr)

save(master_gene_match_df,file="/athena/masonlab/scratch/users/nai2008/master_gene_match_df.rda")
# load("/athena/masonlab/scratch/users/nai2008/master_gene_match_df.rda") # master_gene_match_df [ensembl; entrez; symbol; refseq; chr]

#####################

genes=read.csv('/athena/masonlab/scratch/users/nai2008/mart_export.txt', na.strings="")
genes=genes[which(genes$Chromosome.scaffold.name %in% c(1:22,'X','Y','MT')),]
colnames(genes)=c('Ensembl.ID','chr','Gene.Start','Gene.End','strand','Gene.Name','Gene.Synonym', 'Gene.type','Gene.description')

save(genes,file="/athena/masonlab/scratch/users/nai2008/genes.rda")
# load("/athena/masonlab/scratch/users/nai2008/genes.rda")
























