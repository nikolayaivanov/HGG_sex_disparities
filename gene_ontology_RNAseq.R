#
## Cavatica data (HGG cases 0-21 y.o.); Over Representation Analysis (ORA)
## (only autosomal genes)

library("KEGG.db")
KEGG.PATHID2NAME=unlist(as.list(KEGGPATHID2NAME))
KEGG.PATHID2NAME=data.frame(kegg.path.ID=names(KEGG.PATHID2NAME), kegg.path.name=as.vector(KEGG.PATHID2NAME))

require(goseq)

load('/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/rdas/cases_0to21_DEseq2_results.rda')
rr=rr[which(rr$chr %in% c('X','Y')),]

indicator=rep(0, times=nrow(rr))
indicator[which(rr$padj<=0.1)]=1
aa=indicator
names(aa)=rr$Ensembl

pwf=nullp(aa, "hg19", "ensGene")
GO.KEGG.wall=goseq(pwf,"hg19","ensGene",test.cats=c("GO:CC", "GO:BP", "GO:MF", "KEGG"))
GO.KEGG.wall$over_represented_FDR=p.adjust(GO.KEGG.wall$over_represented_pvalue, method="BH")
GO.KEGG.wall$ontology[which(is.na(GO.KEGG.wall$ontology))]='KEGG'

mm=match(as.vector(GO.KEGG.wall$category),as.vector(KEGG.PATHID2NAME$kegg.path.ID))

indexes=which(!is.na(mm))
mmm=mm[-which(is.na(mm))]

for (i in 1:length(indexes)) {
  GO.KEGG.wall$term[indexes[i]]=as.vector(KEGG.PATHID2NAME$kegg.path.name)[mmm[i]]
}

length(which(GO.KEGG.wall$over_represented_FDR<=0.1)) # 8
length(which(GO.KEGG.wall$over_represented_pvalue<=0.05)) # 82	

# from https://support.bioconductor.org/p/102273/
getGenes <- function(pwf, goterm, genome, ids){
    gene2cat <-  getgo(rownames(pwf), genome, ids, fetch.cats=c("GO:CC","GO:BP","GO:MF", "KEGG"))
    cat2gene <- split(rep(names(gene2cat), sapply(gene2cat, length)),
                      unlist(gene2cat, use.names = FALSE))
    out <- pwf[cat2gene[[goterm]],]
    out <- out[out$DEgenes > 0,]
    out
}

GO.KEGG.wall=GO.KEGG.wall[order(GO.KEGG.wall$over_represented_pvalue),]

nn =  length(which(GO.KEGG.wall$over_represented_pvalue<=0.05))

GO.KEGG.wall.out=GO.KEGG.wall[1:nn,]

GO.KEGG.wall.out$genes.symbols=NA
GO.KEGG.wall.out$genes.Ensembl=NA

library(EnsDb.Hsapiens.v86)
edb = EnsDb.Hsapiens.v86
genes=as.data.frame(genes(edb))

for (i in 1:nrow(GO.KEGG.wall.out)){

  aa=getGenes(pwf, GO.KEGG.wall.out$category[i], "hg19", "ensGene")
  GO.KEGG.wall.out$genes.symbols[i]=paste(genes$gene_name[match(rownames(aa), genes$gene_id)], collapse=", ")
  GO.KEGG.wall.out$genes.Ensembl[i]=paste(rownames(aa), collapse=", ")

}

op=which(GO.KEGG.wall.out$ontology=='KEGG')

for(i in 1:length(op)){
  GO.KEGG.wall.out$category[op[i]]=paste0('KEGG ',GO.KEGG.wall.out$category[op[i]])
}

# GO.KEGG.wall.out=GO.KEGG.wall.out[,c(6,7,1,2,8,9,10)]

write.csv(GO.KEGG.wall.out, file="/athena/masonlab/scratch/users/nai2008/RNAseq/third_round_of_scripts/data_tables/ORA_DE_genes_HGG_samples_0to21.csv", row.names=FALSE)



















