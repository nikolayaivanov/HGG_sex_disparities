##
library(minfi)
library(limma)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")

## remove samples with a a mutant form of IDH

GRset=GRset[,-grep('IDH|MUT|Mutation', pData(GRset)$IDH_status)] # 36/707 samples have an IDH mutation

#####################################
## stratify the dataset
#####################################

# <= 6 yo (cases)
GRset_prepubescent_cases=GRset[,which(pData(GRset)$age <=6 & pData(GRset)$Dx_simplified=='HGG')]

nrow(pData(GRset_prepubescent_cases)) # 25 samples
table(pData(GRset_prepubescent_cases)$sex) # 9F, 16M

# 6-21 (cases)
GRset_6to21_cases=GRset[,which(pData(GRset)$age >6 & pData(GRset)$age<=21 & pData(GRset)$Dx_simplified=='HGG')]

nrow(pData(GRset_6to21_cases)) # 111 samples
table(pData(GRset_6to21_cases)$sex) # 48F  63M

# >21 (cases)
GRset_adult_cases=GRset[,which(pData(GRset)$age >21 & pData(GRset)$Dx_simplified=='HGG')]

nrow(pData(GRset_adult_cases)) # 153 samples
table(pData(GRset_adult_cases)$sex) # 67F  86M

#0-21 (cases)
GRset_allPeds_cases=GRset[,which(pData(GRset)$age <=21 & pData(GRset)$Dx_simplified=='HGG')]

nrow(pData(GRset_allPeds_cases)) # 136 samples
table(pData(GRset_allPeds_cases)$sex) # 57F  79M

##

# <= 6 yo (controls)
GRset_prepubescent_controls=GRset[,which(pData(GRset)$age <=6 & pData(GRset)$Dx=='Control')]

nrow(pData(GRset_prepubescent_controls)) # 51 samples
table(pData(GRset_prepubescent_controls)$sex) # 11F  40M

# 6-21 (controls)
GRset_6to21_controls=GRset[,which(pData(GRset)$age >6 & pData(GRset)$age<=21 & pData(GRset)$Dx=='Control')]

nrow(pData(GRset_6to21_controls)) # 77 samples
table(pData(GRset_6to21_controls)$sex) # 24F  53M

# >21 (controls)
GRset_adult_controls=GRset[,which(pData(GRset)$age >21 & pData(GRset)$Dx=='Control')]

nrow(pData(GRset_adult_controls)) # 254 samples
table(pData(GRset_adult_controls)$sex) # 86F  168M

#0-21 (controls)
GRset_allPeds_controls=GRset[,which(pData(GRset)$age <=21 & pData(GRset)$Dx=='Control')]

nrow(pData(GRset_allPeds_controls)) # 128 samples
table(pData(GRset_allPeds_controls)$sex) # 35F  93M

##

# <= 6 yo (cases & controls)
GRset_prepubescent_cases_and_controls=GRset[,which(pData(GRset)$age <= 6)]

nrow(pData(GRset_prepubescent_cases_and_controls)) # 76 samples
table(pData(GRset_prepubescent_cases_and_controls)$sex) # 20F  56M

# 6-21 (cases & controls)
GRset_6to21_cases_and_controls=GRset[,which(pData(GRset)$age >6 & pData(GRset)$age<=21)]

nrow(pData(GRset_6to21_cases_and_controls)) # 188 samples
table(pData(GRset_6to21_cases_and_controls)$sex) # 72F  116M

# >21 (cases & controls)
GRset_adult_cases_and_controls=GRset[,which(pData(GRset)$age >21)]

nrow(pData(GRset_adult_cases_and_controls)) # 407 samples
table(pData(GRset_adult_cases_and_controls)$sex) # 153F  254M

# 0-21 (cases & controls)
GRset_0to21_cases_and_controls=GRset[,which(pData(GRset)$age <=21)]

nrow(pData(GRset_0to21_cases_and_controls)) # 264 samples
table(pData(GRset_0to21_cases_and_controls)$sex) # 92F  172M

#####################################
## function to carry out male-female DMP analysis
#####################################

DMP_analysis = function(GRset, pdf_filename, table_filename, DMPs_annotated_filename, genes) {

# find DMPs

pd=pData(GRset)
beta=getBeta(GRset)
neg_control_PCs=pd[,12:19]

sex=factor(as.vector(GRset$sex), levels=c('Female','Male')) # Female = 0; Male=1;
mod = model.matrix(~sex)
mod=cbind(mod, neg_control_PCs)

probe_fit=lmFit(object=beta,design=mod)
eb=eBayes(probe_fit)

probeFit=data.frame(probe=rownames(eb$p.value),
	intercept=probe_fit$coefficients[,1],
	slope=probe_fit$coefficients[,2], 
	p.value=eb$p.value[,2],
	fdr=p.adjust(as.vector(eb$p.value[,2]),method='fdr'),
	t=eb$t[,2])

rownames(probeFit)=NULL

o=order(probeFit$fdr)
probeFit=probeFit[o,]
beta=beta[o,]

num_dmps=length(which(probeFit$fdr<=0.05))

if (num_dmps == 0) {

return=c(num_dmps=0,
	dmps_that_overlap_promoters_or_enhancers=0,
	dmps_that_overlap_promoters_or_enhancers_or_genes=0 )

	} else {

percent_dmrs=length(which(probeFit$fdr<=0.05))/nrow(probeFit)

dmps=probeFit[which(probeFit$fdr<=0.05),]

# Map DMPs to genes

anno=getAnnotation(GRset)
mm=match(rownames(beta),rownames(anno))
anno=anno[mm,]

rti=as.data.frame(table(factor(anno$Relation_to_Island)))
num_shore_probes=sum(rti$Freq[grep('Shore', rti$Var1, ignore.case=TRUE)])
num_shelf_probes=sum(rti$Freq[grep('Shelf', rti$Var1, ignore.case=TRUE)])
num_island_probes=rti$Freq[grep('Island', rti$Var1, ignore.case=TRUE)]
num_opensea_probes=rti$Freq[grep('OpenSea', rti$Var1, ignore.case=TRUE)]
num_promoter_probes=length(grep('Promoter',anno$Regulatory_Feature_Group, , ignore.case=TRUE))
num_enhancer_probes=length(which(anno$Enhancer=="TRUE"))

anno_DMPs=anno[1:nrow(dmps),]

ad_df=as.data.frame(anno_DMPs)

gr=GRanges(seqnames=ad_df$chr, 
	ranges=IRanges(ad_df$pos, ad_df$pos),
	strand=ad_df$strand)

match_DMPs_to_genes=matchGenes(gr, genes)
match_DMPs_to_genes=match_DMPs_to_genes[,-c(2,15)]

# ol=findOverlaps(gr, genes, maxgap=5000)

# output table of DMPs with nearest genes

Enhancer=ad_df$Enhancer
Enhancer[which(Enhancer=="")]=FALSE

Promoter=rep('FALSE', times=nrow(dmps))
Promoter[grep('Promoter', ad_df$Regulatory_Feature_Group, ignore.case=TRUE)]='TRUE'

out=cbind(dmps, anno_DMPs[,1:2], match_DMPs_to_genes, Enhancer=Enhancer, Promoter=Promoter)

all_dmps=as.data.frame(out)

save(all_dmps, file=DMPs_annotated_filename) # all DMPs, irrespective of genomic location

out1=out[! as.vector(out$region) %in% c('upstream','downstream') | out$Enhancer == 'TRUE',]

write.csv(out1,file=table_filename, row.names=FALSE) # only the DMPs that overlap genes (uncluding being within 2.5kb downstream of the 3' end), gene promoters, or gene enhancers

return=c(num_dmps=nrow(dmps),
	dmps_that_overlap_promoters_or_enhancers=length(which(out$Enhancer=='TRUE' | out$Promoter=='TRUE')),
	dmps_that_overlap_promoters_or_enhancers_or_genes=nrow(out1) )

# output DMP pdf

pdf(file=pdf_filename)

if (nrow(dmps)>30){
	hist(dmps$slope, xlim = c(-0.4,0.4), xlab='Male versus female\n change in methylation', ylab='Frequency', main='DMPs', col='grey')
		} else {
	hist(dmps$slope, xlab='Male versus female\n change in methylation', ylab='Frequency', main='DMPs', col='grey')
}

rti_dmps=as.data.frame(table(factor(anno_DMPs$Relation_to_Island)))
Island = rti_dmps$Freq[which(rti_dmps$Var1=='Island')]/num_island_probes
Shore= sum(rti_dmps$Freq[grep('Shore', rti_dmps$Var1, ignore.case=TRUE)])/num_shore_probes
Shelf= sum(rti_dmps$Freq[grep('Shelf', rti_dmps$Var1, ignore.case=TRUE)])/num_shelf_probes
OpenSea= rti_dmps$Freq[which(rti_dmps$Var1=='OpenSea')]/num_opensea_probes
Promoter=length(grep('Promoter',anno_DMPs$Regulatory_Feature_Group, , ignore.case=TRUE))/num_promoter_probes
Enhancer=length(which(anno_DMPs$Enhancer=="TRUE"))/num_enhancer_probes

normalized_Relation_to_Island=data.frame(Island=Island, Shore=Shore, Shelf=Shelf, OpenSea=OpenSea, Promoter=Promoter, Enhancer=Enhancer)

bp=barplot(height=as.matrix(normalized_Relation_to_Island), xlab='Genomic Position', ylab='Normalized Frequency', main='DMPs', col='lightgrey')
labs=signif(as.vector(normalized_Relation_to_Island),3)
text(bp, 0, labs, cex=1, pos=3)

par(mar=c(5.1,5.3,4.1,2.1))

if (num_dmps<100) { n= num_dmps} else {n =100 }

mycol=as.vector(pd$sex)
mycol=gsub('Male','blue',mycol)
mycol=gsub('Female','red',mycol)

x=factor(pd$sex, levels=c('Female','Male'))

for (i in 1:n){
	zag1=paste0('probe ', probeFit$probe[i])
	#zag2=paste0('p=',signif(probeFit$p.value[i],3), '; t=',signif(probeFit$t[i],3))
	zag2=paste0('p=',signif(probeFit$p.value[i],3))
	boxplot(as.vector(beta[i,])~x,ylab='DNAm(Beta)', xlab='Sex', main=c(zag1,zag2), cex.main=1, cex.lab=2, cex.axis=2, outline=FALSE, col='lightgrey')
	#stripchart(as.vector(beta[i,])~x, vertical = TRUE, method = "jitter", add = TRUE, pch = 21, col = 'black', bg=factor(mycol,levels=c('red','blue'))) #bg='gray'
	points(	as.vector(beta[i,]) ~ jitter(as.numeric(x), amount=0.2), pch =21, col='black', bg=mycol, cex=1.2)
}

dev.off()

}

return(return)

}

#####################################
## carry out male-female DMP analysis
#####################################

# use hg19 (and not hg38), b/c hg38 was used in the DNAm datasets
# library("TxDb.Hsapiens.UCSC.hg38.knownGene")
# genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg38.knownGene, annotationPackage = 'org.Hs.eg.db', by='gene')
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene, annotationPackage = 'org.Hs.eg.db', by='gene')


# <= 6 yo (cases)
GRset=GRset_prepubescent_cases
pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/pdfs/DMPs_prepubescent_cases.pdf'
table_filename="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/data_tables/DMPs_prepubescent_cases.csv"
DMPs_annotated_filename="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_prepubescent_cases.rda"
DMP_analysis(GRset=GRset, pdf_filename=pdf_filename, table_filename=table_filename, DMPs_annotated_filename=DMPs_annotated_filename, genes=genes)

#                                          num_dmps
#                                                12
#          dmps_that_overlap_promoters_or_enhancers
#                                                 2
# dmps_that_overlap_promoters_or_enhancers_or_genes
#                                                 6

# 6-21 (cases)
GRset=GRset_6to21_cases
pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/pdfs/DMPs_6to21_cases.pdf'
table_filename="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/data_tables/DMPs_6to21_cases.csv"
DMPs_annotated_filename="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_6to21_cases.rda"
DMP_analysis(GRset=GRset, pdf_filename=pdf_filename, table_filename=table_filename, DMPs_annotated_filename=DMPs_annotated_filename, genes=genes)

#                                          num_dmps
#                                                98
#          dmps_that_overlap_promoters_or_enhancers
#                                                31
# dmps_that_overlap_promoters_or_enhancers_or_genes
#                                                45

# >21 (cases)
GRset=GRset_adult_cases
pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/pdfs/DMPs_adult_cases.pdf'
table_filename="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/data_tables/DMPs_adult_cases.csv"
DMPs_annotated_filename="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_adult_cases.rda"
DMP_analysis(GRset=GRset, pdf_filename=pdf_filename, table_filename=table_filename, DMPs_annotated_filename=DMPs_annotated_filename, genes=genes)

#                                          num_dmps
#                                               182
#          dmps_that_overlap_promoters_or_enhancers
#                                                65
# dmps_that_overlap_promoters_or_enhancers_or_genes
#                                                82

# 0-21 (cases)
GRset=GRset_allPeds_cases
pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/pdfs/DMPs_0to21_cases.pdf'
table_filename="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/data_tables/DMPs_0to21_cases.csv"
DMPs_annotated_filename="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_0to21_cases.rda"
DMP_analysis(GRset=GRset, pdf_filename=pdf_filename, table_filename=table_filename, DMPs_annotated_filename=DMPs_annotated_filename, genes=genes)

#                                          num_dmps
#                                               160
#          dmps_that_overlap_promoters_or_enhancers
#                                                54
# dmps_that_overlap_promoters_or_enhancers_or_genes
#                                                71

# <= 6 yo (controls)
GRset=GRset_prepubescent_controls
pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/pdfs/DMPs_prepubescent_controls.pdf'
table_filename="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/data_tables/DMPs_prepubescent_controls.csv"
DMPs_annotated_filename="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_prepubescent_controls.rda"
DMP_analysis(GRset=GRset, pdf_filename=pdf_filename, table_filename=table_filename, DMPs_annotated_filename=DMPs_annotated_filename, genes=genes)

#                                          num_dmps
#                                               169
#          dmps_that_overlap_promoters_or_enhancers
#                                                41
# dmps_that_overlap_promoters_or_enhancers_or_genes
#                                                74

# 6-21 (controls)
GRset=GRset_6to21_controls
pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/pdfs/DMPs_6to21_controls.pdf'
table_filename="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/data_tables/DMPs_6to21_controls.csv"
DMPs_annotated_filename="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_6to21_controls.rda"
DMP_analysis(GRset=GRset, pdf_filename=pdf_filename, table_filename=table_filename, DMPs_annotated_filename=DMPs_annotated_filename, genes=genes)

#                                          num_dmps
#                                               640
#          dmps_that_overlap_promoters_or_enhancers
#                                               180
# dmps_that_overlap_promoters_or_enhancers_or_genes
#                                               286

# >21 (controls)
GRset=GRset_adult_controls
pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/pdfs/DMPs_adult_controls.pdf'
table_filename="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/data_tables/DMPs_adult_controls.csv"
DMPs_annotated_filename="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_adult_controls.rda"
DMP_analysis(GRset=GRset, pdf_filename=pdf_filename, table_filename=table_filename, DMPs_annotated_filename=DMPs_annotated_filename, genes=genes)

#                                          num_dmps
#                                              5739
#          dmps_that_overlap_promoters_or_enhancers
#                                              1720
# dmps_that_overlap_promoters_or_enhancers_or_genes
#                                              2865

# 0-21 (controls)
GRset=GRset_allPeds_controls
pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/pdfs/DMPs_0to21_controls.pdf'
table_filename="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/data_tables/DMPs_0to21_controls.csv"
DMPs_annotated_filename="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_0to21_controls.rda"
DMP_analysis(GRset=GRset, pdf_filename=pdf_filename, table_filename=table_filename, DMPs_annotated_filename=DMPs_annotated_filename, genes=genes)

#                                          num_dmps
#                                              1510
#          dmps_that_overlap_promoters_or_enhancers
#                                               385
# dmps_that_overlap_promoters_or_enhancers_or_genes
#                                               679

#####################################
## function to carry out sex-disease interaction analysis
#####################################

interaction_analysis = function(GRset, pdf_filename, table_filename, genes, interaction_probes_table) {

pd=pData(GRset)
beta=getBeta(GRset)
neg_control_PCs=pd[,12:19]

sex=factor(as.vector(GRset$sex), levels=c('Female','Male')) # Female = 0; Male=1;
Dx=factor(as.vector(GRset$Dx_simplified))
mod = model.matrix(~ sex*Dx)
mod=cbind(mod, neg_control_PCs)

probe_fit=lmFit(object=beta,design=mod)
eb=eBayes(probe_fit)

probeFit=data.frame(probe=rownames(eb$p.value),
	intercept=probe_fit$coefficients[,1],
	coefficient=probe_fit$coefficients[,4], 
	p.value=eb$p.value[,4],
	fdr=p.adjust(as.vector(eb$p.value[,4]),method='fdr'),
	t=eb$t[,4])

rownames(probeFit)=NULL

o=order(probeFit$fdr)
probeFit=probeFit[o,]
beta=beta[o,]

num_intxn_probes=length(which(probeFit$fdr<=0.05))
percent_intxn_probes=length(which(probeFit$fdr<=0.05))/nrow(probeFit)

if (num_intxn_probes !=0 ){

intxn_probes=probeFit[which(probeFit$fdr<=0.05),]

# Map significant interaction probes to genes

anno=getAnnotation(GRset)
mm=match(rownames(beta),rownames(anno))
anno=anno[mm,]
anno_intxn_probes=anno[1:nrow(intxn_probes),]

ad_df=as.data.frame(anno_intxn_probes)

gr=GRanges(seqnames=ad_df$chr, 
	ranges=IRanges(ad_df$pos, ad_df$pos),
	strand=ad_df$strand)

match_intxn_probes_to_genes=matchGenes(gr, genes)

match_intxn_probes_to_genes=match_intxn_probes_to_genes[,-c(2,15)]

# output table of DMPs with nearest genes

Enhancer=ad_df$Enhancer
Enhancer[which(Enhancer=="")]=FALSE

Promoter=rep('FALSE', times=nrow(ad_df))
Promoter[grep('Promoter', ad_df$Regulatory_Feature_Group, ignore.case=TRUE)]='TRUE'

out=cbind(intxn_probes[,c(1,3:6)], anno_intxn_probes[,1:2], match_intxn_probes_to_genes, Enhancer=Enhancer, Promoter=Promoter)

all_intxn_probes=as.data.frame(out)

save(all_intxn_probes, file=interaction_probes_table) # all significant intxn probes, irrespective of genomic location

out1=out[! as.vector(out$region) %in% c('upstream','downstream') | out$Enhancer == 'TRUE',]

write.csv(out1,file=table_filename, row.names=FALSE) # only the significant intx probes that either overlap genes, gene promoters, or gene enhancers

return=c(num_sig_intxn_probes=nrow(intxn_probes),
	num_sig_intxn_probes_that_overlap_promoters_or_enhancers=length(which(out$Enhancer=='TRUE' | out$Promoter=='TRUE')),
	num_sig_intxn_probes_that_overlap_promoters_or_enhancers_or_genes=nrow(out1) )

# make interaction probes pdf

pdf(file=pdf_filename)

sex_short=as.vector(pd$sex)
sex_short[which(sex_short=='Male')]='M'
sex_short[which(sex_short=='Female')]='F'

dx_short=as.vector(pd$Dx_simplified)

xl=paste(sex_short, dx_short,sep=" (")
xl=paste0(xl,')')
x=factor(xl, levels=c('F (Control)', 'M (Control)', 'F (HGG)', 'M (HGG)'))

mycol=xl
mycol[which(mycol=="F (Control)")]='red'
mycol[which(mycol=="F (HGG)")]='red'
mycol[which(mycol=="M (Control)")]='blue'
mycol[which(mycol=="M (HGG)")]='blue'

par(mar=c(5.1,5.3,4.1,2.1))

if (num_intxn_probes<100) { n= num_intxn_probes} else {n =100 }

for (i in 1:n){
	zag1=paste0('probe ', probeFit$probe[i])
	# zag2=paste0('p=',signif(probeFit$p.value[i],3), '; t=',signif(probeFit$t[i],3))
	zag2=paste0('p=',signif(probeFit$p.value[i],3))
	boxplot(as.vector(beta[i,])~x, ylab='DNAm (Beta)', xlab=NULL, main=c(zag1,zag2), cex.main=1, cex.lab=2, cex.axis=1, outline=FALSE, col='lightgrey') # at=c(1,2,4,5)
	#stripchart(as.vector(beta[i,])~x, vertical = TRUE, method = "jitter", add = TRUE, pch = 21, col = 'black', bg='gray') # at=c(1,2,4,5)
	abline(v=2.5, lwd=3, lty=2)
	points(	as.vector(beta[i,]) ~ jitter(as.numeric(x), amount=0.2), pch =21, col='black', bg=mycol, cex=1.1) #at=c(1,2,4,5)
}

dev.off()

} else {return = 'No significant interaction probes.'}

return(return)

}

#####################################
## carry out sex-disease interaction analysis
#####################################

# use hg19 (and not hg38), b/c hg38 was used in the DNAm datasets
# library("TxDb.Hsapiens.UCSC.hg38.knownGene")
# genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg38.knownGene, annotationPackage = 'org.Hs.eg.db', by='gene')
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene, annotationPackage = 'org.Hs.eg.db', by='gene')

# <=6
GRset=GRset_prepubescent_cases_and_controls
table_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/data_tables/interaction_probes_prepubescent.csv'
pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/pdfs/interaction_prepubescent.pdf'
interaction_probes_table='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/intxn_probes_prepubescent.rda'
interaction_analysis(GRset=GRset, pdf_filename=pdf_filename, table_filename=table_filename, genes=genes, interaction_probes_table=interaction_probes_table) 

#                                              num_sig_intxn_probes
#                                                              2913
#          num_sig_intxn_probes_that_overlap_promoters_or_enhancers
#                                                              1201
# num_sig_intxn_probes_that_overlap_promoters_or_enhancers_or_genes
#                                                              1719

# 6-21
GRset=GRset_6to21_cases_and_controls
table_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/data_tables/interaction_probes_6to21.csv'
pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/pdfs/interaction_6to21.pdf'
interaction_probes_table='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/intxn_probes_6to21.rda'
interaction_analysis(GRset=GRset, pdf_filename=pdf_filename, table_filename=table_filename, genes=genes, interaction_probes_table=interaction_probes_table)

#                                              num_sig_intxn_probes
#                                                                 1
#          num_sig_intxn_probes_that_overlap_promoters_or_enhancers
#                                                                 0
# num_sig_intxn_probes_that_overlap_promoters_or_enhancers_or_genes
#                                                                 0

# >21
GRset=GRset_adult_cases_and_controls
table_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/data_tables/interaction_probes_adult.csv'
pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/pdfs/interaction_adult.pdf'
interaction_probes_table='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/intxn_probes_adult.rda'
interaction_analysis(GRset=GRset, pdf_filename=pdf_filename, table_filename=table_filename, genes=genes, interaction_probes_table=interaction_probes_table)

#                                              num_sig_intxn_probes
#                                                               338
#          num_sig_intxn_probes_that_overlap_promoters_or_enhancers
#                                                                92
# num_sig_intxn_probes_that_overlap_promoters_or_enhancers_or_genes
#                                                               176

# 0-21
GRset=GRset_0to21_cases_and_controls
table_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/data_tables/interaction_probes_0to21.csv'
pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/pdfs/interaction_0to21.pdf'
interaction_probes_table='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/intxn_probes_0to21.rda'
interaction_analysis(GRset=GRset, pdf_filename=pdf_filename, table_filename=table_filename, genes=genes, interaction_probes_table=interaction_probes_table)

#                                              num_sig_intxn_probes
#                                                                 9
#          num_sig_intxn_probes_that_overlap_promoters_or_enhancers
#                                                                 2
# num_sig_intxn_probes_that_overlap_promoters_or_enhancers_or_genes
#                                                                 3

#####################################
## Data visualization
#####################################

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")

## remove samples with a a mutant form of IDH

GRset=GRset[,-grep('IDH|MUT|Mutation', pData(GRset)$IDH_status)] # 36/707 samples have an IDH mutation

# 456,513 CpGs total interrogated

promoter_probes=grep('Promoter',anno$Regulatory_Feature_Group, , ignore.case=TRUE)
enhancer_probes=which(anno$Enhancer=="TRUE")

length(unique(c(promoter_probes,enhancer_probes)))

# 185,006 of the interrogated CpGs are in promoters or enhancers

##----------------------------------------
## look at all DMPs
##----------------------------------------

#-------
# [PLOT] compare sex DMPs in HGG patients and controls (stratified by age group)
#-------

dd=data.frame(HGG=c(12,98,160,182),Control=c(169,640,1510,5739))
rownames(dd)=c('Prepubescent (0-6 yo)','6-21 yo','All peds (0-21 yo)','Adults (>21 yo)')

pdf('/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/pdfs/DMPvisualization_num_DMPs_HGG_vs_controls.pdf')

bp=barplot(height=t(dd), beside=T, xlab="Age group", ylab='Frequency', col=c('grey','lightblue'), main='Sex DMPs in HGG and Control Samples')
labs=c(12,169,98,640,160,1510,182,5739)
text(bp, 0, labs, cex=1, pos=3)

legend(x='topleft',legend=c('HGG','Control'), pt.bg=c('grey','lightblue'), col='black', pch=22, cex=1, pt.cex=2, title="Disease state")

dev.off()

#-------
# [PLOT] compare sex DMPs b/w prepubescent and postpubescent subjects (stratified by diesease state)
#-------

dd=data.frame(c(12,168),c(98,640))
colnames(dd)=c('0-6 yo','6-21 yo')
rownames(dd)=c('HGG','Control')

pdf('/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/pdfs/DMPvisualization_all_DMPs_pre_vs_postpubescent.pdf')

bp=barplot(height=t(dd), beside=T, ylab='Frequency', main='Effect of puberty on quantity of sex DMPs', col=c('grey','lightgreen'))
labs=c(12,98,168,640)
text(bp, 0, labs, cex=1, pos=3)

legend(x='topleft',legend=c('Prepubescent (0-6 y.o.)','Postpubescent (6-21 y.o.)'), pt.bg=c('grey','lightgreen'), col='black', pch=22, cex=1, pt.cex=2, title="Age status")

dev.off()

# p-values for the plot above

	# HGG strata:
obs=data.frame(c(12,98),c(456513-12,456513-98))
colnames(obs)=c('DMP','not_DMP')
rownames(obs)=c('prepubescent','post-pubescent')
chisq=chisq.test(obs)
chisq$p.value # 5.278012e-16

	# Control strata:
obs=data.frame(c(168,640),c(456513-168,456513-640))
colnames(obs)=c('DMP','not_DMP')
rownames(obs)=c('prepubescent','post-pubescent')
chisq=chisq.test(obs)
chisq$p.value # 1.021051e-61

#-------
## calculate the percent of DMPs that are lost b/w Controls and HGG (all DMPs)
#-------

# 456,513 CpGs total interrogated

# 0-6 yo

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_prepubescent_cases.rda")
DMPs_prepubescent_cases=all_dmps
load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_prepubescent_controls.rda")
DMPs_prepubescent_controls=all_dmps

mm=match(DMPs_prepubescent_controls$probe,DMPs_prepubescent_cases$probe)
(length(DMPs_prepubescent_controls$probe)-length(which(!is.na(mm))))/length(DMPs_prepubescent_controls$probe) # 0.9349112

obs=data.frame(c(nrow(DMPs_prepubescent_cases),456513-nrow(DMPs_prepubescent_cases)), c(nrow(DMPs_prepubescent_controls), 456513-nrow(DMPs_prepubescent_controls)))
colnames(obs)=c('HGG','Control')
rownames(obs)=c('DMP','not_DMP')
chisq=chisq.test(obs)
chisq$p.value # 4.290864e-31

# 6-21 yo

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_6to21_cases.rda")
DMPs_6to21_cases=all_dmps
load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_6to21_controls.rda")
DMPs_6to21_controls=all_dmps

mm=match(DMPs_6to21_controls$probe, DMPs_6to21_cases$probe)
(length(DMPs_6to21_controls$probe)-length(which(!is.na(mm))))/length(DMPs_6to21_controls$probe) # 0.884375

obs=data.frame(c(nrow(DMPs_6to21_cases),456513-nrow(DMPs_6to21_cases)), c(nrow(DMPs_6to21_controls), 456513-nrow(DMPs_6to21_controls)))
colnames(obs)=c('HGG','Control')
rownames(obs)=c('DMP','not_DMP')
chisq=chisq.test(obs)
chisq$p.value # 2.594958e-88

# 0-21 yo

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_0to21_cases.rda")
DMPs_0to21_cases=all_dmps
load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_0to21_controls.rda")
DMPs_0to21_controls=all_dmps

mm=match(DMPs_0to21_controls$probe, DMPs_0to21_cases$probe)
(length(DMPs_0to21_controls$probe)-length(which(!is.na(mm))))/length(DMPs_0to21_controls$probe) # 0.9211921

obs=data.frame(c(nrow(DMPs_0to21_cases),456513-nrow(DMPs_0to21_cases)), c(nrow(DMPs_0to21_controls), 456513-nrow(DMPs_0to21_controls)))
colnames(obs)=c('HGG','Control')
rownames(obs)=c('DMP','not_DMP')
chisq=chisq.test(obs)
chisq$p.value # 2.105147e-239

# > 21 yo

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_adult_cases.rda")
DMPs_adult_cases=all_dmps
load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_adult_controls.rda")
DMPs_adult_controls=all_dmps

mm=match(DMPs_adult_controls$probe, DMPs_adult_cases$probe)
(length(DMPs_adult_controls$probe)-length(which(!is.na(mm))))/length(DMPs_adult_controls$probe) # 0.9715978

obs=data.frame(c(nrow(DMPs_adult_cases),456513-nrow(DMPs_adult_cases)), c(nrow(DMPs_adult_controls), 456513-nrow(DMPs_adult_controls)))
colnames(obs)=c('HGG','Control')
rownames(obs)=c('DMP','not_DMP')
chisq=chisq.test(obs)
chisq$p.value # 0


##----------------------------------------
## look at DMPs in ENHANCERS & PROMOTERS ONLY
##----------------------------------------

#-------
# [PLOT] compare sex DMPs in HGG patients and controls (stratified by age group)
#-------

dd=data.frame(HGG=c(2,31,54,64),Control=c(41,180,385,1720))
rownames(dd)=c('Prepubescent (0-6 yo)','6-21 yo','All peds (0-21 yo)','>21 yo')
colnames(dd)=c('HGG','Control')

pdf('/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/pdfs/DMPvisualization_num_promoter_and_enhancer_DMPs_HGG_vs_controls.pdf')

bp=barplot(height=t(dd), beside=T, xlab="Age group", ylab='Frequency', col=c('grey','lightblue'), main='Sex DMPs in HGG and Control Samples (promoters/enhancers only)')
labs=c(2,41,31,180,54,385,65,1720)
text(bp, 0, labs, cex=1, pos=3)

legend(x='topleft',legend=c('HGG','Control'), pt.bg=c('grey','lightblue'), col='black', pch=22, cex=1, pt.cex=2, title="Desease state")

dev.off()

#-------
# [PLOT] compare sex DMPs b/w prepubescent and postpubescent subjects (stratified by disease state)
#-------

dd=data.frame(c(2,41),c(31,180))
colnames(dd)=c('0-6 yo','6-21 yo')
rownames(dd)=c('HGG','Control')

pdf('/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/pdfs/DMPvisualization_promoter_enhancer_DMPs_pre_vs_postpubescent.pdf')

bp=barplot(height=t(dd), beside=T, ylab='Frequency', main='Effect of puberty on quantity of sex DMPs (promoters/enhancers only)', col=c('grey','lightgreen'))
labs=c(2,31,41,180)
text(bp, 0, labs, cex=1, pos=3)

legend(x='topleft',legend=c('Prepubescent (0-6 y.o.)','Postpubescent (6-21 y.o.)'), pt.bg=c('grey','lightgreen'), col='black', pch=22, cex=1, pt.cex=2, title="Age status")

dev.off()

# p-values for the plot above

	# HGG strata:
obs=data.frame(c(2,31),c(185006-2,185006-31))
colnames(obs)=c('DMP','not_DMP')
rownames(obs)=c('prepubescent','post-pubescent')
chisq=chisq.test(obs)
chisq$p.value # 1.091442e-06

	# Control strata:
obs=data.frame(c(41,180),c(185006-41,185006-180))
colnames(obs)=c('DMP','not_DMP')
rownames(obs)=c('prepubescent','post-pubescent')
chisq=chisq.test(obs)
chisq$p.value # 1.60711e-20

#-------
# calculate the percent of DMPs that are lost b/w Controls and HGG (promoter/enhancer DMPs only)
#-------

# 185,006 of the interrogated CpGs are in promoters or enhancers

# 0-6 yo

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_prepubescent_cases.rda")
DMPs_prepubescent_cases=all_dmps
DMPs_prepubescent_cases=DMPs_prepubescent_cases[which(DMPs_prepubescent_cases$Promoter=='TRUE' | DMPs_prepubescent_cases$Enhancer=='TRUE'),]

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_prepubescent_controls.rda")
DMPs_prepubescent_controls=all_dmps
DMPs_prepubescent_controls=DMPs_prepubescent_controls[which(DMPs_prepubescent_controls$Promoter=='TRUE' | DMPs_prepubescent_controls$Enhancer=='TRUE'),]

mm=match(DMPs_prepubescent_controls$probe,DMPs_prepubescent_cases$probe)
(length(DMPs_prepubescent_controls$probe)-length(which(!is.na(mm))))/length(DMPs_prepubescent_controls$probe) # 0.9512195

obs=data.frame(c(nrow(DMPs_prepubescent_cases),185006-nrow(DMPs_prepubescent_cases)), c(nrow(DMPs_prepubescent_controls), 185006-nrow(DMPs_prepubescent_controls)))
colnames(obs)=c('HGG','Control')
rownames(obs)=c('DMP','not_DMP')
chisq=chisq.test(obs)
chisq$p.value # 6.82063e-09

# 6-21 yo

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_6to21_cases.rda")
DMPs_6to21_cases=all_dmps
DMPs_6to21_cases=DMPs_6to21_cases[which(DMPs_6to21_cases$Promoter=='TRUE' | DMPs_6to21_cases$Enhancer=='TRUE'),]

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_6to21_controls.rda")
DMPs_6to21_controls=all_dmps
DMPs_6to21_controls=DMPs_6to21_controls[which(DMPs_6to21_controls$Promoter=='TRUE' | DMPs_6to21_controls$Enhancer=='TRUE'),]

mm=match(DMPs_6to21_controls$probe, DMPs_6to21_cases$probe)
(length(DMPs_6to21_controls$probe)-length(which(!is.na(mm))))/length(DMPs_6to21_controls$probe) # 0.8777778

obs=data.frame(c(nrow(DMPs_6to21_cases),185006-nrow(DMPs_6to21_cases)), c(nrow(DMPs_6to21_controls), 185006-nrow(DMPs_6to21_controls)))
colnames(obs)=c('HGG','Control')
rownames(obs)=c('DMP','not_DMP')
chisq=chisq.test(obs)
chisq$p.value # 2.160732e-24

# 0-21 yo

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_0to21_cases.rda")
DMPs_0to21_cases=all_dmps
DMPs_0to21_cases=DMPs_0to21_cases[which(DMPs_0to21_cases$Promoter=='TRUE' | DMPs_0to21_cases$Enhancer=='TRUE'),]

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_0to21_controls.rda")
DMPs_0to21_controls=all_dmps
DMPs_0to21_controls=DMPs_0to21_controls[which(DMPs_0to21_controls$Promoter=='TRUE' | DMPs_0to21_controls$Enhancer=='TRUE'),]

mm=match(DMPs_0to21_controls$probe, DMPs_0to21_cases$probe)
(length(DMPs_0to21_controls$probe)-length(which(!is.na(mm))))/length(DMPs_0to21_controls$probe) # 0.9038961

obs=data.frame(c(nrow(DMPs_0to21_cases),185006-nrow(DMPs_0to21_cases)), c(nrow(DMPs_0to21_controls), 185006-nrow(DMPs_0to21_controls)))
colnames(obs)=c('HGG','Control')
rownames(obs)=c('DMP','not_DMP')
chisq=chisq.test(obs)
chisq$p.value # 5.919949e-56

# > 21 yo

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_adult_cases.rda")
DMPs_adult_cases=all_dmps
DMPs_adult_cases=DMPs_adult_cases[which(DMPs_adult_cases$Promoter=='TRUE' | DMPs_adult_cases$Enhancer=='TRUE'),]

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_adult_controls.rda")
DMPs_adult_controls=all_dmps
DMPs_adult_controls=DMPs_adult_controls[which(DMPs_adult_controls$Promoter=='TRUE' | DMPs_adult_controls$Enhancer=='TRUE'),]

mm=match(DMPs_adult_controls$probe, DMPs_adult_cases$probe)
(length(DMPs_adult_controls$probe)-length(which(!is.na(mm))))/length(DMPs_adult_controls$probe) # 0.9662791

obs=data.frame(c(nrow(DMPs_adult_cases),185006-nrow(DMPs_adult_cases)), c(nrow(DMPs_adult_controls), 185006-nrow(DMPs_adult_controls)))
colnames(obs)=c('HGG','Control')
rownames(obs)=c('DMP','not_DMP')
chisq=chisq.test(obs)
chisq$p.value #0

##----------------------------------------
## visualize the number of significant sex-disease interaction CpGs by age group
##----------------------------------------

intrxn=data.frame(c(2913,1,9,338))
rownames(intrxn)=c('Prepubescent (0-6 yo)','6-21 yo','All peds (0-21 yo)','>21 yo')

pdf('/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/pdfs/num_intrxn_CpGs.pdf')

bp=barplot(height=t(intrxn), xlab="Age group", ylab='Frequency', col='grey', main='Sex-Desease Interaction CpGs')
labs=c(2913,1,9,338)
text(bp, 0, labs, cex=1, pos=3)

dev.off()

#####################################
## Analysis of genes that lose DNAm in HGG
#####################################

## function

lost_sex_DMPs_interrogation = function(cases, controls, lost_DMPs_overlapping_genes_filename) {

# load imprinted genes
ig=read.csv('/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/imprinted_human_genes.csv')

# load human TFs
TFs=read.csv("/athena/masonlab/scratch/users/nai2008/Human_TFs_DatabaseExtract_v_1.01.csv")

mm=match(controls$probe, cases$probe)
loss=controls[which(is.na(mm)),]
loss=loss[,c(1,5,7:ncol(loss))]

loss1=loss[! as.vector(loss$region) %in% c('upstream','downstream'),]

# How much of the lost DMPs overlap human TFs?
mm=match(loss1$name,TFs$HGNC_symbol)
num_TF_DMPs=length(which(!is.na(mm)))

# Regarding the queston above, how many unique TFs are there?
mmm=mm[-which(is.na(mm))]
num_TFs=length(unique(TFs$HGNC_symbol[mmm]))

loss1$TF=FALSE
loss1$TF[!is.na(mm)]=TRUE

# How much of the lost DMPs overlap imprinted genes (genes that are imprinted in humans)?
mm=match(loss1$name,ig$Gene)
num_ig_DMPs=length(which(!is.na(mm)))

# Regarding the queston above, how many unique imprinted genes are there?
mmm=mm[-which(is.na(mm))]
num_ig=length(unique(ig$Gene[mmm]))

loss1$imprinted_gene="-"
loss1$imprinted_gene[!is.na(mm)]=TRUE

loss1$expressed_allele="-"

index_ig=which(loss1$imprinted_gene==TRUE)

if (length(index_ig)>0){
	for (i in 1:length(index_ig)){
		m=match(loss1$name[index_ig[i]],ig$Gene)
		loss1$expressed_allele[index_ig[i]]=ig$ExpressedAllele[m]
	}
}

# Output the data

if(nrow(loss1 != 0)){
	write.csv(loss1, file=lost_DMPs_overlapping_genes_filename, row.names=FALSE) # lost DMPs which overlap genes
}

dat=data.frame(
	num_DMPs_lost=nrow(loss),
	proportion_DMPs_lost=nrow(loss)/nrow(controls),
	how_many_lost_DMPs_overlap_genes_or_enhancers= nrow(loss[! as.vector(loss$region) %in% c('upstream','downstream') | loss$Enhancer == 'TRUE',]),
	how_many_lost_DMPs_overlap_TFs=num_TF_DMPs,
	how_many_TFs=num_TFs,
	how_many_lost_DMPs_overlap_imprinted_genes=num_ig_DMPs,
	how_many_imprinted_genes=num_ig
	)

dat=t(dat)

return(dat)

}

# 0-6 yo -----------------
load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_prepubescent_cases.rda")
DMPs_prepubescent_cases=all_dmps

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_prepubescent_controls.rda")
DMPs_prepubescent_controls=all_dmps

cases=DMPs_prepubescent_cases
controls=DMPs_prepubescent_controls

lost_DMPs_overlapping_genes_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_0to6.csv'

lost_sex_DMPs_interrogation(cases=cases, controls=controls, lost_DMPs_overlapping_genes_filename=lost_DMPs_overlapping_genes_filename)

# num_DMPs_lost                                 158.0000000
# proportion_DMPs_lost                            0.9349112
# how_many_lost_DMPs_overlap_genes_or_enhancers  68.0000000
# how_many_lost_DMPs_overlap_TFs                 11.0000000
# how_many_TFs                                   11.0000000
# how_many_lost_DMPs_overlap_imprinted_genes      0.0000000
# how_many_imprinted_genes                        0.0000000

#6-21 yo -----------------
load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_6to21_cases.rda")
DMPs_6to21_cases=all_dmps

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_6to21_controls.rda")
DMPs_6to21_controls=all_dmps

cases=DMPs_6to21_cases
controls=DMPs_6to21_controls

lost_DMPs_overlapping_genes_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_6to21.csv'

lost_sex_DMPs_interrogation(cases=cases, controls=controls, lost_DMPs_overlapping_genes_filename=lost_DMPs_overlapping_genes_filename)

# num_DMPs_lost                                 566.000000
# proportion_DMPs_lost                            0.884375
# how_many_lost_DMPs_overlap_genes_or_enhancers 254.000000
# how_many_lost_DMPs_overlap_TFs                 39.000000
# how_many_TFs                                   30.000000
# how_many_lost_DMPs_overlap_imprinted_genes      0.000000
# how_many_imprinted_genes                        0.000000

#0-21 yo -----------------
load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_0to21_cases.rda")
DMPs_0to21_cases=all_dmps

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_0to21_controls.rda")
DMPs_0to21_controls=all_dmps

cases=DMPs_0to21_cases
controls=DMPs_0to21_controls

lost_DMPs_overlapping_genes_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_0to21.csv'

lost_sex_DMPs_interrogation(cases=cases, controls=controls, lost_DMPs_overlapping_genes_filename=lost_DMPs_overlapping_genes_filename)

# num_DMPs_lost                                 1391.0000000
# proportion_DMPs_lost                             0.9211921
# how_many_lost_DMPs_overlap_genes_or_enhancers  628.0000000
# how_many_lost_DMPs_overlap_TFs                  97.0000000
# how_many_TFs                                    61.0000000
# how_many_lost_DMPs_overlap_imprinted_genes       1.0000000
# how_many_imprinted_genes                         1.0000000

#>21 yo -----------------
load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_adult_cases.rda")
DMPs_adult_cases=all_dmps

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/annotated_DMPs_adult_controls.rda")
DMPs_adult_controls=all_dmps

cases=DMPs_adult_cases
controls=DMPs_adult_controls

lost_DMPs_overlapping_genes_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_adults.csv'

lost_sex_DMPs_interrogation(cases=cases, controls=controls, lost_DMPs_overlapping_genes_filename=lost_DMPs_overlapping_genes_filename)

# num_DMPs_lost                                 5576.0000000
# proportion_DMPs_lost                             0.9715978
# how_many_lost_DMPs_overlap_genes_or_enhancers 2790.0000000
# how_many_lost_DMPs_overlap_TFs                 396.0000000
# how_many_TFs                                   244.0000000
# how_many_lost_DMPs_overlap_imprinted_genes      23.0000000
# how_many_imprinted_genes                        12.0000000

### NAI
