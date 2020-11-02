##
library(minfi)
library(limma)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")

load('/athena/masonlab/scratch/users/nai2008/gencode_v35_hg19_genes_and_genomicState.rda')

## remove control samples

GRset=GRset[,-which(GRset$Dx=='Control')]
pd=pData(GRset)

# label IDH status

IDH_Status=rep('WT',times=nrow(pData(GRset)))
IDH_Status[grep('IDH|MUT|Mutation', pData(GRset)$IDH_status)]='Mutated'
pData(GRset)$IDH_Status=IDH_Status

table(GRset$IDH_Status)
# Mutated      WT
#      38     300

range(pData(GRset)$age[which(GRset$IDH_Status=='Mutated')]) # 8 71
range(pData(GRset)$age[which(GRset$IDH_Status=='WT')]) # 0.1 85.0

table(pData(GRset)$sex[which(GRset$IDH_Status=='Mutated')])
# Female   Male
#     13     25

table(pData(GRset)$sex[which(GRset$IDH_Status=='WT')])
# Female   Male
#    130    170

##############################################
## find DMPs by IDH status
##############################################

beta=getBeta(GRset)
neg_control_PCs=pd[,13:21]

idh=factor(as.vector(GRset$IDH_Status), levels=c('WT','Mutated')) # WT = 0; Mutated=1;
mod = model.matrix(~idh)
mod=cbind(mod, neg_control_PCs)

probe_fit=lmFit(object=beta,design=mod)
eb=eBayes(probe_fit)

probeFit=data.frame(probes=rownames(eb$p.value),
	intercept=probe_fit$coefficients[,1],
	slope=probe_fit$coefficients[,2], 
	p.value=eb$p.value[,2],
	fdr=p.adjust(as.vector(eb$p.value[,2]),method='fdr'),
	t=eb$t[,2])

rownames(probeFit)=NULL

o=order(probeFit$fdr)
probeFit=probeFit[o,]
beta=beta[o,]

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

# gr=GRanges(seqnames=ad_df$chr, 
# 	ranges=IRanges(ad_df$pos, ad_df$pos),
# 	strand=ad_df$strand)

# match_DMPs_to_genes=matchGenes(gr, genes)
# match_DMPs_to_genes=match_DMPs_to_genes[,-c(2,15)]

# output table of DMPs with nearest genes

Enhancer=ad_df$Enhancer
Enhancer[which(Enhancer=="")]=FALSE

Promoter=rep('FALSE', times=nrow(dmps))
Promoter[grep('Promoter', ad_df$Regulatory_Feature_Group, ignore.case=TRUE)]='TRUE'

out=cbind(dmps, anno_DMPs[,1:2], Enhancer=Enhancer, Promoter=Promoter)

all_dmps=as.data.frame(out)

save(all_dmps, file='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/DMPs_by_IDH_status_allDMPs.rda') # all DMPs, irrespective of genomic location

# out=cbind(dmps, anno_DMPs[,1:2], match_DMPs_to_genes, Enhancer=Enhancer, Promoter=Promoter)
# out1=out[! as.vector(out$region) %in% c('upstream','downstream') | out$Enhancer == 'TRUE',]
# write.csv(out1,file='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/DMPs_by_IDH_status_which_overlap_genes_promoters_or_enhancers.csv', row.names=FALSE) 
# ^ only the DMPs that overlap genes (uncluding being within 2.5kb downstream of the 3' end), gene promoters, or gene enhancers

return=c(num_dmps=nrow(dmps),
	percent_dmps=(length(which(probeFit$fdr<=0.05))/nrow(probeFit))*100,
	dmps_that_overlap_promoters_or_enhancers=length(which(all_dmps$Enhancer=='TRUE' | all_dmps$Promoter=='TRUE')),
	proportion_DMPs_HYPERmethylated_in_IDH_mut=length(which(dmps$slope>0))/nrow(dmps),
	proportion_DMPs_HYPOmethylated_in_IDH_mut=length(which(dmps$slope<0))/nrow(dmps))

#                                   num_dmps
#                               1.913000e+05
#                               percent_dmps
#                               4.186655e+01
#   dmps_that_overlap_promoters_or_enhancers
#                               7.395800e+04
# proportion_DMPs_HYPERmethylated_in_IDH_mut
#                               8.607214e-01
#  proportion_DMPs_HYPOmethylated_in_IDH_mut
#                               1.392786e-01

# output DMP pdf

pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/DMPs_IDH.pdf'
pdf(file=pdf_filename)

#xlab='DNAm difference between IDH WT and Mutant HGG'
#xlab=expression(atop('DNAm difference between IDH WT and Mutant HGG',"[ DNAm in"~bold('IDHmut')~"HGG ] - [ DNAm in"~bold('IDHwt')~"HGG ]"))
#xlab='DNAm difference between IDH WT and Mutant HGG\n [ DNAm in IDHmut HGG ] - [ DNAm in IDHwt HGG ]'

hist(dmps$slope, xlim = c(-0.6, 0.6), xlab=expression("( DNAm in"~bold('IDHmutant')~"HGG ) - ( DNAm in"~bold('IDHwt')~"HGG )"), 
	ylab='Frequency', main='DMPs')
abline(v=0, lwd=3, lty=2, col='red')

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

if (nrow(dmps)<100) { n= nrow(dmps) } else {n =100 }

for (i in 1:n){
	zag1=paste0('probe ', probeFit$probes[i])
	#zag2=paste0('p=',signif(probeFit$p.value[i],3), '; t=',signif(probeFit$t[i],3))
	zag2=paste0('FDR=',signif(probeFit$fdr[i],3))
	boxplot(as.vector(beta[i,])~factor(as.vector(GRset$IDH_Status), levels=c('WT','Mutated')),ylab='DNAm (Beta)', xlab='IDH Status', main=c(zag1,zag2), cex.main=1, cex.lab=2, cex.axis=2, outline=FALSE)
	stripchart(as.vector(beta[i,])~factor(as.vector(GRset$IDH_Status), levels=c('WT','Mutated')), vertical = TRUE, method = "jitter", add = TRUE, pch = 21, col = 'black', bg='lightgrey')
}

dev.off()

##############################################
## find DMPs by IDH status [adjusting for age and sex]
##############################################

beta=getBeta(GRset)
neg_control_PCs=pd[,13:21]

idh=factor(as.vector(GRset$IDH_Status), levels=c('WT','Mutated')) # WT = 0; Mutated=1;
mod = model.matrix(~idh)
mod=cbind(mod, age=GRset$age, sex=GRset$sex, neg_control_PCs)

probe_fit=lmFit(object=beta,design=mod)
eb=eBayes(probe_fit)

probeFit=data.frame(probes=rownames(eb$p.value),
	intercept=probe_fit$coefficients[,1],
	slope=probe_fit$coefficients[,2], 
	p.value=eb$p.value[,2],
	fdr=p.adjust(as.vector(eb$p.value[,2]),method='fdr'),
	t=eb$t[,2])

rownames(probeFit)=NULL

o=order(probeFit$fdr)
probeFit=probeFit[o,]
beta=beta[o,]

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

# gr=GRanges(seqnames=ad_df$chr, 
# 	ranges=IRanges(ad_df$pos, ad_df$pos),
# 	strand=ad_df$strand)

# match_DMPs_to_genes=matchGenes(gr, genes)
# match_DMPs_to_genes=match_DMPs_to_genes[,-c(2,15)]

# output table of DMPs with nearest genes

Enhancer=ad_df$Enhancer
Enhancer[which(Enhancer=="")]=FALSE

Promoter=rep('FALSE', times=nrow(dmps))
Promoter[grep('Promoter', ad_df$Regulatory_Feature_Group, ignore.case=TRUE)]='TRUE'

out=cbind(dmps, anno_DMPs[,1:2], Enhancer=Enhancer, Promoter=Promoter)

all_dmps=as.data.frame(out)

save(all_dmps, file='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/DMPs_by_IDH_status_adjustedForAgeAndSex_allDMPs.rda') # all DMPs, irrespective of genomic location

# out=cbind(dmps, anno_DMPs[,1:2], match_DMPs_to_genes, Enhancer=Enhancer, Promoter=Promoter)
# out1=out[! as.vector(out$region) %in% c('upstream','downstream') | out$Enhancer == 'TRUE',]
# write.csv(out1,file='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/DMPs_by_IDH_status_which_overlap_genes_promoters_or_enhancers.csv', row.names=FALSE) 
# ^ only the DMPs that overlap genes (uncluding being within 2.5kb downstream of the 3' end), gene promoters, or gene enhancers

return=c(num_dmps=nrow(dmps),
	percent_dmps=(length(which(probeFit$fdr<=0.05))/nrow(probeFit))*100,
	dmps_that_overlap_promoters_or_enhancers=length(which(all_dmps$Enhancer=='TRUE' | all_dmps$Promoter=='TRUE')),
	proportion_DMPs_HYPERmethylated_in_IDH_mut=length(which(dmps$slope>0))/nrow(dmps),
	proportion_DMPs_HYPOmethylated_in_IDH_mut=length(which(dmps$slope<0))/nrow(dmps))

#                                   num_dmps
#                               1.922780e+05
#                               percent_dmps
#                               4.208059e+01
#   dmps_that_overlap_promoters_or_enhancers
#                               7.383700e+04
# proportion_DMPs_HYPERmethylated_in_IDH_mut
#                               8.628652e-01
#  proportion_DMPs_HYPOmethylated_in_IDH_mut
#                               1.371348e-01

# output DMP pdf

pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/DMPs_IDH_adjustedForAgeAndSex.pdf'
pdf(file=pdf_filename)

hist(dmps$slope, xlim = c(-0.6, 0.6), xlab=expression("( DNAm in"~bold('IDHmutant')~"HGG ) - ( DNAm in"~bold('IDHwt')~"HGG )"),
 ylab='Frequency', main='DMPs')
abline(v=0, lwd=3, lty=2, col='red')

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

if (nrow(dmps)<100) { n= nrow(dmps) } else {n =100 }

for (i in 1:n){
	zag1=paste0('probe ', probeFit$probes[i])
	#zag2=paste0('p=',signif(probeFit$p.value[i],3), '; t=',signif(probeFit$t[i],3))
	zag2=paste0('FDR=',signif(probeFit$fdr[i],3))
	boxplot(as.vector(beta[i,])~factor(as.vector(GRset$IDH_Status), levels=c('WT','Mutated')),ylab='DNAm (Beta)', xlab='IDH Status', main=c(zag1,zag2), cex.main=1, cex.lab=2, cex.axis=2, outline=FALSE)
	stripchart(as.vector(beta[i,])~factor(as.vector(GRset$IDH_Status), levels=c('WT','Mutated')), vertical = TRUE, method = "jitter", add = TRUE, pch = 21, col = 'black', bg='lightgrey')
}

dev.off()

#NAI
