##
library(minfi)
library(limma)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")

## remove samples with a a mutant form of IDH

GRset=GRset[,-grep('IDH|MUT|Mutation', pData(GRset)$IDH_status)] # 38/638 samples have an IDH mutation

range(pData(GRset)$NeuN_neg_proportion, na.rm=TRUE)
median(range(pData(GRset)$NeuN_neg_proportion, na.rm=TRUE))

# > range(pData(GRset)$NeuN_neg_proportion, na.rm=TRUE)
# [1] 0.3953115 0.9597120
# > median(range(pData(GRset)$NeuN_neg_proportion, na.rm=TRUE))
# [1] 0.6775118

#####################################
## stratify the dataset
#####################################

## CASES

GRset_cases=GRset[,which(pData(GRset)$Dx_simplified=='HGG')]

# prepubescent (cases) [girls <9.5; boys <10.5]
indexes=c(which(pData(GRset_cases)$age < 9.5 & pData(GRset_cases)$sex=='Female'), 
	which(pData(GRset_cases)$age < 10.5 & pData(GRset_cases)$sex=='Male'))
GRset_prepubescent_cases=GRset_cases[,indexes]

nrow(pData(GRset_prepubescent_cases)) # 49 samples
table(pData(GRset_prepubescent_cases)$sex) # 16 F; 33 M
range(pData(GRset_prepubescent_cases)$age)

# postpubescent pediatric (cases) [girls >= 9.5 & <=21; boys >= 10.5 & <=21]
indexes=c(which(pData(GRset_cases)$age >= 9.5 & pData(GRset_cases)$age <= 21 & pData(GRset_cases)$sex=='Female'), 
	which(pData(GRset_cases)$age >= 10.5 & pData(GRset_cases)$age <= 21 & pData(GRset_cases)$sex=='Male'))
GRset_postpubescent_cases=GRset_cases[,indexes]

nrow(pData(GRset_postpubescent_cases)) # 87 samples
table(pData(GRset_postpubescent_cases)$sex) # 40 F; 47 M
range(pData(GRset_postpubescent_cases)$age)

# >21 (cases)
GRset_adult_cases=GRset_cases[,which(pData(GRset_cases)$age >21)]

nrow(pData(GRset_adult_cases)) # 164 samples
table(pData(GRset_adult_cases)$sex) # 74 F; 90 M
range(pData(GRset_adult_cases)$age)

# 0-21 (cases)
GRset_allPeds_cases=GRset_cases[,which(pData(GRset_cases)$age <=21)]

nrow(pData(GRset_allPeds_cases)) # 136 samples
table(pData(GRset_allPeds_cases)$sex) # 56 F; 80 M
range(pData(GRset_allPeds_cases)$age)

## CONTROLS

GRset_controls=GRset[,which(pData(GRset)$Dx_simplified=='Control')]

# prepubescent (controls) [girls <9.5; boys <10.5]
indexes=c(which(pData(GRset_controls)$age < 9.5 & pData(GRset_controls)$sex=='Female'), 
	which(pData(GRset_controls)$age < 10.5 & pData(GRset_controls)$sex=='Male'))
GRset_prepubescent_controls=GRset_controls[,indexes]

nrow(pData(GRset_prepubescent_controls)) # 37 samples
table(pData(GRset_prepubescent_controls)$sex) # 10 F; 27 M
range(pData(GRset_prepubescent_controls)$age)

# postpubescent pediatric (controls) [girls >= 9.5 & <=21; boys >= 10.5 & <=21]
indexes=c(which(pData(GRset_controls)$age >= 9.5 & pData(GRset_controls)$age <= 21 & pData(GRset_controls)$sex=='Female'), 
	which(pData(GRset_controls)$age >= 10.5 & pData(GRset_controls)$age <= 21 & pData(GRset_controls)$sex=='Male'))
GRset_postpubescent_controls=GRset_controls[,indexes]

nrow(pData(GRset_postpubescent_controls)) # 54 samples
table(pData(GRset_postpubescent_controls)$sex) # 17 F; 37 M
range(pData(GRset_postpubescent_controls)$age)

# >21 (controls)
GRset_adult_controls=GRset_controls[,which(pData(GRset_controls)$age >21)]

nrow(pData(GRset_adult_controls)) # 209 samples
table(pData(GRset_adult_controls)$sex) # 65 F; 144 M
range(pData(GRset_adult_controls)$age)

# 0-21 (controls)
GRset_allPeds_controls=GRset_controls[,which(pData(GRset_controls)$age <=21)]

nrow(pData(GRset_allPeds_controls)) # 91 samples
table(pData(GRset_allPeds_controls)$sex) # 27 F; 64 M
range(pData(GRset_allPeds_controls)$age)

## CASES & CONTROLS

# prepubescent (cases & controls)
indexes=c(which(pData(GRset)$age < 9.5 & pData(GRset)$sex=='Female'), 
	which(pData(GRset)$age < 10.5 & pData(GRset)$sex=='Male'))
GRset_prepubescent_cases_and_controls=GRset[,indexes]

nrow(pData(GRset_prepubescent_cases_and_controls)) # 86 samples
table(pData(GRset_prepubescent_cases_and_controls)$sex) # 26 F; 60 M

# postpubescent (cases & controls)
indexes=c(which(pData(GRset)$age >= 9.5 & pData(GRset)$age <= 21 & pData(GRset)$sex=='Female'), 
	which(pData(GRset)$age >= 10.5 & pData(GRset)$age <= 21 & pData(GRset)$sex=='Male'))
GRset_postpubescent_cases_and_controls=GRset[,indexes]

nrow(pData(GRset_postpubescent_cases_and_controls)) # 141 samples
table(pData(GRset_postpubescent_cases_and_controls)$sex) # 57 F; 84 M

# >21 (cases & controls)
GRset_adult_cases_and_controls=GRset[,which(pData(GRset)$age > 21)]

nrow(pData(GRset_adult_cases_and_controls)) # 373 samples
table(pData(GRset_adult_cases_and_controls)$sex) #  139 F; 234 M

# 0-21 (cases & controls)
GRset_allPeds_cases_and_controls=GRset[,which(pData(GRset)$age <= 21)]

nrow(pData(GRset_allPeds_cases_and_controls)) # 227 samples
table(pData(GRset_allPeds_cases_and_controls)$sex) # 83 F; 144 M

#####################################
## function to carry out male-female DMP analysis
#####################################

DMP_analysis = function(GRset, pdf_filename, table_filename, DMPs_annotated_filename, genes, center_hist=TRUE) {

# find DMPs

pd=pData(GRset)
beta=getBeta(GRset)
neg_control_PCs=pd[,13:21]

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

# Determine nearest genes to the DMPs

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

nn=nearest(gr, genes, select='arbitrary', ignore.strand=TRUE)
hits=genes[nn,]
dist=distance(gr,hits,ignore.strand=TRUE)

Enhancer=ad_df$Enhancer
Enhancer[which(Enhancer=="")]=FALSE

Promoter=rep('FALSE', times=nrow(dmps))
Promoter[grep('Promoter', ad_df$Regulatory_Feature_Group, ignore.case=TRUE)]='TRUE'

out=as.data.frame(cbind(dmps, anno_DMPs[,1:2], distance_to_nearest_gene=dist, Enhancer=Enhancer, Promoter=Promoter))
out$distance_to_nearest_gene=as.numeric(out$distance_to_nearest_gene)

all_dmps=as.data.frame(out)

save(all_dmps, file=DMPs_annotated_filename) # all DMPs

out1=out[which(out$distance_to_nearest_gene <= 5000 | out$Enhancer == 'TRUE' | out$Promoter == 'TRUE'),]

# Find genes that overlap or are <= 5kb from genes
oo=as.data.frame(findOverlaps(gr, genes, maxgap=5000, select='all', ignore.strand=TRUE))
hits=genes[oo$subjectHits,]
dist=distance(gr[oo$queryHits],hits,ignore.strand=TRUE)

pt1=dmps[oo$queryHits,]
pt2=anno_DMPs[,1:2]
pt2=pt2[oo$queryHits,]

out2=cbind(pt1, pt2, gene=as.vector(hits$Gene), Ensembl=as.vector(hits$Geneid), distance=dist)

write.csv(out2,file=table_filename, row.names=FALSE) # only the DMPs that overlap or are within 5kb of genes

return=c(num_dmps=nrow(dmps),
	dmps_that_overlap_promoters_or_enhancers=length(which(out$Enhancer=='TRUE' | out$Promoter=='TRUE')),
	dmps_that_overlap_promoters_or_enhancers_or_genes=nrow(out1),
	range_of_abs_values_of_DNAm_differences_bw_sexes_in_sex_DMPs_MIN=signif(min(abs(out$slope)),2),
	range_of_abs_values_of_DNAm_differences_bw_sexes_in_sex_DMPs_MAX=signif(max(abs(out$slope)),2) )

# output DMP pdf

pdf(file=pdf_filename)

if (center_hist=='TRUE'){
	hist(dmps$slope, xlim = c(-max(abs(dmps$slope)),max(abs(dmps$slope))), xlab='Male versus female\n change in methylation', ylab='Frequency', main='DMPs', col='grey')
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
	zag2=paste0('FDR=',signif(probeFit$fdr[i],3))
	boxplot(as.vector(beta[i,])~x,ylab='DNAm (Beta)', xlab='Sex', main=c(zag1,zag2), cex.main=1, cex.lab=2, cex.axis=2, outline=FALSE, col='lightgrey')
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

#load genes
load('/athena/masonlab/scratch/users/nai2008/gencode_v35_hg19_genes_and_genomicState.rda')

# Prepubescent (cases)
DMP_analysis(GRset=GRset_prepubescent_cases, 
	pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/DMPs_prepubescent_cases.pdf', 
	table_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/DMPs_prepubescent_cases.csv", 
	DMPs_annotated_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_prepubescent_cases.rda", 
	genes=genes)

#                                                         num_dmps
#                                                           53.000
#                         dmps_that_overlap_promoters_or_enhancers
#                                                           15.000
#                dmps_that_overlap_promoters_or_enhancers_or_genes
#                                                           50.000
# range_of_abs_values_of_DNAm_differences_bw_sexes_in_sex_DMPs_MIN
#                                                            0.039
# range_of_abs_values_of_DNAm_differences_bw_sexes_in_sex_DMPs_MAX
#                                                            0.400

# Postpubescent (cases)
DMP_analysis(GRset=GRset_postpubescent_cases, 
	pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/DMPs_postpubescent_cases.pdf', 
	table_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/DMPs_postpubescent_cases.csv", 
	DMPs_annotated_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_postpubescent_cases.rda", 
	genes=genes)

#                                                         num_dmps
#                                                           50.000
#                         dmps_that_overlap_promoters_or_enhancers
#                                                           11.000
#                dmps_that_overlap_promoters_or_enhancers_or_genes
#                                                           48.000
# range_of_abs_values_of_DNAm_differences_bw_sexes_in_sex_DMPs_MIN
#                                                            0.037
# range_of_abs_values_of_DNAm_differences_bw_sexes_in_sex_DMPs_MAX
#                                                            0.380

# >21 (cases)
DMP_analysis(GRset=GRset_adult_cases, 
	pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/DMPs_adult_cases.pdf', 
	table_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/DMPs_adult_cases.csv", 
	DMPs_annotated_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_adult_cases.rda", 
	genes=genes)

#                                                         num_dmps
#                                                         199.0000
#                         dmps_that_overlap_promoters_or_enhancers
#                                                          69.0000
#                dmps_that_overlap_promoters_or_enhancers_or_genes
#                                                         190.0000
# range_of_abs_values_of_DNAm_differences_bw_sexes_in_sex_DMPs_MIN
#                                                           0.0031
# range_of_abs_values_of_DNAm_differences_bw_sexes_in_sex_DMPs_MAX
#                                                           0.4000

# 0-21 (cases)
DMP_analysis(GRset=GRset_allPeds_cases, 
	pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/DMPs_allPeds_cases.pdf', 
	table_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/DMPs_allPeds_cases.csv", 
	DMPs_annotated_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_allPeds_cases.rda", 
	genes=genes)

#                                                         num_dmps
#                                                           141.00
#                         dmps_that_overlap_promoters_or_enhancers
#                                                            42.00
#                dmps_that_overlap_promoters_or_enhancers_or_genes
#                                                           136.00
# range_of_abs_values_of_DNAm_differences_bw_sexes_in_sex_DMPs_MIN
#                                                             0.01
# range_of_abs_values_of_DNAm_differences_bw_sexes_in_sex_DMPs_MAX
#                                                             0.39

# Prepubescent (controls)
DMP_analysis(GRset=GRset_prepubescent_controls, 
	pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/DMPs_prepubescent_controls.pdf', 
	table_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/DMPs_prepubescent_controls.csv", 
	DMPs_annotated_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_prepubescent_controls.rda", 
	genes=genes)

#                                                         num_dmps
#                                                           89.000
#                         dmps_that_overlap_promoters_or_enhancers
#                                                           21.000
#                dmps_that_overlap_promoters_or_enhancers_or_genes
#                                                           86.000
# range_of_abs_values_of_DNAm_differences_bw_sexes_in_sex_DMPs_MIN
#                                                            0.012
# range_of_abs_values_of_DNAm_differences_bw_sexes_in_sex_DMPs_MAX
#                                                            0.470

# Postpubescent (controls)
DMP_analysis(GRset=GRset_postpubescent_controls, 
	pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/DMPs_postpubescent_controls.pdf', 
	table_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/DMPs_postpubescent_controls.csv", 
	DMPs_annotated_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_postpubescent_controls.rda", 
	genes=genes)

#                                                         num_dmps
#                                                         276.0000
#                         dmps_that_overlap_promoters_or_enhancers
#                                                          76.0000
#                dmps_that_overlap_promoters_or_enhancers_or_genes
#                                                         264.0000
# range_of_abs_values_of_DNAm_differences_bw_sexes_in_sex_DMPs_MIN
#                                                           0.0036
# range_of_abs_values_of_DNAm_differences_bw_sexes_in_sex_DMPs_MAX
#                                                           0.4800

# >21 (controls)
DMP_analysis(GRset=GRset_adult_controls, 
	pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/DMPs_adult_controls.pdf', 
	table_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/DMPs_adult_controls.csv", 
	DMPs_annotated_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_adult_controls.rda", 
	genes=genes)

#                                                         num_dmps
#                                                        4.702e+03
#                         dmps_that_overlap_promoters_or_enhancers
#                                                        1.494e+03
#                dmps_that_overlap_promoters_or_enhancers_or_genes
#                                                        4.459e+03
# range_of_abs_values_of_DNAm_differences_bw_sexes_in_sex_DMPs_MIN
#                                                        8.900e-04
# range_of_abs_values_of_DNAm_differences_bw_sexes_in_sex_DMPs_MAX
#                                                        4.500e-01

# 0-21 (controls)
DMP_analysis(GRset=GRset_allPeds_controls, 
	pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/DMPs_allPeds_controls.pdf', 
	table_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/DMPs_allPeds_controls.csv", 
	DMPs_annotated_filename="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_allPeds_controls.rda", 
	genes=genes)

#                                                         num_dmps
#                                                         695.0000
#                         dmps_that_overlap_promoters_or_enhancers
#                                                         186.0000
#                dmps_that_overlap_promoters_or_enhancers_or_genes
#                                                         657.0000
# range_of_abs_values_of_DNAm_differences_bw_sexes_in_sex_DMPs_MIN
#                                                           0.0029
# range_of_abs_values_of_DNAm_differences_bw_sexes_in_sex_DMPs_MAX
#                                                           0.4800

###

# Sex-DMP slope ranges

# Prepubescent (cases):	 0.039 - 0.40
# Prepubescent (controls): 0.012 - 0.47

# Postpubescent (cases): 0.037 - 0.38
# Postpubescent (controls): 0.0036 - 0.48

# 0-21 (cases): 0.01 - 0.39
# 0-21 (controls): 0.0029 - 0.48

# >21 (cases): 0.0031 - 0.40
# >21 (controls): 0.00089 - 0.45


#####################################
## function to carry out sex-disease interaction analysis
#####################################

interaction_analysis = function(GRset, pdf_filename, table_filename, genes, interaction_probes_table) {

pd=pData(GRset)
beta=getBeta(GRset)
neg_control_PCs=pd[,13:21]

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

nn=nearest(gr, genes, select='arbitrary', ignore.strand=TRUE)
hits=genes[nn,]
dist=distance(gr,hits,ignore.strand=TRUE)

Enhancer=ad_df$Enhancer
Enhancer[which(Enhancer=="")]=FALSE

Promoter=rep('FALSE', times=nrow(ad_df))
Promoter[grep('Promoter', ad_df$Regulatory_Feature_Group, ignore.case=TRUE)]='TRUE'

out=cbind(intxn_probes[,c(1,3:6)], anno_intxn_probes[,1:2], distance_to_nearest_gene=dist, Enhancer=Enhancer, Promoter=Promoter)
out$distance_to_nearest_gene=as.numeric(out$distance_to_nearest_gene)

all_intxn_probes=as.data.frame(out)

save(all_intxn_probes, file=interaction_probes_table) # all significant intxn probes, irrespective of genomic location

out1=out[which(out$distance_to_nearest_gene <= 5000 | out$Enhancer == 'TRUE' | out$Promoter == 'TRUE'),]

# Find genes that overlap or are <= 5kb from genes
oo=as.data.frame(findOverlaps(gr, genes, maxgap=5000, select='all', ignore.strand=TRUE))
hits=genes[oo$subjectHits,]
dist=distance(gr[oo$queryHits],hits,ignore.strand=TRUE)

pt1=intxn_probes[oo$queryHits,]
pt1=pt1[,c(1,3:6)]
pt2=anno_intxn_probes[,1:2]
pt2=pt2[oo$queryHits,]

out2=cbind(pt1, pt2, gene=as.vector(hits$Gene), Ensembl=as.vector(hits$Geneid), distance=dist)

write.csv(out2,file=table_filename, row.names=FALSE) # only the intxn probes that overlap or are within 5kb of genes

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
	zag2=paste0('FDR=',signif(probeFit$fdr[i],3))
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

#load genes
load('/athena/masonlab/scratch/users/nai2008/gencode_v35_hg19_genes_and_genomicState.rda')

# Prepubescent
interaction_analysis(GRset=GRset_prepubescent_cases_and_controls, 
	pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/interaction_prepubescent.pdf', 
	table_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/interaction_probes_prepubescent.csv', 
	genes=genes, 
	interaction_probes_table='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/intxn_probes_prepubescent.rda') 

# [1] "No significant interaction probes."

# Postpubescent
interaction_analysis(GRset=GRset_postpubescent_cases_and_controls, 
	pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/interaction_postpubescent.pdf', 
	table_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/interaction_probes_postpubescent.csv', 
	genes=genes, 
	interaction_probes_table='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/intxn_probes_postpubescent.rda')

# [1] "No significant interaction probes."

# >21
interaction_analysis(GRset=GRset_adult_cases_and_controls, 
	pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/interaction_adult.pdf', 
	table_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/interaction_probes_adult.csv', 
	genes=genes, 
	interaction_probes_table='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/intxn_probes_adult.rda')

#                                              num_sig_intxn_probes
#                                                                38
#          num_sig_intxn_probes_that_overlap_promoters_or_enhancers
#                                                                 8
# num_sig_intxn_probes_that_overlap_promoters_or_enhancers_or_genes
#                                                                34

# 0-21
interaction_analysis(GRset=GRset_allPeds_cases_and_controls, 
	pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/interaction_allPeds.pdf', 
	table_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/interaction_probes_allPeds.csv', 
	genes=genes, 
	interaction_probes_table='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/intxn_probes_allPeds.rda')

#                                              num_sig_intxn_probes
#                                                                 1
#          num_sig_intxn_probes_that_overlap_promoters_or_enhancers
#                                                                 0
# num_sig_intxn_probes_that_overlap_promoters_or_enhancers_or_genes
#                                                                 1

#####################################
## Data visualization
#####################################

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")

## remove samples with a a mutant form of IDH

GRset=GRset[,-grep('IDH|MUT|Mutation', pData(GRset)$IDH_status)]

# 456928 total probes interrogated

anno=getAnnotation(GRset)
promoter_probes=grep('Promoter',anno$Regulatory_Feature_Group, , ignore.case=TRUE)
enhancer_probes=which(anno$Enhancer=="TRUE")

length(unique(c(promoter_probes,enhancer_probes)))

# 185150 of the interrogated CpGs are in promoters or enhancers

##----------------------------------------
## look at all DMPs
##----------------------------------------

#-------
# [PLOT] compare sex DMPs in HGG patients and controls (stratified by age group)
#-------

dd=data.frame(HGG=c(53,50,141,199),Control=c(89,276,695,4702))
rownames(dd)=c('Prepubescent','Postpubescent','All peds (0-21 yo)','Adults (>21 yo)')

pdf('/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/DMPvisualization_num_DMPs_HGG_vs_controls.pdf')

bp=barplot(height=t(dd), beside=T, xlab="Age group", ylab='Number of sex-DMPs', col=c('grey','lightblue'), main='Sex-DMPs in HGG and Control Samples')
labs=c(53,89,50,276,141,695,199,4702)
text(bp, 0, labs, cex=1, pos=3, offset=1.7)

legend(x='topleft',legend=c('HGG','Control'), pt.bg=c('grey','lightblue'), col='black', pch=22, cex=1, pt.cex=2, title="Disease state")

dev.off()

#-------
# [PLOT] compare sex DMPs b/w prepubescent and postpubescent subjects (stratified by disease state)
#-------

dd=data.frame(c(53,89),c(50,276))
colnames(dd)=c('Prepubescent','Postpubescent')
rownames(dd)=c('HGG','Control')

pdf('/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/DMPvisualization_all_DMPs_pre_vs_postpubescent.pdf')

bp=barplot(height=t(dd), beside=T, ylab='Number of sex-DMPs', main='Effect of puberty on quantity of sex-DMPs', col=c('grey','lightgreen'))
labs=c(53,50,89,276)
text(bp, 0, labs, cex=1, pos=3)

legend(x='topleft',legend=c('Prepubescent','Postpubescent'), pt.bg=c('grey','lightgreen'), col='black', pch=22, cex=1, pt.cex=2, title="Age group")

dev.off()

# p-values for the plot above

	# HGG strata:
obs=data.frame(c(53,50),c(456928-53, 456928-50))
colnames(obs)=c('DMP','not_DMP')
rownames(obs)=c('prepubescent','post-pubescent')
chisq=chisq.test(obs)
chisq$p.value # 0.8437673

	# Control strata:
obs=data.frame(c(89,276),c(456928-89, 456928-276))
colnames(obs)=c('DMP','not_DMP')
rownames(obs)=c('prepubescent','post-pubescent')
chisq=chisq.test(obs)
chisq$p.value # 2.083573e-22

#-------
## calculate the percent of DMPs that are lost b/w Controls and HGG (all DMPs)
#-------

# 456928 CpGs total interrogated

# Prepubescent

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_prepubescent_cases.rda")
DMPs_prepubescent_cases=all_dmps
load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_prepubescent_controls.rda")
DMPs_prepubescent_controls=all_dmps

mm=match(DMPs_prepubescent_controls$probe,DMPs_prepubescent_cases$probe)
(length(DMPs_prepubescent_controls$probe)-length(which(!is.na(mm))))/length(DMPs_prepubescent_controls$probe) # 0.5955056

obs=data.frame(c(nrow(DMPs_prepubescent_cases),456928-nrow(DMPs_prepubescent_cases)), c(nrow(DMPs_prepubescent_controls), 456928-nrow(DMPs_prepubescent_controls)))
colnames(obs)=c('HGG','Control')
rownames(obs)=c('DMP','not_DMP')
chisq=chisq.test(obs)
chisq$p.value # 0.003310164

# Postpubescent

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_postpubescent_cases.rda")
DMPs_postpubescent_cases=all_dmps
load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_postpubescent_controls.rda")
DMPs_postpubescent_controls=all_dmps

mm=match(DMPs_postpubescent_controls$probe, DMPs_postpubescent_cases$probe)
(length(DMPs_postpubescent_controls$probe)-length(which(!is.na(mm))))/length(DMPs_postpubescent_controls$probe) # 0.8369565

obs=data.frame(c(nrow(DMPs_postpubescent_cases),456928-nrow(DMPs_postpubescent_cases)), c(nrow(DMPs_postpubescent_controls), 456928-nrow(DMPs_postpubescent_controls)))
colnames(obs)=c('HGG','Control')
rownames(obs)=c('DMP','not_DMP')
chisq=chisq.test(obs)
chisq$p.value # 1.175985e-35

# 0-21 yo

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_allPeds_cases.rda")
DMPs_allPeds_cases=all_dmps
load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_allPeds_controls.rda")
DMPs_allPeds_controls=all_dmps

mm=match(DMPs_allPeds_controls$probe, DMPs_allPeds_cases$probe)
(length(DMPs_allPeds_controls$probe)-length(which(!is.na(mm))))/length(DMPs_allPeds_controls$probe) # 0.8546763

obs=data.frame(c(nrow(DMPs_allPeds_cases),456928-nrow(DMPs_allPeds_cases)), c(nrow(DMPs_allPeds_controls), 456928-nrow(DMPs_allPeds_controls)))
colnames(obs)=c('HGG','Control')
rownames(obs)=c('DMP','not_DMP')
chisq=chisq.test(obs)
chisq$p.value # 1.299262e-81

# > 21 yo

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_adult_cases.rda")
DMPs_adult_cases=all_dmps
load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_adult_controls.rda")
DMPs_adult_controls=all_dmps

mm=match(DMPs_adult_controls$probe, DMPs_adult_cases$probe)
(length(DMPs_adult_controls$probe)-length(which(!is.na(mm))))/length(DMPs_adult_controls$probe) # 0.9615057

obs=data.frame(c(nrow(DMPs_adult_cases),456928-nrow(DMPs_adult_cases)), c(nrow(DMPs_adult_controls), 456928-nrow(DMPs_adult_controls)))
colnames(obs)=c('HGG','Control')
rownames(obs)=c('DMP','not_DMP')
chisq=chisq.test(obs)
chisq$p.value # 0


##----------------------------------------
## look at DMPs in ENHANCERS & PROMOTERS ONLY
##----------------------------------------

# 185150 of the interrogated CpGs are in promoters or enhancers

#-------
# [PLOT] compare sex DMPs in HGG patients and controls (stratified by age group)
#-------

dd=data.frame(HGG=c(15,11,42,69),Control=c(21,76,186,1494))
rownames(dd)=c('Prepubescent','Postpubescent','All peds (0-21 yo)','>21 yo')
colnames(dd)=c('HGG','Control')

pdf('/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/DMPvisualization_num_promoter_and_enhancer_DMPs_HGG_vs_controls.pdf')

bp=barplot(height=t(dd), beside=T, xlab="Age group", ylab='Number of sex-DMPs', col=c('grey','lightblue'), main='Sex DMPs in HGG and Control Samples (promoters/enhancers only)')
labs=c(15,21,11,76,42,186,69,1494)
text(bp, 0, labs, cex=1, pos=3, offset=1.5)

legend(x='topleft',legend=c('HGG','Control'), pt.bg=c('grey','lightblue'), col='black', pch=22, cex=1, pt.cex=2, title="Disease state")

dev.off()

#-------
# [PLOT] compare sex DMPs b/w prepubescent and postpubescent subjects (stratified by disease state)
#-------

dd=data.frame(c(15,21),c(11,76))
colnames(dd)=c('Prepubescent','Postpubescent')
rownames(dd)=c('HGG','Control')

pdf('/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/DMPvisualization_promoter_enhancer_DMPs_pre_vs_postpubescent.pdf')

bp=barplot(height=t(dd), beside=T, ylab='Number of sex-DMPs', main='Effect of puberty on quantity of sex-DMPs (promoters/enhancers only)', col=c('grey','lightgreen'))
labs=c(15,11,21,76)
text(bp, 0, labs, cex=1, pos=3)

legend(x='topleft',legend=c('Prepubescent','Postpubescent'), pt.bg=c('grey','lightgreen'), col='black', pch=22, cex=1, pt.cex=2, title="Age group")

dev.off()

# p-values for the plot above

	# HGG strata:
obs=data.frame(c(15,11),c(185150-15,185150-11))
colnames(obs)=c('DMP','not_DMP')
rownames(obs)=c('prepubescent','post-pubescent')
chisq=chisq.test(obs)
chisq$p.value # 0.5562846

	# Control strata:
obs=data.frame(c(21,76),c(185150-21,185150-76))
colnames(obs)=c('DMP','not_DMP')
rownames(obs)=c('prepubescent','post-pubescent')
chisq=chisq.test(obs)
chisq$p.value # 4.167859e-08

#-------
# calculate the percent of DMPs that are lost b/w Controls and HGG (promoter/enhancer DMPs only)
#-------

# 185150 of the interrogated CpGs are in promoters or enhancers

# Prepubescent

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_prepubescent_cases.rda")
DMPs_prepubescent_cases=all_dmps
DMPs_prepubescent_cases=DMPs_prepubescent_cases[which(DMPs_prepubescent_cases$Promoter=='TRUE' | DMPs_prepubescent_cases$Enhancer=='TRUE'),]

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_prepubescent_controls.rda")
DMPs_prepubescent_controls=all_dmps
DMPs_prepubescent_controls=DMPs_prepubescent_controls[which(DMPs_prepubescent_controls$Promoter=='TRUE' | DMPs_prepubescent_controls$Enhancer=='TRUE'),]

mm=match(DMPs_prepubescent_controls$probe,DMPs_prepubescent_cases$probe)
(length(DMPs_prepubescent_controls$probe)-length(which(!is.na(mm))))/length(DMPs_prepubescent_controls$probe) # 0.5238095

obs=data.frame(c(nrow(DMPs_prepubescent_cases),185150-nrow(DMPs_prepubescent_cases)), c(nrow(DMPs_prepubescent_controls), 185150-nrow(DMPs_prepubescent_controls)))
colnames(obs)=c('HGG','Control')
rownames(obs)=c('DMP','not_DMP')
chisq=chisq.test(obs)
chisq$p.value # 0.4046339

# Postpubescent

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_postpubescent_cases.rda")
DMPs_postpubescent_cases=all_dmps
DMPs_postpubescent_cases=DMPs_postpubescent_cases[which(DMPs_postpubescent_cases$Promoter=='TRUE' | DMPs_postpubescent_cases$Enhancer=='TRUE'),]

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_postpubescent_controls.rda")
DMPs_postpubescent_controls=all_dmps
DMPs_postpubescent_controls=DMPs_postpubescent_controls[which(DMPs_postpubescent_controls$Promoter=='TRUE' | DMPs_postpubescent_controls$Enhancer=='TRUE'),]

mm=match(DMPs_postpubescent_controls$probe, DMPs_postpubescent_cases$probe)
(length(DMPs_postpubescent_controls$probe)-length(which(!is.na(mm))))/length(DMPs_postpubescent_controls$probe) # 0.8684211

obs=data.frame(c(nrow(DMPs_postpubescent_cases),185150-nrow(DMPs_postpubescent_cases)), c(nrow(DMPs_postpubescent_controls), 185150-nrow(DMPs_postpubescent_controls)))
colnames(obs)=c('HGG','Control')
rownames(obs)=c('DMP','not_DMP')
chisq=chisq.test(obs)
chisq$p.value # 6.774804e-12

# 0-21 yo

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_allPeds_cases.rda")
DMPs_allPeds_cases=all_dmps
DMPs_allPeds_cases=DMPs_allPeds_cases[which(DMPs_allPeds_cases$Promoter=='TRUE' | DMPs_allPeds_cases$Enhancer=='TRUE'),]

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_allPeds_controls.rda")
DMPs_allPeds_controls=all_dmps
DMPs_allPeds_controls=DMPs_allPeds_controls[which(DMPs_allPeds_controls$Promoter=='TRUE' | DMPs_allPeds_controls$Enhancer=='TRUE'),]

mm=match(DMPs_allPeds_controls$probe, DMPs_allPeds_cases$probe)
(length(DMPs_allPeds_controls$probe)-length(which(!is.na(mm))))/length(DMPs_allPeds_controls$probe) # 0.8655914

obs=data.frame(c(nrow(DMPs_allPeds_cases),185150-nrow(DMPs_allPeds_cases)), c(nrow(DMPs_allPeds_controls), 185150-nrow(DMPs_allPeds_controls)))
colnames(obs)=c('HGG','Control')
rownames(obs)=c('DMP','not_DMP')
chisq=chisq.test(obs)
chisq$p.value # 2.710792e-21

# > 21 yo

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_adult_cases.rda")
DMPs_adult_cases=all_dmps
DMPs_adult_cases=DMPs_adult_cases[which(DMPs_adult_cases$Promoter=='TRUE' | DMPs_adult_cases$Enhancer=='TRUE'),]

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_adult_controls.rda")
DMPs_adult_controls=all_dmps
DMPs_adult_controls=DMPs_adult_controls[which(DMPs_adult_controls$Promoter=='TRUE' | DMPs_adult_controls$Enhancer=='TRUE'),]

mm=match(DMPs_adult_controls$probe, DMPs_adult_cases$probe)
(length(DMPs_adult_controls$probe)-length(which(!is.na(mm))))/length(DMPs_adult_controls$probe) # 0.957162

obs=data.frame(c(nrow(DMPs_adult_cases),185150-nrow(DMPs_adult_cases)), c(nrow(DMPs_adult_controls), 185150-nrow(DMPs_adult_controls)))
colnames(obs)=c('HGG','Control')
rownames(obs)=c('DMP','not_DMP')
chisq=chisq.test(obs)
chisq$p.value # 2.701128e-285

##----------------------------------------
## visualize the number of significant sex-disease interaction CpGs by age group
##----------------------------------------

intrxn=data.frame(c(0,0,1,38))
rownames(intrxn)=c('Prepubescent','Postpubescent','All peds (0-21 yo)','>21 yo')

pdf('/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/num_intrxn_CpGs.pdf')

bp=barplot(height=t(intrxn), xlab="Age group", ylab='Number of sex-disease interaction probes', col='lightgrey', main='Sex-Disease Interaction Probes')
labs=c(0,0,1,38)
text(bp, 0, labs, cex=1, pos=3)

dev.off()

#####################################
## Analysis of genes that lose DNAm in HGG
#####################################
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

#load genes
load('/athena/masonlab/scratch/users/nai2008/gencode_v35_hg19_genes_and_genomicState.rda')

# Prepubescent -----------------
load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_prepubescent_cases.rda")
DMPs_prepubescent_cases=all_dmps

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_prepubescent_controls.rda")
DMPs_prepubescent_controls=all_dmps

lost_sex_DMPs_interrogation(cases=DMPs_prepubescent_cases, 
	controls=DMPs_prepubescent_controls, 
	lost_DMPs_overlapping_genes_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_prepubescent.csv', 
	lost_DMPs_overlapping_genes_enrichR_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_prepubescent_for_enrichR.csv',
	lost_DMPs_overlapping_imprinted_genes_filename= '/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_imprinted_genes_prepubescent.csv',
	genes=genes,
	ORA_table_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_prepubescent_ORA.csv')

# num_DMPs_lost                                              "53"
# proportion_DMPs_lost                                       "0.5955056"
# range_of_abs_value_of_DNAm_differenes_between_sexes        "0.012-0.28"
# how_many_lost_DMPs_overlap_genes_or_promoters_or_enhancers "52"
# how_many_lost_DMPs_overlap_TFs                             "9"
# how_many_TFs                                               "10"
# how_many_lost_DMPs_overlap_imprinted_genes                 "0"
# how_many_imprinted_genes                                   "0"
# number_sig_GO_or_KEGG_terms                                "0"
# number_sig_GO_or_KEGG_terms_unadjP                         "479"

#Postpubescent -----------------
load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_postpubescent_cases.rda")
DMPs_postpubescent_cases=all_dmps

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_postpubescent_controls.rda")
DMPs_postpubescent_controls=all_dmps

lost_sex_DMPs_interrogation(cases=DMPs_postpubescent_cases, 
	controls=DMPs_postpubescent_controls, 
	lost_DMPs_overlapping_genes_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_postpubescent.csv', 
	lost_DMPs_overlapping_genes_enrichR_filename= '/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_postpubescent_for_enrichR.csv',
	lost_DMPs_overlapping_imprinted_genes_filename= '/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_imprinted_genes_postpubescent.csv',
	genes=genes,
	ORA_table_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_postpubescent_ORA.csv')

# num_DMPs_lost                                              "231"
# proportion_DMPs_lost                                       "0.8369565"
# range_of_abs_value_of_DNAm_differenes_between_sexes        "0.0036-0.28"
# how_many_lost_DMPs_overlap_genes_or_promoters_or_enhancers "221"
# how_many_lost_DMPs_overlap_TFs                             "36"
# how_many_TFs                                               "31"
# how_many_lost_DMPs_overlap_imprinted_genes                 "0"
# how_many_imprinted_genes                                   "0"
# number_sig_GO_or_KEGG_terms                                "0"
# number_sig_GO_or_KEGG_terms_unadjP                         "1115"

#0-21 yo -----------------
load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_allPeds_cases.rda")
DMPs_allPeds_cases=all_dmps

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_allPeds_controls.rda")
DMPs_allPeds_controls=all_dmps

lost_sex_DMPs_interrogation(cases=DMPs_allPeds_cases, 
	controls=DMPs_allPeds_controls, 
	lost_DMPs_overlapping_genes_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_allPeds.csv', 
	lost_DMPs_overlapping_genes_enrichR_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_allPeds_for_enrichR.csv',
	lost_DMPs_overlapping_imprinted_genes_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_imprinted_genes_allPeds.csv',
	genes=genes,
	ORA_table_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_allPeds_ORA.csv')

# num_DMPs_lost                                              "594"
# proportion_DMPs_lost                                       "0.8546763"
# range_of_abs_value_of_DNAm_differenes_between_sexes        "0.0029-0.24"
# how_many_lost_DMPs_overlap_genes_or_promoters_or_enhancers "560"
# how_many_lost_DMPs_overlap_TFs                             "108"
# how_many_TFs                                               "75"
# how_many_lost_DMPs_overlap_imprinted_genes                 "0"
# how_many_imprinted_genes                                   "0"
# number_sig_GO_or_KEGG_terms                                "0"
# number_sig_GO_or_KEGG_terms_unadjP                         "1142"

#>21 yo -----------------
load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_adult_cases.rda")
DMPs_adult_cases=all_dmps

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_adult_controls.rda")
DMPs_adult_controls=all_dmps

lost_sex_DMPs_interrogation(cases=DMPs_adult_cases, 
	controls=DMPs_adult_controls, 
	lost_DMPs_overlapping_genes_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_adults.csv', 
	lost_DMPs_overlapping_genes_enrichR_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_adults_for_enrichR.csv',
	lost_DMPs_overlapping_imprinted_genes_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_imprinted_genes_adults.csv',
	genes=genes,
	ORA_table_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_adults_ORA.csv')

# num_DMPs_lost                                              "4521"
# proportion_DMPs_lost                                       "0.9615057"
# range_of_abs_value_of_DNAm_differenes_between_sexes        "0.00089-0.26"
# how_many_lost_DMPs_overlap_genes_or_promoters_or_enhancers "4287"
# how_many_lost_DMPs_overlap_TFs                             "844"
# how_many_TFs                                               "471"
# how_many_lost_DMPs_overlap_imprinted_genes                 "48"
# how_many_imprinted_genes                                   "29"
# number_sig_GO_or_KEGG_terms_FDR                            "1"
# number_sig_GO_or_KEGG_terms_unadjP                         "1147"

## Calculate the range of |average DNAm changes b/w sexes| at each age group

# Prepubescent
load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_prepubescent_cases.rda")
DMPs_prepubescent_cases=all_dmps

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_prepubescent_controls.rda")
DMPs_prepubescent_controls=all_dmps

cases=DMPs_prepubescent_cases
controls=DMPs_prepubescent_controls

signif(range(abs(cases$slope)),2)
signif(range(abs(controls$slope)),2)

# > signif(range(abs(cases$slope)),2)
# [1] 0.039 0.400
# > signif(range(abs(controls$slope)),2)
# [1] 0.012 0.470

# Postpubescent
load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_postpubescent_cases.rda")
DMPs_postpubescent_cases=all_dmps

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_postpubescent_controls.rda")
DMPs_postpubescent_controls=all_dmps

cases=DMPs_postpubescent_cases
controls=DMPs_postpubescent_controls

signif(range(abs(cases$slope)),2)
signif(range(abs(controls$slope)),2)

# > signif(range(abs(cases$slope)),2)
# [1] 0.037 0.380
# > signif(range(abs(controls$slope)),2)
# [1] 0.0036 0.4800

# 0-21
load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_allPeds_cases.rda")
DMPs_allPeds_cases=all_dmps

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_allPeds_controls.rda")
DMPs_allPeds_controls=all_dmps

cases=DMPs_allPeds_cases
controls=DMPs_allPeds_controls

signif(range(abs(cases$slope)),2)
signif(range(abs(controls$slope)),2)

# > signif(range(abs(cases$slope)),2)
# [1] 0.01 0.39
# > signif(range(abs(controls$slope)),2)
# [1] 0.0029 0.4800

# > 21
load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_adult_cases.rda")
DMPs_adult_cases=all_dmps

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_adult_controls.rda")
DMPs_adult_controls=all_dmps

cases=DMPs_adult_cases
controls=DMPs_adult_controls

signif(range(abs(cases$slope)),2)
signif(range(abs(controls$slope)),2)

# > signif(range(abs(cases$slope)),2)
# [1] 0.0031 0.4000
# > signif(range(abs(controls$slope)),2)
# [1] 0.00089 0.45000

#####################################
## MGMT differential methylation
#####################################

## How many sex-DMPs overlap MGMT?

# Prepubescent
dmps_overlapping_genes=read.csv('/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/DMPs_prepubescent_cases.csv')
length(unique(dmps_overlapping_genes$probe[which(dmps_overlapping_genes$gene=='MGMT')])) # 0

# Postpubescent
dmps_overlapping_genes=read.csv('/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/DMPs_postpubescent_cases.csv')
length(unique(dmps_overlapping_genes$probe[which(dmps_overlapping_genes$gene=='MGMT')])) # 0

# 0-21 yo
dmps_overlapping_genes=read.csv('/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/DMPs_allPeds_cases.csv')
length(unique(dmps_overlapping_genes$probe[which(dmps_overlapping_genes$gene=='MGMT')])) # 0

# >21 yo
dmps_overlapping_genes=read.csv('/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/DMPs_adult_cases.csv')
length(unique(dmps_overlapping_genes$probe[which(dmps_overlapping_genes$gene=='MGMT')])) # 0


## How many lost sex-DMPs overlap MGMT?

# Prepubescent
lost_dmps_overlapping_genes=read.csv('/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_prepubescent.csv')
length(unique(lost_dmps_overlapping_genes$probe[which(lost_dmps_overlapping_genes$gene=='MGMT')])) # 0

# Postpubescent
lost_dmps_overlapping_genes=read.csv('/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_postpubescent.csv')
length(unique(lost_dmps_overlapping_genes$probe[which(lost_dmps_overlapping_genes$gene=='MGMT')])) # 0

# 0-21 yo
lost_dmps_overlapping_genes=read.csv('/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_allPeds.csv')
length(unique(lost_dmps_overlapping_genes$probe[which(lost_dmps_overlapping_genes$gene=='MGMT')])) # 2
lost_dmps_overlapping_genes[which(lost_dmps_overlapping_genes$gene=='MGMT'),c(1,8,10)]
#          probe gene distance
# 316 cg02941816 MGMT        0
# 544 cg01341123 MGMT      382

# >21 yo
lost_dmps_overlapping_genes=read.csv('/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_adults.csv')
length(unique(lost_dmps_overlapping_genes$probe[which(lost_dmps_overlapping_genes$gene=='MGMT')])) # 6
lost_dmps_overlapping_genes[which(lost_dmps_overlapping_genes$gene=='MGMT'),c(1,8,10)]
#           probe gene distance
# 1397 cg02330106 MGMT      613
# 2236 cg12575438 MGMT      607
# 2392 cg25946389 MGMT      380
# 2603 cg02941816 MGMT        0
# 3809 cg01341123 MGMT      382
# 3839 cg14194875 MGMT      316

