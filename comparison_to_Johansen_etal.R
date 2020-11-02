# compare our results to ones obtaioned by Johansen et al
# Johansen ML, Stetson LC, Vadmal V, Waite K, Berens ME, Connor JR, et al. 
# Gliomas display distinct sex-based differential methylation patterns based on molecular subtype. 
# Neurooncol Adv. 2020;2(1):vdaa002

########### DMPs ###############

# read in our adult HGG DMPs
load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_adult_cases.rda")
our_adult_dmps=all_dmps

# read in Johansen sex-DMP data (IDH wt adult GBM cases; see Table 1 in paper)
johansen_dmps=read.csv("/athena/masonlab/scratch/users/nai2008/DNAm/Johansen_etal_2020/DMPs_adults_IDHwt_GBM.csv")
range(johansen_dmps$adj.P.Val) # 1.440000e-25 4.997335e-02

mm=match(johansen_dmps$probe, our_adult_dmps$probe)
length(which(!is.na(mm))) # 83
nrow(johansen_dmps) # 311
length(which(!is.na(mm)))/nrow(johansen_dmps) # 0.266881

########### DMRs ###############

library(GenomicRanges)

# read in our adult HGG DMRs
load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/DMRs_adult_cases_1000bootstraps.rda")
sigDMRs=bumps$table[bumps$table$fwer<=0.1,]

gr_our=GRanges(seqnames=sigDMRs$chr, 
	ranges=IRanges(sigDMRs$start, sigDMRs$end))

# read in Johansen sex-DMR data (IDH wt adult GBM cases; see Table 1 in paper)
johansen_dmrs=read.csv("/athena/masonlab/scratch/users/nai2008/DNAm/Johansen_etal_2020/DMRs_adults_IDHwt_GBM.csv")

gr_jo=GRanges(seqnames=johansen_dmrs$chr, 
	ranges=IRanges(johansen_dmrs$DMR.start, johansen_dmrs$DMR.end))

oo=findOverlaps(gr_our, gr_jo, maxgap=0, ignore.strand=TRUE)
oo=as.data.frame(oo)

#load genes
load('/athena/masonlab/scratch/users/nai2008/gencode_v35_hg19_genes_and_genomicState.rda')

genes.df[which(genes.df$Geneid=='ENSG00000152467'),]



