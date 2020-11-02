## unzip the samples
# tar -C /athena/masonlab/scratch/users/nai2008/DNAm/Jaffe_NatNeuro_dlpfc_controls/Samples -xvf /athena/masonlab/scratch/users/nai2008/DNAm/Jaffe_NatNeuro_dlpfc_controls/GSE74193_RAW.tar
# gunzip -r /athena/masonlab/scratch/users/nai2008/DNAm/Jaffe_NatNeuro_dlpfc_controls/Samples

##
library(minfi)
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')
library(calibrate)

################################################
#### Read in controls (Jaffe); GEO: GSE74193
################################################

## read in phenotype data
pd=read.csv('https://raw.githubusercontent.com/andrewejaffe/devMeth450k/master/tables/phenotype_data_forGeo.csv')

## remove schizophrenia files
# to_rm=as.vector(pd$Chip[which(pd$Dx=='Schizo')])
# lf=list.files('/athena/masonlab/scratch/users/nai2008/DNAm/Jaffe_NatNeuro_dlpfc_controls/Samples', full.names=TRUE)

# sz_files=list()
# for(i in 1:length(to_rm)){
# gg=grep(to_rm[i], lf)
# sz_files[[i]]=lf[gg]
# }

# files=unlist(sz_files)
# files=paste(files, collapse=" ")
# rm_command=paste('rm', files)
# write.table(rm_command, file='/athena/masonlab/scratch/users/nai2008/DNAm/Jaffe_NatNeuro_dlpfc_controls/rm_Sz_samples.sh', quote=FALSE, col.names=FALSE, row.names=FALSE)
# # run the command above in unix

## read in the idat files
RGset=read.metharray.exp(base = "/athena/masonlab/scratch/users/nai2008/DNAm/Jaffe_NatNeuro_dlpfc_controls/Samples")

## organize pheno data
initial_colnames=colnames(RGset)
colnames(RGset)=paste(ss(initial_colnames,'_',2), ss(initial_colnames,'_',3),sep='_')

mm=match(colnames(RGset),pd$Chip)
pd=pd[mm,]
pd=as(pd,"DataFrame")
pData(RGset)=pd

#check that pd only has controls
which(pd$Dx !='Control') #0

#check that pd has age and sex labeled for all samples
which(is.na(pd$Age)) #0
which(is.na(pd$Gender)) #0

# `bestQC` column below indicates which samples to use to result in one array per subject
RGset=RGset[,pData(RGset)$bestQC]

#drop fetal samples
fetal=which(pData(RGset)$Age<0)
RGset=RGset[,-fetal]

# make sure there's 1 array per subject
pd=pData(RGset)
nrow(pd)==length(unique(pd$BrNum)) #TRUE

#finalize pheno data
pd=pData(RGset)
sex=as.vector(pd$Gender)
sex[which(sex=='M')]='Male'
sex[which(sex=='F')]='Female'
pd=data.frame(age=pd$Age, sex=sex,
	Sentrix_ID=pd$Sentrix_ID,
	Sentrix_Position=pd$Sentrix_Position,
	Dx=pd$Dx,
	Study_ID=rep('LIBD', times=nrow(pd)),
	my_unique_ID=paste0('LIBD_', 1:nrow(pd)),
	IDH_status='control',
	Sampling_time='control',
	NeuN_neg_proportion=pd$NeuN_neg,
	region='DLPFC')

pData(RGset)=as(pd, "DataFrame")
colnames(RGset)=pData(RGset)$my_unique_ID

save(RGset, file="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/RGset_jaffe_controls.rda")
#load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/RGset_jaffe_controls.rda") # annotation: ilmn12.hg19; 385 samples

################################################
## Read in the E_MTAB_5528 450k data
################################################

# read in the idat files
RGset=read.metharray.exp(base = "/athena/masonlab/scratch/users/nai2008/DNAm/DNAm_450k_E_MTAB_5528/RAW_files")

# read in the phenotype data
pd=read.delim('/athena/masonlab/scratch/users/nai2008/DNAm/DNAm_450k_E_MTAB_5528/E-MTAB-5528.sdrf.txt', header=TRUE)

pd_anno=read.delim('/athena/masonlab/scratch/users/nai2008/DNAm/DNAm_450k_E_MTAB_5528/E-MTAB-5528_Annotation.txt', header=TRUE, blank.lines.skip=TRUE, na.strings="")

#organize the pheno data
pd=pd[match(levels(pd$Source.Name),pd$Source.Name),]

pd$Assay.Name=ss(as.vector(pd$Assay.Name),'_Grn')

mm=match(colnames(RGset),pd$Assay.Name)

pd=pd[mm,]

nrow(pd)==length(unique(pd$Source.Name)) #TRUE

mm=match(as.vector(pd$Source.Name), as.vector(pd_anno$Index_ID))
pd_anno=pd_anno[mm,]

all(pd$Source.Name == pd_anno$Index_ID) # check if TRUE; it is
all(colnames(RGset) == pd$Assay.Name) # check if TRUE; it is

sex=as.vector(pd$Characteristics.sex.)
sex[which(sex=='male')]='Male'
sex[which(sex=='female')]='Female'

pd_organized=data.frame( age=as.vector(pd$Characteristics.age.),
	sex=sex,
	Sentrix_ID=ss(as.vector(pd$Assay.Name),'_',1),
	Sentrix_Position=ss(as.vector(pd$Assay.Name),'_',2),
	Dx=as.vector(pd_anno$Diagnosis),
	Study_ID=rep('5528'),
	my_unique_ID=paste0('hgg_5528_', 1:nrow(pd)),
	IDH_status=as.vector(pd_anno$IDH1_all),
	Sampling_time=as.vector(pd_anno$Sampling),
	NeuN_neg_proportion=NA,
	region=pd_anno$Site_2)

pData(RGset)=as(pd_organized,"DataFrame")

colnames(RGset)=pData(RGset)$my_unique_ID

# Some samples don't have sex labeled, drop them
drop=which(as.vector(pData(RGset)$sex)=="  ")
RGset=RGset[,-drop]
pData(RGset)$sex=as.factor(as.vector(pData(RGset)$sex))

# Some samples don't have age labeled, drop them
drop=which(is.na(as.vector(pData(RGset)$age)))
RGset=RGset[,-drop]

# Some samples don't have IDH status labeled, drop them
drop=which(is.na(as.vector(pData(RGset)$IDH_status)))
RGset=RGset[,-drop]

# Remove samples were not gather at diagnosis (but rather at time of second malignancy, relapse, or autopsy)
drop=which(as.vector(pData(RGset)$Sampling_time) != 'Diagnosis')
RGset=RGset[,-drop]

save(RGset, file="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/RGset_E_MTAB_5528.rda")
#load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/RGset_E_MTAB_5528.rda") # annotation: ilmn12.hg19; 91 samples

################################################
## Read in the GSE36278 450k data
################################################

# read in the idat files
RGset=read.metharray.exp(base = "/athena/masonlab/scratch/users/nai2008/DNAm/DNAm_450k_GSE36278/idat", recursive=TRUE)

# read in the phenotype data
pd_main=read.csv('/athena/masonlab/scratch/users/nai2008/DNAm/DNAm_450k_GSE36278/geo_metadata_clean.csv')
nrow(pd_main)==length(unique(pd_main$description)) #TRUE
colnames(pd_main)=c('description','age_yrs','location','sex','IDH_status','disease')
pd_sentrix=read.csv('/athena/masonlab/scratch/users/nai2008/DNAm/DNAm_450k_GSE36278/pGBM_CC_2012_sentrix_2018-06-18_DS.csv')

mm=match(pd_main$description,pd_sentrix$manuscript.ID)
pd_main=pd_main[-which(is.na(mm)),]

mm=match(pd_main$description,pd_sentrix$manuscript.ID)
pd_sentrix=pd_sentrix[mm,]

sex=as.vector(pd_main$sex)
sex[which(sex=='male')]='Male'
sex[which(sex=='female')]='Female'

pd_organized=data.frame( age=as.vector(pd_main$age_yrs),
	sex=sex,
	Sentrix_ID=ss(as.vector(pd_sentrix$sentrix.ID),'_',1),
	Sentrix_Position=ss(as.vector(pd_sentrix$sentrix.ID),'_',2),
	Dx=ss(as.vector(pd_main$disease)," ",2),
	Study_ID=rep('36278'),
	my_unique_ID=paste('hgg_36278_',1:nrow(pd_main),sep=""),
	IDH_status=as.vector(pd_main$IDH_status),
	Sampling_time='Diagnosis',
	NeuN_neg_proportion=NA,
	region=pd_main$location)

mm=match(colnames(RGset),pd_sentrix$sentrix.ID)
pd_organized=pd_organized[mm,]

pd_organized$Dx=as.vector(pd_organized$Dx)

pd_organized$Dx[which(pd_organized$region=='DIPG')]='Diffuse intrinsic pontine glioma'
pd_organized$Dx[which(pd_organized$Dx=='glioblastoma')]='Glioblastoma'

pData(RGset)=as(pd_organized,"DataFrame")

colnames(RGset)=pData(RGset)$my_unique_ID

# check is all samples have age, sex, and IDH_status labeled
which(is.na(as.vector(pData(RGset)$age))) #integer(0)
which(is.na(as.vector(pData(RGset)$sex))) #integer(0)
which(is.na(as.vector(pData(RGset)$IDH_status))) #integer(0)

save(RGset, file="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/RGset_GSE36278.rda")
#load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/RGset_GSE36278.rda") # annotation: ilmn12.hg19; 136 samples

################################################
## Read in the GSE55712 450k data
################################################

# read in the idat files
RGset=read.metharray.exp(base = "/athena/masonlab/scratch/users/nai2008/DNAm/DNAm_450k_GSE55712/idat")

# read in the phenotype data
pd=read.csv('/athena/masonlab/scratch/users/nai2008/DNAm/DNAm_450k_GSE55712/ACVR1_paper_450k_phenoData.csv')

nrow(pd)==length(unique(pd$File)) #TRUE

#organize the pheno data
if (length(which(is.na(as.vector(pd$sex)))) != 0) { pd=pd[-which(is.na(as.vector(pd$sex))),] }
if (length(which(is.na(as.vector(pd$age_yrs)))) != 0) { pd=pd[-which(is.na(as.vector(pd$age_yrs))),] }

sex=as.vector(pd$sex)
sex[which(sex=='F')]='Female'
sex[which(sex=='M')]='Male'

m=match(as.vector(pd$File),colnames(RGset))
RGset=RGset[,m]

pd_organized=data.frame( age=as.vector(pd$age_yrs),
	sex=sex,
	Sentrix_ID=as.vector(pd$Slide),
	Sentrix_Position=as.vector(pd$Array),
	Dx=rep('Glioblastoma',times=nrow(pd)),
	Study_ID=rep('55712'),
	my_unique_ID=paste('hgg_55712_',1:nrow(pd),sep=""),
	IDH_status=as.vector(pd$Phenotype),
	Sampling_time='Diagnosis',
	NeuN_neg_proportion=NA,
	region=pd$location)

pd_organized$Dx=as.vector(pd_organized$Dx)
pd_organized$Dx[which(pd_organized$location=='Pons')]='Diffuse intrinsic pontine glioma'

pData(RGset)=as(pd_organized,"DataFrame")

colnames(RGset)=pData(RGset)$my_unique_ID

# check is all samples have age, sex, and IDH_status labeled
which(is.na(as.vector(pData(RGset)$age))) #integer(0)
which(is.na(as.vector(pData(RGset)$sex))) #integer(0)
which(is.na(as.vector(pData(RGset)$IDH_status))) #integer(0)

save(RGset, file="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/RGset_GSE55712.rda")
#load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/RGset_GSE55712.rda") # annotation: ilmn12.hg19; 60 samples

################################################
## Read in TCGA GBM data
################################################

## download idat files from TCGA using the package TCGAbiolinks

#BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)

projects <- TCGAbiolinks:::getGDCprojects()$project_id
projects <- projects[grepl('^TCGA',projects,perl=T)]
match.file.cases.all <- NULL
for(proj in projects){
    print(proj)
    query <- GDCquery(project = proj,
                      data.category = "Raw microarray data",
                      data.type = "Raw intensities", 
                      experimental.strategy = "Methylation array", 
                      legacy = TRUE,
                      file.type = ".idat",
                      platform = "Illumina Human Methylation 450")
    match.file.cases <- getResults(query,cols=c("cases","file_name"))
    match.file.cases$project <- proj
    match.file.cases.all <- rbind(match.file.cases.all,match.file.cases)
    tryCatch(GDCdownload(query, method = "api", files.per.chunk = 20),
             error = function(e) GDCdownload(query, method = "client"))
}
# This will create a map between idat file name, cases (barcode) and project
readr::write_tsv(match.file.cases.all, path =  "idat_filename_case.txt")

#In unix:
# cd /home/nai2008
# mv /home/nai2008/GDCdata/TCGA-GBM/ /athena/masonlab/scratch/users/nai2008/DNAm/DNAm_TCGA_idat_files/
# mv /home/nai2008/idat_filename_case.txt /athena/masonlab/scratch/users/nai2008/DNAm/DNAm_TCGA_idat_files/
# find /athena/masonlab/scratch/users/nai2008/DNAm/DNAm_TCGA_idat_files/TCGA-GBM/legacy/Raw_microarray_data/Raw_intensities -type f -name "*.idat" -exec cp -f {} /athena/masonlab/scratch/users/nai2008/DNAm/DNAm_TCGA_idat_files/idat_files/ \;
# ^ to move insted of copy, replace '-exec cp' with '-exec mv'

# read in the idat files
RGset=read.metharray.exp(base = "/athena/masonlab/scratch/users/nai2008/DNAm/DNAm_TCGA_idat_files/idat_files/")

## match idat ids to TCGA ids
idat_filename_case=read.table('/athena/masonlab/scratch/users/nai2008/DNAm/DNAm_TCGA_idat_files/idat_filename_case.txt', header=TRUE)
idat_filename_case=idat_filename_case[-grep("Grn",idat_filename_case$file_name),]

# take only the samples from the primary solid tumor
# for more details, see the URLs below:
# https://docs.gdc.cancer.gov/Encyclopedia/pages/images/TCGA-TCGAbarcode-080518-1750-4378.pdf
# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
sample_types=ss(as.vector(idat_filename_case$cases),"-",4)
idat_filename_case=idat_filename_case[grep('01',sample_types),]
idat_filename_case$file_name=ss(as.vector(idat_filename_case$file_name),'_Red.idat',1)
idat_filename_case$submitter_id=paste(ss(as.vector(idat_filename_case$cases),"-",1),ss(as.vector(idat_filename_case$cases),"-",2),ss(as.vector(idat_filename_case$cases),"-",3), sep="-")

nrow(idat_filename_case) #140
length(unique(idat_filename_case$cases)) #140
length(unique(idat_filename_case$submitter_id)) #140
length(unique(idat_filename_case$file_name)) #140

mm=match(colnames(RGset),idat_filename_case$file_name)
RGset=RGset[,-which(is.na(mm))]

mm=match(colnames(RGset),idat_filename_case$file_name)

idat_filename_case=idat_filename_case[mm,]
colnames(RGset)=idat_filename_case$submitter_id

ncol(RGset)==length(unique(colnames(RGset))) #TRUE

# read in clinical data
clin=read.delim('/athena/masonlab/scratch/users/nai2008/DNAm/DNAm_TCGA_idat_files/clinical.tsv', na.strings="--", header=TRUE,)

# drop samples that don't have age or sex labeled
clin=clin[-which(is.na(clin$gender)),]
clin=clin[-which(is.na(clin$age_at_index)),]

mm=match(colnames(RGset),as.vector(clin$submitter_id))
RGset=RGset[,-which(is.na(mm))]

mm=match(colnames(RGset),as.vector(clin$submitter_id))

length(mm) #93
length(unique(mm)) #93

pd=clin[mm,]

all(colnames(RGset)==pd$submitter_id) #TRUE

# read in mutation data from the cBioPortal OncoPrint

cd=read.delim('/athena/masonlab/scratch/users/nai2008/DNAm/DNAm_TCGA_idat_files/PATIENT_DATA_oncoprint.tsv', na.strings="", header=FALSE, row.names=1)
cd=t(cd)
cd=as.data.frame(cd)
rownames(cd)=NULL

# leave only samples which were profiled for mutations
cd=cd[which(cd$Profiled_for_mutations=="Yes"),]

mm=match(colnames(RGset),cd$track_name)
RGset=RGset[,-which(is.na(mm))]
pd=pd[-which(is.na(mm)),]

mm=match(colnames(RGset),cd$track_name)
cd=cd[mm,]

all(colnames(RGset)==cd$track_name) #TRUE
all(colnames(RGset)==pd$submitter_id) #TRUE

pd$IHD1_mutations_status=as.vector(cd$IDH1_MUTATIONS)

# Sampling_time=Primary_Solid_Tumor

sex=as.vector(pd$gender)
sex[which(sex=='female')]='Female'
sex[which(sex=='male')]='Male'

nrow(pd) #88
length(unique(pd$submitter_id)) #88

pd_organized=data.frame(
	age=as.vector(pd$age_at_index),
	sex=sex,
	Sentrix_ID=NA,
	Sentrix_Position=NA,
	Dx=as.vector(pd$primary_diagnosis),
	Study_ID=rep('TCGA'),
	my_unique_ID=paste('TCGA',1:nrow(pd),sep="_"),
	IDH_status=as.vector(pd$IHD1_mutations_status),
	Sampling_time='Primary_tumor',
	NeuN_neg_proportion=NA,
	region=pd$site_of_resection_or_biopsy)

pData(RGset)=as(pd_organized,"DataFrame")

colnames(RGset)=pData(RGset)$my_unique_ID

save(RGset, file="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/RGset_TCGA.rda")
#load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/RGset_TCGA.rda") # annotation: ilmn12.hg19; 88 samples

################################################
## Downstream analysis
################################################

# combine the datsets

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/RGset_jaffe_controls.rda")
RGset_jaffe_controls=RGset
load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/RGset_E_MTAB_5528.rda")
RGset_E_MTAB_5528=RGset
load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/RGset_GSE36278.rda")
RGset_GSE36278=RGset
load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/RGset_GSE55712.rda")
RGset_GSE55712=RGset
load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/RGset_TCGA.rda")
RGset_TCGA=RGset

RGset_combined=combineArrays(RGset_jaffe_controls,RGset_E_MTAB_5528, outType = c("IlluminaHumanMethylation450k"))
RGset_combined=combineArrays(RGset_combined, RGset_GSE36278, outType = c("IlluminaHumanMethylation450k"))
RGset_combined=combineArrays(RGset_combined, RGset_GSE55712, outType = c("IlluminaHumanMethylation450k"))
RGset_combined=combineArrays(RGset_combined, RGset_TCGA, outType = c("IlluminaHumanMethylation450k"))
RGset=RGset_combined

# remove DIPGs from the analysis
pd=as.data.frame(pData(RGset))
drop=which(pd$Dx=='Diffuse intrinsic pontine glioma')
RGset=RGset[,-drop]

save(RGset, file="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/RGset_combined.rda")
#load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/RGset_combined.rda")

#generate a MethylSet
MSet = preprocessRaw(RGset)

# QC
qc = getQC(MSet)
qc.df=as.data.frame(qc)

pd=pData(MSet)

# check if any of the samples are duplicated, and remove them
pd$intensity=paste(qc.df[,1],qc.df[,2],sep='_')
write.csv(pd, file="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/pd.csv")

intensity=paste(qc.df[,1],qc.df[,2],sep='_')
test=duplicated(intensity)

duplicated.sample=vector()
num.test=vector()

for(i in 1: length(test)){

	if(test[i]=='FALSE') {
		duplicated.sample[i]='FALSE'
	} else if (test[i]=='TRUE') {
		num.test[i]=length(which(pd$intensity==pd$intensity[i]))
		sex_check=length(unique(pd$age[which(pd$intensity==pd$intensity[i])]))==1
		age_check=length(unique(pd$sex[which(pd$intensity==pd$intensity[i])]))==1
		if (sex_check=='TRUE' & age_check=='TRUE') {
			duplicated.sample[i]='TRUE'
		} else { duplicated.sample[i]='FALSE'
		}
	}
}

length(which(duplicated.sample=='TRUE')) # 13

RGset=RGset[,-which(duplicated.sample=='TRUE')]

# generate mMed vs uMed plots
MSet = preprocessRaw(RGset)
qc = getQC(MSet)
qc.df=as.data.frame(qc)

pd=pData(MSet)

intensity=paste(qc.df[,1],qc.df[,2],sep='_')
length(which(duplicated(intensity)=='TRUE')) # 2

pdf('/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/QC_mMed_vs_uMed.pdf')

study=as.vector(pd$Study_ID)

study=gsub('LIBD','red',study)
study=gsub('TCGA','blue',study)
study=gsub('36278','darkgreen',study)
study=gsub('5528','gold',study)
study=gsub('55712','cyan',study)

# study=as.factor(as.vector(pd$Study_ID))
# col.rainbow=rainbow(length(levels(study)))
# palette(col.rainbow)

plot(qc.df$mMed, qc.df$uMed, pch=21, bg=study, col='black', xlab='mMed', ylab='uMed', cex=1.4)
legend(x='topleft',legend=c('LIBD','TCGA','36278','5528','55712'), pch=22, col='black', pt.bg=c('red','blue','darkgreen','gold','cyan'), cex=1.2, pt.cex=2, title="Study ID")
#legend(x='topleft',legend=c('LIBD','TCGA','36278','5528','55712'), pch=15, col=c('red','blue','darkgreen','gold','cyan'), cex=1, title="Study ID")
#legend(x='topleft',legend=levels(study), pch=15, col=col.rainbow, cex=.7, title="Study ID")

plot(qc.df$mMed, qc.df$uMed, pch=21, bg='gray', col='black',,xlab='mMed', ylab='uMed', cex=1.3)
textxy(qc.df$mMed,qc.df$uMed,pd$my_unique_ID,cex =.7, offset = .7)

dev.off()

# look at Beta density plot for each sample

beta=getBeta(MSet)

pdf('/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/QC_density_plots.pdf')

densityPlot(MSet)

for (i in 1:ncol(beta)) {

	b=as.vector(beta[,i])
	
	if (length(which(is.na(b)) != 0)){ b=b[-which(is.na(b))] }

	plot(density(b),main=colnames(beta)[i],xlab='Beta')


}

dev.off()

# based on the {mMed vs uMed} & {Beta density} plots, drop the following samples:
	# hgg_5528_50
	# hgg_5528_54
	# hgg_5528_59
	# hgg_5528_98
	# hgg_55712_27

drop=match(c("hgg_5528_50", "hgg_5528_54", "hgg_5528_59","hgg_5528_98","hgg_55712_27"), pData(MSet)$my_unique_ID)
MSet=MSet[,-drop]
RGset=RGset[,-drop]

# check that the reported sex is correct
GMSet=mapToGenome(MSet)
check_sex=getSex(GMSet, cutoff=-2)
ps=as.vector(check_sex$predictedSex)
ps[which(ps=="M")]='Male'
ps[which(ps=="F")]='Female'
length(which((as.vector(GMSet$sex)==ps)=='FALSE')) # in 18 out of 656 samples, the predicted sex disagrees with labeled sex
ss=which((as.vector(GMSet$sex)==ps)=='FALSE')

GMSet=addSex(GMSet)
pd=pData(GMSet)

pdf('/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/disputed_sex_plots.pdf')

study=as.factor(pd$Study_ID)
col.rainbow=rainbow(length(levels(study)))
palette(col.rainbow)
plot(log2(pd$xMed), log2(pd$yMed), pch=21, bg=study, col="black", xlab=bquote(log[2](xMed)), ylab=bquote(log[2](yMed)))
legend(x='topleft',legend=levels(study), col=col.rainbow, pch=15, cex=1, title="Study")

sex=as.factor(pd$predictedSex)
col=c('pink','blue')
palette(col)
plot(log2(pd$xMed), log2(pd$yMed), pch=21, bg=sex, col="black", xlab=bquote(log[2](xMed)), ylab=bquote(log[2](yMed)))
legend(x='topleft',legend=levels(sex), col=col, pch=15, cex=1, title="Predicted sex")

sex=as.factor(pd$predictedSex)
col=c('pink','blue')
palette(col)
plot(log2(pd$xMed), log2(pd$yMed), pch=21, bg=sex, col="black", xlab=bquote(log[2](xMed)), ylab=bquote(log[2](yMed)))
points(log2(pd$xMed[ss]), log2(pd$yMed[ss]), pch=8, col='red')
legend(x='topleft',legend=levels(sex), col=col, pch=15, cex=1, title="Predicted sex")
legend(x='bottomright', legend='Disputed sex', col='red', pch=8)

dev.off()

# remove samples in which the predicted sex disagrees with labeled sex
drop=which((as.vector(GMSet$sex)==ps)=='FALSE')
GMSet=GMSet[,-drop]
MSet=MSet[,-drop]
RGset=RGset[,-drop]

save(RGset, file="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/RGset_HGGandControls_lowQCandSexDisagreement_samples_removed.rda")

# compute neg control PCs
controlProbes = minfi:::.extractFromRGSet450k(RGset)
negControlPCs = prcomp(t(log2(rbind(controlProbes$greenControls$NEGATIVE, controlProbes$redControls$NEGATIVE)+1)))
pca_vars=getPcaVars(negControlPCs) # 7.34e+01 1.02e+01 3.00e+00 9.35e-01 3.98e-01 3.33e-01 1.86e-01 1.73e-01 1.51e-01 1.22e-01
cs=cumsum(pca_vars) # 73.40000 83.60000 86.60000 87.53500 87.93300 88.26600 88.45200 88.62500 88.77600 88.89800

# figure out how many PCs to use
pd=pData(RGset)

Sentrix_ID=paste(pd$Sentrix_ID, pd$Study_ID, sep="_")
Sentrix_Position=paste(pd$Sentrix_Position, pd$Study_ID, sep="_")

pdf('/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/pdfs/neg_control_PCs.pdf')

boxplot(negControlPCs$x[,1]~Sentrix_ID, las=2, xlab='Sentrix_ID', ylab='PC1', main='PC vs. Sentrix_ID')
boxplot(negControlPCs$x[,1]~Sentrix_Position, las=2, xlab='Sentrix_Position', ylab='PC1', main='PC vs. Sentrix_Position')
boxplot(negControlPCs$x[,1]~pd$Study_ID, las=2, xlab='Study_ID', ylab='PC1', main='PC vs. Study_ID')

boxplot(negControlPCs$x[,2]~Sentrix_ID, las=2, xlab='Sentrix_ID', ylab='PC2', main='PC vs. Sentrix_ID')
boxplot(negControlPCs$x[,2]~Sentrix_Position, las=2, xlab='Sentrix_Position', ylab='PC2', main='PC vs. Sentrix_Position')
boxplot(negControlPCs$x[,2]~pd$Study_ID, las=2, xlab='Study_ID', ylab='PC2', main='PC vs. Study_ID')

boxplot(negControlPCs$x[,3]~Sentrix_ID, las=2, xlab='Sentrix_ID', ylab='PC3', main='PC vs. Sentrix_ID')
boxplot(negControlPCs$x[,3]~Sentrix_Position, las=2, xlab='Sentrix_Position', ylab='PC3', main='PC vs. Sentrix_Position')
boxplot(negControlPCs$x[,3]~pd$Study_ID, las=2, xlab='Study_ID', ylab='PC3', main='PC vs. Study_ID')

boxplot(negControlPCs$x[,4]~Sentrix_ID, las=2, xlab='Sentrix_ID', ylab='PC4', main='PC vs. Sentrix_ID')
boxplot(negControlPCs$x[,4]~Sentrix_Position, las=2, xlab='Sentrix_Position', ylab='PC4', main='PC vs. Sentrix_Position')
boxplot(negControlPCs$x[,4]~pd$Study_ID, las=2, xlab='Study_ID', ylab='PC4', main='PC vs. Study_ID')

boxplot(negControlPCs$x[,5]~Sentrix_ID, las=2, xlab='Sentrix_ID', ylab='PC5', main='PC vs. Sentrix_ID')
boxplot(negControlPCs$x[,5]~Sentrix_Position, las=2, xlab='Sentrix_Position', ylab='PC5', main='PC vs. Sentrix_Position')
boxplot(negControlPCs$x[,5]~pd$Study_ID, las=2, xlab='Study_ID', ylab='PC5', main='PC vs. Study_ID')

boxplot(negControlPCs$x[,6]~Sentrix_ID, las=2, xlab='Sentrix_ID', ylab='PC6', main='PC vs. Sentrix_ID')
boxplot(negControlPCs$x[,6]~Sentrix_Position, las=2, xlab='Sentrix_Position', ylab='PC6', main='PC vs. Sentrix_Position')
boxplot(negControlPCs$x[,6]~pd$Study_ID, las=2, xlab='Study_ID', ylab='PC6', main='PC vs. Study_ID')

boxplot(negControlPCs$x[,7]~Sentrix_ID, las=2, xlab='Sentrix_ID', ylab='PC7', main='PC vs. Sentrix_ID')
boxplot(negControlPCs$x[,7]~Sentrix_Position, las=2, xlab='Sentrix_Position', ylab='PC7', main='PC vs. Sentrix_Position')
boxplot(negControlPCs$x[,7]~pd$Study_ID, las=2, xlab='Study_ID', ylab='PC7', main='PC vs. Study_ID')

boxplot(negControlPCs$x[,8]~Sentrix_ID, las=2, xlab='Sentrix_ID', ylab='PC8', main='PC vs. Sentrix_ID')
boxplot(negControlPCs$x[,8]~Sentrix_Position, las=2, xlab='Sentrix_Position', ylab='PC8', main='PC vs. Sentrix_Position')
boxplot(negControlPCs$x[,8]~pd$Study_ID, las=2, xlab='Study_ID', ylab='PC8', main='PC vs. Study_ID')

boxplot(negControlPCs$x[,9]~Sentrix_ID, las=2, xlab='Sentrix_ID', ylab='PC9', main='PC vs. Sentrix_ID')
boxplot(negControlPCs$x[,9]~Sentrix_Position, las=2, xlab='Sentrix_Position', ylab='PC9', main='PC vs. Sentrix_Position')
boxplot(negControlPCs$x[,9]~pd$Study_ID, las=2, xlab='Study_ID', ylab='PC9', main='PC vs. Study_ID')

boxplot(negControlPCs$x[,10]~Sentrix_ID, las=2, xlab='Sentrix_ID', ylab='PC10', main='PC vs. Sentrix_ID')
boxplot(negControlPCs$x[,10]~Sentrix_Position, las=2, xlab='Sentrix_Position', ylab='PC10', main='PC vs. Sentrix_Position')
boxplot(negControlPCs$x[,10]~pd$Study_ID, las=2, xlab='Study_ID', ylab='PC10', main='PC vs. Study_ID')

boxplot(negControlPCs$x[,11]~Sentrix_ID, las=2, xlab='Sentrix_ID', ylab='PC11', main='PC vs. Sentrix_ID')
boxplot(negControlPCs$x[,11]~Sentrix_Position, las=2, xlab='Sentrix_Position', ylab='PC11', main='PC vs. Sentrix_Position')
boxplot(negControlPCs$x[,11]~pd$Study_ID, las=2, xlab='Study_ID', ylab='PC11', main='PC vs. Study_ID')

boxplot(negControlPCs$x[,12]~Sentrix_ID, las=2, xlab='Sentrix_ID', ylab='PC12', main='PC vs. Sentrix_ID')
boxplot(negControlPCs$x[,12]~Sentrix_Position, las=2, xlab='Sentrix_Position', ylab='PC12', main='PC vs. Sentrix_Position')
boxplot(negControlPCs$x[,12]~pd$Study_ID, las=2, xlab='Study_ID', ylab='PC12', main='PC vs. Study_ID')

boxplot(negControlPCs$x[,13]~Sentrix_ID, las=2, xlab='Sentrix_ID', ylab='PC13', main='PC vs. Sentrix_ID')
boxplot(negControlPCs$x[,13]~Sentrix_Position, las=2, xlab='Sentrix_Position', ylab='PC13', main='PC vs. Sentrix_Position')
boxplot(negControlPCs$x[,13]~pd$Study_ID, las=2, xlab='Study_ID', ylab='PC13', main='PC vs. Study_ID')

boxplot(negControlPCs$x[,14]~Sentrix_ID, las=2, xlab='Sentrix_ID', ylab='PC14', main='PC vs. Sentrix_ID')
boxplot(negControlPCs$x[,14]~Sentrix_Position, las=2, xlab='Sentrix_Position', ylab='PC14', main='PC vs. Sentrix_Position')
boxplot(negControlPCs$x[,14]~pd$Study_ID, las=2, xlab='Study_ID', ylab='PC14', main='PC vs. Study_ID')

boxplot(negControlPCs$x[,15]~Sentrix_ID, las=2, xlab='Sentrix_ID', ylab='PC15', main='PC vs. Sentrix_ID')
boxplot(negControlPCs$x[,15]~Sentrix_Position, las=2, xlab='Sentrix_Position', ylab='PC15', main='PC vs. Sentrix_Position')
boxplot(negControlPCs$x[,15]~pd$Study_ID, las=2, xlab='Study_ID', ylab='PC15', main='PC vs. Study_ID')

dev.off()

# use 9 PCs
negControlPCs_use=negControlPCs$x[,1:9]

colnames(negControlPCs_use) = paste0("negControl_", colnames(negControlPCs_use))
pd=pData(RGset)
pd = cbind(pd, negControlPCs_use)
pData(RGset)=pd

# normalize via funcional normalization
GRset=preprocessFunnorm(RGset)

save(GRset, file="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/GRset_funNorm.rda")
# load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/GRset_funNorm.rda")

# remove probes that are common SNPs
GRset # 485512 probes; 638 samples

getSnpInfo_mine = function (object, snpAnno = NULL)
{
    av <- minfi:::.availableAnnotation(object)
    if (is.null(snpAnno)) {
        snpAnno <- grep(pattern = "^SNPs\\.", x = getAnnotationObject(object)@defaults,
            value = TRUE)
    } else {
        snpAnno <- sub("^SNPs\\.", "", snpAnno)
        if (!snpAnno %in% av$classChoices$SNPs) {
            stop(sprintf("snpAnno '%s' is not part of the annotation",
                snpAnno))
        } else {
            snpAnno <- sprintf("SNPs.%s", snpAnno)
        }
    }
    snps <- getAnnotation(object, what = snpAnno)
    snps
}

dropLociWithSnps_mine= function (object, snps = c("CpG", "SBE"), maf = 0, snpAnno = NULL)
{
    minfi:::.isGenomicOrStop(object)
    maf_cols <- paste0(snps, "_maf")
    snpDF <- getSnpInfo_mine(object, snpAnno = snpAnno)
    choices <- c("Probe_maf", "CpG_maf", "SBE_maf")
    if (!all(choices %in% colnames(snpDF))) {
        stop("The specificed 'snpAnno' is not supported by this function")
    }
    if (sum(!(maf_cols %in% choices)) > 0) {
        stop("snps vector argument must be a combination of  \"Probe\", ",
            "\"CpG\" and \"SBE\"")
    }
    if (!is.numeric(maf) || maf < 0 || maf > 1) {
        stop("maf argument must be a numeric value between 0 and 1")
    }
    wh <- Reduce(union, lapply(maf_cols, function(xx) {
        which(snpDF[, xx] >= maf)
    }))
    wh <- sort(wh)
    if (length(wh) == 0)
        return(object)
    object[-wh, ]
}

GRset=dropLociWithSnps_mine(object=GRset, maf=0, snpAnno = "SNPs.147CommonSingle") # 468392 probes (17120 SNP probes removed)

# drop sex chr
anno=getAnnotation(GRset)
GRset=GRset[! anno$chr %in% c("chrX","chrY"),] # 456928 probes, 638 samples

# clean up pheno data
Dx_simplified=as.vector(pData(GRset)$Dx)
Dx_simplified[which(Dx_simplified !='Control')]='HGG'
pData(GRset)$Dx_simplified=as.factor(as.vector(Dx_simplified))

save(GRset, file="/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")
# load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/GRset_normalized_SNPs_and_XY_removed.rda")

# make a csv of phenotype data
pd=pData(GRset)
IDH_Status=rep('WT',times=nrow(pd))
IDH_Status[grep('IDH|MUT|Mutation', pd$IDH_status)]='Mutated'

study=rep(NA, times=nrow(pd))
study[grep('LIBD', pd$my_unique_ID)]='GEO_GSE74193'
study[grep('TCGA', pd$my_unique_ID)]='TCGA'
study[grep('hgg_5528', pd$my_unique_ID)]='ArrayExpress_E-MTAB-5528'
study[grep('hgg_36278', pd$my_unique_ID)]='GEO_GSE36278'
study[grep('hgg_55712', pd$my_unique_ID)]='GEO_GSE55712'

pheno_data=data.frame(age=pd$age,
sex=pd$sex,
Dx=pd$Dx,
region=pd$region,
IDH_status=IDH_Status,
sample_ID=pd$my_unique_ID,
study=study)

write.csv(pheno_data, file='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/pheno_data_for_all_DNAm_samples.csv', row.names=FALSE)

### NAI

