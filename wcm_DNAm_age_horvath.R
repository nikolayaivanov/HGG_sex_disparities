#

#install.packages('RPMM')
library(RPMM)
#install.packages('sqldf')
library(sqldf)
#BiocManager::install("impute")
library(impute)
#install.packages('WGCNA')
library(WGCNA)
library(minfi)

source("/athena/masonlab/scratch/users/nai2008/DNAm/DNAm_age_Horvath/NORMALIZATION.R")

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/GRset_funNorm.rda")

# clean up pheno data
pd=as.data.frame(pData(GRset))
pd$Dx=as.vector(pd$Dx)
pd$Dx[which(pd$Dx=='glioblastoma')]='Glioblastoma'
pd$Dx=as.factor(pd$Dx)

Dx_simplified=pd$Dx
Dx_simplified[which(Dx_simplified !='Control')]='HGG'
pd$Dx_simplified=as.factor(as.vector(Dx_simplified))

pData(GRset)=as(pd,"DataFrame")

trafo= function(x,adult.age=20) { x=(x+1)/(1+adult.age); y=ifelse(x<=1, log( x),x-1);y }
anti.trafo= function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }

probeAnnotation21kdatMethUsed=read.csv("/athena/masonlab/scratch/users/nai2008/DNAm/DNAm_age_Horvath/probeAnnotation21kdatMethUsed.csv")
probeAnnotation27k=read.csv("/athena/masonlab/scratch/users/nai2008/DNAm/DNAm_age_Horvath/datMiniAnnotation27k.csv")
datClock=read.csv("/athena/masonlab/scratch/users/nai2008/DNAm/DNAm_age_Horvath/AdditionalFile3.csv")

p=getBeta(GRset)
pd=pData(GRset)

####################################################################
# run DNAm Age Clock on All Samples
####################################################################

dat0=p
dat0=cbind(Row.Names = rownames(dat0),dat0)

dat0=as.data.frame(dat0)

save(dat0, file="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/Horvath_DNAmAge_allSamples_dat0.rda")

nSamples=dim(dat0)[[2]]-1
nProbes= dim(dat0)[[1]]
# the following command may not be needed. But it is sometimes useful when you use read.csv.sql
dat0[,1]= gsub(x=dat0 [,1],pattern="\"",replacement="")
#Create a log file which will be output into your directory
# The code looks a bit complicated because it serves to create a log file (for error checks etc).
# It will automatically create a log file.
file.remove("LogFile.txt")
file.create("LogFile.txt")
DoNotProceed=FALSE
cat(paste( "The methylation data set contains", nSamples, "samples (e.g. arrays) and ", nProbes, " probes."),file="LogFile.txt")
if (nSamples==0) {DoNotProceed=TRUE; cat(paste( "\n ERROR: There must be a data input error since there seem to be no
samples.\n Make sure that you input a comma delimited file (.csv file)\n that can be read using the R command read.csv.sql .
Samples correspond to columns in that file ."), file="LogFile.txt",append=TRUE) }
if (nProbes==0) {DoNotProceed=TRUE; cat(paste( "\n ERROR: There must be a data input error since there seem to be zero
probes.\n Make sure that you input a comma delimited file (.csv file)\n that can be read using the R command read.csv.sql
CpGs correspond to rows.") , file="LogFile.txt",append=TRUE) }
if ( nSamples > nProbes ) { cat(paste( "\n MAJOR WARNING: It worries me a lot that there are more samples than CpG
probes.\n Make sure that probes correspond to rows and samples to columns.\n I wonder whether you want to first transpose
the data and then resubmit them? In any event, I will proceed with the analysis."),file="LogFile.txt",append=TRUE) }
if ( is.numeric(dat0[,1]) ) { DoNotProceed=TRUE; cat(paste( "\n Error: The first column does not seem to contain probe
identifiers (cg numbers from Illumina) since these entries are numeric values. Make sure that the first column of the file
contains probe identifiers such as cg00000292. Instead it contains ", dat0[1:3,1] ),file="LogFile.txt",append=TRUE) }
if ( !is.character(dat0[,1]) ) { cat(paste( "\n Major Warning: The first column does not seem to contain probe identifiers (cg
numbers from Illumina) since these entries are numeric values. Make sure that the first column of the file contains CpG probe
identifiers such as cg00000292. Instead it contains ", dat0[1:3,1] ),file="LogFile.txt",append=TRUE) }
datout=data.frame(Error=c("Input error. Please check the log file for details","Please read the instructions carefully."),
Comment=c("", "email Steve Horvath."))
if ( ! DoNotProceed ) {
nonNumericColumn=rep(FALSE, dim(dat0)[[2]]-1)
for (i in 2:dim(dat0)[[2]] ){ nonNumericColumn[i-1]=! is.numeric(dat0[,i]) }
if ( sum(nonNumericColumn) >0 ) { cat(paste( "\n MAJOR WARNING: Possible input error. The following samples contain nonnumeric
beta values: ", colnames(dat0)[-1][ nonNumericColumn], "\n Hint: Maybe you use the wrong symbols for missing data. 
4
Make sure to code missing values as NA in the Excel file. To proceed, I will force the entries into numeric values but make sure
this makes sense.\n" ),file="LogFile.txt",append=TRUE) }
XchromosomalCpGs=as.character(probeAnnotation27k$Name[probeAnnotation27k$Chr=="X"])
selectXchromosome=is.element(dat0[,1], XchromosomalCpGs )
selectXchromosome[is.na(selectXchromosome)]=FALSE
meanXchromosome=rep(NA, dim(dat0)[[2]]-1)
if ( sum(selectXchromosome) >=500 ) {
meanXchromosome= as.numeric(apply( as.matrix(dat0[selectXchromosome,-1]),2,mean,na.rm=TRUE)) }
if ( sum(is.na(meanXchromosome)) >0 ) { cat(paste( "\n \n Comment: There are lots of missing values for X chromosomal
probes for some of the samples. This is not a problem when it comes to estimating age but I cannot predict the gender of these
samples.\n " ),file="LogFile.txt",append=TRUE) }
match1=match(probeAnnotation21kdatMethUsed$Name , dat0[,1])
if ( sum( is.na(match1))>0 ) {
missingProbes= probeAnnotation21kdatMethUsed$Name[!is.element( probeAnnotation21kdatMethUsed$Name , dat0[,1])]
DoNotProceed=TRUE; cat(paste( "\n \n Input error: You forgot to include the following ", length(missingProbes), " CpG probes
(or probe names):\n ", paste( missingProbes, sep="",collapse=", ")),file="LogFile.txt",append=TRUE) } 

match1=match(probeAnnotation21kdatMethUsed$Name , dat0[,1])
if ( sum( is.na(match1))>0 ) stop(paste(sum( is.na(match1)), "CpG probes cannot be matched"))
dat1= dat0[match1,]
asnumeric1=function(x) {as.numeric(as.character(x))}
dat1[,-1]=apply(as.matrix(dat1[,-1]),2,asnumeric1)

set.seed(1)
# Do you want to normalize the data (recommended)?
normalizeData=TRUE
source("/athena/masonlab/scratch/users/nai2008/DNAm/DNAm_age_Horvath/StepwiseAnalysis.txt")

if ( sum( datout$Comment != "" ) ==0 ) { cat(paste( "\n The individual samples appear to be fine.
"),file="LogFile.txt",append=TRUE) }
if ( sum( datout$Comment != "" ) >0 ) { cat(paste( "\n Warnings were generated for the following samples.\n",
datout[,1][datout$Comment != ""], "\n Hint: Check the output file for more details."),file="LogFile.txt",append=TRUE) }
}

safety=datout

mm=match(datout$SampleID,pd$my_unique_ID)
datout=cbind(datout,pd[mm,])

all(as.vector(datout$SampleID)==as.vector(datout$my_unique_ID)) #TRUE

# label IDH status
IDH_Status_simplified=rep('WT',times=nrow(datout))
IDH_Status_simplified[grep('IDH|MUT|Mutation', datout$IDH_status)]='Mutated'
IDH_Status_simplified[grep('control', datout$IDH_status, ignore.case=TRUE)]='Control'
datout$IDH_Status_simplified=IDH_Status_simplified

colnames(datout)[which(colnames(datout)=='age')]='RealAge'

save(datout, file="/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/Horvath_DNAmAge_allSamples.rda")


####################################################################
# Plotting
####################################################################

source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

##################
# HGG cases (IDH wild type)
##################
load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/Horvath_DNAmAge_allSamples.rda")

pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/pdfs/Horvath_DNAmAge_HGG_IDHwt.pdf'

datout=datout[which(datout$IDH_Status_simplified == 'WT' & datout$Dx_simplified == 'HGG'),]

females_DNAmAge=datout[which(datout$sex=='Female'),]
males_DNAmAge=datout[which(datout$sex=='Male'),]

DNAmAge_plots(females_DNAmAge=females_DNAmAge, males_DNAmAge=males_DNAmAge, pdf_filename=pdf_filename)

##################
# HGG cases (IDH mutant)
##################

load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/Horvath_DNAmAge_allSamples.rda")

pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/pdfs/Horvath_DNAmAge_HGG_IDHmutant.pdf'

datout=datout[which(datout$IDH_Status_simplified == 'Mutated' & datout$Dx_simplified == 'HGG'),]

females_DNAmAge=datout[which(datout$sex=='Female'),]
males_DNAmAge=datout[which(datout$sex=='Male'),]

DNAmAge_plots(females_DNAmAge=females_DNAmAge, males_DNAmAge=males_DNAmAge, pdf_filename=pdf_filename)

##################
# Controls
##################
load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/Horvath_DNAmAge_allSamples.rda")

pdf_filename='/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/pdfs/Horvath_DNAmAge_controls.pdf'

datout=datout[which(datout$Dx_simplified == 'Control'),]

females_DNAmAge=datout[which(datout$sex=='Female'),]
males_DNAmAge=datout[which(datout$sex=='Male'),]

DNAmAge_plots(females_DNAmAge=females_DNAmAge, males_DNAmAge=males_DNAmAge, pdf_filename=pdf_filename)




