#
#source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

# DMP validation function
validate_dmps = function(GRset.funnorm,studyID) {

GRset.funnorm_subset=GRset.funnorm[,which(pData(GRset.funnorm)$studyID==studyID)]
beta=getBeta(GRset.funnorm_subset)

pd=pData(GRset.funnorm_subset)

mod=model.matrix(~as.factor(pd$sex))
svaobj = sva(beta, mod)

mod = model.matrix(~as.factor(pd$sex))
mod=cbind(mod,svaobj$sv)
probe_fit=lmFit(object=beta,design=mod)
eb=eBayes(probe_fit)

slope=probe_fit$coefficients[,2]
intercept=probe_fit$coefficients[,1]
p.value=eb$p.value[,2]
t=eb$t[,2]
fdr=p.adjust(as.vector(eb$p.value[,2]),method='fdr')

probeFit=data.frame(probes=rownames(eb$p.value), slope=slope, intercept=intercept, 
  p.value=p.value, fdr=fdr, t=t)
rownames(probeFit)=NULL

file=paste('/athena/masonlab/scratch/users/nai2008/DNAm/rdas/probeFit_',studyID,'.rda',sep='')
save(probeFit,file=file)

dmps=which(probeFit$fdr<=.05)

return(length(dmps))

}

# Gene Ontology functions

dogo_entrez <- function(names,universe,species="human", goP = 0.01, #input genes (user set and universe) as entrez
  cond=FALSE, ontology = "BP"){
    if(species=="human"){
    golib="org.Hs.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Hs.egREFSEQ2EG
  } else  if (species == "mouse") {
    golib="org.Mm.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Mm.egREFSEQ2EG
  } else if (species == "rat") {
    golib="org.Rn.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Rn.egREFSEQ2EG
  }
  require(GOstats)
  x=names
  x=x[!is.na(x)]

  params <- new("GOHyperGParams", geneIds = unique(x),
                universeGeneIds = universe,
                annotation = golib,
                ontology = ontology, pvalueCutoff = goP, conditional = cond,
                testDirection="over")
  
  ht=hyperGTest(params)
  tab=summary(ht)
  tmp1=geneIdsByCategory(ht)
  tmp1=tmp1[tab[,1]]
  tab$IDs=sapply(tmp1,function(y) paste(names(x)[x%in%y],collapse=";"))
  return(tab)

}

dogo_RefSeq <- function(names,universe,species="human", goP = 0.01, #input genes (user set and universe) as RefSeq
	cond=FALSE, ontology = "BP"){
    if(species=="human"){
		golib="org.Hs.eg.db"
		library(golib,character.only=TRUE)
		gomap= org.Hs.egREFSEQ2EG
  } else  if (species == "mouse") {
		golib="org.Mm.eg.db"
		library(golib,character.only=TRUE)
		gomap= org.Mm.egREFSEQ2EG
  } else if (species == "rat") {
		golib="org.Rn.eg.db"
		library(golib,character.only=TRUE)
		gomap= org.Rn.egREFSEQ2EG
	}
  require(GOstats)
  x=unlist(mget(as.character(names), gomap,ifnotfound = NA)) # convert inputted RefSeq genes (user set) to entrez 
  x=x[!is.na(x)]
  Universe=unlist(mget(as.character(universe),gomap,ifnotfound = NA)) # convert inputted RefSeq genes (universe) to entrez  
  Universe=unique(c(Universe[!is.na(Universe)],unique(x)))

  params <- new("GOHyperGParams", geneIds = unique(x),
                universeGeneIds = Universe,
                annotation = golib,
                ontology = ontology, pvalueCutoff = goP, conditional = cond,
                testDirection="over")
  ht=hyperGTest(params)
  tab=summary(ht)
  tmp1=geneIdsByCategory(ht)
  tmp1=tmp1[tab[,1]]
  tab$IDs=sapply(tmp1,function(y) paste(names(x)[x%in%y],collapse=";"))
  return(tab)

}

dogo_GeneSymbols <- function(names,universe,species="human", goP = 0.01, #input genes (user set and universe) as Gene Symbol
  cond=FALSE, ontology = "BP"){
    if(species=="human"){
    golib="org.Hs.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Hs.egSYMBOL2EG
  } else  if (species == "mouse") {
    golib="org.Mm.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Mm.egSYMBOL2EG
  } else if (species == "rat") {
    golib="org.Rn.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Rn.egSYMBOL2EG
  }
  require(GOstats)
  x=unlist(mget(as.character(names), gomap,ifnotfound = NA)) # convert inputted Gene Symbols (user set) to entrez 
  x=x[!is.na(x)]
  Universe=unlist(mget(as.character(universe),gomap,ifnotfound = NA)) # convert inputted Gene Symbols (universe) to entrez  
  Universe=unique(c(Universe[!is.na(Universe)],unique(x)))

  params <- new("GOHyperGParams", geneIds = unique(x),
                universeGeneIds = Universe,
                annotation = golib,
                ontology = ontology, pvalueCutoff = goP, conditional = cond,
                testDirection="over")
  ht=hyperGTest(params)
  tab=summary(ht)
  tmp1=geneIdsByCategory(ht)
  tmp1=tmp1[tab[,1]]
  tab$IDs=sapply(tmp1,function(y) paste(names(x)[x%in%y],collapse=";"))
  return(tab)

}

# wrapper for string split and sapply
ss = function(x, pattern, slot=1,...) sapply(strsplit(x,pattern,...), "[", slot)

# command that will return % of variance explained by each PC
getPcaVars = function(pca)  signif(((pca$sdev)^2)/(sum((pca$sdev)^2)),3)*100

# function that takes 2 arguments x1 and x2, and matches all elements of x1 to ALL elements of x2, and returns an index vector corresponding to x2

match_all = function (x1, x2){

indexes=list()

  for (i in 1:length(x1)){ 

    aa = which( x2 %in% x1[i]) 

    if(length(aa)>0){indexes[[i]] = aa 
      } else if (length(aa)==0) { indexes[[i]] = NA }
  }

oo=unlist(indexes)

return(oo)

}

# function to generate a summary file a from fastqc files

fastqcSummary <- function(files){ 
#files=paths to unzipped fastqc folders 

ss = function(x, pattern, slot=1,...) sapply(strsplit(x,pattern,...), "[", slot)

txt_files=paste0(files,'/fastqc_data.txt')

  for (i in 1:length(txt_files)){
  
    scan=as.vector(sapply(txt_files[i], function(x) scan(x,"",sep="\n")))

    ## Total Sequences
    ts=scan[grep('Total Sequences',scan)]
    Total_Sequences=as.numeric(ss(ts,'Sequences\t',2))

    ## Per base sequence quality 
    Per_base_sequence_quality=read.table(txt_files[i],sep='\t',comment.char="",skip=grep('Per base sequence quality',scan)+1, nrows=(grep('Per tile sequence quality',scan)-grep('Per base sequence quality',scan)-3))
    if (ncol(Per_base_sequence_quality)==7){
      colnames(Per_base_sequence_quality)=c('Base', 'Mean', 'Median', 'Lower_Quartile', 'Upper_Quartile', '10th_Percentile', '90th_Percentile') 
    }

    ## Adapter Content
    Adapter_Content=read.table(txt_files[i],sep='\t',comment.char="",skip=grep('Adapter Content',scan)+1, nrows=(grep('Kmer Content',scan)-grep('Adapter Content',scan)-3)) 
    if(ncol(Adapter_Content)==4){
    colnames(Adapter_Content)=c('Position', 'Illumina_Universal_Adapter', 'Illumina_Small_RNA_Adapter', 'Nextera_Transposase_Sequence')
    } else if (ncol(Adapter_Content)==5){
    colnames(Adapter_Content)=c('Position', 'Illumina_Universal_Adapter', 'Illumina_Small_RNA_Adapter', 'Nextera_Transposase_Sequence', 'SOLID_Small_RNA_Adapter')
    }

    fastqcSummary=list(Total_Sequences,Per_base_sequence_quality,Adapter_Content)
    names(fastqcSummary)=c('total_sequences','Per_base_sequence_quality','Adapter_Content')

    output_dir=files[i]
    o=paste0(output_dir,'/fastqcSummary.rda')
    save(fastqcSummary,file=o)

  } 

  #if (length(files)==1) { return(fastqcSummary) }

} #end of function


## finding DMRs

DMRs_v2 = function(input_GRset, DMRs_filename, B) {

# Enable parallelization
require(doParallel)
registerDoParallel(cores = 3)

# Find bumps
library(bumphunter)

pd=pData(input_GRset)
neg_control_PCs=pd[,12:19]
sex=factor(as.vector(pd$sex), levels=c('Female','Male')) # Female = 0; Male=1;

#Arguments for bumphunter
mod=model.matrix(~sex)
mod=cbind(mod, neg_control_PCs)
p=getBeta(input_GRset)

anno=getAnnotation(input_GRset)
chr=as.vector(anno$chr)
pos=as.vector(anno$pos)

bumps = bumphunterEngine(p, mod, chr = chr, 
pos = pos, cutoff= 0.1, nullMethod = "bootstrap",
smooth=TRUE, B=B)

dat=data.frame(chr=bumps$tab$chr,start=bumps$tab$start,end=bumps$tab$end,p.value=bumps$tab$p.value, FWER=bumps$tab$fwer)

save(bumps, dat, file=DMRs_filename)

num_dmrs=length(which(dat$FWER<=0.1))

return(num_dmrs)

}

## finding DNAm blocks

blocks_v2 = function(input_GRset, blocks_filename, B) {

library(minfi)

# Enable parallelization
require(doParallel)
registerDoParallel(cores = 3)

# Find blocks
pd=pData(input_GRset)
neg_control_PCs=pd[,12:19]
sex=factor(as.vector(pd$sex), levels=c('Female','Male')) # Female = 0; Male=1;

mod=model.matrix(~sex)
mod=cbind(mod, neg_control_PCs)

cobj=cpgCollapse(input_GRset, what="Beta")

blocks=blockFinder(cobj$object, design=mod, coef = 2, what = 'Beta', cluster=NULL, nullMethod='bootstrap',
cutoff = 0.1, pickCutoff = FALSE, smooth = TRUE, smoothFunction = locfitByCluster,
B = B, verbose = TRUE, bpSpan = 2.5 * 10^5)

dat=data.frame(chr=blocks$tab$chr,start=blocks$tab$start,end=blocks$tab$end,p.value=blocks$tab$p.value, FWER=blocks$tab$fwer)

save(blocks, dat, file=blocks_filename)

num_blocks=length(which(dat$FWER<=0.1))

return(num_blocks)

}

# function to plot DNAm age (Horvath) wrt actual (chronological) age

DNAmAge_plots = function(females_DNAmAge , males_DNAmAge, pdf_filename) {

pdf(pdf_filename)

## all ages

fit=lm(males_DNAmAge$DNAmAge~males_DNAmAge$RealAge)
eb=summary(fit)
t=eb$coef[2,3]
p=eb$coef[2,4]
intercept=eb$coef[1,1]
slope=eb$coef[2,1]

cc = cor.test(males_DNAmAge$RealAge, males_DNAmAge$DNAmAge, method='pearson', alternative='two.sided')
male_legend_label=paste0('Male: ','cor=',signif(cc$estimate,2),'; p=',signif(cc$p.value,2))

par(mar=c(5.1,5,4.1,2.1))

xlim=c(min(c(males_DNAmAge$RealAge,females_DNAmAge$RealAge)),max(c(males_DNAmAge$RealAge,females_DNAmAge$RealAge)))
ylim=c(min(c(males_DNAmAge$DNAmAge,females_DNAmAge$DNAmAge)),max(c(males_DNAmAge$DNAmAge,females_DNAmAge$DNAmAge)))

plot(males_DNAmAge$RealAge, males_DNAmAge$DNAmAge, xlab='Chronological Age',ylab='DNAm Age', xlim=xlim, ylim=ylim,
  pch=21, col='black', bg='blue', cex=2, cex.axis=2, cex.lab=2, cex.main=1)
abline(intercept,slope,lty='solid',col='blue',lwd=4)

fit=lm(females_DNAmAge$DNAmAge~females_DNAmAge$RealAge)
eb=summary(fit)
t=eb$coef[2,3]
p=eb$coef[2,4]
intercept=eb$coef[1,1]
slope=eb$coef[2,1]

cc = cor.test(females_DNAmAge$RealAge, females_DNAmAge$DNAmAge, method='pearson', alternative='two.sided')
female_legend_label=paste0('Female: ','cor=',signif(cc$estimate,2),'; p=',signif(cc$p.value,2))

points(females_DNAmAge$RealAge, females_DNAmAge$DNAmAge,pch=21, col='black', bg='red', cex=2)
abline(intercept,slope,lty='solid',col='red',lwd=4)

abline(0,1,lty=2,col='black',lwd=4)

legend('topleft', c(male_legend_label, female_legend_label), cex=.8)
legend('bottomright', as.expression(bquote(bold("All ages"))), cex=1.5, bty='n')
#legend("topleft", c('Males', 'Females','Identity'), lwd=4, col=c('blue','red','black'), lty=c(1,1,2), cex=1.3)

## peds only

females_DNAmAge_peds=females_DNAmAge[which(females_DNAmAge$RealAge<=21),]
males_DNAmAge_peds=males_DNAmAge[which(males_DNAmAge$RealAge<=21),]

fit=lm(males_DNAmAge_peds$DNAmAge~males_DNAmAge_peds$RealAge)
eb=summary(fit)
t=eb$coef[2,3]
p=eb$coef[2,4]
intercept=eb$coef[1,1]
slope=eb$coef[2,1]

cc = cor.test(males_DNAmAge_peds$RealAge, males_DNAmAge_peds$DNAmAge, method='pearson', alternative='two.sided')
male_legend_label=paste0('Male: ','cor=',signif(cc$estimate,2),'; p=',signif(cc$p.value,2))

par(mar=c(5.1,5,4.1,2.1))

xlim=c(min(c(males_DNAmAge_peds$RealAge,females_DNAmAge_peds$RealAge)),max(c(males_DNAmAge_peds$RealAge,females_DNAmAge_peds$RealAge)))
ylim=c(min(c(males_DNAmAge_peds$DNAmAge,females_DNAmAge_peds$DNAmAge)),max(c(males_DNAmAge_peds$DNAmAge,females_DNAmAge_peds$DNAmAge)))

plot(males_DNAmAge_peds$RealAge, males_DNAmAge_peds$DNAmAge, xlab='Chronological Age',ylab='DNAm Age', xlim=xlim, ylim=ylim,
  pch=21, col='black', bg='blue', cex=2, cex.axis=2, cex.lab=2, cex.main=1)
abline(intercept,slope,lty='solid',col='blue',lwd=4)

fit=lm(females_DNAmAge_peds$DNAmAge~females_DNAmAge_peds$RealAge)
eb=summary(fit)
t=eb$coef[2,3]
p=eb$coef[2,4]
intercept=eb$coef[1,1]
slope=eb$coef[2,1]

cc = cor.test(females_DNAmAge_peds$RealAge, females_DNAmAge_peds$DNAmAge, method='pearson', alternative='two.sided')
female_legend_label=paste0('Female: ','cor=',signif(cc$estimate,2),'; p=',signif(cc$p.value,2))

points(females_DNAmAge_peds$RealAge, females_DNAmAge_peds$DNAmAge,pch=21, col='black', bg='red', cex=2)
abline(intercept,slope,lty='solid',col='red',lwd=4)

abline(0,1,lty=2,col='black',lwd=4)

legend('topleft', c(male_legend_label, female_legend_label), cex=.8)
legend('bottomright', as.expression(bquote(bold("0-21 y.o."))), cex=1.5, bty='n')
#legend("topleft", c('Males', 'Females','Identity'), lwd=4, col=c('blue','red','black'), lty=c(1,1,2), cex=1.3)

## adults only

females_DNAmAge_adults=females_DNAmAge[which(females_DNAmAge$RealAge>21),]
males_DNAmAge_adults=males_DNAmAge[which(males_DNAmAge$RealAge>21),]

fit=lm(males_DNAmAge_adults$DNAmAge~males_DNAmAge_adults$RealAge)
eb=summary(fit)
t=eb$coef[2,3]
p=eb$coef[2,4]
intercept=eb$coef[1,1]
slope=eb$coef[2,1]

cc = cor.test(males_DNAmAge_adults$RealAge, males_DNAmAge_adults$DNAmAge, method='pearson', alternative='two.sided')
male_legend_label=paste0('Male: ','cor=',signif(cc$estimate,2),'; p=',signif(cc$p.value,2))

par(mar=c(5.1,5,4.1,2.1))

xlim=c(min(c(males_DNAmAge_adults$RealAge,females_DNAmAge_adults$RealAge)),max(c(males_DNAmAge_adults$RealAge,females_DNAmAge_adults$RealAge)))
ylim=c(min(c(males_DNAmAge_adults$DNAmAge,females_DNAmAge_adults$DNAmAge)),max(c(males_DNAmAge_adults$DNAmAge,females_DNAmAge_adults$DNAmAge)))

plot(males_DNAmAge_adults$RealAge, males_DNAmAge_adults$DNAmAge, xlab='Chronological Age',ylab='DNAm Age', xlim=xlim, ylim=ylim,
  pch=21, col='black', bg='blue', cex=2, cex.axis=2, cex.lab=2, cex.main=1)
abline(intercept,slope,lty='solid',col='blue',lwd=4)

fit=lm(females_DNAmAge_adults$DNAmAge~females_DNAmAge_adults$RealAge)
eb=summary(fit)
t=eb$coef[2,3]
p=eb$coef[2,4]
intercept=eb$coef[1,1]
slope=eb$coef[2,1]

cc = cor.test(females_DNAmAge_adults$RealAge, females_DNAmAge_adults$DNAmAge, method='pearson', alternative='two.sided')
female_legend_label=paste0('Female: ','cor=',signif(cc$estimate,2),'; p=',signif(cc$p.value,2))

points(females_DNAmAge_adults$RealAge, females_DNAmAge_adults$DNAmAge,pch=21, col='black', bg='red', cex=2)
abline(intercept,slope,lty='solid',col='red',lwd=4)

abline(0,1,lty=2,col='black',lwd=4)

legend('topleft', c(male_legend_label, female_legend_label), cex=.8)
legend('bottomright', as.expression(bquote(bold(">21 y.o."))), cex=1.5, bty='n')

#plot the legend

# plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
# legend("center", c('Males', 'Females','Identity'), lwd=4, col=c('blue','red','black'), lty=c(1,1,2), cex=1.3)

dev.off()

}



############################
## DRAFTS & ARCHIVE
############################

# ## cset: result of cpgCollapse()$cobj
# ## block450: results of blockFinder
# ## coi = covariate of interest
# ## N: the number of blocks to plot, default=10
# ## blockname = name of block track, default='coi'
# ## filename = where to save plots
# ## scale, in kb. default = 100
# blockPlot = function(cset, blocks450, coi, N=10,
#   blockname = "coi", filename=paste0(blockname,"_blocks.pdf"),scale=100,
#   showMethPanel = TRUE, showGenePanel=TRUE, showDiffPanel=TRUE, 
#   showCancerPanel = TRUE, bty= "o") {
#   panels = c(showMethPanel,showGenePanel, showDiffPanel, showCancerPanel)
  
#   require(GenomicRanges)

#   blocksTable=with(blocks450$table, GRanges(seqnames=chr,ranges=IRanges(start=start,end=end)))
#   colIds=match(c("chr","start","end"),names(blocks450$table))
#   mcols(blocksTable)=blocks450$table[-colIds]

#   plotRegion = blocksTable[1:N]

#   ## annotation based on ensembl
#   cat("Loading Annotation.\n")
#   load("/home/epi/ajaffe/GenomicStates/GenomicState.Hsapiens.ensembl.GRCh37.p12.rda")
#   gs = GenomicState.Hsapiens.ensembl.GRCh37.p12$fullGenome
#   oo = findOverlaps(blocksTable, gs)
#   anno = split(gs[subjectHits(oo)], queryHits(oo))

#   # cancer blocks
#   if(showCancerPanel) {
#     load("/home/epi/ajaffe/Lieber/Projects/450k/devPaper/cancer_blocks_hansen.rda")
#     genome(blocks)="hg19"
#     cancerBlocks = blocks
#   }

#   cat("Ploting.")
#   pdf(filename,height=5,width=10)
#   par(bty=bty)

#   for(i in seq(along=plotRegion)) {
#     cat(".")
#     r = plotRegion[i]
#     tmp=subsetByOverlaps(cset,r)
#     tmp450=sort(subsetByOverlaps(blocksTable,r))
#     if(showCancerPanel) tmpBsmooth=subsetByOverlaps(cancerBlocks,r)
                      
#     beta=getBeta(tmp)
#     x=start(tmp)

#     ii=cset %over% r
#     d=blocks450$coef[ii]
#     sd=blocks450$fitted[ii]
#     ## which rows
#     Index=split(seq_along(coi),coi)
#     mns=sapply(Index,function(ind) rowMeans(beta[,ind]))  
#     smns=apply(mns,2,function(y) limma::loessFit(y,x,span=.2)$fitted)

#     ## paneling, from hector
#     mypar(1,1, brewer.name = "Set1")
#     par(mar=par()$mar+c(0,3,0,0))
#     omar=par()$mar
#     cmar=omar
#     cmar[1]=.5
#     par(mar=cmar)
    
#     layout(cbind(1:sum(panels)),height=c(1+2/3,1, 0.75,0.75)[1:sum(panels)])
#     if(showMethPanel) {
#       matplot(x,beta,col=as.numeric(factor(coi)),type="p",pch=".",
#       xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,1))
#       matplot(x,smns,col=1:2,type="l",lwd=2.5,add=TRUE,lty=1)
#       legend("bottomright",col=seq(along=levels(factor(coi))),
#         lty=1,lwd=2,legend=levels(factor(coi)),cex=.8, bty=bty)
        
#       segments(min(x), .2, min(x)+scale*1000, .2)
#       text(min(x),.05,labels=sprintf("%dkb",scale),pos=4,offset=0)
#       axis(side=2,at=c(.2,.5,.8), cex.axis=1.8)
#       mtext(sprintf("%s:%d-%d", seqnames(r), start(r), end(r)), side=3)
#       mtext("Methylation",side=2, line = 2.5,cex=1.5)

#       legend("topright", paste0("fwer = ",
#         signif(blocks450$tab$fwer[i],3)),bty=bty)
#     }   
#     # annotation
#     if(showGenePanel) {
#       cmar=omar
#       if(!is.na(showDiffPanel)) {
#         cmar[3]=0.5
#       } else {
#         cmar[c(1,3)]=c(0.5,0.5)
#       }
#       par(mar=cmar)

#       plot(x,rep(0,length(x)), type="n",ylim=c(-1.5,1.5),yaxt="n",ylab="",
#          xlab="",cex.axis = 1.5, cex.lab =1.5,xaxt="n")
#       a = as.data.frame(anno[[i]])
#       Strand = ifelse(a$strand == "+", 1, ifelse(a$strand=="-", -1, 0))
#       Col = ifelse(a$theRegion == "exon", "blue", ifelse(a$theRegion == "intron", "lightblue","white"))
#       Lwd = ifelse(a$theRegion == "exon", 1, ifelse(a$theRegion == "intron",0.5,0))
#       axis(2,c(-1,1),c("-","+"),tick=FALSE,las=1,cex.axis = 3)
#       abline(h=0,lty=3)
#       for(k in 1:nrow(a)) {
#         polygon(c(a$start[k],a$end[k],a$end[k],a$start[k]),
#           Strand[k]/2+c(-0.3,-0.3,0.3,0.3)*Lwd[k],col=Col[k])
#       }

#       ## by gene
#       g = split(a, sapply(a$symbol,"[", 1))
#       # g = split(a, a$Symbol)
#       s2 = ifelse(sapply(g, function(x) unique(x$strand))=="+",1,-1)
#       g = sapply(g, function(x) (max(x$end) - min(x$start))/2 + min(x$start) )
      
#       if(length(g) > 0) text(g, y=s2, names(g),font=1,pos=s2+2,cex=0.8)
          
#       mtext("Genes",side=2, line = 2.5,cex=1.5)
#       if(!showDiffPanel) {
#         xtick=pretty(x)
#         axis(side=1,at=xtick,labels=sprintf("%.1fMb",xtick / 1e6))
#       }
#     }

#     ## mean diff
#     if(showDiffPanel) {
#       cmar=omar
#       cmar[c(1,3)]=c(.5,.5)
#       par(mar=cmar)

#       zz=granges(tmp)

#       matplot(x,sd,xaxt="n",ylab="",xlab="",type="n",lty=1,ylim=c(-.6,.6),yaxt="n",pch=21)
#       axis(side=2,at=c(-.3,0,.3),labels=c("-.3","0",".3"))
#       xtick=pretty(x)
#       axis(side=1,at=xtick,labels=sprintf("%.1fMb",xtick / 1e6))
#       mtext("Diff",side=2, line = 2.5,cex=1.5)

#       ii=which(zz$type=="OpenSea")
#       blockgroup=zz$blockgroup[ii]

#       blockIndexes=split(seq(along=blockgroup),blockgroup)
#       for (ind in blockIndexes) {
#         ind=ii[ind]
#         lines(x[ind], sd[ind], lwd=2.5,col="black")
#       }


#       points(x[ii],d[ii],pch=21,cex=1.4,bg="black")
#       axis(side=2,at=c(-2,0,2))
#       abline(h=0,lty=2,col="black")

#       cmar=omar
#       cmar[3]=.5
#       par(mar=cmar)
#       matplot(x,beta,type="n",xaxt="n",yaxt="n",xlab="",
#         ylab="",ylim=c(0,2),bty="n")
#     }
    
#     #  browser()
#     if(showCancerPanel) {
#       col=ifelse(tmp450$value<0 & tmp450$p.value<.05,"blue",ifelse(tmp450$value>0 & tmp450$p.value<.05,"red","black"))
#       rect(start(tmp450),1+1/3,end(tmp450),1+2/3,col=col)
#       if(length(tmpBsmooth) > 0)  rect(start(tmpBsmooth),1/3,end(tmpBsmooth),2/3,col=ifelse(tmpBsmooth$direction=="hypo","blue","red"))
#       axis(side=2,at=c(.5,1.5),labels=c("Hansen et al.",blockname),las=1,lwd=0)
#       legend("bottomleft",pt.bg=c("blue","red"),legend=c("hypo","hyper"),pch=22,cex=.8)
#     }
#   }
#   dev.off()
# }
