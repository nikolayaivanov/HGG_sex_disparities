
## cset: result of minfi::cpgCollapse()$object
## blocks450: results of blockFinder
## coi = covariate of interest
## N: the number of blocks to plot, default=10
## blockname = name of block track, default='coi'
## filename = where to save plots
## scale, in kb. default = 100

blockPlot = function(cset,
  blocks450,
  coi,
  N=10,
  blockname = "coi",
  filename,
  scale=100,
  showMethPanel = TRUE,
  showGenePanel=TRUE,
  showDiffPanel=TRUE, 
  showCancerPanel = FALSE,
  bty= "o") {

  panels = c(showMethPanel,showGenePanel, showDiffPanel, showCancerPanel)
  
  #require(GenomicRanges)

  blocksTable=with(blocks450$table, GRanges(seqnames=chr,ranges=IRanges(start=start,end=end)))
  colIds=match(c("chr","start","end"),names(blocks450$table))
  mcols(blocksTable)=blocks450$table[-colIds]

  plotRegion = blocksTable[1:N]

  ## annotation based on ensembl
  cat("Loading Annotation.\n")
  load("/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/rdas/genomicState.rda")
  gs = genomicState$fullGenome
  oo = findOverlaps(blocksTable, gs)
  anno = split(gs[subjectHits(oo)], queryHits(oo))

  # cancer blocks
  if(showCancerPanel) {
    # read in cancer DNAm bocks (Hansen et al, 2011, PMID: 21706001)
    #load("/home/epi/ajaffe/Lieber/Projects/450k/devPaper/cancer_blocks_hansen.rda")
    blocks=read.csv('/athena/masonlab/scratch/users/nai2008/DNAm/second_round_of_scripts/cancer_blocks_hansen.csv')
    blocks=GRanges(blocks)
    genome(blocks)="hg19"
    cancerBlocks = blocks
  }

  cat("Ploting.")
  pdf(filename,height=5,width=10)
  par(bty=bty)

  for(i in seq(along=plotRegion)) { #//
    cat(".")
    r = plotRegion[i]
    tmp=subsetByOverlaps(cset,r)
    tmp450=sort(subsetByOverlaps(blocksTable,r))
    if(showCancerPanel) tmpBsmooth=subsetByOverlaps(cancerBlocks,r)
                      
    beta=getBeta(tmp)
    x=start(tmp)

    ii=cset %over% r
    d=blocks450$coef[ii]
    sd=blocks450$fitted[ii]
    ## which rows
    Index=split(seq_along(coi),coi)
    mns=sapply(Index,function(ind) rowMeans(beta[,ind]))  
    smns=apply(mns,2,function(y) limma::loessFit(y,x,span=.2)$fitted)

    ## paneling, from hector
    mypar(1,1, brewer.name = "Set1")
    par(mar=par()$mar+c(0,3,0,0))
    omar=par()$mar
    cmar=omar
    cmar[1]=.5
    par(mar=cmar)
    
    layout(cbind(1:sum(panels)),height=c(1+2/3,1, 0.75,0.75)[1:sum(panels)])
    if(showMethPanel) {
      matplot(x,beta,col=as.numeric(factor(coi)),type="p",pch=".",
      xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,1))
      matplot(x,smns,col=1:2,type="l",lwd=2.5,add=TRUE,lty=1)
      legend("bottomright",col=seq(along=levels(factor(coi))),
        lty=1,lwd=2,legend=levels(factor(coi)),cex=.8, bty=bty)
        
      segments(min(x), .2, min(x)+scale*1000, .2)
      text(min(x),.05,labels=sprintf("%dkb",scale),pos=4,offset=0)
      axis(side=2,at=c(.2,.5,.8), cex.axis=1.8)
      mtext(sprintf("%s:%d-%d", seqnames(r), start(r), end(r)), side=3)
      mtext("Methylation",side=2, line = 2.5,cex=1.5)

      legend("topright", paste0("fwer = ",
        signif(blocks450$tab$fwer[i],3)),bty=bty)
    }   
    # annotation
    if(showGenePanel) {
      cmar=omar
      if(!is.na(showDiffPanel)) {
        cmar[3]=0.5
      } else {
        cmar[c(1,3)]=c(0.5,0.5)
      }
      par(mar=cmar)

      plot(x,rep(0,length(x)), type="n",ylim=c(-1.5,1.5),yaxt="n",ylab="",
         xlab="",cex.axis = 1.5, cex.lab =1.5,xaxt="n")
      a = as.data.frame(anno[[i]])
      Strand = ifelse(a$strand == "+", 1, ifelse(a$strand=="-", -1, 0))
      Col = ifelse(a$theRegion == "exon", "blue", ifelse(a$theRegion == "intron", "lightblue","white"))
      Lwd = ifelse(a$theRegion == "exon", 1, ifelse(a$theRegion == "intron",0.5,0))
      axis(2,c(-1,1),c("-","+"),tick=FALSE,las=1,cex.axis = 3)
      abline(h=0,lty=3)
      for(k in 1:nrow(a)) {
        polygon(c(a$start[k],a$end[k],a$end[k],a$start[k]),
          Strand[k]/2+c(-0.3,-0.3,0.3,0.3)*Lwd[k],col=Col[k])
      }

      ## by gene
      g = split(a, sapply(a$symbol,"[", 1))
      # g = split(a, a$Symbol)
      s2 = ifelse(sapply(g, function(x) unique(x$strand))=="+",1,-1)
      g = sapply(g, function(x) (max(x$end) - min(x$start))/2 + min(x$start) )
      
      if(length(g) > 0) text(g, y=s2, names(g),font=1,pos=s2+2,cex=0.8)
          
      mtext("Genes",side=2, line = 2.5,cex=1.5)
      if(!showDiffPanel) {
        xtick=pretty(x)
        axis(side=1,at=xtick,labels=sprintf("%.1fMb",xtick / 1e6))
      }
    }

    ## mean diff
    if(showDiffPanel) {
      cmar=omar
      cmar[c(1,3)]=c(.5,.5)
      par(mar=cmar)

      zz=granges(tmp)

      matplot(x,sd,xaxt="n",ylab="",xlab="",type="n",lty=1,ylim=c(-.6,.6),yaxt="n",pch=21)
      axis(side=2,at=c(-.3,0,.3),labels=c("-.3","0",".3"))
      xtick=pretty(x)
      axis(side=1,at=xtick,labels=sprintf("%.1fMb",xtick / 1e6))
      mtext("Diff",side=2, line = 2.5,cex=1.5)

      ii=which(zz$type=="OpenSea")
      blockgroup=zz$blockgroup[ii]

      blockIndexes=split(seq(along=blockgroup),blockgroup)
      for (ind in blockIndexes) {
        ind=ii[ind]
        lines(x[ind], sd[ind], lwd=2.5,col="black")
      }


      points(x[ii],d[ii],pch=21,cex=1.4,bg="black")
      axis(side=2,at=c(-2,0,2))
      abline(h=0,lty=2,col="black")

      cmar=omar
      cmar[3]=.5
      par(mar=cmar)
      matplot(x,beta,type="n",xaxt="n",yaxt="n",xlab="",
        ylab="",ylim=c(0,2),bty="n")
    }
    
    #  browser()
    if(showCancerPanel) {
      col=ifelse(tmp450$value<0 & tmp450$p.value<.05,"blue",ifelse(tmp450$value>0 & tmp450$p.value<.05,"red","black"))
      rect(start(tmp450),1+1/3,end(tmp450),1+2/3,col=col)
      if(length(tmpBsmooth) > 0)  rect(start(tmpBsmooth),1/3,end(tmpBsmooth),2/3,col=ifelse(tmpBsmooth$Direction.of.Methylation.Change=="hypo","blue","red"))
      axis(side=2,at=c(.5,1.5),labels=c("Hansen et al.",blockname),las=1,lwd=0)
      legend("bottomleft",pt.bg=c("blue","red"),legend=c("hypo","hyper"),pch=22,cex=.8)
    }
  } #//
  
  dev.off()

}
