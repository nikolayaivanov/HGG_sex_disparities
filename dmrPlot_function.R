#
# "http://rafalab.jhsph.edu/CGI/model-based-cpg-islands-hg18.txt"
# "http://rafalab.jhsph.edu/CGI/model-based-cpg-islands-hg19.txt"

# gs stands for genomicState

dmrPlot = function(regions, p, chr, pos, cluster, genes, coi, build="hg19", number=100, gs, species = "human",
  Jitter = FALSE, cols=NULL, lines = FALSE, linesSmooth = TRUE, title = TRUE, Legend = TRUE, colorRamp = FALSE, 
  meanSmooth=TRUE, plotCpG = TRUE, geneAnno = "gene") {
  
  #require(bumphunter)
  #require(derfinder)
  gr = GRanges(regions$chr, IRanges(regions$start, regions$end))
  
  if(build == "hg18") {
    cpg.cur = read.table("http://rafalab.jhsph.edu/CGI/model-based-cpg-islands-hg18.txt",
      header = TRUE, as.is=TRUE)
    library("BSgenome.Hsapiens.UCSC.hg18")
  }
  
  if(build == "hg19") {
    cpg.cur = read.table("http://web.stanford.edu/class/bios221/data/model-based-cpg-islands-hg19.txt",
      header = TRUE, as.is=TRUE)
    library("BSgenome.Hsapiens.UCSC.hg19")
  }
  
  if(build == "mm9") {
    cpg.cur = read.table("http://rafalab.jhsph.edu/CGI/model-based-cpg-islands-mm9.txt",
      header = TRUE, as.is=TRUE)
    library("BSgenome.Mmusculus.UCSC.mm9")
  }
    
  ocpgi=data.frame(chr=I(cpg.cur[,1]), 
    start=as.numeric(cpg.cur[,2]), 
    end=as.numeric(cpg.cur[,3]))
  ADD1 = 1500; PAD = 10
  
  gr2 = GRanges(regions$chr, IRanges(regions$start - ADD1, regions$end + ADD1))
  anno = annotateRegions(gr2, gs)$annotationList
  
  # check regions
  if(is.numeric(coi)) {
    groups=coi
    gNames= sort(unique(coi))
  }

  if(is.character(coi) | is.factor(coi)) {
    groups = factor(coi)
    gNames= levels(groups)
  }
  gIndexes=split(1:length(groups),groups)
  
  brewer.n=max(3,min(length(unique(coi)),11))

  
  if(is.null(cols)) {
    mypar(brewer.n=brewer.n)
  } else if(length(cols) == 1 & cols[1] %in% rownames(brewer.pal.info)) {
    mypar(brewer.name=cols, brewer.n=brewer.n)
  } else {
    mypar()
    palette(cols)
  } 
    
  if(colorRamp) {
    pal = colorRampPalette(palette()[3:brewer.n])
    palette(pal(brewer.n))
  }
  
  cat("Plotting.\n")
  for(j in 1:(min(nrow(regions), number))) {
  
    layout(matrix(1:2,ncol=1),heights=c(0.7,0.3))

    # first plot, region
    par(mar=c(0,4.5,0.25,1.1),oma=c(0,0,2,0))
    
    Index=(regions[j,7]-PAD):(regions[j,8]+PAD)
    Index = Index[Index %in% seq(along=cluster)]
    Index=Index[cluster[Index]==regions[j,6]]
    
    # make DMR plots
    if(Jitter) {
      posx = matrix(pos[Index], nc = ncol(p), nr = length(Index),
        byrow=  FALSE)
      posx = t(apply(posx, 1, function(x) jitter(x,amount = 12)))
      
    } else posx = pos[Index]
          
    matplot(posx, p[Index,], ylim = c(0,1),
      ylab = "", xlab = "",xaxt = "n",cex=0.7,
      cex.axis = 1.7, cex.lab = 1.7, pch=21,
      bg = as.numeric(factor(groups)),col="black",
      xlim = range(pos[Index]), yaxt="n")
    axis(2, at = c(0.2, 0.5, 0.8), cex.axis=1.7)

    xx=pos[Index]
    for(k in seq(along=gIndexes)){
      if(length(gIndexes[[k]]) == 1) yy=p[Index,gIndexes[[k]]]
      if(length(gIndexes[[k]]) > 1)   yy=rowMeans(p[Index,gIndexes[[k]]])
      if(meanSmooth) { 
        fit1=loess(yy~xx,degree=1,span=300/(36*length(xx)),
          family="symmetric")
        lines(xx,fit1$fitted,col=k,lwd=2)
      } else  lines(xx,yy,col=k,lwd=2)

    }
    
    mtext("Methylation",side=2, line = 2.5,cex=1.8)

    if(Legend) {
      if(length(unique(coi)) < 4) {
        legend("topleft",legend=gNames,col=1:length(gNames),
        lty=1, lwd = 4,cex=1,bty="n")
      } else {
        legend("topleft",legend=gNames,col=1:length(gNames),
          pch=19, pt.cex = 2,cex=1.1, nc = 6,bty="n")
      } 
    }
        
    abline(v=(regions$start[j]-15),lty=1)
    abline(v=(regions$end[j]+15),lty=1)

    if(title) mtext(paste0(genes$name[j],"; ", 
      genes$distance[j],"bp from tss:",genes$description[j]), outer=T,cex=1.3)
      
    # add cpgs
    if(plotCpG) {
      thechr=as.character(regions$chr[j])
      chrName = strsplit(thechr, "r")[[1]][2]
      chrName = paste("Chromosome",chrName)
      
      start = pos[Index[1]]
      end = pos[Index[length(Index)]]
      ocpgi2=ocpgi[ocpgi$chr%in%unique(as.character(thechr)),]
      
      ##PLOT CPG ISLANDS
      if(species=="human") seq<-Hsapiens[[as.character(thechr) ]]
      if(species=="mouse") seq<-Mmusculus[[as.character(thechr) ]]
      
      subseq<-subseq(seq,start=start,end=end)
      cpgs=start(matchPattern("CG",subseq))+start-1

      if(plotCpG) rug(cpgs,col="black")  #  previously had 'Rug'

      Index1 = which(ocpgi2[,1]==as.character(thechr) &
           ((ocpgi2[,2] > start & ocpgi2[,2]< end) |
            (ocpgi2[,3] > start & ocpgi2[,3]< end)))
      if(length(Index1)>0) sapply(Index1,function(j) rug(unlist(ocpgi2[j,2:3]),  #  previously had 'Rug'
             col="darkgreen",lwd=3,side=1))
    }
    
    # plot 3
    ##PLOT GENES
    par(mar=c(3.5,4.5,0.25,1.1))

    plot(0,0, type="n", xlim=range(xx),ylim=c(-1.5,1.5),yaxt="n",ylab="",
       xlab="",cex.axis = 1.5, cex.lab =1.5)
    a = as.data.frame(anno[[j]])
    Strand = ifelse(a$strand == "+", 1, ifelse(a$strand=="-", -1, 0))
    Col = ifelse(a$theRegion == "exon", "blue", ifelse(a$theRegion == "intron", "lightblue","white"))
    Lwd = ifelse(a$theRegion == "exon", 1, ifelse(a$theRegion == "intron",0.5,0))
    axis(2,c(-1,1),c("-","+"),tick=FALSE,las=1,cex.axis = 3)
    abline(h=0,lty=3)
    for(k in 1:nrow(a)) {
      polygon(c(a$start[k],a$end[k],a$end[k],a$start[k]),
        Strand[k]/2+c(-0.3,-0.3,0.3,0.3)*Lwd[k],col=Col[k])
    }
    
    if(sum(a$theRegion=="exon") > 0) {
      e = a[a$theRegion=="exon",]
      s2 = Strand[a$theRegion=="exon"]
      g = unlist(e$symbol)
      g[is.na(g)] = ""
      if(length(g) > 0) text(x=e$start + e$width/2,
        y=s2*0.75, g,font=1,pos=s2+2,cex=1.2)
    }
        
    mtext("Genes",side=2, line = 2.5,cex=1.5)
    mtext(chrName,side=1, line = 2,cex=1.4)

    abline(v=(pos[regions[j,7]]-15),lty=1)
    abline(v=(pos[regions[j,8]]+15),lty=1)
    
  }

}


