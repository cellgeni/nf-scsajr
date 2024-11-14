bigWig2Cov = function(bw){
  bw = as.data.frame(bw)
  bw = bw[order(bw$start),]
  start = bw$start[1]
  stop = bw$end[nrow(bw)]
  r = rep(0,stop-start+1)
  for(i in 1:nrow(bw)){
    r[(bw$start[i]:bw$end[i])-start+1] = bw$score[i]
  }
  list(cov=r,x=start:stop)
}

getRecountCov = function(sid,jxn,gene.grange){
  require(recount3)
  require(rtracklayer)

  bigwig = jxn$BigWigURL[colnames(jxn)==sid]

  bw = rtracklayer::import.bw(bigwig,which=gene.grange)
  r = bigWig2Cov(bw)
  j = jxn[,sid]
  jxns = findOverlaps(jxn@rowRanges,gene.grange)
  j = jxn[jxns@from,]
  #introns can be duplicated...
  nms=unique(j@rowRanges@ranges@NAMES)
  j = j[nms,]

  r$juncs =  cbind(as.data.frame(j@rowRanges)[,c('start','end','strand')],score=j@assays@data$counts[,sid])
  r
}

sumCovs = function(l){
  r = l[[1]]
  if(length(l)==1)
    return(r)
  juncs = unique(do.call(rbind,unname(lapply(l,function(c)c$juncs[,1:3]))))
  juncs$score = 0
  juncs[rownames(r$juncs),'score'] = r$juncs$score

  for(i in 2:length(l)){
    if(!all(r$x==l[[i]]$x))
      stop("objects should cover identicall intervals")
    r$cov = r$cov + l[[i]]$cov
    juncs[rownames(l[[i]]$juncs),'score'] = juncs[rownames(l[[i]]$juncs),'score'] + l[[i]]$juncs$score
  }
  r$juncs = juncs
  r
}

loadEnsGTF = function(f,features=NULL){
  r = read.table(f,sep='\t')
  if(!is.null(features ))
    r = r[r$V3 %in% features,]
  a = lapply(strsplit(r$V9,';\\s?',perl=T),function(x){x=strsplit(x,'[ =]');setNames(sapply(x,'[',2),sapply(x,'[',1))})
  names = unique(unlist(lapply(a,names)))
  a = do.call(rbind,lapply(a,'[',names))
  colnames(a) = names
  r = r[,c(1:5,7)]
  colnames(r) = c('chr_id','type','feature','start','stop','strand')
  cbind(r,a)
}


plotTranscripts = function(a,
                           ylim=c(0,length(unique(a$transcript_id))),
                           xlim=c(ifelse(a$strand[1]=='+',min(a$start),max(a$stop)),ifelse(a$strand[1]=='+',max(a$stop),min(a$start))),
                           xlab=a$chr_id[1],
                           new=TRUE,yspace=0.8,exon.col='black',cds.col='black',
                           text.cex = 0.7,
                           ...){
  
  if(!is.na(exon.col))
    a$exon.col = exon.col
  if(!is.na(cds.col))
    a$cds.col = cds.col
  
  transc = split(a,a$transcript_id)
  transc = transc[order(sapply(transc,function(x){max(x$stop)-min(x$start)}))]
  if(new)
    plot(1,t='n',xlim=xlim,ylim=ylim,yaxt='n',ylab='',xlab=xlab,...)
  ystep = (ylim[2]-ylim[1])/length(transc)
  
  for(i in 1:length(transc)){
    y = ylim[1] + ystep*i - ystep/2
    t = transc[[i]]
    lines(c(min(t$start),max(t$stop)),c(y,y))
    f = t$feature == 'exon'
    if(sum(f)>0)
      rect(t$start[f],y-ystep/2*yspace,t$stop[f],y+ystep/2*yspace,col = 'white',border = t$exon.col[f])
    f = t$feature == 'CDS'
    if(sum(f)>0)
      rect(t$start[f],y-ystep/2*yspace,t$stop[f],y+ystep/2*yspace,col = t$cds.col[f],border = t$cds.col[f])
  }
  text(par('usr')[1],seq(ylim[1]+ystep/2,by = ystep,length.out = length(transc)),sapply(transc,function(x)x$transcript_name[1]),
       adj = c(1,0.5),xpd=T,cex=text.cex)
}

#' Extract read coverage from bam files
#'
#' @param bams character vector with paths to bam files
#' @param chr contig name
#' @param start,end coordinates of region
#' @param strand strand, NA for unstranded (default)
#' @param scanBamFlags list of flags to filter reads, see GenomicAlignments::scanBamFlag
#' @param tagFilter list of tags value to filter by. Can be used to filter cells in single cell bams, for example tagFilter=list(CB=barcodes)
#'
#' @return list with three items: x (chr coordinates); cov - number of reads mapped to chr position; juncs - data.frame with introns
#' @export
getReadCoverage = function(bams,chr,start,end,strand=NA,scanBamFlags=list(),tagFilter=list()){
  require(GenomicAlignments)
  if(start>end){
    t = start
    start = end
    end=t
  }
  r = list(x = start:end,
           cov = NULL,
           juncs=NULL,
           chr=chr,
           start=start,
           end=end)
  scanBamFlags$isMinusStrand = strand==-1
  flags = do.call(scanBamFlag,scanBamFlags)
  param = ScanBamParam(flag=flags,which=GRanges(chr, IRanges(start, end)),tagFilter = tagFilter)
  i = 1
  for(b in bams){
    #cat('\r',i,'     ')
    i = i + 1
    bam = readGAlignments(b,param = param)
    cov=coverage(bam)[[chr]][start:end]
    juncs = as.data.frame(summarizeJunctions(bam))
    rownames(juncs)=paste(juncs$seqnames,juncs$start,juncs$end,sep='-')
    if(is.null(r$cov)){
      r$cov=cov
      r$juncs=juncs
    }else{
      r$cov = r$cov + cov
      cmn = intersect(rownames(juncs),rownames(r$juncs))
      r$juncs[cmn,'score'] = r$juncs[cmn,'score'] + juncs[cmn,'score']
      r$juncs = rbind(r$juncs,juncs[setdiff(rownames(juncs),rownames(r$juncs)),])
    }
  }
  invisible(r)
}

#' Plots read coverage
#'
#' @param r read coverage; output of \code{\link{getReadCoverage}}
#' @param min.junc.cov numeric, plots only junctions (introns) with coverage not less than \code{min.junc.cov}
#' @param min.junc.cov.f numeric, plots only junctions (introns) with coverage not less than \code{min.junc.cov.f} of maximal coverage in the region
#' @param plot.junc.only.within logical, plot only junction with both ends within the region, FALSE plots all junctions with at least one end within region. NA plot all junctions overlapping the region.
#' @param ylim,xlim see \code{\link{plot}}
#' @param reverse reverse x coordinates
#' @param junc.col colour for junction line. Individual color could be specified for each junction
#' @param junc.lwd line width for jucntion line
#' @param ... other parameters for plot function
#'
#' @export
plotReadCov = function(r,min.junc.cov=0,min.junc.cov.f=0,plot.junc.only.within=FALSE,ylim=NULL,
                       xlim=range(r$x),reverse=FALSE,junc.col='blue',junc.lwd=3,bottom.mar=0,...){
  f = r$x >= xlim[1] & r$x <=xlim[2]
  r$x = r$x[f]
  r$cov = r$cov[f]
  if(nrow(r$juncs)>0)
    r$juncs$col = junc.col
  r$juncs = r$juncs[r$juncs$start <= xlim[2] & r$juncs$end >=xlim[1] & r$juncs$score >= min.junc.cov,]
  if(!is.na(plot.junc.only.within)){
    if(plot.junc.only.within){
      r$juncs = r$juncs[r$juncs$start > xlim[1] & r$juncs$end < xlim[2],]
    }else{
      r$juncs = r$juncs[(r$juncs$start > xlim[1] & r$juncs$start < xlim[2]) | (r$juncs$end > xlim[1] & r$juncs$end < xlim[2]),]
    }
  }

  start = r$x[1]
  end = r$x[length(r$x)]
  r$cov[c(1,length(r$cov))] = 0
  if(is.null(ylim)){
    ylim = c(0,max(r$cov,ifelse(nrow(r$juncs)>0,max(r$juncs$score),1)))
    ylim = c(-bottom.mar*ylim[2],ylim[2])
  }
  r$juncs = r$juncs[r$juncs$score >= min.junc.cov.f * ylim[2],]
  if(reverse)
    xlim=rev(xlim)
  plot(r$x,r$cov,t='n',ylim=ylim,xlim=xlim,yaxt='n',...)
  axis(2,at=c(0,ylim[2]),labels = c('',ylim[2]))
  polygon(r$x,r$cov,col = 'gray',border=NA)
  if(nrow(r$juncs)>0)
    for(i in 1:nrow(r$juncs)){
      start = r$juncs$start[i]
      stop = r$juncs$end[i]
      # to make junction max height always within region
      if(start<xlim[1] & stop >= xlim[2]){
        start = xlim[1] - mean(xlim)
        stop = xlim[2] + mean(xlim)
      }else if (start < xlim[1]){
        start = max(start,xlim[1] - (stop - xlim[1]))
      }else if (stop > xlim[2]){
        stop = min(stop,xlim[2] + (xlim[2]-start))
      }
      plotArc(start,stop,r$juncs$score[i],col=r$juncs$col[i],lwd=junc.lwd)
    }
  invisible(ylim)
}

#' Plots parabolic arc
#'
#' @param from,to x coordinates of arc
#' @param top highest point of arc
#' @param n number of points
#' @param y.base bottom coordinate of arc
#' @param ... other parameters of lines functoin
plotArc = function(from,to,top,n=100,y.base=0,...){
  len = to - from
  x = seq(from=0,to=len,length.out = n)
  y = x*4*top/len - x^2*(4*top/len^2)
  lines(x+from,y+y.base,...)
}
