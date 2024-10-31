# functions ###############
readNamedMM = function(f){
  require(Matrix)
  if(file.exists(paste0(f,'.mtx'))){
    m = Matrix::readMM(paste0(f,'.mtx'))
    rownames(m) = readLines(paste0(f,'_key1.csv'))
    colnames(m) = readLines(paste0(f,'_key2.csv'))
  }else if(file.exists(paste0(f,'.mtx.gz'))){
    m = Matrix::readMM(paste0(f,'.mtx.gz'))
    rownames(m) = readLines(paste0(f,'_key1.csv.gz'))
    colnames(m) = readLines(paste0(f,'_key2.csv.gz'))
  }else{
    m = NULL
  }
  m
}

loadSC_AS = function(segs,path){
  e = readNamedMM(paste0(path,'.e'))
  i = readNamedMM(paste0(path,'.i'))
  list(e=e[rownames(segs),],
       i=i[rownames(segs),colnames(e)])
}

get_dPSI = function(d,f,mincov=50){
  ii = calcColSums(d$i , f)
  ee = calcColSums(d$e , f)
  ee = ee[rownames(d$i),]
  ii = ii[rownames(d$i),colnames(ee)]
  tt = ii + ee
  psi = ii/tt
  psi[tt<mincov] = NA
  do.call(rbind,apply(psi,1,function(x){
    x = sort(x[!is.na(x)])
    r = data.frame(low_state=NA,high_state=NA,dpsi=NA)
    if(length(x)>1){
      x = x[c(1,length(x))]
      r$low_state = names(x)[1]
      r$high_state = names(x)[2]
      r$dpsi = x[2]-x[1]
    }
    r
  }))
}

findNearestConstantExons = function(seg,sid,psi.thr=0.95){
  gsegs = seg[seg$gene_id %in% seg[sid,'gene_id'],]
  strand = gsegs$strand[1]
  if(strand==0)
    strand = 1
  gsegs = gsegs[order(strand*gsegs$start),]
  up = down = NA
  ucnst = which((is.na(gsegs$mean_psi) | gsegs$mean_psi>=psi.thr) & gsegs$sites %in% c('ad','sd'))
  dcnst = which((is.na(gsegs$mean_psi) | gsegs$mean_psi>=psi.thr) & gsegs$sites %in% c('ad','ae'))
  inx = which(rownames(gsegs) == sid)
  if(length(ucnst)>0 && any(ucnst<inx))
    up = rownames(gsegs)[max(ucnst[ucnst<inx])]
  if(length(dcnst)>0 && any(dcnst>inx))
    down = rownames(gsegs)[min(dcnst[dcnst>inx])]
  if(strand == -1){
    res = c(up=down,down=up)
  }else
    res = c(up=up,down=down)
  res
}

castXYtable = function(x,y,i){
  ys = sort(unique(y))
  xs = sort(unique(x))
  m = matrix(NA,nrow=length(xs),ncol=length(ys),dimnames=list(xs,ys))
  x = as.character(x)
  y = as.character(y)
  for(j in 1:length(x))
    m[x[j],y[j]]= i[j]
  m
}

plotSegmentCoverage = function(sid,usid,dsid,celltypes,
                               seg,barcodes,
                               scanBamFlags=list(isNotPassingQualityControls=FALSE,isDuplicate=FALSE,isSupplementaryAlignment=FALSE,isSecondaryAlignment=FALSE),
                               plot.junc.only.within=NA,
                               min.junc.cov.f = 0.01,
                               min.junc.cov = 3){
  bams = unique(samples[,c('sample_id','bam_path')])
  covs = list()
  start = seg[usid,'start']
  stop = seg[dsid,'stop']
  chr = seg[sid,'chr_id']
  strand = NA
  
  for(ct in celltypes){
    cov = list()
    for(i in seq_len(nrow(bams))){
      tagFilter = list()
      tagFilter$CB = barcodes$barcode[barcodes$sample_id == bams$sample_id[i] & !is.na(barcodes$celltype) & barcodes$celltype==ct] 
      
      if(length(tagFilter$CB)==0) next
      cov[[length(cov)+1]] = getReadCoverage(bams$bam_path[i],
                                             chr,start,stop,strand=strand,scanBamFlags=scanBamFlags,tagFilter = tagFilter)
    }
    covs[[ct]] = sumCovs(cov)
  }
  
  psi = split(pbas$ir[sid,] , pbmeta$celltype)
  psi = lapply(psi,na.omit)
  psi = psi[sapply(psi,length)>0]
  psi = psi[order(sapply(psi,mean))]
  
  l = rbind(c(1,2),
            c(1,3),
            c(1,4))
  layout(l,widths = c(1,3))
  par(bty='n',tcl=-0.2,mgp=c(1.3,0.3,0),mar=c(3,13,1.5,0),oma=c(0,0,3,1))
  boxplot(psi,horizontal=TRUE,las=1,xlab='PSI')
  for(ct in celltypes){
    plotReadCov(covs[[ct]],xlim=c(start,stop),ylab='Coverage',xlab=chr,main=ct,
                plot.junc.only.within=plot.junc.only.within,min.junc.cov = min.junc.cov,
                min.junc.cov.f=min.junc.cov.f)
    #abline(v=c(seg[sid,'start'],seg[sid,'stop']),col='red',lty=2,xpd=FALSE)
  }
  ann = gtf[gtf$gene_id==seg[sid,'gene_id'] & gtf$start < stop & gtf$stop >= start,]
  ann$exon.col = 'black'
  ann$cds.col = 'black'
  f = ann$start >= seg[sid,'start'] & ann$stop <= seg[sid,'stop']
  ann$exon.col[f]=ann$cds.col[f] = 'red'
  
  plotTranscripts(ann,new = T,exon.col = NA,cds.col = NA,xlim=c(start,stop))
  #abline(v=c(seg[sid,'start'],seg[sid,'stop']),col='red',lty=2,xpd=FALSE)
  mtext(paste0(sid,' ', gene.descr[seg[sid,'gene_id'],'name'],'\n',gene.descr[seg[sid,'gene_id'],'descr']),side = 3,outer = TRUE)
}

log_info = function(text,...){
  text = paste0(text,...)
  print(paste0("INFO [",Sys.time(),"]: ",text))
}