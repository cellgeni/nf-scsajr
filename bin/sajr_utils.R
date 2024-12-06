# Conventions:
# 1. parameter to group observations is always named "groupby". If its length is identical to number of observations 
#  than it is assumed to be the factor, otherwise it is a list of column names. See getGroupbyFactor
# 2. 

# TODOs
# 1. plus it depends on ncells that is not propogated by pseudobulk
# 3. change all to SE

library(SummarizedExperiment)
# scsajr ################
DELIMETER = '$'

filterSegmentsAndSamples = function(data,
                                    sample_filter = TRUE,
                                    sites = c('ad','aa','dd'),
                                    min_cov = 10,
                                    sample_min_ncells = 20,
                                    celltype_min_samples = 2,
                                    seg_min_samples=4,
                                    seg_min_samples_fraq=0.2,
                                    seg_min_sd = 0.1,
                                    groupby='celltype'){
  
  sample_filter = sample_filter & data$ncells >= sample_min_ncells
  celltype_tab = table(data[[groupby]][sample_filter])
  sample_filter = sample_filter & data[[groupby]] %in% names(celltype_tab)[celltype_tab>=celltype_min_samples]
  
  
  data = data[rowData(data)$sites %in% sites,sample_filter]
  rowData(data)$nna = rowSums(assay(data,'i') + assay(data,'e') >= min_cov)
  
  seg_filter = rowData(data)$nna >= max(seg_min_samples,ncol(data)*seg_min_samples_fraq)
  data = data[seg_filter,]
  
  assay(data,'i') = as.matrix(assay(data,'i'))
  assay(data,'e') = as.matrix(assay(data,'e'))
  
  assay(data,'psi') = calcPSI(data,min_cov)
  rowData(data)$sd = apply(assay(data,'psi'),1,sd,na.rm=TRUE)
  data = data[rowData(data)$sd > seg_min_sd,]
  data
}

fitASglm = function (data, formula, terms, pseudocount = 0, .parallel = FALSE, 
                     .progress = "none", return.pv = FALSE, overdisp = TRUE, disp.param = NULL) 
{
  require(plyr)
  formula = terms(formula)
  term.names = attr(formula, "term.labels")
  r = alply(1:nrow(data), 1, function(i) {
    t = cbind(assay(data,'i')[i, ], assay(data,'e')[i, ])
    t = t + (t[, 1] + t[, 2]) * pseudocount
    terms$x = t
    t = tryCatch(glm(formula, data = terms, family = "quasibinomial"), 
                 error = function(e) {
                   warning(paste(e$message))
                   return(NA)
                 }, warning = function(w) {
                   warning(paste("Segment ", rownames(data)[i], 
                                 ": ", w$message, sep = ""))
                   return(NA)
                 })
    if (return.pv) 
      t = asqLRT(t, overdisp = overdisp, disp.param = disp.param[i], 
                 term.names = term.names, sid = rownames(data)[i])
    t
  }, .parallel = .parallel, .progress = .progress)
  if (return.pv) {
    r = do.call(rbind, r)
    rownames(r) = rownames(data)
    colnames(r) = c("overdispersion", term.names)
  }
  else {
    names(r) = rownames(data)
    attr(r, "term.labels") = term.names
  }
  r
}

asqLRT = function(g,overdisp,disp.param,term.names,sid){
  r = rep(NA,length(term.names)+1)
  if(is.na(g)[1])
    return(r)
  #set disp parameter
  if(is.null(disp.param)){
    r[1] = summary(g)$dispersion #sum(residuals(g, type = 'pearson')^2)/g$df.residual
    if(g$df.residual == 0){
      r[1] = NA
      if(overdisp)
        warning(paste('Segment ',sid,': cannot account for overdispersion in absense of replicates. Either use replicates or set overdisp parameter to FALSE',sep=''))
    }
  }else
    r[1] = disp.param
  disp = max(1,r[1],na.rm=TRUE)
  if(!overdisp) disp = 1
  a = tryCatch(anova(g,test='Chisq',dispersion=disp),error=function(e){return(NULL)})
  if(!is.null(a)){
    r[-1] = a[2:(1+length(term.names)),5]
  }
  r
}


getGroupbyFactor = function(data,groupby,sep=DELIMETER){
  if('colData' %in% slotNames(data)){ # I'll assume it is either data.frame or SummarizedExperiment
    data = as.data.frame(colData(data)) 
  }
  if(nrow(data) == length(groupby))
    return(groupby)
  if(all(groupby %in% colnames(data))){
    return(do.call(paste,c(data[,groupby,drop=FALSE],sep=sep)))
  }
  stop("Groupby is neither factor itself nor list of names of metadata columns")
}

testAllGroupsAS = function(data,groupby,.parallel=FALSE){
  groupby = getGroupbyFactor(data,groupby)
  pv = fitASglm(data,x ~ group,list(group=groupby),return.pv = TRUE,.parallel=.parallel)
  pv = as.data.frame(pv)
  pv$group_fdr = p.adjust(pv$group,m='BH')
  pv = cbind(pv,get_dPSI(data,data$celltype,min_cov = 50))
}

findMarkerAS = function(data,groupby,.parallel=FALSE,verbose=FALSE){
  groupby = getGroupbyFactor(data,groupby)
  ctpv = ctdpsi = NULL
  ugroups = unique(groupby)
  for(ct in ugroups){
    if(verbose){
      log_info(ct)
    }
    f = groupby == ct
    ctpv = cbind(ctpv,fitASglm(data,x ~ f,list(f=f),return.pv = TRUE,.parallel = .parallel)[,2])
    ctdpsi = cbind(ctdpsi,apply(assay(data,'psi')[,f,drop=F],1,mean,na.rm=T) - apply(assay(data,'psi')[,!f,drop=F],1,mean,na.rm=T))
  }
  colnames(ctpv) = colnames(ctdpsi) = ugroups
  
  list(pv=ctpv,
       fdr = apply(ctpv,2,p.adjust,m='BH'),
       dpsi = ctdpsi)
}


selectMarkers = function(mar,n=5,fdr_thr=0.05,dpsi_thr=0.1,clean_duplicates = TRUE){
  res = NULL
  for(ct in colnames(mar$fdr)){
    f = mar$fdr[,ct] <= fdr_thr & abs(mar$dpsi[,ct]) >= dpsi_thr
    f[is.na(f)] = FALSE
    if(sum(f)==0) next
    t = data.frame(pv  = mar$pv[f,ct],
                   fdr = mar$fdr[f,ct],
                   dpsi = mar$dpsi[f,ct],
                   seg_id = rownames(mar$pv)[f],
                   group=ct)
    t = t[order(abs(t$dpsi),decreasing = TRUE)[1:min(5,nrow(t))],]
    res = rbind(res,t)
  }
  rownames(res) = NULL
  if(clean_duplicates){
    res = lapply(split(res,res$seg_id),function(x)x[order(abs(x$dpsi),decreasing = TRUE)[1],])
    res = do.call(rbind,res)
  }
  res = res[order(abs(res$dpsi),decreasing = TRUE),]
  res
}

pseudobulk = function(data,groupby,cleanDerivatives=c('psi','cpm')){
  groupby = getGroupbyFactor(data,groupby)
  
  for(d in intersect(assayNames(data),cleanDerivatives)){
    assay(data,d) = NULL
  }
  
  res = list()
  for(a in assayNames(data)){
    res[[a]] = calcColSums(assay(data,a),groupby)  
  }
  
  # retain all column that have unique values within each groups
  meta = as.data.frame(colData(data))
  meta = pseudobulk_metadata(meta,groupby)
  meta = meta[colnames(res[[1]]),,drop=FALSE]
  
  res = SummarizedExperiment(assays = res,
                             rowRanges =rowRanges(data),
                             colData=meta)
  res
}

# derivative functions #################
calcPSI = function(data,min_cov=10){
  if(!all(c('i','e') %in% assayNames(data)))
    return(NULL)
  total = (assay(data,'i')+assay(data,'e'))
  psi = assay(data,'i')/total
  psi[total<min_cov] = NA
  as.matrix(psi)
}

calcCPM = function(data){
  if(!('counts' %in% assayNames(data)))
    return(NULL)
  d = assay(data,'counts')
  sweep(d,2,colSums(d),'/')*1e6
}


pseudobulk_metadata = function(meta,groupby,aggregate=list(ncells=sum)){
  groupby = getGroupbyFactor(meta,groupby)
  meta = split(meta,groupby)
  meta = lapply(meta,function(x){
    n = sapply(x,function(x)length(unique(x)))
    res = unique(x[,n==1,drop=FALSE])
    for(col in intersect(names(aggregate),colnames(x))){
      res[,col] = aggregate[[col]](x[,col])
    }
    res
  })
  
  columns = Reduce(intersect,lapply(meta,colnames))
  meta = do.call(rbind,lapply(meta,function(x)x[,columns,drop=FALSE]))
  meta
}

makeSummarizedExperiment = function(data,m){
  seg = data$seg
  data$seg = NULL
  class(data) = 'list'
  rowRanges = GRanges(seg$chr_id,
                      IRanges(start = seg$start,end = seg$stop),
                      strand=ifelse(seg$strand==1,'+',ifelse(seg$strand==-1,'-','*')),
                      feature_id=rownames(seg))
  seg$chr_id = seg$start = seg$stop = seg$strand = NULL
  
  elementMetadata(rowRanges)[names(seg)] = seg
  
  res = SummarizedExperiment(assays=data,
                             rowRanges=rowRanges,
                             colData=m)
  res
}


markerHeatmap = function(data,data_ge,groupby,markers,psi_scale=FALSE,
                         cpm_scale=TRUE,group_col=NULL,col=rev(hcl.colors(100, "RdYlBu")),
                         gene_names_col=NULL,...){
  data = pseudobulk(data,groupby)
  data_ge = pseudobulk(data_ge,groupby)
  
  assay(data,'psi') = calcPSI(data)
  assay(data_ge,'cpm') = calcCPM(data_ge)
  
  # PSI
  psi = t(assay(data,'psi')[markers$seg_id,])
  # choose isoform that is up
  psi[,markers$dpsi<0] = 1 - psi[,markers$dpsi<0]
  
  
  
  if(is.null(group_col))
    group_col = char2col(rownames(psi))
  
  psi.zlim = c(0,1)
  if(psi_scale){
    psi = apply(psi,2,scaleTo)
    psi.zlim = NULL
  }
  
  # group columns
  cor = cor(t(psi),u='pair')
  cor[is.na(cor)] = 0
  mds = cmdscale(1-cor,k = 1)
  groups = rownames(psi)[order(mds[,1])]
  
  markers = markers[order(markers$group),]
  markers = markers[order(match(markers$group,groups)),]
  
  psi = psi[groups,markers$seg_id]
  
  # CPM
  gids = rowData(data)[rownames(markers),'gene_id']
  cpm = t(assay(data_ge,'cpm')[gids,rownames(psi)])
  if(!is.null(gene_names_col)){
    colnames(cpm) = elementMetadata(rowRanges(data_ge))[gids,gene_names_col]
    colnames(psi) = paste0(colnames(cpm),':',colnames(psi))
  }
  cpm.zlim = range(cpm)
  if(cpm_scale){
    cpm = scale(cpm)
    cpm.zlim = max(abs(cpm))
    cpm.zlim = c(-cpm.zlim,cpm.zlim)
  }
  
  imageWithText(psi,'',colAnns = list(ct=rownames(psi)),
                rowAnns = list(ct=markers$group),
                colAnnCols = list(ct=group_col),
                rowAnnCols = list(ct=group_col),
                xaxlab=NULL,
                main='Alternative Splicing',
                col=col,zlim=psi.zlim,...)
  if(!psi_scale){
    plotColorLegend2(grconvertX(1.02, "npc", "nfc"), 1, 
                     grconvertY(0.1, "npc", "nfc"), 
                     grconvertY(0.9,"npc", "nfc"),
                     fullzlim = psi.zlim, zlim = psi.zlim,
                     zfun = identity,
                     z2col = function(x)num2col(x,col),title = 'PSI')
  }
  
  imageWithText(cpm,'',colAnns = list(ct=rownames(psi)),
                rowAnns = list(ct=markers$group),
                colAnnCols = list(ct=group_col),
                rowAnnCols = list(ct=group_col),
                main='Gene Expression',
                col=col,zlim=cpm.zlim,...)
  
  plotColorLegend2(grconvertX(1.02, "npc", "nfc"), 1, 
                   grconvertY(0.1, "npc", "nfc"), 
                   grconvertY(0.9,"npc", "nfc"),
                   fullzlim = cpm.zlim, zlim = cpm.zlim,
                   zfun = identity,
                   z2col = function(x)num2col(x,col),title = ifelse(cpm_scale,'z-score','CPM'))
  
}


# pipeline helpers ######
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


get_dPSI = function(data,groupby,min_cov=50){
  data = pseudobulk(data,groupby)
  psi = calcPSI(data,min_cov = min_cov)
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

findNearestConstantExons = function(data,sid,psi.thr=0.95){
  data = data[rowRanges(data)$gene_id == rowRanges(data)[sid,]@elementMetadata[1,'gene_id'],]
  i = rowSums(assay(data,'i'))
  e = rowSums(assay(data,'e'))
  psi = i/(i+e)
  segs = as.data.frame(rowRanges(data))
  segs$cnst = (is.na(psi) | psi>=psi.thr)
  
  strand = as.character(segs$strand[1])
  strand = ifelse(strand =='-',-1,1)
  segs = segs[order(strand*segs$start),]
  up = down = NA
  ucnst = which(segs$cnst & segs$sites %in% c('ad','sd'))
  dcnst = which(segs$cnst & segs$sites %in% c('ad','ae'))
  inx = which(rownames(segs) == sid)
  if(length(ucnst)>0 && any(ucnst<inx))
    up = rownames(segs)[max(ucnst[ucnst<inx])]
  if(length(dcnst)>0 && any(dcnst>inx))
    down = rownames(segs)[min(dcnst[dcnst>inx])]
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

plotSegmentCoverage = function(sid=NULL,usid=NULL,dsid=NULL,
                               chr=NULL,start = NULL,stop=NULL,
                               celltypes,
                               data,
                               groupby,
                               barcodes,
                               samples,
                               gene.descr,
                               scanBamFlags=list(isNotPassingQualityControls=FALSE,isDuplicate=FALSE,isSupplementaryAlignment=FALSE,isSecondaryAlignment=FALSE),
                               plot.junc.only.within=NA,
                               min.junc.cov.f = 0.01,
                               min.junc.cov = 3,
                               gtf){
  groupby = getGroupbyFactor(data,groupby)
  seg = as.data.frame(rowRanges(data))
  bams = unique(samples[,c('sample_id','bam_path')])
  
  if(is.null(start))
    start = seg[usid,'start']
  if(is.null(stop))
    stop = seg[dsid,'end']
  if(is.null(chr))
    chr = as.character(seg[sid,'seqnames'])
  strand = NA
  
  covs = list()
  for(ct in celltypes){
    cov = list()
    for(i in seq_len(nrow(bams))){
      tagFilter = list()
      tagFilter$CB = barcodes$barcode[barcodes$sample_id == bams$sample_id[i] & !is.na(barcodes$celltype) & barcodes$celltype==ct] 
      
      if(length(tagFilter$CB)==0) next
      cov[[length(cov)+1]] = getReadCoverage(bams$bam_path[i],
                                             chr,start,stop,strand=strand,scanBamFlags=scanBamFlags,tagFilter = tagFilter)
    }
    if(length(cov)>1)
      covs[[ct]] = sumCovs(cov)
    else
      covs[[ct]] = cov[[1]]
  }
  
  l = cbind(1,1+seq_len(1+length(celltypes)))
  
  if(!is.null(sid)){
    data = data[sid,]
    psi = calcPSI(data)[1,]
    psi = split(psi, groupby)
    psi = lapply(psi,na.omit)
    psi = psi[sapply(psi,length)>0]
    psi = psi[order(sapply(psi,mean))]
    
    ann = gtf[gtf$gene_id==seg[sid,'gene_id'],] # in this way it doesn't depends on chr naming
    ann$exon.col = 'black'
    ann$cds.col = 'black'
    f = ann$start >= seg[sid,'start'] & ann$stop <= seg[sid,'end']
    ann$exon.col[f]=ann$cds.col[f] = 'red'
    
  }else{
    l = matrix(seq_len(1+length(celltypes)),ncol=1)
    
    ann = gtf[gtf$start < stop & gtf$stop >= start,] # in this way it doesn't depends on chr naming
    ann$exon.col = 'black'
    ann$cds.col = 'black'
  }
  
  
  layout(l,widths = c(1,3))
  par(bty='n',tcl=-0.2,mgp=c(1.3,0.3,0),mar=c(3,13,1.5,0),oma=c(0,0,3,1))
  if(!is.null(sid))
    boxplot(psi,horizontal=TRUE,las=1,xlab='PSI')
  for(ct in celltypes){
    plotReadCov(covs[[ct]],xlim=c(start,stop),ylab='Coverage',xlab=chr,main=ct,
                plot.junc.only.within=plot.junc.only.within,min.junc.cov = min.junc.cov,
                min.junc.cov.f=min.junc.cov.f)
  }

  
  plotTranscripts(ann,new = T,exon.col = NA,cds.col = NA,xlim=c(start,stop))
  
  if(!is.null(sid))
    mtext(paste0(sid,' ', gene.descr[seg[sid,'gene_id'],'name'],'\n',gene.descr[seg[sid,'gene_id'],'descr']),side = 3,outer = TRUE)
}

readRDSifExists = function(f){
  r = NULL
  if(file.exists(f))
    r = readRDS(f)
  return(r)
}

log_info = function(text,...){
  text = paste0(text,...)
  print(paste0("INFO [",Sys.time(),"]: ",text))
}

# _for summary ##########
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

plotMDS = function(xy,celltypes,ct2col,samples,...){
  pch = as.numeric(factor(samples)) %% 25
  plot(xy,col=ct2col[celltypes],pch=pch,xlab='Dim 1',ylab='Dim 2',...)
  for(ct in unique(celltypes)){
    inx = which(celltypes == ct)
    if(length(inx)>1){
      inx = combn(inx,2)
      segments(xy[inx[1,],1],xy[inx[1,],2],xy[inx[2,],1],xy[inx[2,],2],col=ct2col[ct])
    }
  }
  for(ct in unique(celltypes)){
    inx = which(celltypes == ct)
    text(mean(xy[inx,1]),mean(xy[inx,2]),ct)
  }
}

getDitPlotSize = function(ccr){
  if(is.null(ccr))
    return(c(2,2))
  size = table(ccr@compareClusterResult$Cluster)
  size = size[size>0]
  size = c(height=sum(pmin(5,size)) * 0.15 + 3,
           width=length(size) * 0.6 + 7)
  size
}
