options(error=function(e)quit('no',1))

library(Matrix)
library(visutils)
library(doMC)
library(plyr)

# parse args ##############
# sajr_outs_file, barcodes_file, path2ref (folder), path2bin
args = commandArgs(trailingOnly=TRUE)
writeLines(args,'params.txt')
#args = readLines('params.txt')
path2sajr_outs = args[1]
path2barcodes = args[2]
path2ref = args[3]
path2bin = args[4]
ncores = as.integer(args[5])
DEL='$' # delimeter to separate sample and celltype names
source(paste0(path2bin,'/sajr_utils.R'))
doMC::registerDoMC(ncores)

# load input ############
sajr_outs  = read.table(path2sajr_outs,sep=' ')
colnames(sajr_outs) = c('sample_id','chr','sajr_out','bam_path','strand')
barcodes = read.table(path2barcodes,sep='\t')
colnames(barcodes) = c('sample_id','barcode','celltype')
barcodes$celltype[is.na(barcodes$celltype)] = 'NA'

rownames(barcodes) = paste0(barcodes$sample_id,'|',barcodes$barcode)

if(any(grepl(DEL,c(barcodes$sample_id,barcodes$celltype),fixed = TRUE))){
  stop(paste0("Celltype names or sample IDs contains '",DEL,"' that is not supported, please get rid of it." ))
}

seg = read.csv(paste0(path2ref,"/segments.csv"),row.names = 1)

out.dir = 'rds'
dir.create(out.dir)
log_info('initializing finished')

# load sajr ####################
samples = unique(sajr_outs$sample_id)

pbasl = llply(seq_along(samples),function(i){
  sample_sajr_outs = sajr_outs[sajr_outs$sample_id==samples[i],]
  
  r_ = lapply(seq_len(nrow(sample_sajr_outs)),function(j){
    loadSC_AS(seg[seg$chr_id==sample_sajr_outs$chr[j],],paste0(sample_sajr_outs$sajr_out[j],'/',sample_sajr_outs$chr[j]))
  })
  
  cl = unique(unlist(lapply(r_,function(x)colnames(x$i))))
  
  r = list(i = rbindMatrix(lapply(r_,function(x)x$i),cl),
           e = rbindMatrix(lapply(r_,function(x)x$e),cl))
  r$i = r$i[rownames(seg),]
  r$e = r$e[rownames(seg),]
  
  colnames(r$e) = colnames(r$i) = paste0(sample_sajr_outs$sample_id[i],'|',colnames(r$i))
  saveRDS(r,paste0('rds/',samples[i],'.rds'))
  
  # calc pseudobulks
  cmn = intersect(rownames(barcodes) , colnames(r$i))
  if(length(cmn)==0)
    return(list(i=NULL,e=NULL,cmn = cmn,ncell = 0))
  r$i = r$i[,cmn]
  r$e = r$e[,cmn]
  ncell = rowSums(r$i + r$e > 9)
  f = paste0(barcodes[cmn,'sample_id'],DEL,barcodes[cmn,'celltype'])
  pb = list()
  pb$i = visutils::calcColSums(r$i,f)
  pb$e = visutils::calcColSums(r$e,f)
  pb$cmn = cmn
  pb$ncell = ncell 
  
  rm(r)
  mem = gc()
  mem = round(sum(mem[,2])/2^10,digits = 2)
  pbmem = round(unclass(object.size(pb))/2^30,digits = 2)
  log_info(sample_sajr_outs$sample_id[i]," is loaded. Total mem used: ",mem," Gb. Pb size: ",pbmem," Gb")
  pb
},.parallel = TRUE)

cmnbarcodesl = list()
ncell = rep(0,nrow(seg))

for(i in seq_along(pbasl)){
  cmnbarcodesl[[i]] = pbasl[[i]]$cmn
  ncell = ncell + pbasl[[i]]$ncell
  pbasl[[i]]$cmn = pbasl[[i]]$ncell = NULL
}

gc()
log_info('data loaded')

names(pbasl) = samples
cmnbarcodes = unlist(cmnbarcodesl)

pbas = list(seg=seg,i=NULL,e=NULL)
for(n in names(pbasl)){
  pbas$i = cbind(pbas$i,pbasl[[n]]$i)
  pbas$e = cbind(pbas$e,pbasl[[n]]$e)
  pbasl[[n]] = NULL
}
pbas$seg$ncell = ncell

barcodes = barcodes[cmnbarcodes,]

pbmeta = as.data.frame(do.call(rbind,strsplit(colnames(pbas$i),DEL,fixed = TRUE)))
colnames(pbmeta) = c('sample_id','celltype')
rownames(pbmeta) = colnames(pbas$i)
pbmeta$ncells = as.numeric(table(paste0(barcodes$sample_id,DEL,barcodes$celltype))[rownames(pbmeta)])
s2s = unique(sajr_outs[,c('sample_id','strand')])
rownames(s2s) = s2s$sample_id
pbmeta$strand =  s2s[pbmeta$sample_id,'strand']

pbas_se = makeSummarizedExperiment(pbas,pbmeta)

log_info('pseudobulk is made')

# save ##############

saveRDS(barcodes,paste0(out.dir,'/cell_meta.rds'))
saveRDS(pbas_se,paste0(out.dir,'/pbas.rds'))
