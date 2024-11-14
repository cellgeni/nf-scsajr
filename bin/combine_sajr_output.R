options(error=function(e)quit('no',1))

library(Matrix)
library(SAJR)
library(visutils)
library(doMC)
library(plyr)

# parse args ##############
# samples_file, barcodes_file, path2ref (folder), path2bin
args = commandArgs(trailingOnly=TRUE)
writeLines(args,'params.txt')
#args = readLines('params.txt')
path2samples = args[1]
path2barcodes = args[2]
path2ref = args[3]
path2bin = args[4]
ncores = as.integer(args[5])
DEL='$' # delimeter to separate sample and celltype names
source(paste0(path2bin,'/sajr_utils.R'))
doMC::registerDoMC(ncores)

# load input ############
samples  = read.table(path2samples,sep=' ')
colnames(samples) = c('sajr_out','sample_id','bam_path','bai_path','strand')
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
pbasl = llply(seq_len(nrow(samples)),function(i){
  r = loadSC_AS(seg,paste0(samples$sajr_out[i],'/',samples$sample_id[i]))
  colnames(r$e) = colnames(r$i) = paste0(samples$sample_id[i],'|',colnames(r$i))
  saveRDS(r,paste0('rds/',samples$sample_id[i],'_',samples$strand[i],'.rds'))
  
  # calc pseudobulks
  cmn = intersect(rownames(barcodes) , colnames(r$i))
  r$i = r$i[,cmn]
  r$e = r$e[,cmn]
  ncell = rowSums(r$i + r$e > 9)
  f = paste0(barcodes[cmn,'sample_id'],DEL,barcodes[cmn,'celltype'])
  pb = list()
  pb$i = as.matrix(visutils::calcColSums(r$i,f))
  pb$e = as.matrix(visutils::calcColSums(r$e,f))
  pb$cmn = cmn
  pb$ncell = ncell 
  
  rm(r)
  mem = gc()
  mem = round(sum(mem[,2])/2^10,digits = 2)
  pbmem = round(unclass(object.size(pb))/2^30,digits = 2)
  log_info(samples$sample_id[i]," is loaded. Total mem used: ",mem," Gb. Pb size: ",pbmem," Gb")
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

# chose strand #######
samples$total = sapply(pbasl,function(x)sum(x$i)+sum(x$e))
cov = as.data.frame(castXYtable(samples$sample_id, samples$strand,samples$total))
strand = colnames(cov)[unique(apply(cov,1,which.max))]
if(length(strand)!=1){
  print(cov)
  stop("Cannot autodetermine strand")
}
log_info("autodetermined strand = '",strand,"'")

f = samples$strand == strand
file.remove(paste0('rds/',samples$sample_id[!f],'_',samples$strand[!f],'.rds'))

pbasl = pbasl[f]
samples = samples[f,]
names(pbasl) = samples$sample_id
cmnbarcodes = unlist(cmnbarcodesl[f])

pbas = list(seg=seg,i=NULL,e=NULL)
for(n in names(pbasl)){
  pbas$i = cbind(pbas$i,pbasl[[n]]$i)
  pbas$e = cbind(pbas$e,pbasl[[n]]$e)
  pbasl[[n]] = NULL
}
pbas$seg$ncell = ncell
class(pbas) = c('sajr','list')

# all(colnames(pbas$i)== colnames(pbas$e))

barcodes = barcodes[cmnbarcodes,]

pbmeta = as.data.frame(do.call(rbind,strsplit(colnames(pbas$i),DEL,fixed = TRUE)))
colnames(pbmeta) = c('sample_id','celltype')
rownames(pbmeta) = colnames(pbas$i)
pbmeta$ncells = as.numeric(table(paste0(barcodes$sample_id,DEL,barcodes$celltype))[rownames(pbmeta)])

pbas$ir = pbas$i/(pbas$i+pbas$e)
pbas$ir[pbas$i+pbas$e<10] = NA
pbas$seg$sd = apply(pbas$ir,1,sd,na.rm=TRUE)
pbas$seg$sd[is.na(pbas$seg$sd)] = 0
pbas$seg$nna = apply(!is.na(pbas$ir),1,sum)
log_info('pseudobulk is made')

# save ##############
saveRDS(barcodes,paste0(out.dir,'/cell_meta.rds'))
saveRDS(pbas,paste0(out.dir,'/pb_as.rds'))
saveRDS(pbmeta,paste0(out.dir,'/pb_meta.rds'))