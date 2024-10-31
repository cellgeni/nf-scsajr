options(error=function(e)quit('no',1))

library(Matrix)
library(SAJR)
library(visutils)
library(clusterProfiler)
library(org.Hs.eg.db)

# parse args ##############
# samples_file, barcodes_file, path_to_data, path to ref
args = commandArgs(trailingOnly=TRUE)
writeLines(args,'params.txt')
#args = readLines('params.txt')
path2samples = args[1]
path2barcodes = args[2]
path2data = args[3]
path2ref = args[4]
DEL='$' # delimeter to separate sample and celltype names

source(paste0(path2data,'/../bin/plotCoverage.R'))
source(paste0(path2data,'/../bin/sajr_utils.R'))


samples  = read.table(path2samples,sep=' ')
colnames(samples) = c('sajr_out','sample_id','bam_path','bai_path','strand')
barcodes = read.table(path2barcodes,sep='\t')
colnames(barcodes) = c('sample_id','barcode','celltype')
barcodes$celltype[is.na(barcodes$celltype)] = 'NA'

rownames(barcodes) = paste0(barcodes$sample_id,'|',barcodes$barcode)

seg =  SAJR::loadSAData(path2ref)$seg

gene.descr = readRDS(paste0(path2data,'/gene.descr.rds'))
seg$length = seg$stop-seg$start+1

domain2seg = readRDS(paste0(path2data,'/protein_features/domain2seg.df.rds'))
interpro = readRDS(paste0(path2data,'/protein_features/interpro.rds'))

gtf = readRDS(paste0(path2data,'/gtf.rds'))

if(any(grepl(DEL,c(barcodes$sample_id,barcodes$celltype),fixed = TRUE))){
  stop(paste0("Celltype names or sample IDs contains '",DEL,"' that is not supported, please get rid of it." ))
}
# which segments are identical to exons
ef = gtf$feature=='exon'
seg$is_exon = paste(seg$gene_id,seg$start,seg$stop) %in% paste(gtf$gene_id[ef],gtf$start[ef],gtf$stop[ef])

out.dir = 'rds'
dir.create(out.dir)

# filters ######
# segments
sites = c('ad','aa','dd')
segFilter = function(d,sites,min.sd=0.1){
  # d$seg$type=='ALT' & 
  d$seg$nna >= max(4,nrow(pbmeta)*0.2) & 
    d$seg$sd > min.sd & 
    d$seg$sites %in% sites
}

# samples
mincells = 29
minsamples = 2
log_info('initializing finished')

# load ####################
pbasl = list()
cmnbarcodesl = list()
ncell = rep(0,nrow(seg))

for(i in seq_len(nrow(samples))){
  log_info('load ',samples$sample_id[i])
  r = loadSC_AS(seg,paste0(samples$sajr_out[i],'/',samples$sample_id[i]))
  colnames(r$e) = colnames(r$i) = paste0(samples$sample_id[i],'|',colnames(r$i))
  cmn = intersect(rownames(barcodes) , colnames(r$i))
  r$i = r$i[,cmn]
  r$e = r$e[,cmn]
  saveRDS(r,paste0('rds/',samples$sample_id[i],'_',samples$strand[i],'.rds'))
  ncell = ncell + rowSums(r$i + r$e > 9)
  
  # calc pseudobulks
  f = paste0(barcodes[cmn,'sample_id'],DEL,barcodes[cmn,'celltype'])
  pb = list(seg=seg)
  pb$i = as.matrix(visutils::calcColSums(r$i,f))
  pb$e = as.matrix(visutils::calcColSums(r$e,f))
  pbasl[[length(pbasl)+1]] = pb
  cmnbarcodesl[[length(cmnbarcodesl)+1]] = cmn
}

log_info('data loaded')

# chose strand
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

# save 
saveRDS(barcodes,paste0(out.dir,'/cell_meta.rds'))
saveRDS(pbas,paste0(out.dir,'/pb_as.rds'))
saveRDS(pbmeta,paste0(out.dir,'/pb_meta.rds'))
log_info('pseudobulk is made')

# filter #################
# pbas = readRDS(paste0(out.dir,'/pb_as.rds'))
# pbmeta = readRDS(paste0(out.dir,'/pb_meta.rds'))
pbas_all = pbas 
# pbsamples
sf = pbmeta$ncells>mincells
nsam = table(pbmeta$celltype[sf])
sf = pbmeta$ncells>mincells & pbmeta$celltype %in% names(nsam)[nsam>=minsamples]
# table(pbmeta$celltype,sf)

pbmeta = pbmeta[sf,]
pbas = pbas[,sf]

# segments

f = segFilter(pbas,sites)
# table(f)
# table(pbas$seg$sites[f],pbas$seg$cod[f])

pbas = pbas[f,]

saveRDS(pbas,paste0(out.dir,'/pb_as_filtered.rds'))
saveRDS(pbmeta,paste0(out.dir,'/pb_meta_filtered.rds'))

log_info('samples*celltype filtering: ',length(sf), ' -> ',sum(sf))
log_info('segment filtering: ',length(pbas_all), ' -> ',length(pbas))

# test ###################
# all celltypes together
pv = fitSAGLM(pbas,x ~ celltype,pbmeta,return.pv = TRUE)
pv = as.data.frame(pv)
pv$celltype_fdr = p.adjust(pv$celltype,m='BH')
pv = cbind(pv,get_dPSI(d=pbas,f=pbmeta$celltype,mincov = 50))
saveRDS(pv,paste0(out.dir,'/all_celltypes_test.rds'))
# pv = readRDS(paste0(out.dir,'/all_celltypes_test.rds'))

# per celltype
ctpv = ctdpsi = NULL
cts = unique(pbmeta$celltype)
for(ct in cts){
  #print(ct)
  f = pbmeta$celltype == ct
  ctpv = cbind(ctpv,fitSAGLM(pbas,x ~ f,list(f=f),return.pv = TRUE)[,2])
  ctdpsi = cbind(ctdpsi,apply(pbas$ir[,f,drop=F],1,mean,na.rm=T) - apply(pbas$ir[,!f,drop=F],1,mean,na.rm=T))
}
colnames(ctpv) = colnames(ctdpsi) = cts

perct = list(pv=ctpv,
             fdr = apply(ctpv,2,p.adjust,m='BH'),
             dpsi = ctdpsi)
saveRDS(perct,paste0(out.dir,'/celltypes_marker_test.rds'))

log_info('diff AS finished')
# functional analyses #############
## GO
# gene_uni = unique(pbas$seg$gene_id)
# length(gene_uni)
f = segFilter(pbas_all,sites,min.sd = 0)
gene_uni = unique(pbas_all$seg$gene_id[f])

gids = apply(perct$fdr<0.05 & abs(perct$dpsi)>0.2,2,function(f){f[is.na(f)] = FALSE;unique(seg[rownames(perct$fdr)[f],'gene_id'])})
sgn = pv$celltype_fdr<0.05 & pv$dpsi>0.3
sgn[is.na(sgn)] = FALSE
gids$all = unique(seg[rownames(pv)[sgn],'gene_id'])
#sapply(gids,length)

go  = compareCluster(gids,
                     fun='enrichGO',
                     universe      = gene_uni,
                     pAdjustMethod = "BH",
                     ont='ALL',
                     OrgDb = 'org.Hs.eg.db',
                     keyType = 'ENSEMBL')


saveRDS(go,paste0(out.dir,'/go.rds'))
log_info('GO finished')
# domains
# I'll use only cassette exons
domain_sites2use = 'ad'
seg_uni = rownames(seg)[seg$is_exon & seg$sites %in% domain_sites2use & seg$gene_id %in% gene_uni]
length(seg_uni)

sids = apply(perct$fdr<0.05 & abs(perct$dpsi)>0.2,2,function(f){
  f[is.na(f)] = FALSE
  out = rownames(perct$fdr)[f]
  out[seg[out,'sites'] == 'ad' & seg[out,'is_exon']]
})
sgn = pv$celltype_fdr<0.05 & pv$dpsi>0.3
sgn[is.na(sgn)] = FALSE
sids$all = rownames(pv)[sgn]
sids$all = sids$all[seg[sids$all,'sites'] == 'ad' & seg[sids$all,'is_exon']]
sapply(sids,length)

ipro   = compareCluster(sids,
                        fun='enricher',
                        universe      = seg_uni,
                        pAdjustMethod = "BH",
                        TERM2GENE     = domain2seg$interpro,
                        TERM2NAME = interpro[,-2])
#dotplot(ipro)
saveRDS(ipro,paste0(out.dir,'/ipro.rds'))

log_info('interpro finished')

# Example coverage plots ##########
seg$mean_psi = apply(pbas_all$ir,1,mean,na.rm=TRUE)
N = 100
f = pv$celltype_fdr<0.05 & pv$dpsi > 0.3 & seg[rownames(pv),'sites'] %in% c('aa','ad','dd') & seg[rownames(pv),'is_exon']
f[is.na(f)] = FALSE
toplot = pv[f, ]
toplot = toplot[order(toplot$dpsi,decreasing = TRUE)[1:min(nrow(toplot),N*3)],]
const_exons = t(sapply(rownames(toplot),function(sid){findNearestConstantExons(seg,sid)}))
f = !is.na(const_exons[,1]) & !is.na(const_exons[,2])
toplot = cbind(toplot[f,],const_exons[f,])
toplot = toplot[order(toplot$dpsi,decreasing = TRUE)[1:min(N,nrow(toplot))],]


pdf('examples.pdf',w=12,h=9)
for(i in seq_len(nrow(toplot))){
  cat('\r',i)
  plotSegmentCoverage(sid = rownames(toplot)[i],
                      usid = toplot$up[i],
                      dsid = toplot$down[i],
                      celltypes = c(toplot$high_state[i],toplot$low_state[i]),
                      seg=seg,barcodes=barcodes)  
}
dev.off()
log_info('examples finished')
