options(error=function(e)quit('no',1))

library(Matrix)
library(SAJR)
library(visutils)
library(clusterProfiler)
library(org.Hs.eg.db)

# parse args ###########
# pbsas, pbmeta, mincells, minsamples, path2samples, path2barcodes, path2ref, path2bin
args = commandArgs(trailingOnly=TRUE)
writeLines(args,'params.txt')
#args = readLines('params.txt')

mincells = as.integer(args[3])
minsamples = as.integer(args[4])
path2samples = args[5]
path2barcodes = args[6]
path2ref = args[7]
path2bin = args[8]
source(paste0(path2bin,'/plotCoverage.R'))
source(paste0(path2bin,'/sajr_utils.R'))

# load data #############
pbas = readRDS(args[1])
pbmeta = readRDS(args[2])
out.dir = 'rds' # it is passed from previous processes
samples  = read.table(path2samples,sep='\t')
colnames(samples) = c('sample_id','bam_path')
barcodes = read.table(path2barcodes,sep='\t')
colnames(barcodes) = c('sample_id','barcode','celltype')
# load annotation ###############
gtf = readRDS(paste0(path2ref,'/gtf.rds'))

# this part not necessary exists
gene.descr = readRDS(paste0(path2ref,'/functional_annotation/gene.descr.rds'))
domain2seg = readRDS(paste0(path2ref,'/functional_annotation/domain2seg.df.rds'))
interpro = readRDS(paste0(path2ref,'/functional_annotation/interpro.rds'))

# filters ######
# segments
sites = c('ad','aa','dd')
segFilter = function(d,sites,min.sd=0.1){
  # d$seg$type=='ALT' & 
  d$seg$nna >= max(4,nrow(pbmeta)*0.2) & 
    d$seg$sd > min.sd & 
    d$seg$sites %in% sites
}

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
# perct = readRDS(paste0(out.dir,'/celltypes_marker_test.rds'))

log_info('diff AS finished')
# functional analyses #############
# _GO #############
# gene_uni = unique(pbas$seg$gene_id)
# length(gene_uni)
f = segFilter(pbas_all,sites,min.sd = 0)
gene_uni = unique(pbas_all$seg$gene_id[f])

gids = apply(perct$fdr<0.05 & abs(perct$dpsi)>0.2,2,function(f){f[is.na(f)] = FALSE;unique(pbas_all$seg[rownames(perct$fdr)[f],'gene_id'])})
sgn = pv$celltype_fdr<0.05 & pv$dpsi>0.3
sgn[is.na(sgn)] = FALSE
gids$all = unique(pbas_all$seg[rownames(pv)[sgn],'gene_id'])
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
# _domains ##############
# I'll use only cassette exons
domain_sites2use = 'ad'
seg_uni = rownames(pbas_all$seg)[pbas_all$seg$is_exon & pbas_all$seg$sites %in% domain_sites2use & pbas_all$seg$gene_id %in% gene_uni]
length(seg_uni)

sids = apply(perct$fdr<0.05 & abs(perct$dpsi)>0.2,2,function(f){
  f[is.na(f)] = FALSE
  out = rownames(perct$fdr)[f]
  out[pbas_all$seg[out,'sites'] == 'ad' & pbas_all$seg[out,'is_exon']]
})
sgn = pv$celltype_fdr<0.05 & pv$dpsi>0.3
sgn[is.na(sgn)] = FALSE
sids$all = rownames(pv)[sgn]
sids$all = sids$all[pbas_all$seg[sids$all,'sites'] == 'ad' & pbas_all$seg[sids$all,'is_exon']]
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
pbas_all$seg$mean_psi = apply(pbas_all$ir,1,mean,na.rm=TRUE)
N = 100
f = pv$celltype_fdr<0.05 & pv$dpsi > 0.3 & pbas_all$seg[rownames(pv),'sites'] %in% c('aa','ad','dd') & pbas_all$seg[rownames(pv),'is_exon']
f[is.na(f)] = FALSE
toplot = pv[f, ]
toplot = toplot[order(toplot$dpsi,decreasing = TRUE)[1:min(nrow(toplot),N*3)],]
const_exons = t(sapply(rownames(toplot),function(sid){findNearestConstantExons(pbas_all$seg,sid)}))
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
                      seg=pbas_all$seg,barcodes=barcodes,
                      gtf=gtf)  
}
dev.off()
log_info('examples finished')
