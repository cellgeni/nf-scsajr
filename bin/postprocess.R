options(error=function(e)quit('no',1))

library(Matrix)
library(SAJR)
library(visutils)
library(clusterProfiler)
library(org.Hs.eg.db)
library(plotCoverage)

# parse args ###########
# pbsas, mincells, minsamples, path2samples, path2barcodes, path2ref, path2bin, ncores
args = commandArgs(trailingOnly=TRUE)
writeLines(args,'params.txt')
#args = readLines('params.txt')

mincells = as.integer(args[2])
minsamples = as.integer(args[3])
path2samples = args[4]
path2barcodes = args[5]
path2ref = args[6]
path2bin = args[7]
ncores = as.numeric(args[8])
doMC::registerDoMC(ncores)
source(paste0(path2bin,'/sajr_utils.R'))

# load data #############
out.dir = args[1]
pbas_all = readRDS(paste0(out.dir,'/pbas.rds'))
samples  = read.table(path2samples,sep=' ')
colnames(samples) = c('sample_id','bam_path')
barcodes = read.table(path2barcodes,sep='\t')
colnames(barcodes) = c('sample_id','barcode','celltype')
# load annotation ###############
gtf = readRDS(paste0(path2ref,'/gtf.rds'))
gene.descr = readRDS(paste0(path2ref,'/functional_annotation/gene.descr.rds'))

# this part not necessary exists
domain2seg = readRDSifExists(paste0(path2ref,'/functional_annotation/domain2seg.df.rds'))
interpro   = readRDSifExists(paste0(path2ref,'/functional_annotation/interpro.rds'))

# filter #################
pbas = filterSegmentsAndSamples(pbas_all,celltype_min_samples = minsamples,sample_min_ncells = mincells)
log_info('samples*celltype filtering: ',ncol(pbas_all), ' -> ',ncol(pbas))
log_info('segment filtering: ',nrow(pbas_all), ' -> ',nrow(pbas))

# test ###################
# all celltypes together
pbas@metadata$all_celltype_test = testAllGroupsAS(pbas,'celltype',.parallel=TRUE)

# per celltype
pbas@metadata$markers = findMarkerAS(pbas,'celltype',.parallel=TRUE,verbose=TRUE)
log_info('diff AS finished')

# functional analyses #############
# _GO #############
perct = pbas@metadata$markers
pv = pbas@metadata$all_celltype_test

gene_uni = filterSegmentsAndSamples(pbas_all,seg_min_sd = 0,celltype_min_samples = minsamples,sample_min_ncells = mincells)
gene_uni = unique(rowData(gene_uni)$gene_id)
seg = rowData(pbas_all)

gids = apply(perct$fdr<0.05 & abs(perct$dpsi)>0.2,2,function(f){f[is.na(f)] = FALSE;unique(seg[rownames(perct$fdr)[f],'gene_id'])})
sgn = pv$group_fdr<0.05 & pv$dpsi>0.3
sgn[is.na(sgn)] = FALSE
gids$all = unique(seg[rownames(pv)[sgn],'gene_id'])


tryCatch({
  pbas@metadata$go  = compareCluster(gids,
                       fun='enrichGO',
                       universe      = gene_uni,
                       pAdjustMethod = "BH",
                       ont='ALL',
                       OrgDb = 'org.Hs.eg.db',
                       keyType = 'ENSEMBL')
},error=function(e){
  warning('GO analyses was unsuccesfull')
  print(e)
})
log_info('GO finished')

# _domains ##############
# I'll use only cassette exons
if(!is.null(domain2seg)){
  domain_sites2use = 'ad'
  seg_uni = rownames(pbas_all)[seg$is_exon & seg$sites %in% domain_sites2use & seg$gene_id %in% gene_uni]
  
  sids = apply(perct$fdr<0.05 & abs(perct$dpsi)>0.2,2,function(f){
    f[is.na(f)] = FALSE
    out = rownames(perct$fdr)[f]
    out[seg[out,'sites'] == 'ad' & seg[out,'is_exon']]
  })
  sgn = pv$group_fdr<0.05 & pv$dpsi>0.3
  sgn[is.na(sgn)] = FALSE
  sids$all = rownames(pv)[sgn]
  sids$all = sids$all[seg[sids$all,'sites'] == 'ad' & seg[sids$all,'is_exon']]
  
  pbas@metadata$ipro   = compareCluster(sids,
                          fun='enricher',
                          universe      = seg_uni,
                          pAdjustMethod = "BH",
                          TERM2GENE     = domain2seg$interpro,
                          TERM2NAME = interpro[,-2])
}

saveRDS(pbas,paste0(out.dir,'/pb_as_filtered.rds'))
#pbas = readRDS(paste0(out.dir,'/pb_as_filtered.rds'))

log_info('interpro finished')

# Example coverage plots ##########
N = 100
f = pv$group_fdr<0.05 & pv$dpsi > 0.3 & rowData(pbas_all)[rownames(pv),'sites'] %in% c('aa','ad','dd') & rowData(pbas_all)[rownames(pv),'is_exon']
f[is.na(f)] = FALSE
toplot = pv[f, ]
toplot = toplot[order(toplot$dpsi,decreasing = TRUE)[1:min(nrow(toplot),N*3)],]
const_exons = t(sapply(rownames(toplot),function(sid){findNearestConstantExons(pbas_all,sid)}))
f = !is.na(const_exons[,1]) & !is.na(const_exons[,2])
toplot = cbind(toplot[f,],const_exons[f,])
toplot = toplot[order(toplot$dpsi,decreasing = TRUE)[1:min(N,nrow(toplot))],]


pdf('examples.pdf',w=12,h=9)
for(i in seq_len(nrow(toplot))){
  cat('\r',i)
  plotSegmentCoverage(sid = rownames(toplot)[i],
                      usid = toplot$up[i],
                      dsid = toplot$down[i],
                      data = pbas_all,
                      groupby = pbas_all$celltype,
                      celltypes = c(toplot$high_state[i],toplot$low_state[i]),
                      barcodes=barcodes,
                      samples = samples,
                      gene.descr = gene.descr,
                      gtf=gtf)  
}
dev.off()
log_info('examples finished')
