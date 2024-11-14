library(visutils)
library(Seurat)
library(Matrix)

as  = readRDS('test/limb_visium/scsajr_out/rds/cell_as.rds')
mas = readRDS('test/limb_visium/scsajr_out/rds/cell_meta.rds')
pb  = readRDS('test/limb_visium/scsajr_out/rds/pb_as.rds')
pv  = readRDS('test/limb_visium/scsajr_out/rds/all_celltypes_test.rds')
path2cr = '/nfs/cellgeni/pasham/projects/2302.fetal.skin/data.lustre/visium/spaceranger/'
gene.descr = readRDS('data/gene.descr.rds')

pb$seg$ncell = rowSums(as$i+as$e > 9)
f = pb$seg$sites=='ad' & pb$seg$ncell > 1000 & pb$seg$sd>0.2
t=pb$seg[f,]
t[order(t$ncell),]

sids = unique(mas$sample_id)

vs = lapply(sids,function(s){
  v = myLoad10X_Spatial(paste0(path2cr,'/',s,'/outs'),filter.matrix = T,ens_id = TRUE,slice = s)
  v$barcode = colnames(v)
  v$library_id = s
  colnames(v) = paste0(s,'|',v$barcode)
  v
})
names(vs) = sids

segid = 'ENSG00000138326.s12' # ribosomal protein S24; good
segid = 'ENSG00000092841.s17' # myosin light chain 6; so-so
segid = 'ENSG00000198467.s10' # tropomyosin 2; few spots, so-so
segid = 'ENSG00000130595.s24' # troponin T3, fast skeletal type; few spots, so-so


gid = as$seg[segid,'gene_id']
gene.descr[gid,]

pv[segid,]
i = as$i[segid,]
e = as$e[segid,]
psi = i/(i+e)
#psi[i+e < 10] = NA

par(mfrow=c(2,2),mar=c(0,0,1,6))
for(s in names(vs)){
  plotVisium(vs[[s]],psi[colnames(vs[[s]])],main=s,legend.args = list(title='PSI'),zlim=c(0,0.7),cex=0.)
}


par(mfrow=c(2,2),mar=c(0,0,1,6))
for(s in names(vs)){
  cpm = FetchData(vs[[s]],gid,layer = 'count')[,1]/vs[[s]]$nCount_Spatial*1e4
  plotVisium(vs[[s]],cpm,main=s,zfun = log1p,legend.args = list(title='CP10K'))
}