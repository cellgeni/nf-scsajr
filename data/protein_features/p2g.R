library(AnnotationHub)
inx = commandArgs(trailingOnly=TRUE)[1]
cat("Chunk ",inx,"\n")
ah=AnnotationHub(cache='/lustre/scratch127/cellgen/cellgeni/pasham/rcache')
z = query(ah, c('EnsDb', 'v98','Homo sapiens'))
db = ah[[z$ah_id]]
p = readRDS(paste0('pd',inx,'.rds'))
irs = IRanges(p$prot_dom_start,p$prot_dom_end,names=p$protein_id)
g = proteinToGenome(irs,db)
saveRDS(g,paste0('g',inx,'.rds'))