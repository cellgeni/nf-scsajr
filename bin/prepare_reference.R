options(error=function(e)quit('no',1))

library(SAJR)
library(plotCoverage)

# annotation in gtf, converted annotation in sajr, path2bin
args = commandArgs(trailingOnly=TRUE)

gtf_path = args[1]
sajr_path = args[2]

seg = SAJR::loadSAData(sajr_path)
seg = SAJR::setSplSiteTypes(seg,sajr_path)$seg
seg$length = seg$stop-seg$start+1

gtf = loadEnsGTF(gtf_path)
if(is.null(gtf$transcript_name))
  gtf$transcript_name = gtf$transcript_id


# which segments are identical to exons
ef = gtf$feature=='exon'
seg$is_exon = paste(seg$gene_id,seg$start,seg$stop) %in% paste(gtf$gene_id[ef],gtf$start[ef],gtf$stop[ef])

gene_descr =  SAJR::loadGData(sajr_path)$gene
gene_descr$name = gene_descr$descr = rownames(gene_descr)

# save #############
write.csv(seg,'segments.csv')
saveRDS(gtf,'gtf.rds')
dir.create('functional_annotation')
saveRDS(gene_descr,'functional_annotation/gene.descr.rds')
