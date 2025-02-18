gtf = readRDS('../human_2020A_chr/gtf.rds')
gtf$chr_id = sub('^chr','',gtf$chr_id)
saveRDS(gtf,'gtf.rds')
