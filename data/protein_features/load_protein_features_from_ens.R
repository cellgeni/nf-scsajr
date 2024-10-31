library(AnnotationHub)

# load in protein coordinates
ah=AnnotationHub(cache='/lustre/scratch127/cellgen/cellgeni/pasham/rcache')
z = query(ah, c("EnsDb", "v98",'Homo sapiens'))
db = ah[[z$ah_id]]
pd = proteins(db, columns = c("tx_id", listColumns(db, "protein_domain")))
pd
table(pd$protein_domain_source,useNA='always')
pd = pd[!is.na(pd$protein_domain_source),]
pd = as.data.frame(pd)
saveRDS(pd,'/nfs/cellgeni/pasham/code/sajr-java-sc/protein_features/pd1.rds')
#dir.create('/nfs/cellgeni/pasham/code/sajr-java-sc/protein_features/tmp')

# split in chunks, otherwise too slow
pd$feature_id = seq_len(nrow(pd))
inx = seq_len(nrow(pd))
inx = split(inx,inx %/% 10000)
length(inx)
sum(sapply(inx,length)) 
nrow(pd)
for(i in seq_along(inx)){
  print(i)
  t = pd[inx[[i]],]
  #saveRDS(t,paste0('/nfs/cellgeni/pasham/code/sajr-java-sc/protein_features/tmp/pd',i,'.rds'))
}


# run p2g.sh
# load back
# combine per feature info
gall = list()
for(i in seq_along(inx)){
  cat("\nchunk ",i,'\n')
  p = readRDS(paste0('/nfs/cellgeni/pasham/code/sajr-java-sc/protein_features/tmp/pd',i,'.rds'))
  g = readRDS(paste0('/nfs/cellgeni/pasham/code/sajr-java-sc/protein_features/tmp/g',i,'.rds'))
  colnames(p) = paste0('prot_',colnames(p))
  pcols = paste0('prot_',c("feature_id","tx_id","protein_id","protein_domain_id","protein_domain_source","interpro_accession",'prot_dom_start','prot_dom_end') )
  
  for(j in seq_along(g)){
    #cat('\r',j)
    if(length(g[[j]])>0)
      g[[j]] = cbind(p[j,pcols],as.data.frame(g[[j]]))
  }
  gall = c(gall,g)
}
length(gall)
gall[[1]]

#saveRDS(gall,'/nfs/cellgeni/pasham/code/sajr-java-sc/protein_features/tmp/gall.rds')
gall = readRDS('/nfs/cellgeni/pasham/code/sajr-java-sc/protein_features/tmp/gall.rds')
# concatenate into one table
l = sapply(gall,length)
table(l==0)
gall[l==0]
gall = gall[l > 0]
length(gall)

N = sum(sapply(gall,nrow))
res = list()
for(c in colnames(gall[[1]])){
  res[[c]] = rep(NA,N)
}
print(object.size(res),u='Mb')

t  = Sys.time()
n = 1
for(i in seq_along(gall)){
  if(i %% 1000 == 0){
    r = (Sys.time() - t)/i*length(gall)
    cat('\r',i,r,' ')
  }
  for(c in colnames(gall[[i]])){
    if(is.factor(gall[[i]][,c]))
      gall[[i]][,c] = as.character(gall[[i]][,c])
    res[[c]][n:(n+nrow(gall[[i]])-1)] = gall[[i]][,c]
  }
  n = n + nrow(gall[[i]])
}

sapply(res,function(x)sum(is.na(x)))
resdf = as.data.frame(res)

resdf[1:10,]
table(resdf$seqnames)
table(resdf$strand)

# clean duplicated columns
table(resdf$prot_protein_id==resdf$protein_id)
table(resdf$prot_tx_id==resdf$tx_id)
table(resdf$prot_prot_dom_start==resdf$protein_start)
table(resdf$prot_prot_dom_end==resdf$protein_end)

resdf$prot_protein_id=resdf$prot_tx_id=resdf$prot_prot_dom_start=resdf$prot_prot_dom_end=NULL
resdf[1:10,]
colnames(resdf) = sub('^prot_','',colnames(resdf))
#saveRDS(resdf,'/nfs/cellgeni/pasham/code/sajr-java-sc/protein_features/protein_features_on_genome.rds')
resdf = readRDS('/nfs/cellgeni/pasham/code/sajr-java-sc/protein_features/protein_features_on_genome.rds')

table(resdf$protein_domain_source)
table(is.na(resdf$interpro_accession))
table(resdf$interpro_accession=='',resdf$protein_domain_source)
length(unique(resdf$interpro_accession))

# add genes
#source('/nfs/cellgeni/pasham/code/plotCoverage/plotCoverage.R')
# gtf = loadEnsGTF('/nfs/cellgeni/STAR/human/2020A/GRCh38_v32_modified.gtf')
#saveRDS(gtf,'/nfs/cellgeni/pasham/code/sajr-java-sc/2020A_unfiltered.gtf.rds')
gtf = readRDS('/nfs/cellgeni/pasham/code/sajr-java-sc/2020A_unfiltered.gtf.rds')
gtf[1:2,]
gtf = gtf[gtf$feature=='transcript',]
resdf$gene_id = setNames(gtf$gene_id,gtf$transcript_id)[resdf$tx_id]
table(resdf$seqnames,is.na(resdf$gene_id))
# so missed transcripts are just from chrs that are not here (patches)

# links with segments
segs = read.csv('/nfs/cellgeni/pasham/code/sajr-java-sc/refdata-gex-GRCh38-2020-A.segments.csv',row.names = 1)
table(resdf$seqnames %in% sub('^chr','',sub('^chrM$','MT',segs$chr_id)),is.na(resdf$gene_id))
# make per source (+all) seg2id lists
library(GenomicRanges)
segr = GRanges(sub('^chr','',segs$chr_id),IRanges(segs$start,segs$stop),strand = ifelse(segs$strand=='1','+','-'))
pfet = GRanges(resdf$seqnames,IRanges(resdf$start,resdf$end),strand = resdf$strand)
overlap = findOverlaps(segr,pfet,type='any')
overlap = data.frame(from=overlap@from,to=overlap@to)
overlap$same_gene = segs$gene_id[overlap$from] == resdf$gene_id[overlap$to]
table(overlap$same_gene)
overlap$seg_id = rownames(segs)[overlap$from]
overlap$interpro_id = resdf$interpro_accession[overlap$to]
overlap$domain_id = resdf$protein_domain_id[overlap$to]
overlap$domain_source = resdf$protein_domain_source[overlap$to]
# lets remove domain from another genes
overlap = overlap[overlap$same_gene,]

seg2ann = list()
seg2ann$interpro = lapply(split(overlap$interpro_id,overlap$seg_id),function(x)unique(x[x!='']))
seg2ann$interpro = seg2ann$interpro[sapply(seg2ann$interpro,length)>0]
for(s in unique(overlap$domain_source)){
  print(s)
  f = overlap$domain_source == s
  seg2ann[[s]] = lapply(split(overlap$domain_id[f],overlap$seg_id[f]),function(x)unique(x[x!='']))
}
seg2ann$all= lapply(split(c(overlap$interpro_id,overlap$domain_id),overlap$seg_id),function(x)unique(x[x!='']))

sapply(seg2ann,function(a)sum(sapply(a,length)==0))

#saveRDS(seg2ann,'/nfs/cellgeni/pasham/code/sajr-java-sc/protein_features/seg2domain.rds')

# domain2seg
doamin2seg = list()
f = overlap$interpro_id!=''
doamin2seg$interpro = lapply(split(overlap$seg_id[f],overlap$interpro_id[f]),function(x)unique(x[x!='']))

for(s in unique(overlap$domain_source)){
  print(s)
  f = overlap$domain_source == s
  doamin2seg[[s]] = lapply(split(overlap$seg_id[f],overlap$domain_id[f]),function(x)unique(x[x!='']))
}
sids = c(overlap$seg_id,overlap$seg_id)
f = sids != ''
doamin2seg$all= lapply(split(sids[f],c(overlap$interpro_id,overlap$domain_id)[f]),function(x)unique(x[x!='']))

sapply(doamin2seg,function(a)sum(sapply(a,length)==0))
#saveRDS(doamin2seg,'/nfs/cellgeni/pasham/code/sajr-java-sc/protein_features/domain2seg.rds')

# as data.frame
doamin2seg = list()
f = overlap$interpro_id!=''
doamin2seg$interpro = data.frame('domain'=overlap$interpro_id[f],'seg'=overlap$seg_id[f])

for(s in unique(overlap$domain_source)){
  print(s)
  f = overlap$domain_source == s
  doamin2seg[[s]] = data.frame('domain'=overlap$domain_id[f],'seg'=overlap$seg_id[f])
}

doamin2seg$all = do.call(rbind,doamin2seg)

sapply(doamin2seg,dim)
#saveRDS(doamin2seg,'/nfs/cellgeni/pasham/code/sajr-java-sc/protein_features/domain2seg.df.rds')


# load entry descriptions
interpro = read.table('https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list',sep='\t',header=TRUE,quote = "")
rownames(interpro) = interpro$ENTRY_AC
#saveRDS(interpro,'/nfs/cellgeni/pasham/code/sajr-java-sc/protein_features/interpro.rds')

dim(interpro)
interpro[1:2,]
accs = unique(resdf$interpro_accession)
table(accs %in% interpro$ENTRY_AC)
setdiff(accs , interpro$ENTRY_AC)
