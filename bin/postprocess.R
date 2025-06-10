options(error = function(e) quit("no", 1))

library(Matrix)
library(SAJR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(plotCoverage)
library(scsajr)
library(visutils)


## Parse arguments
# Rscript postprocess.R pbsas, mincells, minsamples, path2samples, path2barcodes, path2ref, path2bin, ncores
args <- commandArgs(trailingOnly = TRUE)
writeLines(args, "params.txt")

# args = readLines('params.txt')
out_dir <- args[1] # output directory???
mincells <- as.integer(args[2])
minsamples <- as.integer(args[3])
path2samples <- args[4]
path2barcodes <- args[5]
path2ref <- args[6]
path2bin <- args[7]
ncores <- as.numeric(args[8])


## Import utility functions
doMC::registerDoMC(ncores)


## Load inputs
pbas_all <- readRDS(file.path(out_dir, "pbas.rds"))
samples <- read.table(path2samples, sep = " ", col.names = c("sample_id", "bam_path"))
barcodes <- read.table(path2barcodes, sep = "\t", col.names = c("sample_id", "barcode", "celltype"))

gtf <- readRDS(file.path(path2ref, "gtf.rds"))
gene.descr <- readRDS(file.path(path2ref, "functional_annotation/gene.descr.rds"))

# optional domain annotation
domain2seg <- readRDSifExists(file.path(path2ref, "functional_annotation/domain2seg.df.rds"))
domain_descr <- readRDSifExists(file.path(path2ref, "functional_annotation/all_domain_descr.rds"))


## Filtering
# Removes pseudbulk groups with fewer than `minsamples` samples
#  and segments with fewer than `mincells` cells
pbas <- scsajr::filter_segments_and_samples(
  pbas_all,
  celltype_min_samples = minsamples,
  sample_min_ncells = mincells
)
scsajr::log_info("samples*celltype filtering: ", ncol(pbas_all), " -> ", ncol(pbas))
scsajr::log_info("segment filtering: ", nrow(pbas_all), " -> ", nrow(pbas))


## Differential-splicing tests
# Test all celltypes together
pbas@metadata$all_celltype_test <- scsajr::test_all_groups_as(pbas, "celltype", .parallel = TRUE)

# Find marker segments per celltype
pbas@metadata$markers <- scsajr::find_marker_as(pbas, "celltype", .parallel = TRUE, verbose = TRUE)

scsajr::log_info("diff AS finished")


## Gene Ontology (GO) enrichment analysis
perct <- pbas@metadata$markers
pv <- pbas@metadata$all_celltype_test

# Build a universe of all genes
gene_uni <- unique(rowData(
  scsajr::filter_segments_and_samples(
    pbas_all,
    seg_min_sd = 0,
    celltype_min_samples = minsamples,
    sample_min_ncells = mincells
  )
)$gene_id)

segs <- SummarizedExperiment::as.data.frame(pbas_all@rowRanges)

# For each celltype, select genes with |ΔPSI|>0.2 & FDR<0.05
gids <- apply(perct$fdr < 0.05 & abs(perct$dpsi) > 0.2, 2, function(f) {
  f[is.na(f)] <- FALSE
  unique(segs[rownames(perct$fdr)[f], "gene_id"])
})

# Add an “all” category from the global test (dpsi>0.3 & FDR<0.05)
sgn <- pv$group_fdr < 0.05 & pv$dpsi > 0.3
sgn[is.na(sgn)] <- FALSE
gids$all <- unique(segs[rownames(pv)[sgn], "gene_id"])

# Run GO enrichment analysis
pbas@metadata$go <- tryCatch(
  {
    clusterProfiler::compareCluster(
      gids,
      fun = "enrichGO",
      universe = gene_uni,
      OrgDb = "org.Hs.eg.db",
      keyType = "ENSEMBL",
      pAdjustMethod = "BH",
      ont = "ALL"
    )
  },
  error = function(e) {
    warning("GO analyses was unsuccesfull")
    print(e)
  }
)
scsajr::log_info("GO finished")


## InterPro domain enrichment analysis
# Use only cassette exons
if (!is.null(domain2seg)) {
  # take only domains with description
  domain2seg <- domain2seg$all
  domain2seg <- domain2seg[domain2seg$domain %in% domain_descr$ENTRY_AC, ]

  domain_sites2use <- "ad"
  seg_uni <- rownames(pbas_all)[segs$is_exon & segs$sites %in% domain_sites2use & segs$gene_id %in% gene_uni]

  sids <- apply(perct$fdr < 0.05 & abs(perct$dpsi) > 0.2, 2, function(f) {
    f[is.na(f)] <- FALSE
    out <- rownames(perct$fdr)[f]
    out[segs[out, "sites"] == "ad" & segs[out, "is_exon"]]
  })
  sgn <- pv$group_fdr < 0.05 & pv$dpsi > 0.3
  sgn[is.na(sgn)] <- FALSE
  sids$all <- rownames(pv)[sgn]
  sids$all <- sids$all[segs[sids$all, "sites"] == "ad" & segs[sids$all, "is_exon"]]

  pbas@metadata$ipro <- clusterProfiler::compareCluster(
    sids,
    fun = "enricher",
    universe = seg_uni,
    TERM2GENE = domain2seg,
    TERM2NAME = domain_descr[, c("ENTRY_AC", "ENTRY_NAME")],
    pAdjustMethod = "BH"
  )
}


## Save filtered pseudobulk object
saveRDS(pbas, file.path(out_dir, "pb_as_filtered.rds"))
scsajr::log_info("interpro finished")


## Example coverage plots
# pbas = readRDS(paste0(out.dir,'/pb_as_filtered.rds'))
markers <- scsajr::select_all_markers(
  pbas@metadata$markers,
  pbas@metadata$all_celltype_test,
  dpsi_thr = 0.2,
  n = Inf
)

# Select top 200 by |dpsi| and p-value
N <- min(200, nrow(markers), max(sum(abs(markers$dpsi) > 0.5), 100))
markers <- markers[order(abs(markers$dpsi), decreasing = TRUE)[seq_len(N)], ]
# an attempt to sort examples by both pv and dpsi
markers <- markers[order(-log10(markers$pv) / 100 + abs(markers$dpsi), decreasing = TRUE), ]

# to save RAM we'll keep only data we need for plotting
pbas_mar <- pbas_all[markers$seg_id, ]
pbas_all$all <- "all"
pb_all <- scsajr::pseudobulk(pbas_all, "all")
gc()


# For each marker segment, plot coverage
dir.create(paste0(out_dir, "/examples_coverage"))
dir.create("examples")

l_ply(seq_along(markers$seg_id), function(i) {
  sid <- markers$seg_id[i]
  gid <- segs[sid, "gene_id"]

  # get plot coords,a dn previous coverage if exists
  coors <- scsajr::get_plot_coords_for_seg(sid, pb_all, gene.descr)
  pdf(paste0("examples/", substr(10000 + i, 2, 100), "_", gene.descr[gid, "name"], "_", sid, ".pdf"), w = 12, h = 12)

  covs <- NULL
  rdsf <- paste0(out_dir, "/examples_coverage/", sid, ".rds")
  if (file.exists(rdsf)) {
    covs <- readRDS(rdsf)
  }

  covs <- scsajr::plot_segment_coverage(
    chr = segs[sid, "seqnames"],
    start = coors$start, stop = coors$stop,
    covs = covs,
    sid = sid,
    data_as = pbas_mar,
    groupby = "celltype",
    celltypes = NULL,
    barcodes = barcodes,
    samples = samples,
    gene.descr = gene.descr,
    plot.junc.only.within = NA,
    min.junc.cov.f = 0.02,
    min.junc.cov = 3,
    ylim_by_junc = TRUE,
    gtf = gtf,
    oma = c(6, 14, 3, 1)
  )
  dev.off()
  if (!file.exists(rdsf)) {
    saveRDS(covs, rdsf)
  }
}, .parallel = TRUE)

scsajr::log_info("examples finished")
