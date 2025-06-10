options(error = function(e) quit("no", 1))

# devtools::install_github("iaaka/sajr")
# devtools::install_github("cellgeni/scsajr")
library(SAJR)
library(plotCoverage)
library(scsajr)

# Read command-line args: GTF path and SAJR segments file (converted annotation)
args <- commandArgs(trailingOnly = TRUE)
gtf_path <- args[1]
sajr_path <- args[2]

# Load SAJRsegments object
seg <- SAJR::loadSAData(sajr_path)

# Annotate which segments are coding exons
segs <- scsajr::add_is_coding_by_ens_gtf(gtf_path = gtf_path, segs) # input for func: segs or seg?

# Assign segment types (e.g. ALT, EXN, INT)
seg <- SAJR::setSplSiteTypes(seg, sajr_path)$seg

# Compute segment lengths
seg$length <- seg$stop - seg$start + 1

# Load original GTF file
gtf <- plotCoverage::loadEnsGTF(gtf_path)
if (is.null(gtf$transcript_name)) {
  gtf$transcript_name <- gtf$transcript_id
}

# Flag which segments are identical to exons
# Compares gene, start, and stop
ef <- gtf$feature == "exon"
seg$is_exon <- paste(seg$gene_id, seg$start, seg$stop) %in%
  paste(gtf$gene_id[ef], gtf$start[ef], gtf$stop[ef])

# Build gene description data frame
gene_descr <- SAJR::loadGData(sajr_path)$gene
gene_descr$name <- rownames(gene_descr)
gene_descr$descr <- rownames(gene_descr)

# Save the processed data
utils::write.csv(seg, "segments.csv")
saveRDS(gtf, "gtf.rds")
dir.create("functional_annotation")
saveRDS(gene_descr, "functional_annotation/gene.descr.rds")
