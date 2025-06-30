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
segments <- SAJR::loadSAData(sajr_path)

# Annotate which segments are coding exons
segments <- scsajr::add_is_coding_by_ens_gtf(gtf_path = gtf_path, segments)

# Assign segment types (e.g. ALT, EXN, INT)
segments <- SAJR::setSplSiteTypes(segments, sajr_path)
segments <- segments$seg

# Compute segment lengths
segments$length <- segments$stop - segments$start + 1

# Load original GTF file
gtf <- plotCoverage::loadEnsGTF(gtf_path)
if (is.null(gtf$transcript_name)) {
  gtf$transcript_name <- gtf$transcript_id
}

# Flag which segments are identical to exons
# Compares gene, start, and stop
ef <- gtf$feature == "exon"
segments$is_exon <- paste(segments$gene_id, segments$start, segments$stop) %in%
  paste(gtf$gene_id[ef], gtf$start[ef], gtf$stop[ef])

# Build gene description data frame
gene_descr <- SAJR::loadGData(sajr_path)$gene
gene_descr$descr <- rownames(gene_descr)
gene_descr$name <- rownames(gene_descr)

# Save the processed data
utils::write.csv(segments, "segments.csv")
saveRDS(gtf, "gtf.rds")
dir.create("functional_annotation")
saveRDS(gene_descr, "functional_annotation/gene.descr.rds")
