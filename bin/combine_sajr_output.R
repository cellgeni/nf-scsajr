options(error = function(e) quit("no", 1))

library(Matrix)
library(doMC)
library(plyr)
library(scsajr)
library(visutils)


## Parse arguments
# Rscript combine_sajr_output.R sajr_outs_file, barcodes_file, path2ref (folder), path2bin, ncores
args <- commandArgs(trailingOnly = TRUE)
writeLines(args, "params.txt") # Store parsed parameters

# args = readLines('params.txt')
path2sajr_outs <- args[1]
path2barcodes <- args[2]
path2ref <- args[3]
path2bin <- args[4]
ncores <- as.integer(args[5])
delimiter <- "$" # delimeter to separate sample and celltype


## Import utility functions
doMC::registerDoMC(ncores)


## Load inputs
sajr_outs <- utils::read.table(path2sajr_outs, sep = " ")
colnames(sajr_outs) <- c("sample_id", "chr", "sajr_out", "bam_path", "strand")

barcodes <- utils::read.table(path2barcodes,
  sep = "\t",
  col.names = c("sample_id", "barcode", "celltype")
)
barcodes$celltype[is.na(barcodes$celltype)] <- "NA"
rownames(barcodes) <- paste0(barcodes$sample_id, "|", barcodes$barcode)

# Check for unsupported characters (val: delimiter)
if (any(grepl(delimiter, c(barcodes$sample_id, barcodes$celltype), fixed = TRUE))) {
  stop(paste0("Celltype names or sample IDs contains '", delimiter, "' that is not supported, please get rid of it."))
}

segments_path <- file.path(path2ref, "segments.csv")
seg <- utils::read.csv(segments_path, row.names = 1)


## Initialize output directory
out_dir <- "rds"
dir.create(out_dir)
scsajr::log_info("initializing output directory finished")


## Load SAJR outputs and make pseudobulks per sample in parallel
samples <- unique(sajr_outs$sample_id)

pbasl <- llply(seq_along(samples), function(i) {
  # For each sample, extract SAJR output lists
  sample <- samples[i]
  sample_sajr_outs <- subset(sajr_outs, sample_id == sample)

  # Load exon,intron matrices for each chromosome of the sample
  r_list <- lapply(seq_len(nrow(sample_sajr_outs)), function(j) {
    chr_prefix <- paste0(sample_sajr_outs$sajr_out[j], "/", sample_sajr_outs$chr[j])

    # Returns a list with $e (exon matrix) and $i (intron matrix)
    scsajr::load_sc_as(
      seg[seg$chr_id == sample_sajr_outs$chr[j], ],
      chr_prefix
    )
  })

  # Identify the union of all cell barcodes across the chromosomes
  cl <- unique(unlist(lapply(r_list, function(x) colnames(x$i))))

  # Row-bind all intron matrices and exon matrices into a single matrix
  i_mat <- scsajr::rbind_matrix(lapply(r_list, function(x) x$i), cl)
  e_mat <- scsajr::rbind_matrix(lapply(r_list, function(x) x$e), cl)

  # Ensure rows follow the master 'seg' ordering
  i_mat <- i_mat[rownames(seg), , drop = FALSE]
  e_mat <- e_mat[rownames(seg), , drop = FALSE]

  # Prefix each column with "sample|cell"
  colnames(e_mat) <- colnames(i_mat) <- paste0(sample, "|", colnames(i_mat))

  # Save matrices to RDS files
  saveRDS(
    list(i = i_mat, e = e_mat),
    file.path("rds", paste0(sample, ".rds"))
  )


  # Load intron only matrices per chromosome
  intron_list <- lapply(seq_len(nrow(sample_sajr_outs)), function(j) {
    scsajr::read_named_mm(paste0(sample_sajr_outs$sajr_out[j], "/", sample_sajr_outs$chr[j], ".intron"))
  })
  intron_mat <- scsajr::rbind_matrix(intron_list)
  colnames(intron_mat) <- paste0(sample, "|", colnames(intron_mat))
  saveRDS(
    intron_mat,
    file.path("rds", paste0(sample, "intron.rds"))
  )


  # Determine which cells are present in the sample
  cmn <- intersect(rownames(barcodes), colnames(i_mat))
  if (length(cmn) == 0) {
    # no overlapping cells → return an empty result for this sample
    return(list(i = NULL, e = NULL, cmn = cmn, ncell = 0))
  }
  # Filter matrices to only include common cells
  sub_i <- i_mat[, cmn]
  sub_e <- e_mat[, cmn]

  # For each segment, count how many cells have ≥10 total reads
  ncell <- rowSums(sub_i + sub_e > 9)

  # Build the grouping factor "sample|celltype" for pseudobulk
  f <- paste0(
    barcodes[cmn, "sample_id"],
    delimiter,
    barcodes[cmn, "celltype"]
  )

  # Sum counts within each (sample, celltype) group
  pb_i <- visutils::calcColSums(sub_i, f)
  pb_e <- visutils::calcColSums(sub_e, f)

  # Create a pseudobulk object
  pb <- list(
    i     = pb_i,
    e     = pb_e,
    cmn   = cmn,
    ncell = ncell
  )

  # Clean up the loaded data to free memory
  rm(i_mat, e_mat, sub_i, sub_e, r_list, introns)
  mem <- sum(gc()[, 2]) / 2^10 # MB used
  pbmem <- object.size(pb) / 2^20 # MB of pseudobulk
  scsajr::log_info(sample_sajr_outs$sample_id[i], " is loaded. Total mem used: ", mem, " Gb. Pb size: ", pbmem, " Gb")

  # Return the pseudobulk object
  pb
}, .parallel = TRUE)


## Combine pseudobulks from all samples into a single list
# Collect common barcodes and their counts across all samples
cmnbarcodesl <- list()
ncell <- rep(0, nrow(seg))

for (i in seq_along(pbasl)) {
  cmnbarcodesl[[i]] <- pbasl[[i]]$cmn
  ncell <- ncell + pbasl[[i]]$ncell
  pbasl[[i]]$cmn <- pbasl[[i]]$ncell <- NULL
}

gc()
scsajr::log_info("data loaded")

names(pbasl) <- samples
cmnbarcodes <- unlist(cmnbarcodesl)

# Assemble final segment annotation
pbas <- list(seg = seg, i = NULL, e = NULL)
for (n in names(pbasl)) {
  pbas$i <- cbind(pbas$i, pbasl[[n]]$i)
  pbas$e <- cbind(pbas$e, pbasl[[n]]$e)
  pbasl[[n]] <- NULL
}

# Add total cell count to annotations
pbas$seg$ncell <- ncell


## Build pseudobulk metadata
barcodes <- barcodes[cmnbarcodes, ]
pbmeta <- as.data.frame(do.call(rbind, strsplit(colnames(pbas$i), delimiter, fixed = TRUE)))
colnames(pbmeta) <- c("sample_id", "celltype")
rownames(pbmeta) <- colnames(pbas$i)

# Number of cells in each pseudobulk
pbmeta$ncells <- as.numeric(
  table(paste0(barcodes$sample_id, delimiter, barcodes$celltype))
)[rownames(pbmeta)]

# Add strand info from SAJR outputs
s2s <- unique(sajr_outs[, c("sample_id", "strand")])
rownames(s2s) <- s2s$sample_id
pbmeta$strand <- s2s[pbmeta$sample_id, "strand"]

## Wrap into a SummarizedExperiment object and save
pbas_se <- scsajr::make_summarized_experiment(pbas, pbmeta)
saveRDS(barcodes, paste0(out_dir, "/cell_meta.rds"))
saveRDS(pbas_se, paste0(out_dir, "/pbas.rds"))
scsajr::log_info("pseudobulk is made")
