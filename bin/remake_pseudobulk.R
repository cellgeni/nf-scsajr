options(error = function(e) quit("no", 1))

library(Matrix)
library(SAJR)
library(clusterProfiler)
library(doMC)
library(org.Hs.eg.db)
library(plyr)
library(scsajr)
library(visutils)

## Parse arguments
# Rscript remake_pseudobulk.R samples_file, barcodes_file, rds, path2ref (folder), path2bin, ncores
args <- commandArgs(trailingOnly = TRUE)
writeLines(args, "params.txt")

# args = readLines('params.txt')
path2samples <- args[1]
path2barcodes <- args[2]
path2rds <- args[3]
path2ref <- args[4]
path2bin <- args[5]
ncores <- as.integer(args[6])
delimiter <- "$" # delimeter to separate sample and celltype names

doMC::registerDoMC(ncores)


## Load inputs
samples <- utils::read.table(path2samples, sep = "\t", col.names = c("sample_id", "bam_path"))
barcodes <- utils::read.table(
  path2barcodes,
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


## Load SAJR
pbasl <- plyr::llply(seq_len(nrow(samples)), function(i) {
  sample <- samples$sample_id[i]

  rds_file <- list.files(path2rds, pattern = sample, full.names = TRUE)
  if (length(rds_file) != 1) {
    stop(paste0("Sample '", sample, "' has more than one rds!"))
  }

  r <- readRDS(rds_file)
  i_mat <- r$i
  e_mat <- r$e

  # calc pseudobulks
  cmn <- intersect(rownames(barcodes), colnames(i_mat))
  sub_i <- i_mat[, cmn]
  sub_e <- e_mat[, cmn]

  ncell <- rowSums(sub_i + sub_e > 9)

  f <- paste0(
    barcodes[cmn, "sample_id"],
    delimiter,
    barcodes[cmn, "celltype"]
  )

  pb_i <- visutils::calcColSums(sub_i, f)
  pb_e <- visutils::calcColSums(sub_e, f)

  pb <- list(
    i     = pb_i,
    e     = pb_e,
    cmn   = cmn,
    ncell = ncell
  )

  rm(i_mat, e_mat, sub_i, sub_e, pb_i, pb_e)
  mem <- gc()
  mem <- round(sum(mem[, 2]) / 2^10, digits = 2)
  pbmem <- round(unclass(object.size(pb)) / 2^30, digits = 2)
  scsajr::log_info(sample, " is loaded. Total mem used: ", mem, " Gb. Pb size: ", pbmem, " Gb")
  pb
}, .parallel = TRUE)

cmnbarcodesl <- list()
ncell <- rep(0, nrow(seg))

for (i in seq_along(pbasl)) {
  cmnbarcodesl[[i]] <- pbasl[[i]]$cmn
  ncell <- ncell + pbasl[[i]]$ncell
  pbasl[[i]]$cmn <- pbasl[[i]]$ncell <- NULL
}
gc()
scsajr::log_info("data loaded")

# Combine
names(pbasl) <- samples$sample_id
cmnbarcodes <- unlist(cmnbarcodesl)

pbas <- list(seg = seg, i = NULL, e = NULL)
for (n in names(pbasl)) {
  pbas$i <- cbind(pbas$i, pbasl[[n]]$i)
  pbas$e <- cbind(pbas$e, pbasl[[n]]$e)
  pbasl[[n]] <- NULL
}
pbas$seg$ncell <- ncell
class(pbas) <- c("sajr", "list")

# all(colnames(pbas$i)== colnames(pbas$e))

barcodes <- barcodes[cmnbarcodes, ]

pbmeta <- as.data.frame(do.call(rbind, strsplit(colnames(pbas$i), delimiter, fixed = TRUE)))
colnames(pbmeta) <- c("sample_id", "celltype")
rownames(pbmeta) <- colnames(pbas$i)
pbmeta$ncells <- as.numeric(
  table(paste0(barcodes$sample_id, delimiter, barcodes$celltype))[rownames(pbmeta)]
)

pbas_se <- scsajr::make_summarized_experiment(pbas, pbmeta)
scsajr::log_info("pseudobulk is made")

# save ##############
saveRDS(barcodes, paste0(out_dir, "/cell_meta.rds"))
saveRDS(pbas_se, paste0(out_dir, "/pbas.rds"))
