# 10k Human PBMCs, 5' v2.0, Chromium X (with intronic reads)
# https://www.10xgenomics.com/datasets/10k-human-pbmcs-5-v2-0-chromium-x-with-intronic-reads-2-standard
mkdir 10x.PBMC
cd 10x.PBMC

path="https://cf.10xgenomics.com/samples/cell-vdj/6.1.2/10k_PBMC_5pv2_nextgem_Chromium_X_intron_10k_PBMC_5pv2_nextgem_Chromium_X_intron/10k_PBMC_5pv2_nextgem_Chromium_X_intron_10k_PBMC_5pv2_nextgem_Chromium_X_intron_count"
curl ${path}_filtered_feature_bc_matrix.h5 -o filtered_feature_bc_matrix.h5
curl ${path}_sample_alignments.bam -o alignment.bam
curl ${path}_sample_alignments.bam.bai -o alignment.bam.bai
curl ${path}_analysis.tar.gz -o clusters.tar.gz

# prepare input files
tar -xzf clusters.tar.gz
echo -e "sample\t$( pwd )/alignment.bam" > samples.tsv
awk -F, 'BEGIN { OFS="\t" }; NR>1 {print  "sample",$1,$2}' analysis/clustering/graphclust/clusters.csv > barcodes.tsv

# clone repo
git clone https://github.com/cellgeni/nf-scsajr

# run pipeline
nextflow run nf-scsajr \
 --SAMPLEFILE samples.tsv \
 --BARCODEFILE barcodes.tsv \
 --minsamples 1 \
 -resume

# it should generate output similar to one in "output" dir (plus rds files)