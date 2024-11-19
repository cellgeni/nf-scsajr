# 10k Human PBMCs, 5' v2.0, Chromium X (with intronic reads)
# https://www.10xgenomics.com/datasets/human-glioblastoma-multiforme-5-v-1-targeted-neuroscience-panel-1-standard-4-0-0
mkdir 10x.GBM
cd 10x.GBM

path="https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-vdj/4.0.0/Parent_SC5v1_Human_Glioblastoma/Parent_SC5v1_Human_Glioblastoma"
wget ${path}_possorted_genome_bam.bam -O alignment.bam
wget ${path}_possorted_genome_bam.bam.bai -O alignment.bam.bai
wget ${path}_analysis.tar.gz -O clusters.tar.gz

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