# Human Glioblastoma Multiforme: 5’v1 Targeted, Neuroscience Panel
# https://www.10xgenomics.com/datasets/human-glioblastoma-multiforme-5-v-1-targeted-neuroscience-panel-1-standard-4-0-0
mkdir 10x.GBM
cd 10x.GBM

path="https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/Targeted_SC5v1_Human_Glioblastoma_Neuroscience/Targeted_SC5v1_Human_Glioblastoma_Neuroscience"
curl ${path}_filtered_feature_bc_matrix.h5 -o filtered_feature_bc_matrix.h5
curl ${path}_possorted_genome_bam.bam -o alignment.bam
curl ${path}_possorted_genome_bam.bam.bai -o alignment.bam.bai
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