# Human Glioblastoma Multiforme: 5’v1 Whole Transcriptome Analysis
# https://www.10xgenomics.com/datasets/human-glioblastoma-multiforme-5-v-1-whole-transcriptome-analysis-1-standard-4-0-0
mkdir 10x.GBM
cd 10x.GBM

url="https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-vdj/4.0.0/Parent_SC5v1_Human_Glioblastoma/Parent_SC5v1_Human_Glioblastoma"
curl ${url}_raw_feature_bc_matrix.h5 -o raw_feature_bc_matrix.h5
curl ${url}_filtered_feature_bc_matrix.h5 -o filtered_feature_bc_matrix.h5
curl ${url}_possorted_genome_bam.bam -o alignment.bam
curl ${url}_possorted_genome_bam.bam.bai -o alignment.bam.bai
curl ${url}_analysis.tar.gz -o clusters.tar.gz

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
 --ref ./nf-scsajr/ref/human_2020A_chr

# it should generate output similar to one in "output" dir (plus rds files)
