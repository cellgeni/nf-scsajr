# nf-scsajr

Pipeline to quantify alternative splicing in single cell data. It uses segments - part of gene between two adjacent splicing sites - as a feature. For each segment it quantify number of inclusion (confirming inclusion of the segment into RNA) and exclusion UMIs (confirming exclusion of the segment from RNA) in each barcode. Resultant matrices are extremely sparse, so it calcualtes per sample*celltype pseudobulks and uses GLM with quasibinomial distribution to look for differences between celltype. Afterwards it performs GO and interpro domain enrichment analyses and generates html output, coverage plot for up to 100 best events and rds files with counts and pv.


It uses [java read counter](https://github.com/iaaka/sajr-java-sc) for read counting and [sajr](https://github.com/iaaka/sajr) for statistical analyses of alternative. splicing


# Input data
Pipeline requires bam files and celltype annotation. Bam files are provided as tsv file with two columns and no header:
```
sample1 /path/to/bam1
sample2 /path/to/bam2
```
Celltype annotation specified by another tsv file:
```
sample barcode celltype
```

# Run
```
nextflow run main.nf \
 --SAMPLEFILE samples.tsv \
 --BARCODEFILE barcodes.tsv \
 --bam_on_irods true \
 -resume
```
For now pipeline works only with human 2020A reference. It can be relatively easily applied for other species/annotation withou GO/domain enrichment part. However it is not implemented yet, please contact me if you want to run pipeline on any other reference.

The pipeline relies on junction reads. So, while technically it can work on any bam with cell barcode (CB bam tag) set, it will hardly detect anything in 3' data while single-nuclei data can be very noisy. It seems to work reasonably well on single cell 5' short reads or on 3' long reads.  

# TODO
0. Change to SummarizedExperiment to store splicing data
1. Make workflow for reference creation
2. Make code work without go and/or interpro (other species)
3. Autodetect strand using part of data (not whole as now)
4. Autodetect chr/not-chr annotation?
5. Fail smartly if filtering leaves no segments
6. Speed up quantification (rewrite it entirely?). Make per-gene? paralelaze?
