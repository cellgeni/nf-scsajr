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
By default pipeline uses 2020-A human reference, if you want to use another reference you'll need first convert gtf file by `java -jar sajr.ss.jar gff2sajr sajr.config -ann_foreign=input.gtf -ann_out=output.sajr` and then add `--ref output.sajr` to pipeline command.

# TODO
1. Make workflow for reference creation
2. Make code work without go and/or interpro (other psecies)
3. Autodetect strand using part of data (not whole as now)
4. Autodetect chr/not-chr annotation?
5. Fail smartly if filtering leaves no segments
6. Speed up quantification (rewrite it entirely?). Make per-gene? paralelaze?