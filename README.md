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
# Create reference
nf-scSAJR is distributed with pre-build human 2020A reference. It includes gene descriptions and protein domain annotation.
If you are working with non-human species or if you want to use other annotation version you can build th reference from gtf file using:
```
nextflow main.nf \
 -entry reference \
 --gtf annotation.gtf \
 --outdir <path2reference> \
 -resume
```
Please keep in mind that this reference will not have domain annotation and gene descriptions so there will be no domain enrichment analyses in pipeline output obtained using this such reference. 

# Run
```
nextflow main.nf \
 --SAMPLEFILE samples.tsv \
 --BARCODEFILE barcodes.tsv \
 --outdir sajr_out \
 --ref ref/human_2020A_chr \
 -resume
```

The pipeline relies on junction reads. So, while technically it can work on any bam with cell barcode (CB bam tag) set, it will hardly detect anything in 3' data while single-nuclei data can be very noisy. It seems to work reasonably well on single cell 5' short reads or on 3' long reads.  

# TODO
1. Autodetect strand using part of data (not whole as now)
2. Autodetect chr/not-chr annotation?
3. Fail smartly if filtering leaves no segments
4. Speed up quantification (rewrite it entirely?). Make per-gene? paralelaze?
