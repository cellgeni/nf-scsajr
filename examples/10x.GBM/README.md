Here is a code to run nf-scSAJR on [public 10x GBM data](https://www.10xgenomics.com/datasets/human-glioblastoma-multiforme-5-v-1-targeted-neuroscience-panel-1-standard-4-0-0).

`alignment.bam`, `alignment.bam.bai`, and `clusters.tar.gz` are placeholders in `input` directory.

Please prepare the data as written in `run.sh`

Expected outputs are in `output` directory:
1. [summary.html](https://html-preview.github.io/?url=https://github.com/cellgeni/nf-scsajr/blob/main/examples/10x.GBM/output/summary.html).
2. [coverage plots for top examples](output/examples)

Caution: `sample.rds` and `sampleintron.rds` are not included in the repo due to file size. 

Please run `run.sh` you run `10x_gbm.Rmd`.
