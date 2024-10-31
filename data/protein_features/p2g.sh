#! /bin/bash -e
#BSUB -J prots[1-73]
#BSUB -o logs/%J.%I.prot.log
#BSUB -e logs/%J.%I.prot.err
#BSUB -n 1
#BSUB -M3000
#BSUB -R "span[hosts=1] select[mem>3000] rusage[mem=3000]"

cd /nfs/cellgeni/pasham/code/sajr-java-sc/protein_features/tmp/

Rscript ../p2g.R ${LSB_JOBINDEX}
