#!/bin/bash

## activate nf-core conda environment
source /home/preskaa/miniforge3/bin/activate nf-core

module load java/20.0.1
## specify params
outdir=/data1/shahs3/users/preskaa/SarcAtlas/data/APS033_ont_transcript_assembly/bambu_nf_test/results
pipelinedir=$HOME/bambu-nf
samplesheet=${pipelinedir}/assets/samplesheet.csv

mkdir -p ${outdir}
cd ${outdir}

nextflow run ${pipelinedir}/main.nf \
    -profile singularity \
    -work-dir ${outdir}/work \
    --outdir ${outdir} \
    --input ${samplesheet} \
    -resume