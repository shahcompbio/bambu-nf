#!/bin/bash

## activate nf-core conda environment
source /home/preskaa/miniforge3/bin/activate nf-core

module load java/20.0.1
## specify params
outdir=/data1/shahs3/users/preskaa/SarcAtlas/data/APS033_ont_transcript_assembly/bambu_nf_test/results
pipelinedir=$HOME/bambu-nf
samplesheet=${pipelinedir}/assets/samplesheet.csv
fasta=/data1/shahs3/isabl_data_lake/assemblies/GRCh38-P14/GRCh38.primary_assembly.genome.fa
gtf=/data1/shahs3/isabl_data_lake/assemblies/GRCh38-P14/gencode.v45.primary_assembly.annotation.gtf

mkdir -p ${outdir}
cd ${outdir}

nextflow run ${pipelinedir}/main.nf \
    -profile singularity \
    -work-dir ${outdir}/work \
    --outdir ${outdir} \
    --input ${samplesheet} \
    --filter_acc_reads \
    --fasta ${fasta} \
    --gtf ${gtf} \
    -resume