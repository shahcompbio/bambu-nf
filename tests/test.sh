#!/bin/bash
env_nf=$HOME/miniforge3/envs/m1/env_nf
source $HOME/miniforge3/bin/activate ${env_nf}
## specify params
outdir=$HOME/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/SarcAtlas/APS033_bambu_nf_testing/results
pipelinedir=$HOME/VSCodeProjects/bambu-nf
mkdir -p ${outdir}
cd ${outdir}

nextflow run ${pipelinedir}/main.nf \
    -profile arm,docker,test \
    -work-dir ${outdir}/work \
    --outdir ${outdir} \
    --single_sample \
    --NDR 0.1 \
    -resume
