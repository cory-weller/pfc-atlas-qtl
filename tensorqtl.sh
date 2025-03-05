#!/bin/bash
module load singularity
. /usr/local/current/singularity/app_conf/sing_binds

export SINGULARITY_BINDPATH="/usr/local/apps,/usr/local/CUDA,/usr/local/cuDNN,/usr/local/Qt,$SINGULARITY_BINDPATH"

export SINGULARITYENV_PYTHONNOUSERSITE=1
export SINGULARITYENV_PYTHONPATH="/usr/local/lib/python3.10/dist-packages:$SINGULARITYENV_PYTHONPATH"
export SINGULARITYENV_R_LIBS_SITE="/usr/local/lib/R/site-library" #:/opt/conda/envs/tensorqtl/lib/R/library"
my_command="shell"
if [[ $# > 0 ]]; then
    my_command="exec"
fi

selfdir="$(dirname $(readlink -f ${BASH_SOURCE[0]}))"

singularity --silent "$my_command" --home '/data/CARD_singlecell/users/wellerca/pfc-atlas-qtl' --nv /usr/local/apps/tensorqtl/1.0.9/libexec/TensorQTL-1.0.9_from_docker.sif "$@"
