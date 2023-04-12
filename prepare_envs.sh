#!/usr/bin/env bash
usage(){
    cat << HEREDOC
This script prepares all necessary environments and external libraries.

usage: `basename $0`
HEREDOC
}
[[ ! -z "${@}" ]] && usage && exit 0

here="`dirname ${0}`"
here=/projects/enviro_lab/scripts/covid_analysis # temp
cd "$here"

# lab-specific module
module load anaconda3

# Uncomment a block below to install the environment

# # To install nextflow (if java version is inadequate)
# conda create -y -p ./conda/env-nextflow python=3.9
# conda activate ./conda/env-nextflow
# conda install -y -c bioconda nextflow
# conda deactivate

# # install artic
# (
# git clone https://github.com/artic-network/fieldbioinformatics
# cd fieldbioinformatics
# conda env create -f environment.yml -p ../conda/env-artic
# conda activate ../conda/env-artic
# python setup.py install
# conda deactivate
# )

# # install kraken2
# conda create -y -p ./conda/env-kraken2
# conda activate ./conda/env-kraken2
# conda install -y -c bioconda kraken2
# conda deactivate

# # install porechop
# (
# conda create -p ./conda/env-porechop python=3.9 -y --mkdir
# git clone https://github.com/rrwick/Porechop.git
# cd Porechop
# python3 setup.py install
# conda deactivate
# )

# # install homopolish
# conda config --set channel_priority flexible
# conda create -y -p ./conda/env-homopolish -c conda-forge -c bioconda homopolish=0.4.1=pyhdfd78af_1

# # # install pangolin (this ran out of memory and couldn't finish, for me.... See Below.)
# # conda create --prefix ./conda/env-pangolin
# # conda activate ./conda/env-pangolin
# # conda install -c bioconda -c conda-forge -c defaults pangolin
# # conda deactivate

# # install pangolin alternative
# ( git clone https://github.com/cov-lineages/pangolin.git
# cd pangolin
# conda env create -f environment.yml -p ../conda/env-pangolin
# conda activate ../conda/env-pangolin
# pip install .
# conda deactivate
# )
# # get dataset - not required
# conda activate ./conda/env-pangolin
# pangolin
# conda deactivate

# # install nextclade
# conda create -c bioconda -c conda-forge --yes --prefix ./conda/env-nextclade nextclade
# conda activate ./conda/env-nextclade
# nextclade dataset get --name 'sars-cov-2' --output-dir ./nextclade_data
