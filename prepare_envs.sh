#!/usr/bin/env bash
#SBATCH --mem-per-cpu=8GB
set -eu

usage(){
    cat << HEREDOC
This script prepares all necessary environments and external libraries.

usage: `basename $0`
HEREDOC
}
[[ ! -z "${@}" ]] && usage && exit 0

here="`dirname ${0}`"
# here=/projects/enviro_lab/scripts/covid_analysis # temp
# alt
# get dirname of this script
if [[ ${0} = *"slurm"* ]]; then
COMMAND=$(scontrol show job ${SLURM_JOBID} | awk -F= '/Command=/{print $2}')
else
COMMAND="${0} ${@}"
fi
# echo "CMD: ${COMMAND}"
here=`echo "$COMMAND" | while read cmd args; do dirname $cmd; done;`
cd "$here"

# lab-specific module
module load anaconda3 &>/dev/null

# Uncomment a block below to install the environment

# To install nextflow (if java version is inadequate)
if [[ ! -d conda/env-nextflow ]]; then
echo "Creating env-nextflow"
conda create -y -p ./conda/env-nextflow python=3.9
conda activate ./conda/env-nextflow
conda install -y -c bioconda nextflow
conda deactivate
fi

# install artic
if [[ ! -d conda/env-artic ]]; then
echo "Creating env-artic"
git clone https://github.com/artic-network/fieldbioinformatics # TODO
cd fieldbioinformatics
echo "in fieldbioinformatics"
ls
conda create -c bioconda -c conda-forge -c defaults -p ../conda/env-artic -y
conda activate ../conda/env-artic
conda install -f environment.yml
# conda env create -f environment.yml -p ../conda/env-artic
python setup.py install
conda deactivate
cd ..
fi

# # install kraken2
# if [[ ! -d conda/env-kraken ]]; then
# echo "Creating env-kraken"
# conda create -y -p ./conda/env-kraken2
conda activate ./conda/env-kraken2
conda install -y -c bioconda kraken2
conda deactivate
# fi

# install porechop # TODO: uncomment stuff
if [[ ! -d conda/env-porechop ]]; then
echo "Creating env-porechop"
conda create -p ./conda/env-porechop python=3.9 -y --mkdir
conda activate ./conda/env-porechop
[[ ! -d Porechop ]] && git clone https://github.com/rrwick/Porechop.git
cd Porechop
python3 setup.py install
conda deactivate
cd ..
fi

# install homopolish
if [[ ! -d conda/env-homopolish ]]; then
echo "Creating env-homopolish"
conda config --set channel_priority flexible
conda create -y -p ./conda/env-homopolish -c conda-forge -c bioconda homopolish=0.4.1=pyhdfd78af_1
fi

# # install pangolin (this ran out of memory and couldn't finish, for me.... See Below.)
# conda create --prefix ./conda/env-pangolin
# conda activate ./conda/env-pangolin
# conda install -c bioconda -c conda-forge -c defaults pangolin
# conda deactivate

# install pangolin alternative
if [[ ! -d conda/env-pangolin ]]; then
echo "Creating env-pangolin"
[[ ! -d pangolin ]] && git clone https://github.com/cov-lineages/pangolin.git
cd pangolin
conda env create -f environment.yml -p ../conda/env-pangolin
conda activate ../conda/env-pangolin
pip install .
conda deactivate
cd ..
fi
# # get dataset - not required
# conda activate ./conda/env-pangolin
# pangolin
# conda deactivate

# install nextclade
if [[ ! -d conda/env-nextclade ]]; then
echo "Creating env-nextclade"
conda create -c bioconda -c conda-forge --yes --prefix ./conda/env-nextclade nextclade
conda activate ./conda/env-nextclade
nextclade dataset get --name 'sars-cov-2' --output-dir ./nextclade_data
conda deactivate
fi