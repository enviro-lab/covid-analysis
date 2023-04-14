#!/usr/bin/env bash
set -eu

data_type=$1
modify_file=last_${data_type}_update.time
one_week_ago=`date +%s --date='-1 week'`

doUpdate=false

if [[ -f $modify_file ]]; then
    last_modified=`stat -c %y`
    [[ $last_modified -lt $one_week_ago ]] && doUpdate=true
else doUpdate=true
fi

# fail if conda env not active
which $data_type | grep -q "conda/env-$data_type"
if [[ $? -ne 0 ]]; then
# if [[ ! `conda env list | grep env-$data_type | grep -q '*'` ]]; then
    echo "conda env 'conda/env-$data_type' must be activated when running this script"
    exit 1
fi

if [[ $doUpdate = false ]]; then exit 0; fi

# do the update
if [[ $data_type = pangolin ]]; then
    pangolin --update
    # # if the above doesn't work, these is the recommended steps for updating the development version
    # ( cd pangolin
    # git pull
    # conda env update -f environment.yml
    # pip install .
    # pangolin --update-data
    # )

elif [[ $data_type = nextclade ]]; then
    conda update nextclade
    nextclade dataset get --name=sars-cov-2 --output-dir='nextclade_data'
fi

# write update file
printf '' > $modify_file