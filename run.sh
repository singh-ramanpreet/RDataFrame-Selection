#!/bin/bash

redhat_version=$(lsb_release -rs | cut -f1 -d.)
if [[ $(echo "${redhat_version} == 6" | bc) -eq 1 ]]; then
    echo "Running on slc6 machine"
    source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-slc6-gcc8-opt/setup.sh

elif [[ $(echo "${redhat_version} == 7" | bc) -eq 1 ]]; then
    echo "Running on centos7 machine"
    source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-centos7-gcc8-opt/setup.sh

else
    echo "Run on SLC 6 or Cent OS 7"
    exit
fi

sample_list_txt=${1}
isMC=${2}
year=${3}
root -q -n -b vbs_flat_ntupler.cc+\(\"${sample_list_txt}\",${isMC},${year}\)
