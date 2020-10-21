#!/bin/bash

redhat_version=$(lsb_release -rs | cut -f1 -d.)
if [[ $(echo "${redhat_version} == 6" | bc) -eq 1 ]]; then
    echo "Running on slc6 machine"
    source /cvmfs/sft.cern.ch/lcg/views/LCG_98python3/x86_64-slc6-gcc8-opt/setup.sh

elif [[ $(echo "${redhat_version} == 7" | bc) -eq 1 ]]; then
    echo "Running on centos7 machine"
    source /cvmfs/sft.cern.ch/lcg/views/LCG_98python3/x86_64-centos7-gcc8-opt/setup.sh

else
    echo "Run on SLC 6 or Cent OS 7"
    exit
fi

sample_list_txt=${1}
year=${2}
eos_output_dir=${3}

root -q -n -b vbs_flat_ntupler.cc+\(\"${sample_list_txt}\",${year}\)
rm *.pcm
rm *.d
rm *.so

if [[ ! -z "${eos_output_dir}" ]]; then
    echo ${eos_output_dir}
    xrdfs root://cmseos.fnal.gov/ mkdir -p ${eos_output_dir}
    xrdcp -f *.root root://cmseos.fnal.gov/${eos_output_dir}/
    rm *.root
    rm *docker*
fi
