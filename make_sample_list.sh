#!/bin/bash
shopt -s expand_aliases

directory=${1}
alias eosfind="eos root://cmseos.fnal.gov find"
listOfSamples=$(eosfind -d --maxdepth 1 ${directory})

for d in ${listOfSamples}; do
    fixed_d=${d//path=/}
    write=$(eosfind -name *.root ${fixed_d} > $(basename ${fixed_d}).txt)
done
