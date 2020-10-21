#!/bin/bash
shopt -s expand_aliases

# args
directory=${1}
year=${2}

alias eosfind="eos root://cmseos.fnal.gov find"
listOfSamples=$(eosfind -d ${directory} | cut -d"=" -f2- | grep -E "[0-9]{6}_[0-9]{6}/$")

echo "Making \"${year}\" directory for txt files"
mkdir -p ${year}

for sample in ${listOfSamples}
do
  sampleName=$(echo $sample | sed "s|$directory||" | sed -E 's/(\/[0-9]{6}_[0-9]{6}\/)$//' | sed 's/\//_/g')
  sampleName=${sampleName//\//_}
  nSubDirs=$(eosfind -d ${sample} | grep -E "[0-9]{6}_[0-9]{6}/$" | wc -l)

  if [[ $nSubDirs != "1" ]]
  then
    echo "This $sampleName has multiple sub-directories"
  else
    if [[ -f ${year}/${sampleName}.txt ]]
    then
      echo "Replacing ${year}/${sampleName}.txt"
    fi
    write=$(eosfind -name *.root ${sample} > ${year}/${sampleName}.txt)
  fi

done
