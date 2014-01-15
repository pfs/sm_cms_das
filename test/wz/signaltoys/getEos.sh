#!/bin/bash

source ~/setup.sh
echo $CMSSW_BASE
dir=$1

for x in ` /afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select ls -l  eos/cms/store/cmst3/user/psilva/CMSDAS_v2/summary   |  awk '{print $9}'`; do
cmsStageIn  /store/cmst3/user/psilva/CMSDAS_v2/summary/${x} .
done

for x in `ls *.root | sed "s@_@ @g" | awk '{print $1}' | uniq`; do 
    hadd ${x}.root $x*.root
done