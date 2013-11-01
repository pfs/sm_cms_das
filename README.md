SM exercise for CMS DAS
=======================

Before starting
---------------

register to git at https://github.com/ 

git config --global user.github <your github username>

Installation
------------

export SCRAM_ARCH=slc5_amd64_gcc462

cmsrel CMSSW_5_3_12_patch2

cd CMSSW_5_3_12_patch2/src/

cmsenv

wget -q -O - --no-check-certificate https://raw.github.com/pfs/sm_cms_das/master/TAGS.txt | sh

scram b -j 9

Creating ntuples
----------------

cmsRun test/runSManalyzer_data_cfg.py

cmsRun test/runSManalyzer_mc_cfg.py

or use the crab/multicrab files under test/grid

Analysis scripts
----------------

python test/wz/runEventSelection.py SManalysis.root

Authors
-------