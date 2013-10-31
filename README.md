SM exercise for CMS DAS
=======================

Register to git at https://github.com/ 
-------------------------------------

git config --global user.github <your github username>

Installation
------------

export SCRAM_ARCH=slc5_amd64_gcc462

cmsrel CMSSW_5_3_12_patch2

cd CMSSW_5_3_12_patch2/src/

cmsenv

wget -q -O - --no-check-certificate https://raw.github.com/pfs/sm_cms_das/master/TAGS.txt | sh

Creating ntuples
----------------

Analysis scripts
----------------


Authors
-------