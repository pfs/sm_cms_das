#!/bin/bash

eval `scramv1 runtime -sh`
export RIVET_ANALYSIS_PATH=./
export RIVET_REF_PATH=./

rm RivetCMS_2012_I941555.so
g++ -o "RivetCMS_2012_I941555.so" -shared -fPIC -I/afs/cern.ch/cms/slc5_amd64_gcc472/external/rivet/1.8.2-cms8/include -I/afs/cern.ch/cms/slc5_amd64_gcc472/external/hepmc/2.06.07/include -I/afs/cern.ch/cms/slc5_amd64_gcc472/external/fastjet/3.0.1/include -I/afs/cern.ch/cms/slc5_amd64_gcc472/external/gsl/1.10/include -I/afs/cern.ch/cms/slc5_amd64_gcc472/external/boost/1.51.0-cms2/include -pedantic -ansi -Wall -Wno-long-long -O2 -std=c++0x CMS_2012_I941555.cc


# rivet -q --analysis=CMS_2012_I941555 /nfsdisk/perrozzi/CMSSW_4_4_5/src/WMassNNLO/resbos/a1_scan_kc1/resbos_all.hep

# mv Rivet.aida rivet-CMS-EWK-10-010-2.aida

# rivet-mkhtml rivet-CMS-EWK-10-010-2.aida


# eval `scramv1 runtime -sh`; export RIVET_ANALYSIS_PATH=./; export RIVET_REF_PATH=./; rm RivetCMS_2012_I941555.so; g++ -o "RivetCMS_2012_I941555.so" -shared -fPIC -I/afs/cern.ch/cms/slc5_amd64_gcc472/external/rivet/1.8.2-cms8/include -I/afs/cern.ch/cms/slc5_amd64_gcc472/external/hepmc/2.06.07/include -I/afs/cern.ch/cms/slc5_amd64_gcc472/external/fastjet/3.0.1/include -I/afs/cern.ch/cms/slc5_amd64_gcc472/external/gsl/1.10/include -I/afs/cern.ch/cms/slc5_amd64_gcc472/external/boost/1.51.0-cms2/include -pedantic -ansi -Wall -Wno-long-long -O2 -std=c++0x CMS_2012_I941555.cc; 

# eval `scramv1 runtime -sh`; export RIVET_ANALYSIS_PATH=./; export RIVET_REF_PATH=./; rm RivetCMS_2012_I1107658.so; g++ -o "RivetCMS_2012_I1107658.so" -shared -fPIC -I/afs/cern.ch/cms/slc5_amd64_gcc472/external/rivet/1.8.2-cms8/include -I/afs/cern.ch/cms/slc5_amd64_gcc472/external/hepmc/2.06.07/include -I/afs/cern.ch/cms/slc5_amd64_gcc472/external/fastjet/3.0.1/include -I/afs/cern.ch/cms/slc5_amd64_gcc472/external/gsl/1.10/include -I/afs/cern.ch/cms/slc5_amd64_gcc472/external/boost/1.51.0-cms2/include -pedantic -ansi -Wall -Wno-long-long -O2 -std=c++0x CMS_2012_I1107658.cc; 
# rivet --analysis=CMS_2012_I941555 /afs/cern.ch/work/p/perrozzi/private/WMassMC/resbosPN/test/resbos.hep -q; rivet-mkhtml Rivet.aida

root2flat ../vpt_unfolded.root
mv dataunfolded.dat d02-x01-y03.dat
flat2aida d02-x01-y03.dat
rm *.dat
mv d02-x01-y03.aida CMS_2012_I941555.aida

sed -i 's/dataunfolded/d02-x01-y03/g' CMS_2012_I941555.aida
sed -i 's/name\=\"d02\-x01\-y03\"\ dimension\=\"2\"/name\=\"d02\-x01\-y03\"\ dimension\=\"2\"\ path\=\"\/REF\/CMS_2012_I941555\"\ title\=\"(1\/SIG)*D(SIG)\/DQT\ IN\ GEV**\-1\ \"\>/g' CMS_2012_I941555.aida
sed -i 's/.*title="d02-x01-y03".*/\ /' CMS_2012_I941555.aida
  
# exit 1  


cmsRun rivetPythiaZ2NoFSRM50To130_cfg.py
 
rivet-mkhtml DYtoMuMu-Z2-pythia6.aida
