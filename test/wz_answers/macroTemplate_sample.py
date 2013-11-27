#!/usr/bin/env python
import math
import ROOT
from ROOT import TFile

#wrapper to help the analysis of multiple files
from UserCode.sm_cms_das.Tools import getFilesForProcess

## Run by calling in the shell:
## python macroTemplate_sample.py

#get the files for a sample
files=getFilesForProcess(jsonUrl="../wz/wz_samples.json",tag="WplusToMuNu",inDir="../../results") #/store/cmst3/user/psilva/CMSDAS_v2/summary/")

#init variables
lumi=1
xsec=0
totalSelected=0
totalGenerated=0

#loop over files
for f in files:
    inF=TFile.Open(f)
    tree=inF.Get("data/data")
    totalSelected=totalSelected+tree.GetEntries("cat==-13")
    totalGenerated=totalGenerated+inF.Get("iniEvents")[0]
    xsec=inF.Get("crossSection")[0]
    inF.Close()

#compute expected number of events and show it
nExpected=xsec*lumi*totalSelected/totalGenerated
nExpectedErr=xsec*lumi*math.sqrt(totalSelected)/totalGenerated
print "N=%3.3f +/- %3.3f"%(nExpected,nExpectedErr)
