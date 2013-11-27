import math,ROOT
from UserCode.sm_cms_das.Tools import getFilesForProcess
files=getFilesForProcess(jsonUrl="test/wz/wz_samples.json",tag="WplusToMuNu",inDir="/store/cmst3/user/psilva/CMSDAS_v2/summary/")
lumi=1
xsec=0
totalSelected=0
totalGenerated=0
for f in files:
    inF=ROOT.TFile.Open(f)
    tree=inF.Get("data/data")
    totalSelected=totalSelected+tree.GetEntries("cat==-13")
    totalGenerated=totalGenerated+inF.Get("iniEvents")[0]
    xsec=inF.Get("crossSection")[0]
    inF.Close()
    nExpected=xsec*lumi*totalSelected/totalGenerated
    nExpectedErr=xsec*lumi*math.sqrt(totalSelected)/totalGenerated
print "N=%3.3f +/- %3.3f"%(nExpected,nExpectedErr)
