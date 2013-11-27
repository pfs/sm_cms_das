import math,ROOT
from UserCode.sm_cms_das.Tools import *

h=ROOT.TH1F('leg1pt',';Transverse momentum [GeV];Events',50,0,250)
h.Sumw2()
fOut=ROOT.TFile('plotter.root','RECREATE')

processes=getProcesses(jsonUrl="test/wz/wz_samples.json",getMC=True,getData=True)

for idx,p in enumerate(processes):
    files=getFilesForProcess(jsonUrl="test/wz/wz_samples.json",tag=p,inDir="/store/cmst3/user/psilva/CMSDAS_v2/summary/")
    lumi=19
    xsec=0
    totalSelected=0
    totalGenerated=0

    hproc=h.Clone(p+'_leg1pt')
    hproc.SetTitle(p)
    hproc.SetFillStyle(1001)
    hproc.SetFillColor(idx+1)
    hproc.SetDirectory(0)
    
    for f in files:
        inF=ROOT.TFile.Open(f)
        tree=inF.Get("data/data")
        totalSelected=totalSelected+tree.GetEntries("abs(cat)==13")
        totalGenerated=totalGenerated+inF.Get("iniEvents")[0]
        xsec=inF.Get("crossSection")[0]
        h.SetDirectory(inF)
        tree.Draw('leg1_pt>>leg1pt','abs(cat)==13','goff')
        hproc.Add(h)
        h.SetDirectory(0)
        inF.Close()
    if xsec>=0 :
        nExpected=xsec*lumi*totalSelected/totalGenerated
        nExpectedErr=xsec*lumi*math.sqrt(totalSelected)/totalGenerated
        hproc.Scale(xsec*lumi/totalGenerated)
    else :
        nExpected=totalSelected
        nExpectedErr=0
    print "N(%s)=%3.3f +/- %3.3f"%(p,nExpected,nExpectedErr)

    fOut.mkdir(p).cd()
    hproc.Write()

fOut.Close()
