import math,ROOT
from UserCode.sm_cms_das.Tools import *

h=ROOT.TH1F('leg1pt',';Transverse momentum [GeV];Events',50,0,250)
h.Sumw2()
mcstack=ROOT.THStack('mcstack','mcstack')

processes=getProcesses(jsonUrl="test/wz/wz_samples.json")

for idx,p in enumerate(processes):
    files=getFilesForProcess(jsonUrl="test/wz/wz_samples.json",tag=p,inDir="/store/cmst3/user/psilva/CMSDAS_v2/summary/")
    lumi=1
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
        print h.Integral(),hproc.Integral()
        h.SetDirectory(0)
        inF.Close()
        nExpected=xsec*lumi*totalSelected/totalGenerated
        nExpectedErr=xsec*lumi*math.sqrt(totalSelected)/totalGenerated
    print "N(%s)=%3.3f +/- %3.3f"%(p,nExpected,nExpectedErr)

    hproc.Scale(xsec*lumi/totalGenerated)
    mcstack.Add(hproc,"hist")

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
c=ROOT.TCanvas('c','c',500,500)
mcstack.Draw()
mcstack.GetXaxis().SetTitle(h.GetXaxis().GetTitle())
mcstack.GetYaxis().SetTitle(h.GetYaxis().GetTitle())
c.BuildLegend()
c.SaveAs('leg1pt.png')
